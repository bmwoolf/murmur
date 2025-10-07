import os, sys, subprocess, tempfile, shutil, gzip, csv
from pathlib import Path

"""
Purpose of this file:
- take a raw (or gz) VCF and annotate with VEP and SnpEff (GRCh37/GRCh38)
- standardize to a tidy TSV with: CHROM, POS, REF, ALT, SYMBOL (gene),
  Consequence, IMPACT, Transcript, Protein_change, dbSNP, AF fields (if available)
- optionally emit a simple gene→dose CSV for the CPA step

Prereqs (choose one path):
1) Conda-installed tools + local caches:
   - ensembl-vep, htslib, perl, snpeff
   - VEP cache at ~/.vep or --dir_cache
   - SnpEff data at $SNPEFF/data (e.g., GRCh38.105)
2) Docker images:
   - ensemblorg/ensembl-vep
   - pcingola/snpeff
"""

class Config:
    reference = "GRCh38"
    threads = 20  # CPU threads available- there isnt any GPU-based VEP or SnpEff alternatives
    vep_cache_dir = str(Path.home() / ".vep")
    snpeff_data_dir = os.environ.get("SNPEFF_DATA_DIR", "/usr/local/share/snpeff/data")
    
    # common VEP plugins you might enable later: LoF, CADD, dbNSFP, gnomAD
    vep_plugins = []               # ie ["LoF,loftee_path:/path/to/loftee, ..."]
    
    # output columns we will keep in the final TSV:
    out_cols = [
        "CHROM","POS","ID","REF","ALT",
        "SYMBOL","Consequence","IMPACT","Feature","HGVSp","HGVSc",
        "BIOTYPE","SIFT","PolyPhen",
        "gnomADg_AF","MAX_AF","CLIN_SIG"
    ]

# mapping heuristic for later dose estimation (used elsewhere)
IMPACT_TO_DOSE = {"HIGH":0.8, "MODERATE":0.5, "LOW":0.2, "MODIFIER":0.1}

# utils
def ensure_tools_available():
    """Check VEP and SnpEff are callable; raise with hints if not."""
    for tool in ["vep", "snpeff"]:
        try:
            subprocess.run([tool, "-h"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
        except Exception:
            raise RuntimeError(f"{tool} not found. Install via conda or use docker wrappers.")

def detect_reference_build(vcf_path: str) -> str:
    """
    Peek VCF header for '##reference=' or '##contig=' hints.
    Fallback to Config.reference.
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if not line.startswith("##"):
                break
            if "GRCh37" in line or "b37" in line:
                return "GRCh37"
            if "GRCh38" in line or "hg38" in line:
                return "GRCh38"
    return Config.reference

def bgzip_and_index(vcf_in: str) -> str:
    """
    Ensure input is bgzipped & indexed (tabix).
    Returns path to .vcf.gz
    """
    if vcf_in.endswith(".vcf.gz"):
        return vcf_in
    gz = vcf_in + ".gz"
    subprocess.run(["bgzip","-c",vcf_in], check=True, stdout=open(gz,"wb"))
    subprocess.run(["tabix","-p","vcf",gz], check=True)
    return gz


# VEP
def run_vep(vcf_gz: str, out_tsv: str, reference: str):
    """
    Run VEP to TSV with useful fields. Uses cache; offline mode preferred.
    """
    assembly_flag = "--assembly GRCh38" if reference=="GRCh38" else "--assembly GRCh37"
    fields = ",".join([
        "Uploaded_variation","Location","Allele","Gene","SYMBOL","Feature","BIOTYPE",
        "Consequence","IMPACT","HGVSc","HGVSp","SIFT","PolyPhen",
        "MAX_AF","gnomADg_AF","CLIN_SIG"
    ])
    plugin_args = []
    for p in Config.vep_plugins:
        plugin_args += ["--plugin", p]

    cmd = [
        "vep",
        "--input_file", vcf_gz,
        "--format", "vcf",
        "--vcf",
        "--output_file", "STDOUT",
        "--cache", "--offline", assembly_flag,
        "--dir_cache", Config.vep_cache_dir,
        "--symbol", "--canonical", "--nearest symbol",
        "--fork", str(Config.threads),
        "--force_overwrite",
        "--tab", "--fields", fields
    ] + plugin_args

    # VEP with --tab prints header + rows (tab-separated)
    with open(out_tsv, "w") as fout:
        subprocess.run(cmd, check=True, stdout=fout)


# SnpEff
def snpeff_genome_name(reference: str) -> str:
    """
    Choose appropriate SnpEff genome database name.
    Examples: 'GRCh38.105' or 'GRCh37.75' depending on what you've downloaded.
    """
    # adjust to your local installation
    return "GRCh38.105" if reference=="GRCh38" else "GRCh37.75"

def run_snpeff(vcf_gz: str, out_vcf: str, reference: str):
    """
    Run SnpEff to annotate VCF with ANN field.
    """
    genome = snpeff_genome_name(reference)
    cmd = [
        "snpEff",
        "-Xmx8g",
        genome,
        vcf_gz
    ]
    with open(out_vcf, "w") as fout:
        subprocess.run(cmd, check=True, stdout=fout)

def snpeff_to_tsv(ann_vcf: str, out_tsv: str):
    """
    Parse ANN field to a tidy TSV with core columns (SYMBOL, Consequence, IMPACT, etc.).
    In practice, use a VCF parser (cyvcf2) and split ANN pipe-delimited fields.
    """
    # Pseudocode:
    # with cyvcf2.VCF(ann_vcf) as v:
    #   write header
    #   for rec in v:
    #       for ann in rec.INFO.get("ANN","").split(","):
    #           fields = ann.split("|")
    #           SYMBOL = fields[3]; Consequence = fields[1]; IMPACT = fields[2]; Feature = fields[6]
    #           HGVSc = fields[9]; HGVSp = fields[10]; BIOTYPE = fields[7]
    #           write row with CHROM, POS, ID, REF, ALT + these fields
    pass

# merge + standardize
def merge_vep_snpeff(vep_tsv: str, snpeff_tsv: str, out_tsv: str):
    """
    Join on (CHROM, POS, REF, ALT) and prefer enriched fields when duplicated.
    Strategy:
      - Start from VEP rows as base.
      - If SnpEff has a stronger IMPACT for same variant, keep that.
      - Fill missing SYMBOL/HGVSp from the other tool.
    """
    # Use pandas merge, custom resolver for conflicts, write out `out_tsv` with Config.out_cols
    pass

def write_gene_dose_csv(annot_tsv: str, out_csv: str):
    """
    Collapse variant-level rows to gene-level dose estimates for CPA.
    - For each gene: dose = max(IMPACT_TO_DOSE[IMPACT]) across its variants.
    """
    # pd.read_csv(annot_tsv, sep="\t") -> groupby("SYMBOL") -> map IMPACT->dose -> max -> write CSV
    pass

# orchestration
def annotate_vcf(vcf_path: str, out_dir: str) -> dict:
    """
    Full flow:
      1) Ensure bgzip+tabix.
      2) Detect reference build (or use config).
      3) Run VEP -> vep.tsv
      4) Run SnpEff -> snpeff.vcf; parse -> snpeff.tsv
      5) Merge -> annotated.tsv (canonical tidy table)
      6) Optional: gene_dose.csv for CPA
    Returns paths.
    """
    ensure_tools_available()
    os.makedirs(out_dir, exist_ok=True)
    vcf_gz = bgzip_and_index(vcf_path)
    ref = detect_reference_build(vcf_gz)

    vep_tsv = str(Path(out_dir) / "vep.tsv")
    run_vep(vcf_gz, vep_tsv, ref)

    snp_vcf = str(Path(out_dir) / "snpeff.vcf")
    run_snpeff(vcf_gz, snp_vcf, ref)

    snp_tsv = str(Path(out_dir) / "snpeff.tsv")
    snpeff_to_tsv(snp_vcf, snp_tsv)

    merged = str(Path(out_dir) / "annotated.tsv")
    merge_vep_snpeff(vep_tsv, snp_tsv, merged)

    gene_dose = str(Path(out_dir) / "gene_dose.csv")
    write_gene_dose_csv(merged, gene_dose)

    return {"vep_tsv": vep_tsv, "snpeff_vcf": snp_vcf, "snpeff_tsv": snp_tsv, "annotated_tsv": merged, "gene_dose_csv": gene_dose}
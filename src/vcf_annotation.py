import os, sys, subprocess, tempfile, shutil, gzip, csv, yaml
from pathlib import Path
import cyvcf2

"""
Purpose of this file:
- take a raw (or gz) VCF and annotate with VEP and SnpEff (GRCh37/GRCh38)
- standardize to a tidy TSV with: CHROM, POS, REF, ALT, SYMBOL (gene),
  Consequence, IMPACT, Transcript, Protein_change, dbSNP, AF fields (if available)
- optionally emit a simple gene→dose CSV for the CPA step

Prereqs:
- Docker engine running locally
- Docker images:
   - ensemblorg/ensembl-vep
   - staphb/snpeff
- Local cache directories for VEP and SnpEff databases
"""

def load_config():
    """Load configuration from config.yml"""
    config_path = Path(__file__).parent.parent / "config.yml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

class Config:
    _config = load_config()
    
    reference = _config['system']['reference_genome']
    threads = _config['pipeline']['vcf_annotation']['threads']
    vep_cache_dir = _config['pipeline']['vcf_annotation']['vep_cache_dir']
    snpeff_data_dir = _config['pipeline']['vcf_annotation']['snpeff_data_dir']
    
    # Docker image names
    vep_docker_image = _config['pipeline']['vcf_annotation']['vep_docker_image']
    snpeff_docker_image = _config['pipeline']['vcf_annotation']['snpeff_docker_image']
    
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
    """Check Docker is available and pull required images if needed."""
    # Check Docker is running
    try:
        subprocess.run(["docker", "info"], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    except Exception:
        raise RuntimeError("Docker is not running or not available. Please start Docker and try again.")
    
    # Pull Docker images if they don't exist locally
    images_to_pull = [
        Config.vep_docker_image,
        Config.snpeff_docker_image
    ]
    
    for image in images_to_pull:
        try:
            # Check if image exists locally
            result = subprocess.run(
                ["docker", "images", "-q", image], 
                capture_output=True, 
                text=True, 
                check=True
            )
            if not result.stdout.strip():
                print(f"Pulling Docker image: {image}")
                subprocess.run(["docker", "pull", image], check=True)
                print(f"Successfully pulled {image}")
            else:
                print(f"Using existing image: {image}")
        except Exception as e:
            raise RuntimeError(f"Failed to pull Docker image {image}: {e}")

def setup_databases():
    """Setup local cache directories and download databases using Docker."""
    config = load_config()
    reference = config['system']['reference_genome']
    
    # Setup VEP cache directory
    vep_cache_dir = os.path.expanduser(Config.vep_cache_dir)
    os.makedirs(vep_cache_dir, exist_ok=True)
    
    # Setup SnpEff data directory
    snpeff_data_dir = os.path.expanduser(Config.snpeff_data_dir)
    os.makedirs(snpeff_data_dir, exist_ok=True)
    
    # Check if VEP cache exists, if not download using Docker
    if not os.path.exists(os.path.join(vep_cache_dir, "homo_sapiens")):
        print("Downloading VEP cache using Docker...")
        try:
            subprocess.run([
                "docker", "run", "--rm",
                "-v", f"{vep_cache_dir}:/opt/vep/.vep",
                Config.vep_docker_image,
                "bash", "-c", "vep_install -a cf -s homo_sapiens -y GRCh38 -c /opt/vep/.vep"
            ], check=True)
            print("VEP cache downloaded successfully!")
        except Exception as e:
            print(f"Warning: Failed to download VEP cache: {e}")
            print("VEP will run in online mode (slower)")
    else:
        print("VEP cache already exists")
    
    # Check if SnpEff database exists, if not download using Docker
    genome_name = snpeff_genome_name(reference)
    db_path = os.path.join(snpeff_data_dir, genome_name)
    
    if not os.path.exists(db_path):
        print(f"Downloading SnpEff database {genome_name} using Docker...")
        try:
            subprocess.run([
                "docker", "run", "--rm",
                "-v", f"{snpeff_data_dir}:/data",
                Config.snpeff_docker_image,
                "snpEff", "download", "-v", genome_name
            ], check=True)
            print("SnpEff database downloaded successfully!")
        except Exception as e:
            print(f"Warning: Failed to download SnpEff database: {e}")
            print("SnpEff may not work properly without the database")
    else:
        print("SnpEff database already exists")


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
    Run VEP using Docker container with useful fields. Uses cache; offline mode preferred.
    """
    assembly_flag = "--assembly=GRCh38" if reference=="GRCh38" else "--assembly=GRCh37"
    fields = ",".join([
        "Uploaded_variation","Location","Allele","Gene","SYMBOL","Feature","BIOTYPE",
        "Consequence","IMPACT","HGVSc","HGVSp","SIFT","PolyPhen",
        "MAX_AF","gnomADg_AF","CLIN_SIG"
    ])
    plugin_args = []
    for p in Config.vep_plugins:
        plugin_args += ["--plugin", p]

    # Get absolute paths for Docker volume mounting
    vcf_abs_path = os.path.abspath(vcf_gz)
    vcf_dir = os.path.dirname(vcf_abs_path)
    vcf_filename = os.path.basename(vcf_abs_path)
    out_abs_path = os.path.abspath(out_tsv)
    out_dir = os.path.dirname(out_abs_path)
    out_filename = os.path.basename(out_abs_path)
    vep_cache_abs_path = os.path.expanduser(Config.vep_cache_dir)

    # VEP command within Docker container (using database since cache setup failed)
    vep_cmd = [
        "vep",
        "--input_file", f"/data/input/{vcf_filename}",
        "--format", "vcf",
        "--database", assembly_flag,  # --assembly GRCh38
        "--symbol", "--canonical", "--nearest", "symbol",
        "--fork", str(Config.threads),
        "--force_overwrite",
        "--no_stats",
        "--tab", "--fields", fields
    ] + plugin_args

    # Docker run command
    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{vcf_dir}:/data/input",
        Config.vep_docker_image
    ] + vep_cmd

    print(f"Running VEP with Docker: {Config.vep_docker_image}")
    with open(out_tsv, "w") as fout:
        subprocess.run(docker_cmd, check=True, stdout=fout)


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
    Run SnpEff using Docker container to annotate VCF with ANN field.
    """
    genome = snpeff_genome_name(reference)
    
    # Get absolute paths for Docker volume mounting
    vcf_abs_path = os.path.abspath(vcf_gz)
    vcf_dir = os.path.dirname(vcf_abs_path)
    vcf_filename = os.path.basename(vcf_abs_path)
    out_abs_path = os.path.abspath(out_vcf)
    out_dir = os.path.dirname(out_abs_path)
    out_filename = os.path.basename(out_abs_path)
    snpeff_data_abs_path = os.path.expanduser(Config.snpeff_data_dir)

    # SnpEff command within Docker container
    snpeff_cmd = [
        "snpEff",
        "-Xmx8g",
        "-dataDir", "/data/snpeff_data",
        genome,
        f"/data/input/{vcf_filename}"
    ]

    # Docker run command
    docker_cmd = [
        "docker", "run", "--rm",
        "-v", f"{vcf_dir}:/data/input",
        "-v", f"{out_dir}:/data/output",
        "-v", f"{snpeff_data_abs_path}:/data/snpeff_data",
        Config.snpeff_docker_image
    ] + snpeff_cmd

    print(f"Running SnpEff with Docker: {Config.snpeff_docker_image}")
    with open(out_vcf, "w") as fout:
        subprocess.run(docker_cmd, check=True, stdout=fout)

def snpeff_to_tsv(ann_vcf: str, out_tsv: str):
    """
    Parse ANN field to a tidy TSV with core columns (SYMBOL, Consequence, IMPACT, etc.).
    Uses cyvcf2 to parse VCF and extract SnpEff annotations.
    """

    try:
        with open(out_tsv, 'w') as fout:
            # write header
            header = ["CHROM", "POS", "ID", "REF", "ALT", "SYMBOL", "Consequence", "IMPACT", 
                     "Feature", "HGVSc", "HGVSp", "BIOTYPE", "SIFT", "PolyPhen"]
            fout.write("\t".join(header) + "\n")
            
            # parse VCF and extract annotations
            with cyvcf2.VCF(ann_vcf) as v:
                for record in v:
                    chrom = record.CHROM
                    pos = record.POS
                    vid = record.ID if record.ID else "."
                    ref = record.REF
                    
                    # process each alternate allele
                    for i, alt in enumerate(record.ALT):
                        alt_str = str(alt)
                        
                        # get ANN field (SnpEff annotations)
                        ann_field = record.INFO.get("ANN", "")
                        if not ann_field:
                            # no annotation, write basic info
                            row = [chrom, str(pos), vid, ref, alt_str, ".", ".", ".", 
                                   ".", ".", ".", ".", ".", "."]
                            fout.write("\t".join(row) + "\n")
                            continue
                        
                        # parse each annotation (comma-separated)
                        for ann in ann_field.split(","):
                            if not ann.strip():
                                continue
                                
                            fields = ann.split("|")
                            if len(fields) < 11:
                                continue
                                
                            # SnpEff ANN field structure:
                            # 0: Allele, 1: Annotation, 2: Impact, 3: Gene Name, 4: Gene ID,
                            # 5: Feature Type, 6: Feature ID, 7: Transcript Biotype,
                            # 8: Rank/Total, 9: HGVS.c, 10: HGVS.p, 11: cDNA_position,
                            # 12: CDS_position, 13: Protein_position, 14: Distance_to_feature,
                            # 15: Errors/Warnings/Info
                            
                            symbol = fields[3] if fields[3] else "."
                            consequence = fields[1] if fields[1] else "."
                            impact = fields[2] if fields[2] else "."
                            feature = fields[6] if fields[6] else "."
                            hgvsc = fields[9] if fields[9] else "."
                            hgvsp = fields[10] if fields[10] else "."
                            biotype = fields[7] if fields[7] else "."
                            
                            # for SIFT and PolyPhen, we'd need additional SnpEff databases
                            sift = "."
                            polyphen = "."
                            
                            row = [chrom, str(pos), vid, ref, alt_str, symbol, consequence, 
                                   impact, feature, hgvsc, hgvsp, biotype, sift, polyphen]
                            fout.write("\t".join(row) + "\n")
                            
    except Exception as e:
        raise RuntimeError(f"Failed to parse SnpEff VCF {ann_vcf}: {e}")

# merge + standardize
def merge_vep_snpeff(vep_tsv: str, snpeff_tsv: str, out_tsv: str):
    """
    Join on (CHROM, POS, REF, ALT) and prefer enriched fields when duplicated.
    Strategy:
      - Start from VEP rows as base.
      - If SnpEff has a stronger IMPACT for same variant, keep that.
      - Fill missing SYMBOL/HGVSp from the other tool.
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for merging annotations. Install with: pip install pandas")
    
    # Read both annotation files
    vep_df = pd.read_csv(vep_tsv, sep="\t", dtype=str)
    snpeff_df = pd.read_csv(snpeff_tsv, sep="\t", dtype=str)
    
    # Ensure consistent column names
    common_cols = ["CHROM", "POS", "REF", "ALT"]
    merge_cols = common_cols + ["SYMBOL", "Consequence", "IMPACT", "Feature", "HGVSc", "HGVSp", "BIOTYPE"]
    
    # Filter to common columns that exist in both
    vep_cols = [col for col in merge_cols if col in vep_df.columns]
    snpeff_cols = [col for col in merge_cols if col in snpeff_df.columns]
    
    vep_subset = vep_df[vep_cols].copy()
    snpeff_subset = snpeff_df[snpeff_cols].copy()
    
    # merge on variant coordinates
    merged = pd.merge(vep_subset, snpeff_subset, 
                     on=common_cols, 
                     how='outer', 
                     suffixes=('_vep', '_snpeff'))
    
    # resolve conflicts - prefer higher impact
    impact_order = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, ".": 0}

def resolve_impact(row):
    vep_impact = row.get('IMPACT_vep', '.')
    snpeff_impact = row.get('IMPACT_snpeff', '.')
    
    if impact_order.get(vep_impact, 0) >= impact_order.get(snpeff_impact, 0):
        return row.get('SYMBOL_vep', ''), row.get('Consequence_vep', ''), vep_impact, \
                row.get('Feature_vep', ''), row.get('HGVSc_vep', ''), row.get('HGVSp_vep', ''), \
                row.get('BIOTYPE_vep', '')
    else:
        return row.get('SYMBOL_snpeff', ''), row.get('Consequence_snpeff', ''), snpeff_impact, \
                row.get('Feature_snpeff', ''), row.get('HGVSc_snpeff', ''), row.get('HGVSp_snpeff', ''), \
                row.get('BIOTYPE_snpeff', '')

    try:
        # apply resolution
        resolved = merged.apply(resolve_impact, axis=1, result_type='expand')
        resolved.columns = ['SYMBOL', 'Consequence', 'IMPACT', 'Feature', 'HGVSc', 'HGVSp', 'BIOTYPE']
        
        # combine with variant coordinates
        final_df = pd.concat([
            merged[common_cols].reset_index(drop=True),
            resolved.reset_index(drop=True)
        ], axis=1)
        
        # fill missing values
        final_df = final_df.fillna('.')
        
        # select output columns as defined in Config
        output_cols = [col for col in Config.out_cols if col in final_df.columns]
        final_output = final_df[output_cols]
        
        # write to output file
        final_output.to_csv(out_tsv, sep="\t", index=False)
        
    except Exception as e:
        raise RuntimeError(f"Failed to merge VEP and SnpEff annotations: {e}")

def write_gene_dose_csv(annot_tsv: str, out_csv: str):
    """
    Collapse variant-level rows to gene-level dose estimates for CPA.
    - For each gene: dose = max(IMPACT_TO_DOSE[IMPACT]) across its variants.
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas is required for gene dose calculation. Install with: pip install pandas")
    
    try:
        # read annotated variants
        df = pd.read_csv(annot_tsv, sep="\t", dtype=str)
        
        # filter to variants with gene symbols
        df = df[df['SYMBOL'] != '.']
        
        if df.empty:
            # create empty file with headers
            pd.DataFrame(columns=['gene', 'dose']).to_csv(out_csv, index=False)
            return
        
        # map impact to dose
        def impact_to_dose(impact):
            return IMPACT_TO_DOSE.get(impact, 0.1)
        
        df['dose'] = df['IMPACT'].apply(impact_to_dose)
        
        # group by gene and take maximum dose
        gene_doses = df.groupby('SYMBOL')['dose'].max().reset_index()
        gene_doses.columns = ['gene', 'dose']
        
        # sort by dose descending
        gene_doses = gene_doses.sort_values('dose', ascending=False)
        
        # write to CSV
        gene_doses.to_csv(out_csv, index=False)
        
    except Exception as e:
        raise RuntimeError(f"Failed to generate gene dose CSV: {e}")

# orchestration
def annotate_vcf(vcf_path: str, out_dir: str) -> dict:
    """
    Full flow:
      1) Ensure tools and databases are available
      2) Ensure bgzip+tabix.
      3) Run VEP -> vep.tsv
      4) Run SnpEff -> snpeff.vcf; parse -> snpeff.tsv
      5) Merge -> annotated.tsv (canonical tidy table)
      6) Optional: gene_dose.csv for CPA
    Returns paths.
    """
    ensure_tools_available()
    setup_databases() # VEP and SnpEff
    os.makedirs(out_dir, exist_ok=True)
    vcf_gz = bgzip_and_index(vcf_path)
    ref = Config.reference

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
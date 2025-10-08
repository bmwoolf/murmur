import os, sys, subprocess, tempfile, shutil, gzip, csv, yaml
from pathlib import Path
import cyvcf2

"""
Purpose of this file:
- take a raw (or gz) VCF and annotate with VEP (GRCh37/GRCh38)
- standardize to a tidy TSV with: CHROM, POS, REF, ALT, SYMBOL (gene),
  Consequence, IMPACT, Transcript, Protein_change, dbSNP, AF fields (if available)
- optionally emit a simple gene→dose CSV for the CPA step

Prereqs:
- Docker engine running locally
- Docker image: ensemblorg/ensembl-vep
- Local cache directories for VEP databases
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
        Config.vep_docker_image
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
      3) Run VEP -> annotated.tsv (canonical tidy table)
      4) Optional: gene_dose.csv for CPA
    Returns paths.
    """
    ensure_tools_available()
    setup_databases() # VEP only
    os.makedirs(out_dir, exist_ok=True)
    vcf_gz = bgzip_and_index(vcf_path)
    ref = Config.reference

    # Run VEP and use its output as the final annotated file
    annotated_tsv = str(Path(out_dir) / "annotated.tsv")
    run_vep(vcf_gz, annotated_tsv, ref)

    gene_dose = str(Path(out_dir) / "gene_dose.csv")
    write_gene_dose_csv(annotated_tsv, gene_dose)

    return {"annotated_tsv": annotated_tsv, "gene_dose_csv": gene_dose}
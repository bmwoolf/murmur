# from .ingest_xatlas import load_subset
# from .train_cpa import CPA, prepare_labels
# from .vcf_map import vcf_to_gene_dose
# from .run_cpa import simulate_transcriptome
# from .pathdnn import load_gene_sets, pathway_scores
# from .phenotype import trait_scores_from_gwas
from .vcf_annotation import annotate_vcf
import yaml
import os
from pathlib import Path

def load_config():
    """Load configuration from config.yml"""
    config_path = Path(__file__).parent.parent / "config.yml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

# run full pipeline
def run():
    """
    Run the complete pipeline:
    
    VCF → annotation → perturbation → pathways → phenotype
    
    Uses configuration from config.yml for all file paths and settings.
    
    Returns:
        dict: Results from the annotation step (and eventually full pipeline)
    """
    config = load_config()
    
    # Get paths from config
    vcf_path = config['pipeline']['input']['vcf_path']
    out_dir = config['pipeline']['output']['base_dir']
    xatlas_h5ad = config['pipeline']['input']['xatlas_h5ad']
    msigdb_gmt = config['pipeline']['input']['msigdb_gmt']
    gwas_map = config['pipeline']['input']['gwas_map']
    
    # Validate input files exist
    if not os.path.exists(vcf_path):
        raise FileNotFoundError(f"VCF file not found: {vcf_path}")
    
    # step 1: annotate variants (VEP + SnpEff)
    print("~~~step 1: annotating variants with VEP and SnpEff~~~")
    print(f"   Input VCF: {vcf_path}")
    print(f"   Output directory: {out_dir}")
    annotation_results = annotate_vcf(vcf_path, out_dir)
    consequence_tsv = annotation_results["annotated_tsv"]
    
    # # step 2: load X-atlas data and prepare for training
    # print("~~~step 2: loading X-atlas data~~~")
    # adata = load_subset(xatlas_h5ad)
    # X, g, d, c, le_gene, le_cell = prepare_labels(adata)

    # # step 3: train/load CPA model
    # print("~~~step 3: setting up CPA model~~~")
    # model = CPA(n_genes=X.shape[1],
    #             n_gene_ids=len(le_gene.classes_),
    #             n_celltypes=len(le_cell.classes_))
    # # TODO: train model here (or load checkpoint)

    # # step 4: simulate perturbations from VCF variants
    # print("~~~step 4: simulating transcriptome perturbations~~~")
    # gene_dose = vcf_to_gene_dose(vcf_path, consequence_tsv)
    # expr_pred = simulate_transcriptome(model, gene_dose, le_gene, le_cell, cell_type="HEK293T")

    # # step 5: model pathway disruptions
    # print("~~~step 5: analyzing pathway disruptions~~~")
    # masks, pathway_names = load_gene_sets(msigdb_gmt, gene_index=adata.var_names)
    # p_scores = pathway_scores(expr_pred, masks)

    # # step 6: predict phenotype effects
    # print("~~~step 6: predicting phenotype effects~~~")
    # trait_scores = trait_scores_from_gwas(expr_pred, adata.var_names, gwas_map)
    
    # print("Pipeline completed successfully!")
    # return expr_pred, dict(zip(pathway_names, p_scores)), trait_scores
    
    print("Pipeline completed successfully!")
    print(f"Annotated variants saved to: {consequence_tsv}")
    return annotation_results

if __name__ == "__main__":
    import os
    
    print("Running murmur VCF annotation pipeline...")
    print("Using configuration from config.yml")
    
    try:
        results = run()
        
        print("\nPipeline completed!")
        print("Results:")
        for key, path in results.items():
            if os.path.exists(path):
                size = os.path.getsize(path)
                print(f"   {key}: {path} ({size} bytes)")
            else:
                print(f"   {key}: {path} (not found)")
                
    except Exception as e:
        print(f"\nPipeline failed: {e}")
        print("\nMake sure you have:")
        print("   1. Generated your config.yml: ./config/generate_config.sh")
        print("   2. VEP installed and configured")
        print("   3. SnpEff installed with GRCh38 database")
        print("   4. Required Python packages: pandas, cyvcf2")
        print("   5. Your VCF file path is correct in config.yml")

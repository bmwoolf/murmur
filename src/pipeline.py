from .ingest_xatlas import load_subset
from .train_cpa import CPA, prepare_labels
from .vcf_map import vcf_to_gene_dose
from .run_cpa import simulate_transcriptome
from .pathdnn import load_gene_sets, pathway_scores
from .phenotype import trait_scores_from_gwas
from .vcf_annotation import annotate_vcf

# run full pipeline
def run(vcf_path, out_dir, xatlas_h5ad, msigdb_gmt, gwas_map):
    """
    Run the complete pipeline:
    
    VCF → annotation → perturbation → pathways → phenotype
    
    Args:
        vcf_path: Path to input VCF file
        out_dir: Directory for intermediate and output files
        xatlas_h5ad: Path to X-atlas data file
        msigdb_gmt: Path to MSigDB gene sets file
        gwas_map: Path to GWAS trait mapping file
    
    Returns:
        tuple: (expr_pred, pathway_scores_dict, trait_scores)
    """
    
    # step 1: annotate variants (VEP + SnpEff)
    print("~~~step 1: annotating variants with VEP and SnpEff~~~")
    annotation_results = annotate_vcf(vcf_path, out_dir)
    consequence_tsv = annotation_results["annotated_tsv"]
    
    # step 2: load X-atlas data and prepare for training
    print("~~~step 2: loading X-atlas data~~~")
    adata = load_subset(xatlas_h5ad)
    X, g, d, c, le_gene, le_cell = prepare_labels(adata)

    # step 3: train/load CPA model
    print("~~~step 3: setting up CPA model~~~")
    model = CPA(n_genes=X.shape[1],
                n_gene_ids=len(le_gene.classes_),
                n_celltypes=len(le_cell.classes_))
    # TODO: train model here (or load checkpoint)

    # step 4: simulate perturbations from VCF variants
    print("~~~step 4: simulating transcriptome perturbations~~~")
    gene_dose = vcf_to_gene_dose(vcf_path, consequence_tsv)
    expr_pred = simulate_transcriptome(model, gene_dose, le_gene, le_cell, cell_type="HEK293T")

    # step 5: model pathway disruptions
    print("~~~step 5: analyzing pathway disruptions~~~")
    masks, pathway_names = load_gene_sets(msigdb_gmt, gene_index=adata.var_names)
    p_scores = pathway_scores(expr_pred, masks)

    # step 6: predict phenotype effects
    print("~~~step 6: predicting phenotype effects~~~")
    trait_scores = trait_scores_from_gwas(expr_pred, adata.var_names, gwas_map)
    
    print("Pipeline completed successfully!")
    return expr_pred, dict(zip(pathway_names, p_scores)), trait_scores

from .ingest_xatlas import load_subset
from .train_cpa import CPA, prepare_labels
from .vcf_map import vcf_to_gene_dose
from .run_cpa import simulate_transcriptome
from .pathdnn import load_gene_sets, pathway_scores
from .phenotype import trait_scores_from_gwas

def run(vcf_path, consequence_tsv, xatlas_h5ad, msigdb_gmt, gwas_map):
    adata = load_subset(xatlas_h5ad)
    X, g, d, c, le_gene, le_cell = prepare_labels(adata)

    model = CPA(n_genes=X.shape[1],
                n_gene_ids=len(le_gene.classes_),
                n_celltypes=len(le_cell.classes_))
    # TODO: train model here (or load checkpoint)

    gene_dose = vcf_to_gene_dose(vcf_path, consequence_tsv)
    expr_pred = simulate_transcriptome(model, gene_dose, le_gene, le_cell, cell_type="HEK293T")

    masks, pathway_names = load_gene_sets(msigdb_gmt, gene_index=adata.var_names)
    p_scores = pathway_scores(expr_pred, masks)

    trait_scores = trait_scores_from_gwas(expr_pred, adata.var_names, gwas_map)
    return expr_pred, dict(zip(pathway_names, p_scores)), trait_scores

import scanpy as sc
import anndata as ad

def load_subset(h5ad_path, n_genes=1000, n_cells=200000):
    adata = sc.read_h5ad(h5ad_path, backed="r")
    # pick highly variable genes (or fixed list)
    X = adata[:, :n_genes][:n_cells].to_memory()  # quick subset
    sc.pp.normalize_total(X, target_sum=1e4); sc.pp.log1p(X)
    # keep needed obs columns: guide, target_gene, dose, cell_type
    keep = ["guide", "target_gene", "dose", "cell_type"]
    X.obs = X.obs[keep]
    return X

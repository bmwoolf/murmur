"""
Ingest and preprocess X-Atlas/Orion dataset.

The X-Atlas/Orion dataset contains genome-wide Perturb-seq data from:
- HCT116 (colorectal cancer cell line)
- HEK293T (kidney embryonic cell line)

Reference: https://huggingface.co/datasets/Xaira-Therapeutics/X-Atlas-Orion
"""

import scanpy as sc
import anndata as ad
import yaml
import pandas as pd
import numpy as np
from pathlib import Path
from datasets import load_dataset

def load_config():
    """Load configuration from config.yml"""
    config_path = Path(__file__).parent.parent / "config.yml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def load_xatlas_from_hf(cell_line="HCT116", max_cells=None):
    """
    Load X-Atlas/Orion dataset from Hugging Face in streaming mode.
    
    Args:
        cell_line: "HCT116" or "HEK293T"
        max_cells: Maximum number of cells to load (None for all)
    
    Returns:
        AnnData object with perturbation data
    """
    print(f"Loading X-Atlas/Orion dataset: {cell_line}")
    
    # load from Hugging Face
    ds = load_dataset(
        "Xaira-Therapeutics/X-Atlas-Orion", 
        streaming=True, 
        split=cell_line
    )
    
    # load gene metadata
    gene_metadata = load_dataset(
        "Xaira-Therapeutics/X-Atlas-Orion",
        "gene_metadata"
    )['train'].to_pandas()
    
    print(f"✓ Loaded gene metadata: {len(gene_metadata)} genes")
    
    # collect cells
    cells = []
    cell_count = 0
    
    for cell in ds:
        cells.append(cell)
        cell_count += 1
        
        if max_cells and cell_count >= max_cells:
            break
            
        if cell_count % 10000 == 0:
            print(f"  Loaded {cell_count:,} cells...")
    
    print(f"✓ Loaded {len(cells):,} cells")
    
    # convert to AnnData
    # the dataset has gene_token_id and gene_expression as sparse representation
    # we need to convert this to a proper count matrix
    
    return cells, gene_metadata

def load_xatlas_streaming(cell_line="HCT116", n_cells=10000, n_genes=1000):
    """
    Load X-Atlas data in streaming mode for CPA training.
    
    Args:
        cell_line: "HCT116" or "HEK293T"
        n_cells: Number of cells to stream
        n_genes: Number of genes to use
    
    Returns:
        AnnData object ready for CPA
    """
    print(f"streaming {n_cells:,} cells from {cell_line}...")
    
    # load gene metadata
    gene_metadata = load_dataset(
        "Xaira-Therapeutics/X-Atlas-Orion",
        "gene_metadata"
    )['train'].to_pandas()
    
    # stream cells
    ds = load_dataset(
        "Xaira-Therapeutics/X-Atlas-Orion", 
        streaming=True, 
        split=cell_line
    )
    
    # collect cells and build sparse matrix
    cells_data = []
    cell_count = 0
    
    print("~~~collecting cells~~~")
    for cell in ds:
        if cell_count >= n_cells:
            break
            
        # convert sparse representation to dense for this cell
        gene_indices = cell['gene_token_id']
        gene_counts = cell['gene_expression']
        
        # create dense vector for this cell
        dense_vector = np.zeros(len(gene_metadata))
        for idx, count in zip(gene_indices, gene_counts):
            if idx < len(gene_metadata):  # safety check
                dense_vector[idx] = count
        
        cells_data.append({
            'expression': dense_vector,
            'gene_target': cell['gene_target'],
            'guide_target': cell['guide_target'],
            'sample': cell['sample'],
            'cell_barcode': cell['cell_barcode'],
            'total_counts': cell['total_counts'],
            'pct_counts_mt': cell['pct_counts_mt'],
            'pass_guide_filter': cell['pass_guide_filter']
        })
        
        cell_count += 1
        if cell_count % 1000 == 0:
            print(f"  Loaded {cell_count:,}/{n_cells:,} cells...")
    
    print(f"✓ Loaded {len(cells_data):,} cells")
    
    # build AnnData object
    print("~~~building AnnData object~~~")
    
    # expression matrix (cells × genes)
    X = np.array([cell['expression'] for cell in cells_data])
    
    # cell metadata
    obs_data = {
        'perturbation': [cell['gene_target'] for cell in cells_data],
        'batch': [cell['sample'] for cell in cells_data],
        'cell_barcode': [cell['cell_barcode'] for cell in cells_data],
        'guide_target': [cell['guide_target'] for cell in cells_data],
        'total_counts': [cell['total_counts'] for cell in cells_data],
        'pct_counts_mt': [cell['pct_counts_mt'] for cell in cells_data],
        'pass_guide_filter': [cell['pass_guide_filter'] for cell in cells_data]
    }
    
    # gene metadata
    var_data = {
        'gene_name': gene_metadata['gene_name'].values,
        'ensembl_id': gene_metadata['ensembl_id'].values
    }
    
    # create AnnData
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(obs_data),
        var=pd.DataFrame(var_data)
    )
    
    # preprocess
    print("~~~preprocessing~~~")
    
    # normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # select highly variable genes
    if n_genes and n_genes < adata.n_vars:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_genes)
        adata = adata[:, adata.var['highly_variable']]
    
    print(f"✓ final streaming dataset: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"✓ perturbations: {adata.obs['perturbation'].nunique()}")
    
    return adata

def load_subset(h5ad_path=None, n_genes=1000, n_cells=200000, cell_line="HCT116"):
    """
    Load subset of X-Atlas data for CPA training.
    
    Args:
        h5ad_path: Path to local h5ad file (if available)
        n_genes: Number of highly variable genes to keep
        n_cells: Maximum number of cells to load
        cell_line: "HCT116" or "HEK293T"
    
    Returns:
        AnnData object preprocessed for CPA
    """
    config = load_config()
    max_ram_gb = config['pipeline']['memory']['max_ram_usage_gb']
    
    if h5ad_path and Path(h5ad_path).exists():
        print(f"Loading from local h5ad: {h5ad_path}")
        adata = sc.read_h5ad(h5ad_path, backed="r")
        
        # subset cells and genes
        if n_cells:
            adata = adata[:n_cells]
        
        adata = adata.to_memory()
        
    else:
        print(f"Loading from Hugging Face (streaming mode)")
        print("Note: Streaming mode loads data on-demand without caching")
        return load_xatlas_streaming(cell_line=cell_line, n_cells=n_cells, n_genes=n_genes)
    
    # preprocess for CPA
    print("~~~preprocessing for CPA~~~")
    
    # normalize and log-transform
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # select highly variable genes
    if n_genes:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_genes)
        adata = adata[:, adata.var['highly_variable']]
    
    # the X-Atlas dataset has these columns in obs:
    # - guide_target: guide identity
    # - gene_target: gene targeted by guide
    # - sample: GEM batch
    # - cell_barcode: cell identifier
    
    # map to CPA expected format
    if 'gene_target' in adata.obs.columns:
        adata.obs['perturbation'] = adata.obs['gene_target']
    
    if 'sample' in adata.obs.columns:
        adata.obs['batch'] = adata.obs['sample']
    
    # keep essential columns
    essential_cols = ['perturbation', 'batch', 'cell_barcode', 'guide_target']
    keep_cols = [c for c in essential_cols if c in adata.obs.columns]
    adata.obs = adata.obs[keep_cols]
    
    print(f"✓ final dataset: {adata.n_obs} cells × {adata.n_vars} genes")
    
    return adata

if __name__ == "__main__":
    # test loading
    config = load_config()
    h5ad_path = config['pipeline']['input']['xatlas_h5ad']
    
    print("~~~testing X-Atlas data loading~~~")
    adata = load_subset(h5ad_path=h5ad_path, n_cells=1000, n_genes=500)
    print(adata)

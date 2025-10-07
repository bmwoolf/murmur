# src/pathways.py
import numpy as np, pandas as pd

def load_gene_sets(msigdb_gmt, gene_index):
    # returns a list of binary masks per pathway over your gene order
    ...

def pathway_scores(expr_vec, masks):
    # simple average as baseline
    return np.array([expr_vec[m.astype(bool)].mean() if m.sum() else 0.0 for m in masks])

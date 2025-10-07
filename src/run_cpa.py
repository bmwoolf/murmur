import torch, numpy as np
from .train_cpa import CPA
from sklearn.preprocessing import LabelEncoder

def simulate_transcriptome(model, gene_dose_df, le_gene, le_cell, cell_type="HEK293T"):
    genes = gene_dose_df["gene"].tolist()
    doses = torch.tensor(gene_dose_df["dose"].values, dtype=torch.float32)
    gidx = torch.tensor(le_gene.transform([g if g in le_gene.classes_ else le_gene.classes_[0] for g in genes]))
    cidx = torch.tensor([le_cell.transform([cell_type])[0]] * len(genes))
    preds = model(gidx, doses, cidx)  # (N, n_genes)
    # combine multiple perturbations: sum or mean; start with mean
    return preds.mean(0).detach().cpu().numpy()  # vector[n_genes]

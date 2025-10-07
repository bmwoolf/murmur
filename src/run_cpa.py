import torch, numpy as np
from .train_cpa import CPA
from sklearn.preprocessing import LabelEncoder
import yaml
from pathlib import Path

def load_config():
    """Load configuration from config.yml"""
    config_path = Path(__file__).parent.parent / "config.yml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def simulate_transcriptome(model, gene_dose_df, le_gene, le_cell, cell_type="HEK293T"):
    config = load_config()
    gpu_memory_fraction = config['pipeline']['memory']['gpu_memory_fraction']
    
    genes = gene_dose_df["gene"].tolist()
    doses = torch.tensor(gene_dose_df["dose"].values, dtype=torch.float32)
    gidx = torch.tensor(le_gene.transform([g if g in le_gene.classes_ else le_gene.classes_[0] for g in genes]))
    cidx = torch.tensor([le_cell.transform([cell_type])[0]] * len(genes))
    preds = model(gidx, doses, cidx)  # (N, n_genes)
    # combine multiple perturbations: sum or mean; start with mean
    return preds.mean(0).detach().cpu().numpy()  # vector[n_genes]

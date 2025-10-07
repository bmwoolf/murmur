import torch, torch.nn as nn
from sklearn.preprocessing import LabelEncoder
import yaml
from pathlib import Path

class CPA(nn.Module):
    def __init__(self, n_genes, n_gene_ids, n_celltypes, d_lat=128):
        super().__init__()
        self.gene_emb = nn.Embedding(n_gene_ids, d_lat)
        self.cell_emb = nn.Embedding(n_celltypes, d_lat)
        self.dose_mlp = nn.Sequential(nn.Linear(1, d_lat), nn.SiLU(), nn.Linear(d_lat, d_lat))
        self.decoder = nn.Sequential(nn.Linear(3*d_lat, 512), nn.SiLU(),
                                     nn.Linear(512, n_genes))
    def forward(self, gene_idx, dose, cell_idx):
        z = torch.cat([self.gene_emb(gene_idx),
                       self.dose_mlp(dose.view(-1,1)),
                       self.cell_emb(cell_idx)], dim=1)
        return self.decoder(z)  # predicted Δexpr

def load_config():
    """Load configuration from config.yml"""
    config_path = Path(__file__).parent.parent / "config.yml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def prepare_labels(adata):
    config = load_config()
    batch_size = config['pipeline']['cpa_model']['batch_size']
    
    le_gene = LabelEncoder().fit(adata.obs["target_gene"])
    le_cell = LabelEncoder().fit(adata.obs["cell_type"])
    g = torch.tensor(le_gene.transform(adata.obs["target_gene"]), dtype=torch.long)
    c = torch.tensor(le_cell.transform(adata.obs["cell_type"]), dtype=torch.long)
    d = torch.tensor(adata.obs["dose"].values, dtype=torch.float32)
    X = torch.tensor(adata.X.A if hasattr(adata.X, "A") else adata.X, dtype=torch.float32)
    return X, g, d, c, le_gene, le_cell

# Training loop (MSE on Δexpr or expr; add gene-wise weighting if desired)

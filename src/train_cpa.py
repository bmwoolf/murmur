import torch, torch.nn as nn
import torch.optim as optim
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
import yaml
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import os

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
    """Prepare X-Atlas data for CPA training."""
    print("~~~preparing labels for CPA training~~~")
    
    # Use perturbation column (gene_target) as the target gene
    perturbation_col = 'perturbation' if 'perturbation' in adata.obs.columns else 'gene_target'
    
    # Encode categorical variables
    le_gene = LabelEncoder().fit(adata.obs[perturbation_col])
    le_cell = LabelEncoder().fit(adata.obs['batch'])  # Use batch as cell type for now
    
    # Create tensors
    g = torch.tensor(le_gene.transform(adata.obs[perturbation_col]), dtype=torch.long)
    c = torch.tensor(le_cell.transform(adata.obs['batch']), dtype=torch.long)
    
    # For X-Atlas, we don't have explicit dose values, so we'll use 1.0 for all perturbations
    # In real CPA, this would be the perturbation strength (0.0 = no perturbation, 1.0 = full)
    d = torch.ones(len(adata), dtype=torch.float32)
    
    # Expression matrix
    X = torch.tensor(adata.X, dtype=torch.float32)
    
    print(f"✓ Prepared {len(adata)} cells with {len(le_gene.classes_)} perturbations")
    print(f"✓ Gene expression shape: {X.shape}")
    
    return X, g, d, c, le_gene, le_cell

def train_cpa_model(adata, model_save_path=None, epochs=100):
    """
    Train CPA model on X-Atlas data.
    
    Args:
        adata: AnnData object from X-Atlas
        model_save_path: Path to save trained model
        epochs: Number of training epochs
    
    Returns:
        Trained CPA model and label encoders
    """
    config = load_config()
    batch_size = config['pipeline']['cpa_model']['batch_size']
    learning_rate = config['pipeline']['cpa_model']['learning_rate']
    max_epochs = config['pipeline']['cpa_model']['max_epochs']
    
    print("=" * 60)
    print("CPA Model Training")
    print("=" * 60)
    
    # Prepare data
    X, g, d, c, le_gene, le_cell = prepare_labels(adata)
    
    # Split data
    print("~~~splitting data~~~")
    indices = np.arange(len(X))
    train_idx, val_idx = train_test_split(indices, test_size=0.2, random_state=42)
    
    train_dataset = torch.utils.data.TensorDataset(
        X[train_idx], g[train_idx], d[train_idx], c[train_idx]
    )
    val_dataset = torch.utils.data.TensorDataset(
        X[val_idx], g[val_idx], d[val_idx], c[val_idx]
    )
    
    train_loader = torch.utils.data.DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = torch.utils.data.DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
    
    print(f"✓ Training samples: {len(train_idx)}")
    print(f"✓ Validation samples: {len(val_idx)}")
    
    # Initialize model
    print("~~~initializing CPA model~~~")
    model = CPA(
        n_genes=X.shape[1],
        n_gene_ids=len(le_gene.classes_),
        n_celltypes=len(le_cell.classes_)
    )
    
    # Move to GPU if available
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = model.to(device)
    print(f"✓ Using device: {device}")
    
    # Optimizer and loss
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    criterion = nn.MSELoss()
    
    # Training loop
    print(f"~~~training for {max_epochs} epochs~~~")
    train_losses = []
    val_losses = []
    
    for epoch in range(max_epochs):
        # Training
        model.train()
        train_loss = 0.0
        for batch_X, batch_g, batch_d, batch_c in tqdm(train_loader, desc=f"Epoch {epoch+1}/{max_epochs}"):
            batch_X, batch_g, batch_d, batch_c = batch_X.to(device), batch_g.to(device), batch_d.to(device), batch_c.to(device)
            
            optimizer.zero_grad()
            
            # Forward pass
            pred_delta = model(batch_g, batch_d, batch_c)
            
            # Loss: predict the difference from baseline (non-perturbed) expression
            # For now, we'll use the actual expression as target
            # In real CPA, this would be Δexpr = expr_perturbed - expr_baseline
            target = batch_X
            loss = criterion(pred_delta, target)
            
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
        
        # Validation
        model.eval()
        val_loss = 0.0
        with torch.no_grad():
            for batch_X, batch_g, batch_d, batch_c in val_loader:
                batch_X, batch_g, batch_d, batch_c = batch_X.to(device), batch_g.to(device), batch_d.to(device), batch_c.to(device)
                
                pred_delta = model(batch_g, batch_d, batch_c)
                target = batch_X
                loss = criterion(pred_delta, target)
                val_loss += loss.item()
        
        avg_train_loss = train_loss / len(train_loader)
        avg_val_loss = val_loss / len(val_loader)
        
        train_losses.append(avg_train_loss)
        val_losses.append(avg_val_loss)
        
        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}: Train Loss = {avg_train_loss:.4f}, Val Loss = {avg_val_loss:.4f}")
    
    print("✓ Training completed!")
    
    # Save model
    if model_save_path:
        os.makedirs(os.path.dirname(model_save_path), exist_ok=True)
        torch.save({
            'model_state_dict': model.state_dict(),
            'le_gene': le_gene,
            'le_cell': le_cell,
            'n_genes': X.shape[1],
            'train_losses': train_losses,
            'val_losses': val_losses
        }, model_save_path)
        print(f"✓ Model saved to: {model_save_path}")
    
    return model, le_gene, le_cell

def load_trained_cpa(model_path):
    """Load a trained CPA model."""
    checkpoint = torch.load(model_path, map_location='cpu')
    
    model = CPA(
        n_genes=checkpoint['n_genes'],
        n_gene_ids=len(checkpoint['le_gene'].classes_),
        n_celltypes=len(checkpoint['le_cell'].classes_)
    )
    
    model.load_state_dict(checkpoint['model_state_dict'])
    model.eval()
    
    return model, checkpoint['le_gene'], checkpoint['le_cell']

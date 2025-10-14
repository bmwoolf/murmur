import torch, numpy as np
import pandas as pd
import os
from .train_cpa import CPA, load_trained_cpa
from sklearn.preprocessing import LabelEncoder
import yaml
from pathlib import Path

def load_config():
    """Load configuration from config.yml"""
    config_path = Path(__file__).parent.parent / "config.yml"
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def simulate_transcriptome(model, gene_dose_df, le_gene, le_cell, cell_type="HCT116_Batch1", device=None):
    """
    Simulate transcriptome perturbations from personal variants.
    
    Args:
        model: Trained CPA model
        gene_dose_df: DataFrame with columns ['gene', 'dose']
        le_gene: Label encoder for genes
        le_cell: Label encoder for cell types
        cell_type: Cell type to simulate (should match training data)
        device: PyTorch device
    
    Returns:
        Predicted expression changes as numpy array
    """
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    model.eval()
    model = model.to(device)
    
    print(f"~~~simulating perturbations for {len(gene_dose_df)} genes~~~")
    print(f"   Cell type: {cell_type}")
    print(f"   Device: {device}")
    
    # Filter genes that exist in the model's training data
    valid_genes = gene_dose_df[gene_dose_df['gene'].isin(le_gene.classes_)]
    if len(valid_genes) < len(gene_dose_df):
        print(f"   Warning: {len(gene_dose_df) - len(valid_genes)} genes not found in training data")
    
    if len(valid_genes) == 0:
        print("   Warning: No valid genes found, returning zeros")
        return np.zeros(model.decoder[-1].out_features)
    
    genes = valid_genes["gene"].tolist()
    doses = torch.tensor(valid_genes["dose"].values, dtype=torch.float32)
    
    # Encode genes
    gidx = torch.tensor(le_gene.transform(genes), dtype=torch.long)
    
    # Encode cell type
    if cell_type not in le_cell.classes_:
        print(f"   Warning: Cell type '{cell_type}' not found, using first available")
        cell_type = le_cell.classes_[0]
    cidx = torch.tensor([le_cell.transform([cell_type])[0]] * len(genes), dtype=torch.long)
    
    # Move to device
    gidx, doses, cidx = gidx.to(device), doses.to(device), cidx.to(device)
    
    # Predict perturbations
    with torch.no_grad():
        preds = model(gidx, doses, cidx)  # (N, n_genes)
        
        # Combine multiple perturbations: mean of all predicted changes
        combined_pred = preds.mean(0).detach().cpu().numpy()
    
    print(f"✓ Predicted expression changes for {len(genes)} genes")
    return combined_pred

def predict_personal_variants(vcf_annotations_path, model_path, output_path=None):
    """
    Predict transcriptome effects of personal variants using trained CPA model.
    
    Args:
        vcf_annotations_path: Path to VEP annotations (TSV file)
        model_path: Path to trained CPA model
        output_path: Path to save predictions
    
    Returns:
        Dictionary with prediction results
    """
    print("=" * 60)
    print("Personal Variant Prediction with CPA")
    print("=" * 60)
    
    # Load trained model
    print("~~~loading trained CPA model~~~")
    model, le_gene, le_cell = load_trained_cpa(model_path)
    print(f"✓ Model loaded: {len(le_gene.classes_)} genes, {len(le_cell.classes_)} cell types")
    
    # Load gene doses from VCF annotations
    print("~~~loading gene doses from VCF annotations~~~")
    gene_dose_path = str(Path(vcf_annotations_path).parent / "gene_dose.csv")
    
    if not Path(gene_dose_path).exists():
        raise FileNotFoundError(f"Gene dose file not found: {gene_dose_path}")
    
    gene_dose_df = pd.read_csv(gene_dose_path)
    print(f"✓ Loaded {len(gene_dose_df)} genes with variants")
    
    # Simulate transcriptome
    predictions = {}
    
    # Predict for different cell types if available
    for cell_type in le_cell.classes_:
        print(f"\n~~~predicting for {cell_type}~~~")
        pred_expr = simulate_transcriptome(model, gene_dose_df, le_gene, le_cell, cell_type)
        predictions[cell_type] = pred_expr
    
    # Save results
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Create results DataFrame
        results_df = pd.DataFrame(predictions)
        results_df.to_csv(output_path, index=False)
        print(f"✓ Predictions saved to: {output_path}")
    
    return {
        'predictions': predictions,
        'gene_dose_df': gene_dose_df,
        'model_info': {
            'n_genes': len(le_gene.classes_),
            'n_celltypes': len(le_cell.classes_),
            'trained_genes': le_gene.classes_.tolist()
        }
    }

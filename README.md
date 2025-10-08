# murmur

VCF → annotation → perturbation → pathways → phenotype

## training pipeline 
1. download VCF
2. annotate variants (VEP + SnpEff)
3. measure perturbations (CPA trained/fine-tuned on X-atlas, factoring in personal variants)
4. model disregulated pathways via chain reactions (PathDNN + Reactome)
5. measure phenotype differences (DeepDR + GWAS Catalog)
6. return visualized results (make it beautiful)

## Frontend user experience 
1. upload VCF file (Nucleus, Nebula Genomics)
2. user sees progress bars for:
    - variant annotation (VEP/SnpEff)
    - functional simulation (CPA)
    - pathway mapping (PathDNN)
    - phenotype inference (DeepDR)
    - report generation
3. a report shows how their unique variants ripple through genes, pathways, and traits

## Backend computations
1. VCF is uploaded
2. annotate variants (VEP + SnpEff)
3. measure perturbations (CPA trained/fine-tuned on X-atlas, factoring in personal variants)
4. model disregulated pathways via chain reactions (PathDNN + Reactome)
5. measure phenotype differences (DeepDR + GWAS Catalog)
6. return results to user in pretty viualization

## Environment setup
Packages:
- Python 3.10
- PyTorch 2.6.0 with CUDA 12.4
- scanpy & anndata (single-cell genomics analysis)
- pandas & numpy
- scikit-learn
- cyvcf2 (VCF file processing)
- networkx (graph analysis)
- matplotlib & seaborn
- scvi-tools (single-cell variational inference)

## Quick Start

```bash
# 1. activate the environment
conda activate murmur

# 2. fill out machine config
./config/generate_config.sh

# 3. run pipeline
./src/pipeline.py
```

## Example usage
```python
expr_pred, pathway_scores, trait_scores = run(
    vcf_path="path/to/sample.vcf",
    out_dir="results/",
    xatlas_h5ad="data/xatlas.h5ad",
    msigdb_gmt="data/msigdb.gmt",
    gwas_map="data/gwas_traits.tsv"
)
```
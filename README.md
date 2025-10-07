# murmur

DNA → perturbation → phenotype

## workflow

```
input: VCF variants
    ↓
CPA + X-atlas (perturbation)
    ↓
PathDNN + Reactome (pathway)
    ↓
DeepDR + GWAS Catalog (phenotype)
    ↓
output: phenotype probabilities
```
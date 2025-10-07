# murmur

VCF → annotation → perturbation → pathways → phenotype

## training pipeline 
1. download VCF
2. annotate variants (VEP + SnpEff)
3. measure perturbations (CPA trained/fine-tuned on X-atlas, factoring in personal variants)
4. model disregulated pathways via chain reactions (PathDNN + Reactome)
5. measure phenotype differences (DeepDR + GWAS Catalog)
6. return visualized results (make it beautiful)


## user inference 
1. upload VCF file (Nucleus, Nebula Genomics)
2. user sees progress bars for:
    - variant annotation (VEP/SnpEff)
    - functional simulation (CPA)
    - pathway mapping (PathDNN)
    - phenotype inference (DeepDR)
    - report generation
3. a report shows how their unique variants ripple through genes, pathways, and traits



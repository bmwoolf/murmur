# murmur

VCF → annotation → perturbation → pathways → phenotype


## user steps 
1. upload VCF file (Nucleus, Nebula Genomics)
2. user sees progress bars for:
    - variant annotation (VEP/SnpEff)
    - functional simulation (CPA)
    - pathway mapping (PathDNN)
    - phenotype inference (DeepDR)
    - report generation
3. a report shows how their unique variants ripple through genes, pathways, and traits


## pipeline 
1. VCF is uploaded
2. annotate variants (VEP + SnpEff)
3. measure perturbations (CPA trained/fine-tuned on X-atlas, factoring in personal variants)
4. model disregulated pathways via chain reactions (PathDNN + Reactome)
5. measure phenotype differences (DeepDR + GWAS Catalog)
6. return results to user in pretty viualization
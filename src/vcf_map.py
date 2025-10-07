import cyvcf2, pandas as pd

def vcf_to_gene_dose(vcf_path, consequence_tsv):
    # consequence_tsv: from VEP/SnpEff with gene + impact
    cons = pd.read_csv(consequence_tsv, sep="\t")
    rows = []
    for _, r in cons.iterrows():
        impact = r.get("IMPACT","MODIFIER")
        gene = r["SYMBOL"]
        if impact == "HIGH": dose = 0.8
        elif impact == "MODERATE": dose = 0.5
        elif impact == "LOW": dose = 0.2
        else: dose = 0.1
        # Noncoding/regulatory could use motif/Enformer scores to set dose later
        rows.append({"gene": gene, "dose": dose})
    return pd.DataFrame(rows).groupby("gene")["dose"].max().reset_index()

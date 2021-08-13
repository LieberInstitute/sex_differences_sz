"""
This script is used to generate the supplementary data for
the sex differences manuscript. Specifically, it labels
and combines functional enrichment analysis for stringent
male-specific analysis separated by direction of effect.
"""
import pandas as pd

def tissue_map(tissue):
    return {"caudate": "Caudate", "dlpfc": "DLPFC",
            "hippocampus": "Hippocampus"}[tissue]


def get_enrichment(tissue, label):
    config = {
        "All_DEGs": "DEGs_functional_enrichment.tsv",
        "Upregulated": "upreg_DEGs_functional_enrichment.tsv",
        "Downregulated": "downreg_DEGs_functional_enrichment.tsv",
    }
    enrich_file = "../../../%s/male_analysis/" % tissue +\
        "gprofiler_analysis/_m/%s" % config[label]
    df = pd.read_csv(enrich_file, sep='\t')
    df["Direction"] = label
    df["Tissue"] = tissue_map(tissue)
    return df


def main():
    df = pd.DataFrame()
    for tissue in ["caudate", "dlpfc", "hippocampus"]:
        for labelx in ["All_DEGs", "Upregulated", "Downregulated"]:
            df = pd.concat([df, get_enrichment(tissue, labelx)], axis=0)
    df.to_csv("functional_enrichment_analysis_maleSZ_3brain_regions.txt",
              index=False, sep='\t')


if __name__ == '__main__':
    main()

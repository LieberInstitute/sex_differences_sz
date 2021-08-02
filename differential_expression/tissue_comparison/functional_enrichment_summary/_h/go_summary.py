## Combine functional enrichment analysis
import numpy as np
import pandas as pd

def map_tissue(tissue):
    return {"caudate": "Caudate", "dlpfc": "DLPFC",
            "hippocampus": "Hippocampus"}[tissue]


def annotate_GO(tissue, fn, label):
    df = pd.read_csv(fn, sep='\t').sort_values('p_value')
    df["Log10"] = -np.log10(df['p_value'])
    df["Tissue"] = map_tissue(tissue)
    df["Bias"] = label
    return df


def extract_GO(tissue):
    config = {
        'All': '../../../%s/gprofiler_analysis/' % tissue +\
        '_m/allDEGs_functional_enrichment.txt',
        'Female': '../../../%s/gprofiler_analysis/_m/' % tissue +\
        'female_bias_DEGs_functional_enrichment.txt',
        'Male': '../../../%s/gprofiler_analysis/_m/' % tissue +\
        'male_bias_DEGs_functional_enrichment.txt',
    }
    go_df = []
    for bias in ["All", "Female", "Male"]:
        go_df.append(annotate_GO(tissue, config[bias], bias))
    df = pd.concat(go_df, axis=0)
    return df


def main():
    bigdf = pd.DataFrame()
    for tissue in ["caudate", "dlpfc", "hippocampus"]:
        bigdf = pd.concat([bigdf, extract_GO(tissue)], axis=0)
    bigdf.to_csv("functional_enrichment_analysis_3brain_regions_sex.txt",
                 sep='\t', index=False)


if __name__ == '__main__':
    main()

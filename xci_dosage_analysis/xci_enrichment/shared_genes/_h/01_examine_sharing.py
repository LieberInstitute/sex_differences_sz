"""
This script conducts enrichment analysis for XCI and sharing.
"""
import numpy as np
import pandas as pd
import session_info
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

@lru_cache()
def get_deg():
    fn = here("differential_expression/tissue_comparison/summary_table",
              "_m/differential_expression_analysis_4features_sex.txt.gz")
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == "Gene")].copy()


@lru_cache()
def get_shared_genes():
    df = get_deg()[(get_deg()["adj.P.Val"] <= 0.05)].copy()
    shared_genes = list(df.groupby("Feature").size()[(df.groupby("Feature").size() == 3)].index)
    return shared_genes


@lru_cache()
def get_xci():
    fn = "../../_h/xci_status_hg19.txt"
    df = pd.read_csv(fn, sep='\t')
    df["ensemblID"] = df['Gene ID'].str.replace("\\..*", "", regex=True)
    return df.loc[:, ["ensemblID", "Combined XCI status"]]


@lru_cache()
def merge_data():
    return get_deg().merge(get_xci(), on="ensemblID", how="left")\
        .drop_duplicates(subset="ensemblID")


def cal_fishers(status, direction):
    df = merge_data(); shared_genes = get_shared_genes()
    df["Shared"] = ["Yes" if x in shared_genes else "No" for x in df.Feature]
    if direction == "Up":
        df = df[(df["Direction"] > 0)].copy()
    elif direction == "Down":
        df = df[(df["Direction"] < 0)].copy()
    else:
        df = df
    table = [[np.sum((df["Shared"] == "Yes") &
                     (df["Combined XCI status"] == status)),
              np.sum((df["Shared"] == "Yes") &
                     (df["Combined XCI status"] != status))],
             [np.sum((df["Shared"] == "No") &
                     (df["Combined XCI status"] == status)),
              np.sum((df["Shared"] == "No") &
                     (df["Combined XCI status"] != status))]]
    print(table)
    return fisher_exact(table)
    

def calculate_enrichment():
    dir_lt = []; fdr_lt = []; xci_lt = []; oddratio_lt = []; pvals = []
    for direction in ['Up', 'Down', 'All']:
        for status in get_xci().loc[:, 'Combined XCI status'].unique():
            odd_ratio, pval = cal_fishers(status, direction)
            xci_lt.append(status); pvals.append(pval);
            oddratio_lt.append(odd_ratio); dir_lt.append(direction)
    _, fdr, _, _ = multipletests(pvals, method='bonferroni')
    # Generate dataframe
    return pd.DataFrame({'XCI_status': xci_lt, 'OR': oddratio_lt,
                         'PValue': pvals, "Bonferroni": fdr,
                         'Direction': dir_lt})
                

def main():
    ## Run enrichment
    calculate_enrichment()\
        .to_csv('xci_enrichment_analysis.txt', sep='\t', index=False)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()

"""
This script conducts enrichment analysis for XCI and sharing.
"""
import numpy as np
import pandas as pd
import session_info
from pyhere import here
from functools import lru_cache
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import fdrcorrection

@lru_cache()
def get_deg(tissue):
    fn = here("differential_expression/tissue_comparison/summary_table",
              "_m/differential_expression_analysis_4features_sex.txt.gz")
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == "Gene") &
              (df["Tissue"] == tissue)].copy()


@lru_cache()
def get_xci():
    fn = "../../_h/xci_status_hg19.txt"
    df = pd.read_csv(fn, sep='\t')
    df["ensemblID"] = df['Gene ID'].str.replace("\\..*", "", regex=True)
    return df.loc[:, ["ensemblID", "Combined XCI status"]]


@lru_cache()
def merge_data(tissue):
    return get_deg(tissue).merge(get_xci(), on="ensemblID", how="left")


def cal_fishers(status, tissue, direction):
    df = merge_data(tissue)
    if direction == "Up":
        df = df[(df["Direction"] > 0)].copy()
    elif direction == "Down":
        df = df[(df["Direction"] < 0)].copy()
    else:
        df = df
    table = [[np.sum((df["adj.P.Val"] <= 0.05) &
                     (df["Combined XCI status"] == status)),
              np.sum((df["adj.P.Val"] <= 0.05) &
                     (df["Combined XCI status"] != status))],
             [np.sum((df["adj.P.Val"] >  0.05) &
                     (df["Combined XCI status"] == status)),
              np.sum((df["adj.P.Val"] >  0.05) &
                     (df["Combined XCI status"] != status))]]
    print(table)
    return fisher_exact(table)
    

def calculate_enrichment():
    region_lt = []; dir_lt = []; fdr_lt = []; xci_lt = []; pval_lt = [];
    oddratio_lt = []
    for tissue in ["Caudate", "Dentate.Gyrus", "DLPFC", "Hippocampus"]:
        pvals = []
        for direction in ['Up', 'Down', 'All']:
            for status in get_xci().loc[:, 'Combined XCI status'].unique():
                odd_ratio, pval = cal_fishers(status, tissue, direction)
                xci_lt.append(status); pvals.append(pval); region_lt.append(tissue)
                oddratio_lt.append(odd_ratio); dir_lt.append(direction)
        _, fdr = fdrcorrection(pvals) # FDR correction per comparison and version
        pval_lt = np.concatenate((pval_lt, pvals))
        fdr_lt = np.concatenate((fdr_lt, fdr))
    # Generate dataframe
    return pd.DataFrame({'Tissue': region_lt, 'XCI_status': xci_lt,
                         'OR': oddratio_lt, 'PValue': pval_lt,
                         "FDR": fdr_lt, 'Direction': dir_lt})
                

def main():
    ## Run enrichment
    calculate_enrichment()\
        .to_csv('xci_enrichment_analysis.txt', sep='\t', index=False)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()

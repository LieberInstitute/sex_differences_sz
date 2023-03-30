"""
This script extracts the tissue-specific features.
"""
import numpy as np
import pandas as pd
import session_info
from functools import lru_cache

@lru_cache()
def load_sig():
    fn = "../../_m/differential_expression_analysis_4features_sex.txt.gz"
    return pd.read_csv(fn, sep='\t')


@lru_cache()
def get_tissues(feature):
    cc = load_sig()[(load_sig()["Type"] == feature) &
                    (load_sig()["adj.P.Val"] <= 0.05) &
                    (load_sig()["Tissue"] == "Caudate")].copy()
    dd = load_sig()[(load_sig()["Type"] == feature) &
                    (load_sig()["adj.P.Val"] <= 0.05) &
                    (load_sig()["Tissue"] == "DLPFC")].copy()
    hh = load_sig()[(load_sig()["Type"] == feature) &
                    (load_sig()["adj.P.Val"] <= 0.05) &
                    (load_sig()["Tissue"] == "Hippocampus")].copy()
    return cc, dd, hh


@lru_cache()
def region_specific(feature):
    cc, dd, hh = get_tissues(feature)
    caud8 = set(cc.Feature) - set(dd.Feature) - set(hh.Feature)
    dlpfc = set(dd.Feature) - set(cc.Feature) - set(hh.Feature)
    hippo = set(hh.Feature) - set(cc.Feature) - set(dd.Feature)
    return pd.concat([cc[(cc["Feature"].isin(caud8))],
                      dd[(dd["Feature"].isin(dlpfc))],
                      hh[(hh["Feature"].isin(hippo))]], axis=0)


@lru_cache()
def get_discordant(feature):
    cc, dd, hh = get_tissues(feature)
    shared = set(cc.Feature) & set(dd.Feature) & set(hh.Feature)
    shared_df0 = load_sig()[(load_sig()["Feature"].isin(shared))]\
        .loc[:, ["Tissue", "Feature", "logFC", "Chrom_Type"]]
    print(shared_df0.groupby(["Tissue", "Chrom_Type"]).size())
    shared_df = shared_df0.pivot_table(
        values="logFC",columns=["Tissue"],index=["Feature"])\
                          .apply(np.sign)
    discordant = shared_df.eq(shared_df.loc[:, "Caudate"], axis=0).all(axis=1)
    return cc.merge(pd.DataFrame({"Concordant":discordant}).reset_index(),
                    on="Feature")\
             .drop(["Tissue", "adj.P.Val", "logFC"], axis=1)


def main():
    df1 = pd.DataFrame(); df2= pd.DataFrame()
    for feature in ["Gene", "Transcript", "Exon", "Junction"]:
        df1 = pd.concat([df1, region_specific(feature)], axis=0)
        df2 = pd.concat([df2, get_discordant(feature)], axis=0)
    df1.to_csv("BrainSeq_sex_region_specific.tsv", sep='\t',index=False)
    df2.to_csv("BrainSeq_sex_shared_disconcordant.tsv", sep='\t',
               index=False)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()

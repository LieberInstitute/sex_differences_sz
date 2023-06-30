"""
This script summarized DE results by chromosome
"""
import session_info
import pandas as pd
from pyhere import here
from scipy.stats import kstest
from functools import lru_cache

@lru_cache()
def load_data():
    fn = "BrainSeq_sexGenotypes_4features_3regions.txt.gz"
    return pd.read_csv(fn, sep='\t')


@lru_cache()
def filter_data():
    df = load_data()
    return df[(df["feature_type"] == "Gene")].copy()


@lru_cache()
def get_annot():
    fn = here("input/counts/text_files_counts/_m/",
              "caudate/gene_annotation.txt")
    return pd.read_csv(fn, sep='\t', usecols=[1,5,6])


@lru_cache()
def get_shat():
    fn = here("eqtl_analysis/interaction_analysis/genes/_m",
              "shat_interaction_3regions.txt.gz")
    return pd.read_csv(fn, sep='\t', usecols=[0]).drop_duplicates()


@lru_cache()
def background_genes():
    return get_annot().merge(get_shat(), left_on="gencode_id",
                             right_on="phenotype_id")


def main():
    # Load annotation
    annot = background_genes().groupby("seqnames")\
                              .size().reset_index()
    annot["percentage"] = annot.loc[:, 0] / annot.loc[:, 0].sum()
    annot.rename(columns={0: "n"}, inplace=True)
    
    # Load si-eQTL eGenes
    df = filter_data()\
        .loc[:, ["gene_id", "seqnames"]]\
        .drop_duplicates()
    print(df.shape)
    dx = df.groupby("seqnames").size().reset_index()
    dx["percentage"] = dx.loc[:, 0] / dx.loc[:, 0].sum()
    dx.rename(columns={0: "n"}, inplace=True)

    # Compare with kstest
    dt = pd.merge(annot, dx, on="seqnames",
                  suffixes=["_total", "_eGenes"])
    print(dt)
    print(kstest(annot.percentage, dx.percentage))
    dt.to_csv("by_chrom.summary.tsv", sep="\t", index=False)
    ## Reproducibilty information
    session_info.show()


if __name__ == "__main__":
    main()
    
    

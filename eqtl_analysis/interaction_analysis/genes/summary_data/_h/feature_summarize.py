"""
This script summarized DE results using mash model.
"""
import session_info
import pandas as pd
from pyhere import here
from functools import lru_cache

@lru_cache()
def get_mash_eqtl(tissue, fdr):
    # lfsr < 0.05 for significant eQTL by features
    df = pd.read_csv("../../_m/lfsr.sex_interaction.txt.gz", sep='\t')\
           .loc[:, ["effect", "gene_id", "variant_id", tissue]]\
           .rename(columns={tissue: "lfsr"})
    return df[(df["lfsr"] < fdr)].copy()


@lru_cache()
def get_mash_es(tissue):
    # get effect size
    return pd.read_csv("../../_m/posterior_mean.sex_interaction.txt.gz",
                     sep='\t')\
             .loc[:, ["gene_id", "variant_id", tissue]]\
             .rename(columns={tissue: "posterior_mean"})


@lru_cache()
def get_annotation(feature="genes"):
    config = {
        "genes": here("input/counts/text_files_counts/_m",
                      "caudate/gene_annotation.txt"),
    }
    return pd.read_csv(config[feature], sep='\t')\
             .loc[:, ["name", "seqnames", "start", "end",
                      "gene_name", "gencode_id"]]


@lru_cache()
def annotate_eqtl(tissue, fdr):
    df = get_mash_eqtl(tissue, fdr)\
        .merge(get_mash_es(tissue),
               on=["gene_id", "variant_id"])
    return df.merge(get_annotation(), left_on="gene_id",
                    right_on="gencode_id")\
             .rename(columns={"name": "feature_id"})


def get_eqtl_result_by_tissue(tissue, fdr=0.05):
    df = annotate_eqtl(tissue, fdr)
    df["feature_type"] = "Gene"
    df["region"] = tissue
    return df


def main():
    bigdata1 = []
    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        data1 = get_eqtl_result_by_tissue(tissue, 1)
        bigdata1.append(data1)
    df1 = pd.concat(bigdata1)
    cols = ["region", "gene_id", "variant_id", "gencode_id", "gene_name",
            "seqnames", "start", "end", "lfsr", "posterior_mean", "feature_type"]
    df1.sort_values(["region", "feature_type", "lfsr", "posterior_mean"])\
       .loc[:, cols]\
       .to_csv("BrainSeq_sexGenotypes_genes_3regions.txt.gz",
               sep='\t', index=False)
    ## Reproducibilty information
    session_info.show()


if __name__ == "__main__":
    main()

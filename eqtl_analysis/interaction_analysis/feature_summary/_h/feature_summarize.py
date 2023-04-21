"""
This script summarized DE results using mash model.
"""
import session_info
import pandas as pd
from pyhere import here
from functools import lru_cache

@lru_cache()
def get_mash_eqtl(feature, tissue, fdr):
    # lfsr < 0.05 for significant eQTL by features
    df = pd.read_csv(f"../../{feature}/_m/lfsr.sex_interaction.txt.gz", sep='\t')\
           .loc[:, ["effect", "gene_id", "variant_id", tissue]]\
           .rename(columns={tissue: "lfsr"})
    return df[(df["lfsr"] < fdr)].copy()


@lru_cache()
def get_mash_es(feature, tissue):
    # get effect size
    return pd.read_csv(f"../../{feature}/_m/posterior_mean.sex_interaction.txt.gz",
                     sep='\t')\
             .loc[:, ["gene_id", "variant_id", tissue]]\
             .rename(columns={tissue: "posterior_mean"})


@lru_cache()
def get_annotation(feature):
    config = {
        "genes": here("input/counts/text_files_counts/_m",
                      "caudate/gene_annotation.txt"),
        "transcripts": here("input/counts/text_files_counts/_m",
                            "caudate/tx_annotation.txt"),
        "exons": here("input/counts/text_files_counts/_m",
                      "caudate/exon_annotation.txt"),
        "junctions": here("input/counts/text_files_counts/_m",
                          "caudate/jxn_annotation.txt"),
    }
    return pd.read_csv(config[feature], sep='\t')\
             .loc[:, ["name", "seqnames", "start", "end",
                      "gene_name", "gencode_id"]]


@lru_cache()
def annotate_eqtl(feature, tissue, fdr):
    df = get_mash_eqtl(feature, tissue, fdr)\
        .merge(get_mash_es(feature, tissue),
               on=["gene_id", "variant_id"])
    if feature == "genes":
        return df.merge(get_annotation(feature), left_on="gene_id",
                        right_on="gencode_id")\
                 .rename(columns={"name": "feature_id"})
    else:
        return df.merge(get_annotation(feature), left_on="gene_id",
                        right_on="name").drop(["name"], axis=1)


@lru_cache()
def extract_features(tissue, fdr):
    # Extract DE from mash model
    genes = annotate_eqtl("genes", tissue, fdr)
    trans = annotate_eqtl("transcripts", tissue, fdr)
    exons = annotate_eqtl("exons", tissue, fdr)
    juncs = annotate_eqtl("junctions", tissue, fdr)
    return genes, trans, exons, juncs


def print_summary(tissue, fdr=0.05):
    genes, trans, exons, juncs = extract_features(tissue, fdr)
    if tissue == "Caudate":
        w_mode = "w"
    else:
        w_mode = "a"
    statement = "Significant DE (lfsr < 0.05) in %s" % tissue
    with open("summarize_results.log", mode=w_mode) as f:
        print(statement, file=f)
        for variable in ["effect", "gene_id", "gencode_id"]:
            print(variable, file=f)
            gg = len(set(genes[variable]))
            tt = len(set(trans[variable]))
            ee = len(set(exons[variable]))
            jj = len(set(juncs[variable]))
            print("\nGene:\t\t%d\nTranscript:\t%d\nExon:\t\t%d\nJunction:\t%d\n" %
                  (gg, tt, ee, jj), file=f)


def get_eqtl_result_by_tissue(tissue, fdr=0.05):
    genes, trans, exons, juncs = extract_features(tissue, fdr)
    genes["feature_type"] = "Gene"
    trans["feature_type"] = "Transcript"
    exons["feature_type"] = "Exon"
    juncs["feature_type"] = "Junction"
    df = pd.concat([genes, trans, exons, juncs])
    df["feature_type"] = df.feature_type.astype("category")\
                           .cat.reorder_categories(["Gene", "Transcript", "Exon", "Junction"])
    df["region"] = tissue
    return df


def main():
    bigdata1 = []
    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        print_summary(tissue)
        data1 = get_eqtl_result_by_tissue(tissue)
        bigdata1.append(data1)
    df1 = pd.concat(bigdata1)
    cols = ["region", "gene_id", "variant_id", "gencode_id", "gene_name",
            "seqnames", "start", "end", "lfsr", "posterior_mean", "feature_type"]
    df1.sort_values(["region", "feature_type", "lfsr", "posterior_mean"])\
       .loc[:, cols]\
       .to_csv("BrainSeq_sexGenotypes_4features_3regions.txt.gz",
               sep='\t', index=False)
    ## Reproducibilty information
    session_info.show()


if __name__ == "__main__":
    main()

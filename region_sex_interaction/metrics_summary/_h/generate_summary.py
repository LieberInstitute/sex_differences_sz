"""
This script summarized DE results.
"""
import numpy as np
import pandas as pd
import session_info
from pyhere import here
from functools import lru_cache

@lru_cache()
def get_annotation(feature):
    """
    Get the annotation file based on feature.
    """
    feat_lt = {"genes": "gene", "transcripts": "tx", 
               "exons": "exon", "junctions": "jxn"}
    new_feature = feat_lt[feature]
    fn = here("input/counts/text_files_counts/_m",
              f"caudate/{new_feature}_annotation.txt")
    return pd.read_csv(fn, sep="\t")

    
@lru_cache()
def get_de(feature, comp, fdr):
    """
    Using feature and comp, select the correct DE results
    file and filter based on FDR (input).
    """
    df = pd.read_csv(f"../../_m/{feature}/diffExpr_{comp}_sex_full.txt",
                     sep="\t", index_col=0)\
           .rename(columns={"gene_id":"gencodeID",
                            "gencodeGeneID": "gencodeID",
                            "gene_name": "Symbol"})
    df["ensemblID"] = df.gencodeID.str.replace("\\..*", "", regex=True)
    return df[(df["adj.P.Val"] < fdr)].copy()


@lru_cache()
def annotate_de(feature, comp, fdr):
    """
    Annotate DE features
    """
    cols = ["Feature", "seqnames", "start", "end", "width", "gencodeID",
            "ensemblID", "Symbol", "logFC", "AveExpr", "t", "P.Value",
            "adj.P.Val", "SE"]
    return get_de(feature, comp, fdr)\
        .merge(get_annotation(feature),
               left_index=True, right_on="name")\
        .rename(columns={"name": "Feature"})\
        .loc[:, cols].sort_values("adj.P.Val")



@lru_cache()
def extract_features(comp, fdr):
    # Extract DE from mash model
    genes = annotate_de("genes", comp, fdr)
    trans = annotate_de("transcripts", comp, fdr)
    exons = annotate_de("exons", comp, fdr)
    juncs = annotate_de("junctions", comp, fdr)
    genes["Type"] = "Gene"; trans["Type"] = "Transcript"
    exons["Type"] = "Exon"; juncs["Type"] = "Junction"
    return genes, trans, exons, juncs


def annotate_chrom(df):
    df.loc[:, "Chrom_Type"] = "Other"
    df.loc[df["seqnames"].isin(["chrX", "chrY"]), "Chrom_Type"] = "Allosome"
    df.loc[df["seqnames"].str.contains("chr\d+"), "Chrom_Type"] = "Autosome"
    df.loc[df["seqnames"] == "chrM", "Chrom_Type"] = "Mitochondria"
    return df


def get_DEGs_result_by_comparison(comp, fdr=0.05):
    # Male = 1; Female = -1
    genes, trans, exons, juncs = extract_features(comp, fdr)
    df = pd.concat([genes, trans, exons, juncs])
    df["Type"] = df.Type.astype("category")\
                   .cat.reorder_categories(["Gene", "Transcript",
                                            "Exon", "Junction"])
    df["Comp"] = comp
    return annotate_chrom(df)


def print_summary(comp, fdr=0.05):
    df = get_DEGs_result_by_comparison(comp, fdr)
    if comp == "CvD":
        w_mode = "w"
    else:
        w_mode = "a"
    statement = f"Significant DE (FDR < 0.05) in {comp}"
    with open("summarize_results.log", mode=w_mode) as f:
        print(statement, file=f)
        genes = df[(df["Type"] == "Gene")].copy()
        trans = df[(df["Type"] == "Transcript")].copy()
        exons = df[(df["Type"] == "Exon")].copy()
        juncs = df[(df["Type"] == "Junction")].copy()
        for variable in ["Feature", "gencodeID"]:
            print(variable, file=f)
            gg = len(set(genes[variable]))
            tt = len(set(trans[variable]))
            ee = len(set(exons[variable]))
            jj = len(set(juncs[variable]))
            print(f"\nGene:\t\t{gg}\nTranscript:\t{tt}"+\
                  f"\nExon:\t\t{ee}\nJunction:\t{jj}\n", file=f)


def merge_data():
    bigdata = []
    for comp in ["CvD", "CvH", "DvH"]:
        print_summary(comp)
        data = get_DEGs_result_by_comparison(comp, 1)
        bigdata.append(data)
    return pd.concat(bigdata)


def main():
    df = merge_data()
    # Summary
    print("\nSummary:")
    gene = df[(df["Type"] == "Gene") & (df["adj.P.Val"] < 0.05)]
    print("Unique DEGs: %d" %
          gene.drop_duplicates(subset="gencodeID").shape[0])
    print(gene.groupby(["Comp"]).size())
    # Output
    df.sort_values(["Comp", "Type", "P.Value"])\
      .to_csv("differential_expression_region_interaction_sex_4features_full.txt.gz",
              sep='\t', index=False)
    df[(df["adj.P.Val"] < 0.05)]\
        .sort_values(["Comp", "Type", "P.Value"])\
        .to_csv("differential_expression_region_interaction_sex_4features.txt",
                sep='\t', index=False)
    # Session infomation
    session_info.show()


if __name__ == '__main__':
    main()

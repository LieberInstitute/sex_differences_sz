"""
This script summarized DE results.
"""
import numpy as np
import pandas as pd
import session_info
from pyhere import here
from functools import lru_cache

@lru_cache()
def get_annotation(feature, tissue):
    """
    Get the annotation file based on feature and tissue.
    """
    feat_lt = {"genes": "gene", "transcripts": "tx", 
               "exons": "exon", "junctions": "jxn"}
    new_feature = feat_lt[feature]
    fn = here("input/counts/text_files_counts/_m",
              f"{tissue.lower()}/{new_feature}_annotation.txt")
    return pd.read_csv(fn, sep="\t")

    
@lru_cache()
def get_de(feature, tissue, fdr):
    """
    Using feature and tissue, select the correct DE results
    file and filter based on FDR (input).
    """
    df = pd.read_csv(here(f"differential_expression/{tissue.lower()}/_m",
                          f"{feature}/diffExpr_maleVfemale_full.txt"),
                     sep="\t", index_col=0)\
           .rename(columns={"gene_id":"gencodeID",
                            "gencodeGeneID": "gencodeID",
                            "gene_name": "Symbol"})
    df["ensemblID"] = df.gencodeID.str.replace("\\..*", "", regex=True)
    return df[(df["adj.P.Val"] < fdr)].copy()


@lru_cache()
def annotate_de(feature, tissue, fdr):
    """
    Annotate DE features
    """
    cols = ["Feature", "seqnames", "start", "end", "width", "gencodeID",
            "ensemblID", "Symbol", "logFC", "AveExpr", "t", "P.Value",
            "adj.P.Val", "SE", "B"]
    return get_de(feature, tissue, fdr)\
        .merge(get_annotation(feature, tissue),
               left_index=True, right_on="name")\
        .rename(columns={"name": "Feature"})\
        .loc[:, cols].sort_values("adj.P.Val")



@lru_cache()
def extract_features(tissue, fdr):
    # Extract DE from mash model
    genes = annotate_de("genes", tissue, fdr)
    trans = annotate_de("transcripts", tissue, fdr)
    exons = annotate_de("exons", tissue, fdr)
    juncs = annotate_de("junctions", tissue, fdr)
    genes["Type"] = "Gene"; trans["Type"] = "Transcript"
    exons["Type"] = "Exon"; juncs["Type"] = "Junction"
    return genes, trans, exons, juncs


def print_summary(tissue, fdr=0.05):
    genes, trans, exons, juncs = extract_features(tissue, fdr)
    if tissue == "Caudate":
        w_mode = "w"
    else:
        w_mode = "a"
    statement = f"Significant DE (FDR < 0.05) in {tissue}"
    with open("summarize_results.log", mode=w_mode) as f:
        print(statement, file=f)
        for variable in ["Feature", "gencodeID", "Symbol"]:
            print(variable, file=f)
            gg = len(set(genes[variable]))
            tt = len(set(trans[variable]))
            ee = len(set(exons[variable]))
            jj = len(set(juncs[variable]))
            print(f"\nGene:\t\t{gg}\nTranscript:\t{tt}"+\
                  f"\nExon:\t\t{ee}\nJunction:\t{jj}\n", file=f)


def get_DEGs_result_by_tissue(tissue, fdr=0.05):
    genes, trans, exons, juncs = extract_features(tissue, fdr)
    df = pd.concat([genes, trans, exons, juncs])
    df["Type"] = df.Type.astype("category")\
                   .cat.reorder_categories(["Gene", "Transcript",
                                            "Exon", "Junction"])
    df["Tissue"] = tissue
    return df


def main():
    bigdata = []
    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        print_summary(tissue)
        data = get_DEGs_result_by_tissue(tissue, 1)
        bigdata.append(data)
    df = pd.concat(bigdata)
    # Summary
    print("\nSummary:")
    gene = df[(df["Type"] == "Gene") & (df["adj.P.Val"] < 0.05)]
    gene.loc[:, "Direction"] = np.sign(gene.logFC) # Male = 1; Female = -1
    print("Unique DEGs: %d" %
          gene.drop_duplicates(subset="gencodeID").shape[0])
    print(gene.groupby(["Tissue", "Direction"]).size())
    # Output
    df.sort_values(["Tissue", "Type", "P.Value"])\
      .to_csv("differential_expression_analysis_4features_sex.txt.gz",
              sep='\t', index=False)
    # Session infomation
    session_info.show()


if __name__ == '__main__':
    main()

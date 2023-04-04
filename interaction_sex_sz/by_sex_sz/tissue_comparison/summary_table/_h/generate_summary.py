"""
This script summarized DE results.
"""
import numpy as np
import pandas as pd
import session_info
from pyhere import here
from functools import lru_cache
from scipy.stats import mannwhitneyu

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
def get_resid(feature, tissue):
    """
    Using feature and tissue to select residualized expression.
    """
    return pd.read_csv(here("interaction_sex_sz/interaction_model",
                            f"{tissue.lower()}/_m/{feature}",
                            "residualized_expression.tsv"),
                       sep="\t", index_col=0)\
             .transpose()


@lru_cache()
def get_pheno():
    return pd.read_csv(here("input/phenotypes/_m/phenotypes.csv"))\
             .loc[:, ["RNum", "Sex", "Dx", "Race", "Age"]]


@lru_cache()
def merge_resid(feature, tissue):
    return pd.merge(get_pheno(), get_resid(feature, tissue),
                    left_on="RNum", right_index=True)


@lru_cache()
def subset_resid(tissue, feature, sex):
    sex_dict = {"Male": "F", "Female": "M"}
    df = merge_resid(feature, tissue)
    ctl = df[(df["Dx"] == "Control") &
             (df["Sex"] == sex_dict[sex])].copy()
    sz = df[(df["Dx"] == "SCZD") &
            (df["Sex"] == sex_dict[sex])].copy()
    return ctl, sz


@lru_cache()
def get_de(feature, tissue, sex, fdr):
    """
    Using feature and tissue, select the correct DE results
    file and filter based on FDR (input).
    """
    df = pd.read_csv(here("interaction_sex_sz/by_sex_sz",
                          f"{tissue.lower()}/{sex.lower()}_analysis/_m",
                          f"{feature}/diffExpr_szVctl_full.txt"),
                     sep="\t", index_col=0)\
           .rename(columns={"gene_id":"gencodeID",
                            "gencodeGeneID": "gencodeID",
                            "gene_name": "Symbol"})
    df["ensemblID"] = df.gencodeID.str.replace("\\..*", "", regex=True)
    return df[(df["adj.P.Val"] < fdr)].copy()


@lru_cache()
def annotate_de(feature, tissue, sex, fdr):
    """
    Annotate DE features
    """
    cols = ["Feature", "seqnames", "gencodeID", "ensemblID",
            "Symbol", "logFC", "P.Value", "adj.P.Val", "SE"]
    return get_de(feature, tissue, sex, fdr)\
        .merge(get_annotation(feature, tissue),
               left_index=True, right_on="name")\
        .rename(columns={"name": "Feature"})\
        .loc[:, cols].sort_values("adj.P.Val")


@lru_cache()
def annotate_pvals(feature, tissue, sex, fdr):
    if fdr == 1:
        return annotate_de(feature, tissue, sex, fdr)
    else:
        df = annotate_de(feature, tissue, sex, fdr)
        ctl, sz = subset_resid(tissue, feature, sex)
        pval_df = []
        for feature_id in df.Feature:
            stat, pval = mannwhitneyu(ctl[feature_id], sz[feature_id])
            pval_df.append(pval)
        df["other_pval"] = pval_df
        return df


@lru_cache()
def extract_de(feature, tissue, fdr):
    ff = annotate_pvals(feature, tissue, "Female", fdr)
    mm = annotate_pvals(feature, tissue, "Male", fdr)
    ff["Sex"] = "Female"; mm["Sex"] = "Male"
    return ff, mm


@lru_cache()
def get_unique(feature, tissue, fdr):
    ff, mm = extract_de(feature, tissue, fdr)
    ff_genes = list(set(ff.Feature) - set(mm.Feature))
    mm_genes = list(set(mm.Feature) - set(ff.Feature))
    ff = ff[(ff["Feature"].isin(ff_genes))].copy()
    mm = mm[(mm["Feature"].isin(mm_genes))].copy()
    return ff, mm


@lru_cache()
def filter_pvals(feature, tissue, fdr):
    if fdr == 0.05:
        ff, mm = get_unique(feature, tissue, fdr)
        ff = ff[(ff["other_pval"] < fdr)].copy()
        mm = mm[(mm["other_pval"] < fdr)].copy()
    else:
        ff, mm = extract_de(feature, tissue, fdr)
    return pd.concat([ff, mm], axis=0)


@lru_cache()
def extract_features(tissue, fdr):
    # Extract DE from mash model
    genes = filter_pvals("genes", tissue, fdr)
    trans = filter_pvals("transcripts", tissue, fdr)
    exons = filter_pvals("exons", tissue, fdr)
    juncs = filter_pvals("junctions", tissue, fdr)
    genes["Type"] = "Gene"; trans["Type"] = "Transcript"
    exons["Type"] = "Exon"; juncs["Type"] = "Junction"
    return genes, trans, exons, juncs


def annotate_chrom(df):
    df.loc[:, "Chrom_Type"] = "Other"
    df.loc[df["seqnames"].isin(["chrX", "chrY"]), "Chrom_Type"] = "Allosome"
    df.loc[df["seqnames"].str.contains("chr\d+"), "Chrom_Type"] = "Autosome"
    df.loc[df["seqnames"] == "chrM", "Chrom_Type"] = "Mitochondria"
    return df


def get_DEGs_result_by_tissue(tissue, fdr=0.05):
    # Schizophrenia = 1; Control = -1
    genes, trans, exons, juncs = extract_features(tissue, fdr)
    df = pd.concat([genes, trans, exons, juncs])
    df.loc[:, "Direction"] = np.sign(df.logFC)
    df["Type"] = df.Type.astype("category")\
                   .cat.reorder_categories(["Gene", "Transcript",
                                            "Exon", "Junction"])
    df["Tissue"] = tissue
    return annotate_chrom(df)


def print_summary(tissue, fdr=0.05):
    df = get_DEGs_result_by_tissue(tissue, fdr)
    if tissue == "Caudate":
        w_mode = "w"
    else:
        w_mode = "a"
    statement = f"Significant DE (FDR < 0.05) in {tissue}"
    with open("summarize_results.log", mode=w_mode) as f:
        print(statement, file=f)
        for sex in ["Female", "Male"]:
            print(sex, file=f)
            genes = df[(df["Sex"] == sex) &
                       (df["Type"] == "Gene")].copy()
            trans = df[(df["Sex"] == sex) &
                       (df["Type"] == "Transcript")].copy()
            exons = df[(df["Sex"] == sex) &
                       (df["Type"] == "Exon")].copy()
            juncs = df[(df["Sex"] == sex) &
                       (df["Type"] == "Junction")].copy()
            for variable in ["Feature", "gencodeID"]:
                print(variable, file=f)
                gg = len(set(genes[variable]))
                tt = len(set(trans[variable]))
                ee = len(set(exons[variable]))
                jj = len(set(juncs[variable]))
                print(f"\nGene:\t\t{gg}\nTranscript:\t{tt}"+\
                      f"\nExon:\t\t{ee}\nJunction:\t{jj}\n", file=f)


def merge_data():
    bigdata1 = []; bigdata2 = [];
    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        print_summary(tissue)
        data2 = get_DEGs_result_by_tissue(tissue, 0.05)
        data1 = get_DEGs_result_by_tissue(tissue, 1)
        bigdata1.append(data1); bigdata2.append(data2)
    return pd.concat(bigdata1), pd.concat(bigdata2)


def main():
    df1, df2 = merge_data()
    # Summary
    print("\nSummary:")
    gene = df2[(df2["Type"] == "Gene")].copy()
    ff = gene[(gene["Sex"] == "Female")].copy()
    mm = gene[(gene["Sex"] == "Male")].copy()
    print("Unique DEGs (Female): %d" %
          ff.drop_duplicates(subset="gencodeID").shape[0])
    print("Unique DEGs (Male): %d" %
          mm.drop_duplicates(subset="gencodeID").shape[0])
    print(gene.groupby(["Tissue", "Sex", "Direction"]).size())
    # Output
    ## Full data
    df1.sort_values(["Tissue", "Sex", "Type", "P.Value"])\
       .to_csv("differential_expression_schizophrenia_by_sex_4features.txt.gz",
               sep='\t', index=False)
    ## FDR significant
    df2.sort_values(["Tissue", "Sex", "Type", "P.Value"])\
       .to_csv("differential_expression_schizophrenia_by_sex_4features.sig.txt.gz",
               sep='\t', index=False)
    # Session infomation
    session_info.show()


if __name__ == '__main__':
    main()

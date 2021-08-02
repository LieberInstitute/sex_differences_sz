## Generate supplementary data for DEGs
import pandas as pd

def map_tissue(tissue):
    return {"Caudate": "caudate", "DLPFC": "dlpfc",
            "Hippocampus": "hippocampus"}[tissue]


def annotate_leafcutter(tissue):
    fn1 = "../../../%s/localsplicing/_m/" % map_tissue(tissue) +\
        "leafcutter_ds_cluster_significance.txt"
    fn2 = "../../../%s/localsplicing/visualization/_m/" % map_tissue(tissue) +\
        "cluster_ds_results_annotated.txt"
    clusters = pd.read_csv(fn1, sep='\t')
    clusters[["chr", "clusterID"]] = \
        clusters.cluster.str.split(":", expand=True)
    annot = pd.read_csv(fn2, sep='\t').drop("gene", axis=1)
    df = pd.merge(clusters[["clusterID", "cluster", "genes", "loglr", "p",
                            "p.adjust"]],
                  annot, on="clusterID", how='left')
    df["Tissue"] = tissue
    return df


def get_tissues_DEG(tissue):
    cols = ["Feature", "gencodeID", "ensemblID", "Symbol", "logFC",
            "AveExpr", "t", "P.Value", "adj.P.Val", "Type"]
    gg = pd.read_csv("../../../%s/_m/genes/diffExpr_maleVfemale_full.txt" %
                     map_tissue(tissue), sep='\t', index_col=0)
    gg["Feature"] = gg.index; gg["Type"] = "Gene"
    tt = pd.read_csv("../../../%s/_m/transcripts/diffExpr_maleVfemale_full.txt" %
                     map_tissue(tissue), sep='\t', index_col=0)\
           .rename(columns={"gene_id": "gencodeID", "gene_name": "Symbol"})
    tt["ensemblID"] = tt.gencodeID.str.replace("\\..*", "", regex=True)
    tt["Feature"] = tt.index; tt["Type"] = "Transcript"
    ee = pd.read_csv("../../../%s/_m/exons/diffExpr_maleVfemale_full.txt" %
                     map_tissue(tissue), sep='\t', index_col=0)
    ee["Feature"] = ee.index; ee["Type"] = "Exon"
    jj = pd.read_csv("../../../%s/_m/junctions/diffExpr_maleVfemale_full.txt" %
                     map_tissue(tissue), sep='\t', index_col=0)\
           .drop(["Symbol"], axis=1)\
           .rename(columns={"newGeneID": "gencodeID", "newGeneSymbol": "Symbol"})
    jj["ensemblID"] = jj.gencodeID.str.replace("\\..*", "", regex=True)
    jj["Feature"] = jj.index; jj["Type"] = "Junction"
    df = pd.concat([gg.reset_index().loc[:, cols],
                    tt.reset_index().loc[:, cols],
                    ee.reset_index().loc[:, cols],
                    jj.reset_index().loc[:, cols]], axis=0)
    df["Tissue"] = tissue
    return df


def main():
    # Get data
    feature_df = []; leafcutter_df = []
    for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
        df1 = get_tissues_DEG(tissue)
        df2 = annotate_leafcutter(tissue)
        feature_df.append(df1)
        leafcutter_df.append(df2)
    sex_df = pd.concat(feature_df, axis=0)
    splicing_df = pd.concat(leafcutter_df, axis=0)
    print(sex_df.shape)
    print(sex_df.groupby(["Tissue", "Type"]).size())
    print("Significant (FDR < 0.05):")
    print(sex_df[(sex_df["adj.P.Val"] < 0.05)].groupby(["Tissue", "Type"]).size())
    # Save file
    sex_df.to_csv("differential_expression_analysis_4features_sex.txt.gz",
                  sep='\t', index=False)
    splicing_df.to_csv("differential_splicing_analysis_3brainregions_sex.txt.gz",
                       sep='\t', index=False)


if __name__ == '__main__':
    main()

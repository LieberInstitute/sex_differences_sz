## This script extracts XCI and generates a heatmap.

import session_info
import polars as pl
import seaborn as sns
from pyhere import here

def get_degs(region):
    fn = here(f"differential_expression/{region.lower()}",
              "_m/genes/diffExpr_maleVfemale_full.txt")
    return pl.read_csv(fn, separator="\t")\
             .filter(pl.col("adj.P.Val") < 0.05)\
             .rename({"": "feature_id", "gencodeID": "gencode_id",
                      "Symbol": "gene_name"})\
             .select([pl.col("feature_id"), pl.col("gencode_id"),
                      pl.col("gene_name"), pl.col("logFC"),
                      pl.col("adj.P.Val"), pl.col("SE")])


def get_xci(region):
    fn = here("xci_dosage_analysis/xci_enrichment/_h",
              "xci_status_hg19.txt")
    return pl.read_csv(fn, separator="\t")\
        .select([pl.col("Gene name"),pl.col("Combined XCI status")])\
        .rename({"Gene name": "gene_name",
                 "Combined XCI status": "xci_status"})\
        .join(get_degs(region), on=["gene_name"])


def get_res_df(region):
    fn = here(f"differential_expression/{region.lower()}",
              "_m/genes/residualized_expression.tsv")
    return pl.read_csv(fn, separator="\t")\
             .join(get_xci(region), on=["feature_id"],
                   how="semi")


def get_pheno(region):
    fn = here("input/phenotypes/_m/phenotypes.csv")
    if region == "Hippocampus":
        region = "HIPPO"
    return pl.read_csv(fn)\
             .filter((pl.col("Region") == region) &
                     (pl.col("Age") > 17) &
                     (pl.col("Dx").is_in(["SCZD", "Control"])))\
             .select(["RNum", "Dx", "Sex", "Age"]).unique()


def plotting_heatmap(region, width, height):
    res_df = get_res_df(region).to_pandas().set_index("feature_id")
    pheno_df = get_pheno(region)\
        .to_pandas().set_index("RNum").loc[res_df.columns, :]
    xci_df = get_xci(region).to_pandas().set_index("feature_id")
    xci = xci_df.pop("xci_status"); xmap = xci_df.pop("gene_name")
    sex = pheno_df.pop("Sex");
    lut = dict(zip(["F", "M"], "rbg"))
    lux = dict(zip(['escape', 'variable', 'inactive'],
                   sns.color_palette()))
    row_colors = sex.map(lut); col_colors = xci.map(lux)
    g = sns.clustermap(res_df.T, row_colors=row_colors,
                       col_colors=col_colors, cmap="mako",
                       yticklabels=False, figsize=(width,height),
                       cbar_pos=(0.03,0.5,0.03,0.2),
                       vmin=-5, vmax=5)
    g.ax_heatmap.axes.set_xticklabels(xmap)
    g.savefig(f"heatmap.xci_status.{region.lower()}.pdf")


def main():
    ## Plotting the heatmap
    plotting_heatmap("Caudate", 16, 10)
    plotting_heatmap("DLPFC", 10, 8)
    plotting_heatmap("Hippocampus", 14, 10)
    ## Reproducibility
    session_info.show()


if __name__ == "__main__":
    main()

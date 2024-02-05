## Examine RXE of male individuals in the hippocampus.

import session_info
import polars as pl
from pyhere import here
import statsmodels.api as sm
from scipy.stats import mannwhitneyu
from statsmodels.formula.api import ols

def get_pheno():
    fn = here("input/phenotypes/_m/phenotypes.csv")
    return pl.scan_csv(fn)\
             .filter((pl.col("Region") == "HIPPO") &
                     (pl.col("Sex") == "M"))\
             .unique()\
             .select(["RNum", "Sex", "Dx", "Age"])


def get_degs(region="Hippocampus"):
    fn = here(f"differential_expression/{region.lower()}",
              "_m/genes/diffExpr_maleVfemale_full.txt")
    return pl.scan_csv(fn, separator="\t")\
             .filter(pl.col("adj.P.Val") < 1)\
             .rename({"": "feature_id", "gencodeID": "gencode_id",
                      "Symbol": "gene_name"})\
             .select([pl.col("feature_id"), pl.col("gencode_id"),
                      pl.col("gene_name"), pl.col("logFC"),
                      pl.col("adj.P.Val"), pl.col("SE")])


def get_xci(region="Hippocampus"):
    fn = here("xci_dosage_analysis/xci_enrichment/_h",
              "xci_status_hg19.txt")
    return pl.scan_csv(fn, separator="\t")\
        .select([pl.col("Gene name"),pl.col("Combined XCI status")])\
        .rename({"Gene name": "gene_name",
                 "Combined XCI status": "xci_status"})\
        .join(get_degs(region), on=["gene_name"])


def xci_annotation():
    fn = here("input/counts/text_files_counts/_m",
              "hippocampus/gene_annotation.txt")
    return pl.scan_csv(fn, separator="\t")\
             .select(["name", "seqnames"])\
             .rename({"name": "feature_id",
                      "seqnames": "chrom"})\
             .join(get_xci(), on=["feature_id"])


def get_annotation():
    fn = here("input/counts/text_files_counts/_m",
              "hippocampus/gene_annotation.txt")
    return pl.scan_csv(fn, separator="\t")\
             .select(["name", "gene_name", "seqnames"])\
             .rename({"name": "feature_id",
                      "seqnames": "chrom"})


def get_logTPM():
    fn = here("input/counts/text_files_counts/tpm/_m",
              "hippocampus/gene.log2tpm.csv")
    return pl.scan_csv(fn)\
             .rename({"name": "feature_id"})\
             .join(get_annotation(), on=["feature_id"])\
             .with_columns(
                 pl.when(pl.col("chrom") == "chrX")\
                 .then(pl.lit("X"))\
                 .when(pl.col("chrom").str.contains("chr\d+"))\
                 .then(pl.lit("Autosome"))\
                 .otherwise(pl.lit("Other")).alias("chrom_type")
             )


def subset_data():
    samples = get_logTPM()\
        .select(pl.col("^.*(R).*$")).collect().columns
    return get_logTPM()\
        .with_columns(sum=pl.sum_horizontal("^.*(R).*$"))\
        .filter((pl.col("sum") >= (0.2 * len(samples))) &
                (pl.col("chrom_type").is_in(["X", "Autosome"])))\
        .select(pl.col("*")\
                .exclude("gene_name","chrom","sum"))


def cal_rxe(XCI, STATUS):
    if XCI:
        if STATUS is not None:
            xci = xci_annotation()\
                .filter(pl.col("xci_status") == STATUS)
            df = subset_data().join(xci, on="feature_id", how="anti")
        else:
            df = subset_data().join(xci_annotation(), on="feature_id", how="anti")
    else:
        df = subset_data()
    dx = df.collect()\
           .to_pandas()\
           .set_index("feature_id")\
           .groupby("chrom_type")\
           .mean(numeric_only=True)\
           .transpose().reset_index()
    return pl.DataFrame(dx)\
             .rename({"index": "RNum"})\
             .with_columns((pl.col("X") - pl.col("Autosome")).alias("RXE"))\
             .lazy()


def annot_samples(XCI=False, STATUS=None):
    return cal_rxe(XCI, STATUS).join(get_pheno(), on="RNum")


def test_diagnosis(dt):
    ctl = dt.filter(pl.col("Dx") == "Control").select("RXE").collect()
    sz  = dt.filter(pl.col("Dx") == "SCZD").select("RXE").collect()
    _, pval = mannwhitneyu(ctl, sz)
    return print(f'Mann-WhitneyU for ctl vs sz (RXE): {pval[0]:.6}')


def examine_inactive():
    xci = xci_annotation()\
        .filter(pl.col("xci_status").is_in(["escape", "variable"]))
    dx = subset_data().join(xci, on="feature_id", how="anti")\
                      .collect()\
                      .to_pandas()\
                      .set_index("feature_id")\
                      .groupby("chrom_type")\
                      .mean(numeric_only=True)\
                      .transpose().reset_index()
    return pl.DataFrame(dx)\
             .rename({"index": "RNum"})\
             .with_columns((pl.col("X") - pl.col("Autosome")).alias("RXE"))\
             .lazy().join(get_pheno(), on="RNum").collect()


def main():
    ## Full dataset
    print("Full dataset")
    test_diagnosis(annot_samples())
    ## All XCI
    print("Removing all XCI genes")
    test_diagnosis(annot_samples(True))
    ## XCI Status
    for STATUS in ["escape", "variable", "inactive"]:
        print(STATUS.upper())
        test_diagnosis(annot_samples(True, STATUS))
    ## Examine inactive
    print("RXE")
    examine_inactive().group_by("Dx")\
                      .agg(median=pl.median("RXE"),
                           mean=pl.mean("RXE"),
                           sd=pl.col("RXE").std(),
                           n=pl.count("RXE"))
    print("X")
    examine_inactive().group_by("Dx")\
                      .agg(median=pl.median("X"),
                           mean=pl.mean("X"),
                           sd=pl.col("X").std(),
                           n=pl.count("X"))
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()

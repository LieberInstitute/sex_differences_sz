#### This script compares GENCODE version 25 with version 41.

import session_info
import polars as pl
import seaborn as sns
from pyhere import here
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

def load_annotation(fn, region, sex):
    return pl.read_csv(fn, separator="\t")\
             .filter((pl.col("Type").is_in(["Gene", "gene"])) &
                     (pl.col("Tissue") == region) &
                     (pl.col("Sex") == sex))\
             .select(pl.col(["gencodeID", "ensemblID",
                             "logFC", "adj.P.Val"]))

    
def get_signif(fn, region, sex):
    return load_annotation(fn, region, sex)\
        .filter(pl.col("adj.P.Val") < 0.05)


def merge_data(fnc, fn1, fn2, region, sex, variable):
    df1 = fnc(fn1, region, sex)
    df2 = fnc(fn2, region, sex)
    return df1.join(df2, on=variable, suffix="_v41")


def corr_effect(df):
    return spearmanr(df["logFC"], df["logFC_v41"])


def scatter_plot(df, region, sex):
    fn = f"{region.lower()}.effect_size."+\
        f"{sex.lower()}.between_annotations.pdf"
    sns.set_theme()
    sns.set_style("ticks")
    sns.set_context("paper")
    sns.set_palette("colorblind")
    ax = sns.regplot(x="logFC", y="logFC_v41", data=df)
    ax.set(xlabel="GENCODE (v25)", ylabel="GENCODE (v41)")
    plt.tight_layout()
    plt.savefig(fn)
    plt.close()


def main():
    v25_file = here("supp_data/old_version",
                    "differential_expression_schizophrenia_by_sex_4features.txt.gz")
    v41_file = "../../summary_table/_m/"+\
        "differential_expression_schizophrenia_by_sex_4features.txt.gz"
    with open("annotation_comparison_statistics.log", "w") as f:
        for region in ["Caudate", "DLPFC", "Hippocampus"]:
            print(region, file=f)
            for sex in ["Female", "Male"]:
                print(sex, file=f)
                ## Prepare data
                df1 = merge_data(load_annotation, v25_file,
                                 v41_file, region, sex, "gencodeID")
                df2 = merge_data(load_annotation, v25_file,
                                 v41_file, region, sex, "ensemblID")
                df3 = merge_data(get_signif, v25_file,
                                 v41_file, region, sex, "ensemblID")
                ## Calculate correlation
                rho1, pval1 = corr_effect(df1)
                print(f"GENE_ID:\t rho > {rho1:.3}\tR2 ~ {rho1**2:.3}\t"+\
                      f"p-value < {pval1:.2}", file=f)
                rho2, pval2 = corr_effect(df2)
                print(f"ENSEMBL_ID:\t rho > {rho2:.3}\tR2 ~ {rho2**2:.3}\t"+\
                      f"p-value < {pval2:.2}", file=f)
                rho3, pval3 = corr_effect(df3)
                print(f"ENSEMBL_ID (FDR < 0.05):\t rho > {rho3:.3}\tR2 ~ {rho3**2:.3}\t"+\
                      f"p-value < {pval3:.2}", file=f)
                ## Plot all features
                scatter_plot(df2, region, sex)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()

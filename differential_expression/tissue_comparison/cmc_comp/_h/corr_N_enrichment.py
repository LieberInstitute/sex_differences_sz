"""
This script examines the Spearman correlation of effect sizes
between brain regions.
"""
import session_info
import pandas as pd
from pyhere import here
from functools import lru_cache
from scipy.stats import spearmanr
from scipy.stats import fisher_exact
from rpy2.robjects import r, globalenv, pandas2ri

@lru_cache()
def get_de():
    fn = "../../summary_table/_m/"+\
        "differential_expression_analysis_4features_sex.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == "Gene")].set_index("ensemblID").copy()


@lru_cache()
def get_sig():
    return get_de()[(get_de()["adj.P.Val"] < 0.05)].copy()


@lru_cache()
def get_cmc(label):
    cmc_lt = {"MSSM":here("differential_expression/cmc_dlpfc/_m",
                          "mssm_penn_pitt_maleVfemale.tsv"),
              "NIMH":here("differential_expression/cmc_dlpfc/_m",
                          "nimh_hbcc_maleVfemale.tsv")}
    df = pd.read_csv(cmc_lt[label], sep='\t')
    df["ensemblID"] = df.Geneid.str.replace("\\..*", "",
                                            regex=True)
    return df.set_index("ensemblID")


@lru_cache()
def get_cmc_sig(label):
    return get_cmc(label)[(get_cmc(label)["adj.P.Val"] < 0.05)].copy()


@lru_cache()
def prepare_data(tissue, label):
    bs_df = get_de().loc[(get_de()["Tissue"] == tissue),
                         ["logFC", "adj.P.Val"]]
    cmc_df = get_cmc(label).loc[:, ["logFC", "adj.P.Val"]]
    return pd.merge(bs_df, cmc_df, left_index=True,
                    right_index=True, suffixes=("_bs", "_cmc"))


def corr_beta(tissue, label, SIG):
    df = prepare_data(tissue, label)
    if SIG:
        df = df[(df["adj.P.Val_bs"] < 0.05) &
                (df["adj.P.Val_cmc"] < 0.05)].copy()
    return spearmanr(df["logFC_bs"], df["logFC_cmc"])


def cal_fishers(tissue, label):
    df = prepare_data(tissue, label)
    table = [[sum((df["adj.P.Val_bs"]< 0.05) & ((df["adj.P.Val_cmc"]< 0.05))),
              sum((df["adj.P.Val_bs"]< 0.05) & ((df["adj.P.Val_cmc"]>=0.05)))],
             [sum((df["adj.P.Val_bs"]>=0.05) & ((df["adj.P.Val_cmc"]< 0.05))),
              sum((df["adj.P.Val_bs"]>=0.05) & ((df["adj.P.Val_cmc"]>=0.05)))]]
    # print(table)
    return fisher_exact(table)


def plotNsave_corr(tissue, label):
    pandas2ri.activate()
    df = prepare_data(tissue, label)
    df = df[(df["adj.P.Val_bs"] < 0.05) &
            (df["adj.P.Val_cmc"] < 0.05)].copy()
    globalenv['df'] = df
    globalenv['tissue'] = tissue
    globalenv['label'] = label
    r('''
    save_plot <- function(p, fn, w, h){
        for(ext in c('.png', '.pdf')){
            ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    ## Plotting
    xlab = "Effect Size (CMC DLPFC)"
    ylab = paste0("Effect Size (", tissue, ")")
    fn = tolower(paste("effectsize_scatter", tissue, label, sep="_"))
    pp = ggpubr::ggscatter(df, x="logFC_cmc", y="logFC_bs", add="reg.line",
                           xlab=xlab, ylab=ylab, size=1,
                           add.params=list(color="blue", fill="lightgray"),
                           conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                           cor.coeff.args=list(method="spearman", label.sep="\n"),
                           ggtheme=ggpubr::theme_pubr(base_size=15)) +
         ggpubr::font("xy.title", face="bold", size=16)
    save_plot(pp, fn, 6, 6)
    ''')


def main():
    ## Calculate rho
    with open("correlation_statistics.log", "w") as f:
        for label in ["MSSM", "NIMH"]:
            for tissue in ["Caudate", "DLPFC", "Hippocampus"]:
                print(f"{tissue} vs {label}:\n", file=f)
                plotNsave_corr(tissue, label)
                odds, pval2 = cal_fishers(tissue, label)
                print(f"Enrichment:\t Odds Ratio > {odds:.3}\t"+\
                      f"p-value < {pval2:.2}", file=f)
                for SIG in [True, False]:
                    print(f"Significant DEGs: {SIG}", file=f)
                    rho, pval1 = corr_beta(tissue, label, SIG)
                    print(f"Correlation:\t rho > {rho:.3}\tR2 ~ {rho**2:.3}\t"+\
                          f"p-value < {pval1:.2}", file=f)
                
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()

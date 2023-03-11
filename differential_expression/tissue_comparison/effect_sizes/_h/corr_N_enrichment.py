"""
This script examines the Spearman correlation of effect sizes
between brain regions.
"""
import session_info
import pandas as pd
from functools import lru_cache
from scipy.stats import spearmanr
from scipy.stats import fisher_exact
from rpy2.robjects.packages import importr
from rpy2.robjects import r, globalenv, pandas2ri

@lru_cache()
def get_de():
    fn = "../../summary_table/_m/"+\
        "differential_expression_analysis_4features_sex.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == "Gene")].set_index("Feature").copy()


@lru_cache()
def get_sig():
    return get_de()[(get_de()["adj.P.Val"] < 0.05)].copy()


@lru_cache()
def prepare_data(tissue1, tissue2):
    df1 = get_sig().loc[(get_sig()["Tissue"] == tissue1), "logFC"]
    df2 = get_sig().loc[(get_sig()["Tissue"] == tissue2), "logFC"]
    sig_genes = list(set(df1.index) | set(df2.index))
    df = get_de().loc[sig_genes, ["Tissue", "logFC"]]
    return df[(df["Tissue"].isin([tissue1, tissue2]))]\
        .reset_index()\
        .pivot(index="Feature", columns="Tissue", values="logFC")\
        .dropna()


@lru_cache()
def prepare_data(tissue1, tissue2):
    df1 = get_sig().loc[(get_sig()["Tissue"] == tissue1), "logFC"]
    df2 = get_sig().loc[(get_sig()["Tissue"] == tissue2), "logFC"]
    sig_genes = list(set(df1.index) | set(df2.index))
    df = get_de().loc[sig_genes, ["Tissue", "logFC"]]
    return df[(df["Tissue"].isin([tissue1, tissue2]))]\
        .reset_index()\
        .pivot(index="Feature", columns="Tissue", values="logFC")\
        .dropna()


@lru_cache()
def prepare_data_pvals(tissue1, tissue2):
    return get_de()[(get_de()["Tissue"].isin([tissue1, tissue2]))]\
        .reset_index()\
        .pivot(index="Feature", columns="Tissue", values="adj.P.Val")\
        .dropna()


def corr_beta(tissue1, tissue2):
    df = prepare_data(tissue1, tissue2)
    return spearmanr(df[tissue1], df[tissue2])


def cal_fishers(tissue1, tissue2):
    df = prepare_data_pvals(tissue1, tissue2)
    table = [[sum((df[tissue1]< 0.05) & ((df[tissue2]< 0.05))),
              sum((df[tissue1]< 0.05) & ((df[tissue2]>=0.05)))],
             [sum((df[tissue1]>=0.05) & ((df[tissue2]< 0.05))),
              sum((df[tissue1]>=0.05) & ((df[tissue2]>=0.05)))]]
    # print(table)
    return fisher_exact(table)


def plotNsave_corr(tissue1, tissue2):
    pandas2ri.activate()
    globalenv['df'] = prepare_data(tissue1, tissue2)
    globalenv['tissue1'] = tissue1
    globalenv['tissue2'] = tissue2
    r('''
    save_plot <- function(p, fn, w, h){
        for(ext in c('.png', '.pdf')){
            ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    ## Plotting
    xlab = paste0("Effect Size\n(", tissue1, ")")
    ylab = paste0("Effect Size\n(", tissue2, ")")
    fn = tolower(paste("effectsize_scatter", tissue1, tissue2, sep="_"))
    pp = ggpubr::ggscatter(df, x=tissue1, y=tissue2, add="reg.line", size=1,
                           xlab=xlab, ylab=ylab,panel.labs.font=list(face="bold"),
                           add.params=list(color="blue", fill="lightgray"),
                           conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                           cor.method="spearman", cor.coeff.args=list(label.sep="\n"),
                           ggtheme=ggpubr::theme_pubr(base_size=15), ncol=3)
    save_plot(pp, fn, 6, 6)
    ''')


def main():
    ## Calculate rho
    with open("rho_statistics.log", "w") as f:
        ## Caudate vs DLPFC
        plotNsave_corr("Caudate", "DLPFC")
        rho, pval1 = corr_beta("Caudate", "DLPFC")
        odds, pval2 = cal_fishers("Caudate", "DLPFC")
        print("Caudate vs DLPFC:\n", file=f)
        print(f"Correlation:\t rho > {rho:.3}\tR2 ~ {rho**2:.3}\t"+\
              f"p-value < {pval1:.2}", file=f)
        print(f"Enrichment:\t Odds Ratio > {odds:.3}\t"+\
              f"p-value < {pval2:.2}", file=f)
        ## Caudate vs Hippocampus
        plotNsave_corr("Caudate", "Hippocampus")
        rho, pval1 = corr_beta("Caudate", "Hippocampus")
        odds, pval2 = cal_fishers("Caudate", "Hippocampus")
        print("Caudate vs Hippocampus:\n", file=f)
        print(f"Correlation:\t rho > {rho:.3}\tR2 ~ {rho**2:.3}\t"+\
              f"p-value < {pval1:.2}", file=f)
        print(f"Enrichment:\t Odds Ratio > {odds:.3}\t"+\
              f"p-value < {pval2:.2}", file=f)
        ## DLPFC vs Hippocampus
        plotNsave_corr("DLPFC", "Hippocampus")
        rho, pval1 = corr_beta("DLPFC", "Hippocampus")
        odds, pval2 = cal_fishers("DLPFC", "Hippocampus")
        print("DLPFC vs Hippocampus:\n", file=f)
        print(f"Correlation:\t rho > {rho:.3}\tR2 ~ {rho**2:.3}\t"+\
              f"p-value < {pval1:.2}", file=f)
        print(f"Enrichment:\t Odds Ratio > {odds:.3}\t"+\
              f"p-value < {pval2:.2}", file=f)
        ## Session information
        session_info.show()


if __name__ == "__main__":
    main()

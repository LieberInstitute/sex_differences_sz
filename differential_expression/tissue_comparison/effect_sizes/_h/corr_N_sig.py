"""
This script examines the Spearman correlation of effect sizes
between brain regions.
"""
import session_info
import pandas as pd
from functools import lru_cache
from scipy.stats import spearmanr
from rpy2.robjects.packages import importr
from rpy2.robjects import r, globalenv, pandas2ri

@lru_cache()
def get_de(feature):
    fn = "../../summary_table/_m/"+\
        "differential_expression_analysis_4features_sex.txt.gz"
    df = pd.read_csv(fn, sep='\t')
    return df[(df["Type"] == feature)].set_index("Feature").copy()


@lru_cache()
def get_sig(feature):
    return get_de(feature)[(get_de(feature)["adj.P.Val"] < 0.05)].copy()


@lru_cache()
def prepare_data(feature, tissue1, tissue2):
    df1 = get_sig(feature).loc[(get_sig(feature)["Tissue"] == tissue1),
                               ["logFC", "adj.P.Val"]]
    df2 = get_sig(feature).loc[(get_sig(feature)["Tissue"] == tissue2),
                               ["logFC", "adj.P.Val"]]
    return pd.merge(df1, df2, left_index=True, right_index=True,
                    suffixes=(f"_{tissue1}", f"_{tissue2}"))


def corr_beta(feature, tissue1, tissue2):
    df = prepare_data(feature, tissue1, tissue2)
    df = df[(df[f"adj.P.Val_{tissue1}"] < 0.05) &
            (df[f"adj.P.Val_{tissue2}"] < 0.05)].copy()
    return spearmanr(df[f"logFC_{tissue1}"],
                     df[f"logFC_{tissue2}"])


def plotNsave_corr(feature, tissue1, tissue2):
    pandas2ri.activate()
    df = prepare_data(feature, tissue1, tissue2)
    df = df[(df[f"adj.P.Val_{tissue1}"] < 0.05) &
            (df[f"adj.P.Val_{tissue2}"] < 0.05)].copy()
    globalenv['df'] = df
    globalenv['tissue1'] = tissue1
    globalenv['tissue2'] = tissue2
    globalenv['feature'] = feature
    r('''
    save_plot <- function(p, fn, w, h){
        for(ext in c('.png', '.pdf')){
            ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    ## Plotting
    xlab = paste0("Effect Size (", tissue1, ")")
    ylab = paste0("Effect Size (", tissue2, ")")
    if(!dir.exists("sig")){dir.create("sig")}
    outdir = paste0("sig/", feature)
    fn = tolower(paste(outdir, "effectsize_scatter", tissue1, tissue2, sep="."))
    pp = ggpubr::ggscatter(df, x=paste0("logFC_",tissue1), y=paste0("logFC_",tissue2),
                           add="reg.line", xlab=xlab, ylab=ylab, size=1,
                           add.params=list(color="blue", fill="lightgray"),
                           conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                           cor.coeff.args=list(method="spearman", label.sep="\n"),
                           ggtheme=ggpubr::theme_pubr(base_size=15)) +
         ggpubr::font("xy.title", face="bold", size=16)
    save_plot(pp, fn, 6, 6)
    ''')


def main():
    ## Calculate rho
    with open("sig/correlation_statistics.log", "w") as f:
        for feature in ["Gene", "Transcript", "Exon", "Junction"]:
            print(feature, file=f)
            ## Caudate vs DLPFC
            plotNsave_corr(feature, "Caudate", "DLPFC")
            rho, pval1 = corr_beta(feature, "Caudate", "DLPFC")
            print("Caudate vs DLPFC:\n", file=f)
            print(f"Correlation:\t rho > {rho:.3}\tR2 ~ {rho**2:.3}\t"+\
                  f"p-value < {pval1:.2}", file=f)
            ## Caudate vs Hippocampus
            plotNsave_corr(feature, "Caudate", "Hippocampus")
            rho, pval1 = corr_beta(feature, "Caudate", "Hippocampus")
            print("Caudate vs Hippocampus:\n", file=f)
            print(f"Correlation:\t rho > {rho:.3}\tR2 ~ {rho**2:.3}\t"+\
                  f"p-value < {pval1:.2}", file=f)
            ## DLPFC vs Hippocampus
            plotNsave_corr(feature, "DLPFC", "Hippocampus")
            rho, pval1 = corr_beta(feature, "DLPFC", "Hippocampus")
            print("DLPFC vs Hippocampus:\n", file=f)
            print(f"Correlation:\t rho > {rho:.3}\tR2 ~ {rho**2:.3}\t"+\
                  f"p-value < {pval1:.2}", file=f)
    ## Session information
    session_info.show()


if __name__ == "__main__":
    main()

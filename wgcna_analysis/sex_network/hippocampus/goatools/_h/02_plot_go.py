"""
This script is used to plot GO analysis produced by GOATOOLS.
It communicates with R, so rpy2 needs to be installed.
"""
import numpy as np
import pandas as pd
from glob import glob
import session_info, re
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def limit_go_name(s, maxlen=45):
    if len(s) > maxlen:
        return s[:maxlen - 3] + "..."
    else:
        return s


def go_df_for_plotting(df, name):
    df = df[df['enrichment']=='e'].copy()
    df["Log10"] = -np.log10(df['p_fdr_bh'])
    df['Feature'] = name
    df['prettyname'] = df['name'].apply(limit_go_name)
    fac = []
    for ii in range(df.shape[0]):
        xx, yy = df[['ratio_in_study']].iloc[ii, 0].split('/')
        zz, tt = df[['ratio_in_pop']].iloc[ii, 0].split('/')
        fac.append((int(xx) / int(yy)) / (int(zz) / int(tt)))
    df['geneRatio'] = fac
    return df.sort_values('p_uncorrected')


def plot_go(df, name, filename):
    godf = go_df_for_plotting(df, name).sort_values('p_uncorrected').head(15)
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_godf = ro.conversion.py2rpy(godf)
    ro.globalenv['r_godf'] = r_godf
    ro.globalenv['r_filename'] = filename
    ro.r("""
    library(dplyr)
    library(ggplot2)
    save_plot <- function(p, fn, w, h){
        for(ext in c('.svg', '.pdf')){
            ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    plot_GO <- function(df){
        cbPalette <- ggpubr::get_palette(palette = "npg", 4)
        gg1 = ggplot(df, aes(x=Log10, y=prettyname, size=geneRatio)) +
            geom_point(shape=18, alpha=0.8, col='#f8766d') +
            labs(y='', x='-Log10 (FDR)') + facet_grid('.~Feature') +
            geom_vline(xintercept = -log10(0.05), linetype = "dotted") +
            scale_size_continuous(range = c(2.5, 12)) +
            theme_bw(base_size=15) +
            theme(axis.title=element_text(size=18, face='bold'),
                  strip.text=element_text(size=18, face='bold'))
        return(gg1)
    }
    ## Plotting
    df = r_godf
    df$prettyname <- factor(df$prettyname,
                            levels=unique(df$prettyname[order(df$Log10, df$p_uncorrected, df$name, decreasing=FALSE)]))
    gg = plot_GO(df); save_plot(gg, r_filename, 8, 6)
    """)


def main():
    for fn in glob('GO_analysis_mash_*.xlsx'):
        m = re.search(r'(\w+)', re.split(r'_', fn)[3])
        module_number = m.group(0)
        name = "Module %s" % module_number
        filename = 'module_%s_go_enrichment' % module_number
        df = pd.read_excel(fn)
        plot_go(df, name, filename)
        print(filename)
    ## Session information
    session_info.show()


if __name__ == '__main__':
    main()

"""
This script is used to plot GO analysis produced by GOATOOLS.
It communicates with R, so rpy2 needs to be installed.
"""

import re, glob
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def limit_go_name(s, maxlen=45):
    if len(s) > maxlen:
        return s[:maxlen - 3] + '...'
    else:
        return s


def go_df_for_plotting(df, name):
    df = df[df['enrichment']=='e'].copy()
    df['Log10'] = -np.log10(df['p_fdr_bh'])
    df['Feature'] = name
    df['prettyname'] = df['name'].apply(limit_go_name)
    fac = []
    for ii in range(df.shape[0]):
        xx, yy = df[['ratio_in_study']].iloc[ii, 0].split('/')
        zz, tt = df[['ratio_in_pop']].iloc[ii, 0].split('/')
        fac.append((int(xx) / int(yy)) / (int(zz) / int(tt)))
    df['geneRatio'] = fac
    return df.drop(columns=['study_items']).sort_values('p_uncorrected')



def plot_go(df, name, filename):
    godf = go_df_for_plotting(df, name).sort_values('p_uncorrected').head(15)
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_godf = ro.conversion.py2rpy(godf)
    ro.globalenv['r_godf'] = r_godf
    ro.globalenv['r_filename'] = filename
    ro.r("""
    library(ggplot2)
    df1 = r_godf
    #df1$fac1 = -log2(df1[, 'geneRatio'])
    df1$prettyname <- factor(df1$prettyname,
                             levels=unique(df1$prettyname[order(df1$Log10, df1$p_uncorrected,
                                                                df1$name, decreasing=FALSE)]))
    gg1 = (ggplot(df1, aes(x=Log10, y=prettyname, size=geneRatio)) +
           geom_point(shape=18, col='#f8766d') +
           labs(y='', x='-log10(FDR)') + theme_bw() +
           facet_grid('.~Feature') +
           geom_vline(xintercept = -log10(0.05), linetype = "dotted") +
           theme(axis.text=element_text(size=14),
                 axis.title=element_text(size=18, face='bold'),
                 strip.text=element_text(size=18, face='bold'),
                 ))
    #print(r_filename)
    ggsave(file=paste(sep='', r_filename, '.pdf'), plot=gg1, width=8, height=6)
    ggsave(file=paste(sep='', r_filename, '.svg'), plot=gg1, width=8, height=6)
    ggsave(file=paste(sep='', r_filename, '.png'), plot=gg1, width=8, height=6)
    """)


def main():
    for fn in glob.glob('../../_m/GO_analysis_module_*.xlsx'):
        m = re.search('module(\w+)', fn)
        module_number = m.groups(1)
        name = "Module %s" % module_number
        filename = 'module%s_go_enrichment' % module_number
        df = pd.read_excel(fn)
        plot_go(df, name, filename)
        print(filename)


if __name__ == '__main__':
    main()

## Calculate enrichment for DE genes

import functools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

@functools.lru_cache()
def get_wgcna_modules():
    return pd.read_csv("../../_m/modules.csv", index_col=0)


@functools.lru_cache()
def get_degs():
    return set(pd.read_csv('../../../../../differential_expression/dlpfc/'+\
                           '_m/genes/diffExpr_maleVfemale_FDR05.txt',
                           sep='\t', usecols=[0], index_col=0).index)


@functools.lru_cache()
def get_mhc_genes():
    return set(pd.read_csv('../../../../../input/counts/mhc_region_genes/'+\
                           '_m/mhc_genes.csv')['gene_id'])


def fet(a, b, u):
    # a, b, u are sets
    # u is the universe
    yes_a = u.intersection(a)
    yes_b = u.intersection(b)
    no_a = u - a
    no_b = u - b
    m = [[len(yes_a.intersection(yes_b)), len(no_a.intersection(yes_b)) ],
         [len(yes_a.intersection(no_b)), len(no_a.intersection(no_b))]]
    return fisher_exact(m)


def enrichment_rows():
    mod = get_wgcna_modules().module.unique()
    u = set(get_wgcna_modules().index)
    for ii in range(len(mod)): # for each module
        a = set(get_wgcna_modules()[(get_wgcna_modules().module) == mod[ii]].index)
        b = set(get_wgcna_modules()[(get_wgcna_modules().module) == mod[ii]].index) - get_mhc_genes()
        yield (mod[ii],
               len(a),
               *fet(a, get_degs(), u),
               *fet(b, get_degs() - get_mhc_genes(), u))


def plot_heatmap(edf, label):
    df = edf.sort_values("N_Genes", ascending=False)
    df2 = np.log(df.loc[:, ['%s_OR' % label]]).replace([np.inf, -np.inf], 0)
    df2.columns = ['%s' % label]
    df2.index = ["Module %s (%d genes)" % (x,y) for x,y in zip(df2.index,
                                                               df['N_Genes'])]
    df3 = df.loc[:, ['%s_FDR' % label]]
    fig, ax = plt.subplots(figsize=(6,10))
    p = sns.heatmap(df2, cmap='coolwarm', annot=df3, yticklabels=df2.index,
                    center=0, cbar_kws={'label': 'Log10(Enrichment Ratio)'},
                    vmin=-2, vmax=2)
    p.set_title("Enrichment/depletion DE genes in WGCNA modules\n(FDR values)")
    p.get_figure().savefig('wgcna_module_enrichment_%s.pdf' % label,
                           bbox_inches='tight')


def main():
    edf = pd.DataFrame\
            .from_records(enrichment_rows(),
                          columns=['Module_ID', 'N_Genes', 'DEG_OR', 'DEG_P',
                                   'DEG_noMHC_OR', 'DEG_noMHC_P'],
                          index='Module_ID')
    edf['DEG_FDR'] = multipletests(edf['DEG_P'], method='fdr_bh')[1]
    edf['DEG_noMHC_FDR'] = multipletests(edf['DEG_noMHC_P'], method='fdr_bh')[1]
    edf = edf.loc[:, ['N_Genes', 'DEG_OR', 'DEG_P', 'DEG_FDR', 'DEG_noMHC_OR',
                      'DEG_noMHC_P', 'DEG_noMHC_FDR']]
    edf.to_csv('wgcna_module_enrichment.csv')
    plot_heatmap(edf, "DEG")
    plot_heatmap(edf, "DEG_noMHC")


if __name__ == "__main__":
    main()

"""
Examining the variation explained of gene expression for sex-specific
expression differences by:
  1. All DEGs
  2. Allosomes only
  3. Autosomes only

Also examine correlation for all DEGs accumentively.
"""

import numpy as np
import pandas as pd
from pyhere import here
import errno, os, argparse
from functools import lru_cache
from scipy.stats import linregress
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from rpy2.robjects import r, globalenv, pandas2ri

def mkdir_p(directory):
    """
    Make a directory if it does not already exist.

    Input: Directory name
    """
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


@lru_cache()
def get_pheno_data():
    """
    Load phenotype data.
    """
    fn = here('input/phenotypes/_m/phenotypes.csv')
    cols = ["SAMPLE_ID", "BrNum", "Dx", "Race", "Sex", "Age"]
    df = pd.read_csv(fn).loc[:, cols]
    df["RNum"] = df.SAMPLE_ID.str.replace("\\_.*", "", regex=True)
    return df.drop_duplicates()


@lru_cache()
def get_residualized(tissue):
    '''
    Load residualization file.
    '''
    new_tissue = tissue.lower()
    fn = here(f'differential_expression/{new_tissue}',
              '_m/genes/residualized_expression.tsv')
    return pd.read_csv(fn, sep='\t', index_col=0).transpose()


@lru_cache()
def get_deg(tissue):
    """
    Get DE features from limma-voom model
    """
    fn = here('differential_expression/tissue_comparison/summary_table/_m',
              'differential_expression_analysis_4features_sex.txt.gz')
    df = pd.read_csv(fn, sep='\t')
    return df[(df["adj.P.Val"] < 0.05) &
              (df["Type"] == "Gene") &
              (df["Tissue"] == tissue)].copy()


@lru_cache()
def get_autosomes(tissue):
    """
    Get annotated results and select autosomes
    """
    df = get_deg(tissue)
    return df[(df["Chrom_Type"] == "Autosome")].copy()


@lru_cache()
def get_allosomes(tissue):
    """
    Get annotated results and select allosomes
    """
    df = get_deg(tissue)
    return df[(df["Chrom_Type"] == "Allosome")].copy()


@lru_cache()
def get_deg_res_df(tissue, num, fnc, FILTER):
    """
    Merge DE and residualized data after selecting N features.
    """
    if FILTER:
        newList = list(fnc(tissue)\
                       .sort_values("P.Value").head(num).Feature)
    else:
        newList = list(fnc(tissue).Feature)
    return get_residualized(tissue)[newList]


def get_explained_variance(df):
    x = StandardScaler().fit_transform(df)
    pca = PCA(n_components=2).fit(x)
    pc1 = pca.explained_variance_ratio_[0]
    pc2 = pca.explained_variance_ratio_[1]
    print("Explained Variance\nPC1:\t%0.5f\nPC2:\t%0.5f" % (pc1, pc2))


def cal_pca(df):
    x = StandardScaler().fit_transform(df)
    pca = PCA(n_components=2).fit_transform(x)
    return pd.DataFrame(data=pca, columns=['PC1', 'PC2'], index=df.index)


def get_corr(dft):
    '''This calculates R^2 correlation via linear regression:
         - used to calculate relationship between 2 arrays
         - the arrays are principal components 1 or 2 (PC1, PC2) AND ancestry
         - calculated on a scale of 0 to 1 (with 0 being no correlation)
        Inputs:
            dft: Data frame with continuous variable and PCs
        Outputs:
          1. r2
          2. p-value, two-sided test
            - whose null hypothesis is that two sets of data are uncorrelated
          3. slope (beta): directory of correlations
    '''
    xx = dft.Sex.astype("category").cat.codes; yy = dft.PC1
    slope, intercept, r_value, p_value, std_err = linregress(xx, yy)
    return r_value**2, p_value


def get_pca_df(tissue, num, fnc, FILTER):
    '''
    new_pheno: This gets the correct size of samples using the the first two
               columns of residualized expression
      - the residualized expression data frame, has the correct samples
      - output new_pheno shape row numbers should be the same as res_df row numbers
    '''
    expr_res = get_deg_res_df(tissue, num, fnc, FILTER)
    pheno_df = get_pheno_data()
    # Generate pheno data frame with correct samples
    new_pheno = pheno_df.merge(expr_res.iloc[:, 0:1], left_on="RNum",
                               right_index=True)\
                        .drop(expr_res.iloc[:, 0:1].columns, axis=1)\
                        .set_index("RNum")
    principalDf = cal_pca(expr_res)
    ##get_explained_variance(expr_res)
    return pd.concat([principalDf, new_pheno], axis = 1)


def plotNsave_corr(tissue, fn, num, fnc, FILTER):
    df = get_pca_df(tissue, num, fnc, FILTER)
    rho, pval = get_corr(df)
    label = 'PC1 R2: %.2f\nP-value: %.2e' % (rho, pval)
    pandas2ri.activate()
    globalenv['df'] = df
    globalenv['fn'] = fn
    globalenv['ntitle'] = '\n'.join([label])
    r('''
    library(ggpubr)
    save_plot <- function(p, fn, w, h){
        for(ext in c('.png', '.pdf')){
            ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
        }
    }
    ## Plotting
    pp <- ggscatter(df, x="PC1", y="PC2", color="Sex", size=2,
                    alpha=0.75, legend="bottom", palette="npg",
                    ggtheme=theme_pubr(base_size=15, border=TRUE))
    pp <- ggpar(pp, title=ntitle) + \
        theme(plot.title = element_text(hjust = 0.5))
    save_plot(pp, fn, 6, 6)
    ''')


def pc_recursive_autosomes(tissue):
    new_tissue = tissue.lower()
    geneList = list(get_autosomes(tissue).Feature)
    pheno_df = get_pheno_data(); pvals = []; rsq = []; nums = []
    for num in range(2, len(geneList)+1):
        expr_res = get_deg_res_df(tissue, num, get_autosomes, True)
        new_pheno = pheno_df.merge(expr_res.iloc[:, 0:1],
                                   right_index=True, left_on="RNum")\
                            .drop(expr_res.iloc[:, 0:1].columns, axis=1)\
                            .set_index("RNum")
        dft = pd.concat([cal_pca(expr_res), new_pheno], axis=1)
        r2, pval = get_corr(dft)
        nums.append(num); pvals.append(pval); rsq.append(r2)
    pd.DataFrame({"DEGs":nums, "PValue":pvals, "Rsq": rsq})\
        .sort_values("DEGs", ascending=True)\
        .to_csv(f"{new_tissue}/variance_explained_rsq_recursive_autosomes.csv",
                index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tissue', type=str)
    args=parser.parse_args()
    # Script calling
    tissue = args.tissue; new_tissue = tissue.lower(); mkdir_p(new_tissue)
    with open(f"{new_tissue}/summarize_explained_variance.log", mode="w") as f:
        # Allosomes
        print("Explained Variance (%%) of PC1: Allosomes", file=f)
        r2, pval = get_corr(get_pca_df(tissue, 0, get_allosomes, False))
        print(f"All DEGs\tr2: {r2:.1%}\tp-value: {pval:.2}", file=f)
        plotNsave_corr(tissue, f"{new_tissue}/pca_allosomes_all", 0,
                       get_allosomes, False)
        r2, pval = get_corr(get_pca_df(tissue, 100, get_allosomes, True))
        print(f"Top 100 DEGs\tr2: {r2:.1%}\tp-value: {pval:.2}", file=f)
        plotNsave_corr(tissue, f"{new_tissue}/pca_allosomes_top100", 100,
                       get_allosomes, True)
        r2, pval = get_corr(get_pca_df(tissue, 10, get_allosomes, True))
        print(f"Top 10 DEGs\tr2: {r2:.1%}\tp-value: {pval:.2}", file=f)
        plotNsave_corr(tissue, f"{new_tissue}/pca_allosomes_top10", 10,
                       get_allosomes, True)
        # Autosomes
        print("Explained Variance (%%) of PC1: Autosomes", file=f)
        r2, pval = get_corr(get_pca_df(tissue, 0, get_autosomes, False))
        print(f"All DEGs\tr2: {r2:.1%}\tp-value: {pval:.2}", file=f)
        plotNsave_corr(tissue, f"{new_tissue}/pca_autosomes_all", 0,
                       get_autosomes, False)
        r2, pval = get_corr(get_pca_df(tissue, 100, get_autosomes, True))
        print(f"Top 100 DEGs\tr2: {r2:.1%}\tp-value: {pval:.2}", file=f)
        plotNsave_corr(tissue, f"{new_tissue}/pca_autosomes_top100", 100,
                       get_autosomes, True)
        r2, pval = get_corr(get_pca_df(tissue, 10, get_autosomes, True))
        print(f"Top 10 DEGs\tr2: {r2:.1%}\tp-value: {pval:.2}", file=f)
        plotNsave_corr(tissue, f"{new_tissue}/pca_autosomes_top10", 10,
                       get_autosomes, True)
        # Recursive run
        pc_recursive_autosomes(tissue)


if __name__ == '__main__':
    main()

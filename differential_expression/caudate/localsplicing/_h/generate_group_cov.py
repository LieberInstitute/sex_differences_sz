#!/usr/bin/python

import argparse
import functools
import pandas as pd
from rpy2.robjects import r, pandas2ri

@functools.lru_cache()
def get_pheno(fn):
    pheno = pd.read_csv(fn, index_col=0)
    return pheno[(pheno['Dx'].isin(['Control', 'Schizo'])) &
                 (pheno['Age'] > 17) &
                 (pheno['Race'].isin(['AA', 'CAUC']))]


@functools.lru_cache()
def get_counts(fn):
    return pd.read_csv(fn, sep='\t', nrows=500, index_col=0).T


@functools.lru_cache()
def merge_dataframe(fn1, fn2):
    return pd.merge(get_pheno(fn1), get_counts(fn2),
                    left_index=True, right_index=True)\
             .drop(get_counts(fn2).columns, axis=1)


def generate_covs(fn1, fn2):
    pandas2ri.activate()
    r('''
    voomDat = '../../_m/genes/voomSVA.RData'
    load(voomDat)
    dft = data.frame(v$design)
    dft['RNum'] = rownames(v$design)
    ''')
    sva_plus = r['dft'].set_index('RNum')
    pheno = merge_dataframe(fn1, fn2)
    new_pheno = pd.merge(pheno[['Sex', 'Dx', 'Age']],
                         sva_plus.iloc[:, 4:],
                         left_index=True, right_index=True)
    return new_pheno.sort_values('Sex')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pheno_file', type=str)
    parser.add_argument('--counts_file', type=str)
    args=parser.parse_args()
    covs = generate_covs(args.pheno_file,
                         args.counts_file)
    covs.to_csv('groups_cov_file.txt', sep=' ',
                header=False, index=True)


if __name__ == '__main__':
    main()

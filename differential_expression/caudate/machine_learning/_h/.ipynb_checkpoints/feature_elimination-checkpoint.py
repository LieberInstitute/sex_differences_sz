#!/usr/bin/env python
""" 
This package has several function to run feature elimination for random forest classifier and
logistic regression. For random forest classification, it assumes Out-of-Bag (OOB) is set to 
True. For both random forest and logistic regression, three measurements are calculated for 
feature selection:

1. Normalized mutual information
2. Accuracy
3. Area under the curve (AUC) ROC curve

The package has been split in to three additional scripts for:
1. Random forest feature elimination (AP)
2. Logistic regression feature elimination (KJB)
3. Rank features function (TK)

Original author Apuã Paquola (AP).
Edits and package management by Kynon Jade Benjamin (KJB)
Feature ranking modified from Tarun Katipalli (TK) ranking function.
"""

__author__ = 'Apuã Paquola'

import numpy as np
import pandas as pd
from plotnine import *
from random_forest import *
from logistic_regression import *
from warnings import filterwarnings
from matplotlib.cbook import mplDeprecation
from sklearn.metrics import balanced_accuracy_score
filterwarnings("ignore", category=mplDeprecation)
filterwarnings('ignore', category=UserWarning, module='plotnine.*')
filterwarnings('ignore', category=DeprecationWarning, module='plotnine.*')


def n_features_iter(nf, keep_rate):
    """
    Determines the features to keep.

    Args:
    nf: current number of features
    keep_rate: percentage of features to keep

    Yields:
    int: number of features to keep
    """    
    while nf != 1:
        nf = max(1, int(nf * keep_rate))
        yield nf


def feature_elimination(estimator, X, Y, features, fold, out_dir='.',
                        elimination_rate=0.2, RANK=True):
    """
    Runs random forest feature elimination step over iterator process.

    Args:
    estimator: Random forest classifier object
    X: a data frame of training data
    Y: a vector of sample labels from training data set
    features: a vector of feature names
    fold: current fold
    out_dir: output directory. default '.'
    elimination_rate: percent rate to reduce feature list. default .2

    Yields:
    dict: a dictionary with number of features, normalized mutual
          information score, accuracy score, auc roc curve and array of the indexes
          for features to keep
    """
    d = dict()
    pfirst = None
    keep_rate = 1-elimination_rate
    for p in rf_fe(estimator, X, Y, n_features_iter(X.shape[1], keep_rate), 
                   features, fold, out_dir, RANK):
        if pfirst is None:
            pfirst = p
        d[p[0]] = p

    return d, pfirst


def feature_elimination_lr(estimator, X, Y, features, fold, out_dir='.',
                           elimination_rate=0.2, dev_size=0.2, RANK=True,
                           SEED=True):
    """
    Runs random forest feature elimination step over iterator process.

    Args:
    estimator: Random forest classifier object
    X: a data frame of training data
    Y: a vector of sample labels from training data set
    features: a vector of feature names
    fold: current fold
    out_dir: output directory. default '.'
    elimination_rate: percent rate to reduce feature list. default .2

    Yields:
    dict: a dictionary with number of features, normalized mutual
          information score, accuracy score, auc roc curve and array of the indexes
          for features to keep
    """
    d = dict()
    pfirst = None
    keep_rate = 1-elimination_rate
    for p in lr_fe(estimator, X, Y, n_features_iter(X.shape[1], keep_rate), 
                   features, fold, out_dir, dev_size, SEED, RANK):
        if pfirst is None:
            pfirst = p
        d[p[0]] = p
    return d, pfirst


def save_plot(p, fn, width=7, height=7):
    '''Save plot as svg, png, and pdf with specific label and dimension.'''
    for ext in ['.svg', '.png', '.pdf']:
        p.save(fn+ext, width=width, height=height)
        

def plot_nmi(d, fold, output_dir):
    """
    Plot feature elimination results for normalized mutual information.

    Args:
    d: feature elimination class dictionary
    fold: current fold
    out_dir: output directory. default '.'

    Yields:
    graph: plot of feature by NMI, automatically saves files as png and svg
    """
    df_elim = pd.DataFrame([{'n features':k,
                             'normalized mutual information':d[k][1]} for k in d.keys()])
    gg = ggplot(df_elim, aes(x='n features', y='normalized mutual information'))\
        + geom_point() + scale_x_log10() + theme_light()
    save_plot(gg, output_dir+"/nmi_fold_%d" % (fold))
    print(gg)


def plot_roc(d, fold, output_dir):
    """
    Plot feature elimination results for AUC ROC curve.

    Args:
    d: feature elimination class dictionary
    fold: current fold
    out_dir: output directory. default '.'

    Yields:
    graph: plot of feature by AUC, automatically saves files as png and svg
    """
    df_elim = pd.DataFrame([{'n features':k,
                             'ROC AUC':d[k][3]} for k in d.keys()])
    gg = ggplot(df_elim, aes(x='n features', y='ROC AUC'))\
        + geom_point() + scale_x_log10() + theme_light()
    save_plot(gg, output_dir+"/roc_fold_%d" % (fold))
    print(gg)


def plot_acc(d, fold, output_dir):
    """
    Plot feature elimination results for accuracy.

    Args:
    d: feature elimination class dictionary
    fold: current fold
    out_dir: output directory. default '.'

    Yields:
    graph: plot of feature by accuracy, automatically saves files as png and svg
    """
    df_elim = pd.DataFrame([{'n features':k,
                             'Accuracy':d[k][3]} for k in d.keys()])
    gg = ggplot(df_elim, aes(x='n features', y='Accuracy'))\
        + geom_point() + scale_x_log10() + theme_light()
    save_plot(gg, output_dir+"/acc_fold_%d" % (fold))
    print(gg)


# def plot_scores(d, alpha, output_dir):
#     df_nmi = pd.DataFrame([{'n features':k, 'Score':d[k][1]} for k in d.keys()])
#     df_nmi['Type'] = 'NMI'
#     df_acc = pd.DataFrame([{'n features':k, 'Score':d[k][2]} for k in d.keys()])
#     df_acc['Type'] = 'Acc'
#     df_roc = pd.DataFrame([{'n features':k, 'Score':d[k][3]} for k in d.keys()])
#     df_roc['Type'] = 'ROC'
#     df_elim = pd.concat([df_nmi, df_acc, df_roc], axis=0)
#     gg = ggplot(df_elim, aes(x='n features', y='Score', color='Type'))\
#         + geom_point() + scale_x_log10() + theme_light()
#     gg.save(output_dir+"/scores_wgt_%.2f.png" % (alpha))
#     gg.save(output_dir+"/scores_wgt_%.2f.svg" % (alpha))
#     print(gg)

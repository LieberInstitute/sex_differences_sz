#!/usr/bin/env python
"""
This script contains the logistic regression modification of the original
random forest feature elimination package. Instead of Out-of-Bag, it creates
a developmental test set from the training data. 

Developed by Kynon Jade Benjamin.
Edits by Apuã Paquola.
"""

__author__ = 'Kynon J Benjamin'

import numpy as np
import pandas as pd
from itertools import chain
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
from rank_function import features_rank_fnc
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import normalized_mutual_info_score


class LogisticRegression_FI(LogisticRegression):
    """
    Add feature importance to Logistic Regression class similar to 
    random forest output.
    """
    def fit(self, *args, **kwargs):
        super(LogisticRegression_FI, self).fit(*args, **kwargs)
        self.feature_importances_ = np.abs(self.coef_).flatten()

        
def dev_predictions(estimator, X):
    """
    Extracts predictions using a development fold for logistic
    regression classifier classes.
    
    Args:
    estimator: Logistic regression classifier object
    X: a data frame of normalized values from developmental dataset
    
    Yields:
    vector: Development set predicted labels
    """
    return estimator.predict(X)
        
    
def dev_score_roc(estimator, X, Y):
    """
    Calculates the area under the ROC curve score
    for the develomental dataset predictions.

    Args:
    estimator: Logistic regression classifier object
    X: a data frame of normalized values from developmental dataset
    Y: a vector of sample labels from developmental data set

    Yields:
    float: AUC ROC score    
    """
    if len(np.unique(Y)) > 2:
        labels_pred = estimator.predict_proba(X) 
        kwargs = {'multi_class': 'ovr'}
    else:
        labels_pred = dev_predictions(estimator, X)
        kwargs = {}
    return roc_auc_score(Y, labels_pred, **kwargs)


def dev_score_nmi(estimator, X, Y):
    """
    Calculates the normalized mutual information score
    from the developmental dataset predictions.

    Args:
    estimator: Random forest classifier object
    X: a data frame of normalized values from developmental dataset
    Y: a vector of sample labels from developmental dataset

    Yields:
    float: normalized mutual information score
    """
    labels_pred = dev_predictions(estimator, X)
    return normalized_mutual_info_score(Y, labels_pred,
                                        average_method='arithmetic')


def dev_score_accuracy(estimator, X, Y):
    """
    Calculates the accuracy score from the developmental dataset
    predictions.

    Args:
    estimator: Random forest classifier object
    X: a data frame of normalized values from developmental dataset
    Y: a vector of sample labels from developmental dataset

    Yields:
    float: accuracy score
    """
    labels_pred = dev_predictions(estimator, X)
    return accuracy_score(Y, labels_pred)


def lr_fe_step(estimator, X, Y, n_features_to_keep, features, 
               fold, out_dir, dev_size, SEED, RANK):
    """
    Split training data into developmental dataset and apply 
    logistic regression to developmental dataset, rank features, 
    and conduct feature elimination, single steps.

    Args:
    estimator: Logistic regression classifier object
    X: a data frame of training data
    Y: a vector of sample labels from training data set
    n_features_to_keep: number of features to keep
    features: a vector of feature names
    fold: current fold
    out_dir: output directory. default '.'
    dev_size: developmental size. default '0.20'
    SEED: random state. default 'True'
    RANK: run feature ranking. default 'True'

    Yields:
    dict: a dictionary with number of features, normalized mutual
          information score, accuracy score, auc roc score and 
          selected features
    """
    kwargs = {'random_state': 13, 
              'test_size': dev_size} if SEED else {'test_size': dev_size}
    X1, X2, Y1, Y2 = train_test_split(X, Y, **kwargs)
    # print(X.shape[1], n_features_to_keep)
    assert n_features_to_keep <= X1.shape[1]
    estimator.fit(X1, Y1)
    rank = np.argsort(estimator.feature_importances_)
    selected = rank[-n_features_to_keep:]
    features_rank_fnc(features, rank, n_features_to_keep, fold, 
                      out_dir, RANK)
    return {'n_features': X1.shape[1],
            'nmi_score': dev_score_nmi(estimator, X2, Y2),
            'accuracy_score': dev_score_accuracy(estimator, X2, Y2),
            'roc_auc_score': dev_score_roc(estimator, X2, Y2),
            'selected': selected}


def lr_fe(estimator, X, Y, n_features_iter, features, fold, out_dir,
          dev_size, SEED, RANK):
    """
    Iterate over features to by eliminated by step.

    Args:
    estimator: Logistic regression classifier object
    X: a data frame of training data
    Y: a vector of sample labels from training data set
    n_features_iter: iterator for number of features to keep loop
    features: a vector of feature names
    fold: current fold
    out_dir: output directory. default '.'
    dev_size: developmental size. default '0.20'
    SEED: random state. default 'True'
    RANK: run feature ranking. default 'True'

    Yields:
    list: a list with number of features, normalized mutual
          information score, accuracy score, auc roc curve and 
          array of the indices for features to keep
    """
    indices = np.array(range(X.shape[1]))
    for nf in chain(n_features_iter, [1]):
        p = lr_fe_step(estimator, X, Y, nf, features, fold, 
                       out_dir, dev_size, SEED, RANK)
        yield p['n_features'], p['nmi_score'], p['accuracy_score'], p['roc_auc_score'], indices
        indices = indices[p['selected']]
        X = X[:, p['selected']]
        features = features[p['selected']]

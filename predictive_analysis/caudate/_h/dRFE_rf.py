"""
This script uses dRFEtools to classify female versus male
in BrainSeq Phase 2 DLPFC, Hippocampus, and BrainSeq Phase 3
Caudate Nucleus. This is a RiboZero total RNA datasets.
"""
import numpy as np
import polars as pl
import pandas as pd
from pyhere import here
import re, errno, os, argparse
import dRFEtools, session_info
import matplotlib.pyplot as plt
from matplotlib import rcParams
from rpy2.robjects import r, pandas2ri, globalenv
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, accuracy_score
from sklearn.metrics import normalized_mutual_info_score as nmi

def mkdir_p(directory):
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def R_function():
    pandas2ri.activate()
    r('''
    library(dplyr)
    ml_residuals <- function(train_indices, test_indices)
    {
                    # Subset for training data
    expression_train = expr[, train_indices]
                    # Null model
    null_model_train = v$design[train_indices, ] %>% as.data.frame %>%
        select(-c(starts_with("Male"))) %>% as.matrix
                    # Fit model
    fit_train = limma::lmFit(expression_train, design=null_model_train)
                    # Calculate residuals from training data
    residuals_train = expression_train -(fit_train$coefficients %*% t(null_model_train))
    residuals_train_sd = apply(residuals_train, 1, sd)
    residuals_train_mean = apply(residuals_train, 1, mean)
                    # Normalize residuals
    residuals_train_norm = (residuals_train - residuals_train_mean) / residuals_train_sd
                    # Subset for test data
    expression_test = expr[, test_indices]
                    # Null model
    null_model_test = v$design[test_indices, ] %>% as.data.frame %>%
        select(-c(starts_with("Male"))) %>% as.matrix
                    # Apply training to test data and normalize
    residuals_test = expression_test - (fit_train$coefficients %*% t(null_model_test))
    residuals_test_norm = (residuals_test - residuals_train_mean) / residuals_train_sd
    return (list(residuals_train_norm = residuals_train_norm,
                 residuals_test_norm = residuals_test_norm))
    }
    ''')


def residualize(train_index, test_index):
    pandas2ri.activate()
    R_function()
    globalenv['train_index'] = train_index+1 # Adjust for python
    globalenv['test_index'] = test_index+1 # Adjust for python
    r('''
    res_mx      <- ml_residuals(train_index, test_index)
    X_test      <- t(res_mx$residuals_test_norm)
    X_train     <- t(res_mx$residuals_train_norm)
    ''')
    x_train = pl.DataFrame(data=r['X_train'], schema=r['gnames'])
    x_test = pl.DataFrame(data=r['X_test'], schema=r['gnames'])
    return x_train, x_test


def optimize_rf(X, Y, features, estimator, cv, outdir, feature):
    fold_num = 0
    for train_index, test_index in cv.split(X, Y):
        Y_train, Y_test = Y[train_index], Y[test_index]
        fold_num += 1
    fold_num -= 1
    X_train, X_test = residualize(train_index, test_index)
    X_train = X_train.select(pl.col(features))
    X_test = X_test.select(pl.col(features))
    ## Edit to remove "|"
    if feature in ["genes", "transcripts"]:
        new_features = [re.sub("\\|", "_", x) for x in features]
    else:
        new_features = [re.sub("\\-$", "minus",
                            re.sub("\\+$", "plus",
                                   re.sub("\:", "_", x))) for x in features]
    d, pfirst = dRFEtools.rf_rfe(estimator, X_train.to_numpy(),
                                 Y_train.ravel(), np.array(new_features),
                                 fold_num, outdir, elimination_rate=0.1,
                                 RANK=False)
    print("Optimize fraction")
    for frac in [0.2, 0.25, 0.3, 0.35]:
        plt.clf()
        dRFEtools.optimize_lowess_plot(d, fold_num, outdir, frac=frac,
                                       step_size=0.05,classify=True,
                                       save_plot=True)
    ## Default 0.35 looks good for frac
    print("Optimize step size")
    for step in [0.01, 0.02, 0.03, 0.04]:
        plt.clf()
        dRFEtools.optimize_lowess_plot(d, fold_num, outdir, frac=0.30,
                                       step_size=step,classify=True,
                                       save_plot=True)


def annotation(tissue, feature):
    feat_lt = {"genes": "gene", "transcripts": "tx",
               "exons": "exon", "junctions": "jxn"}
    annot_file = here("input/counts/text_files_counts/_m",
                      f"{tissue}/{feat_lt[feature]}_annotation.txt")
    return pl.read_csv(annot_file, sep='\t')\
             .filter(pl.col("seqnames").str.contains(r"chr\d+"))


def extract_feature_annotation(pred_feat, path, fold, tissue, feature):
    # Get important features
    dft = pd.DataFrame.from_records(pred_feat,
                                    columns=['feature_importance',
                                             'Geneid'])
    dft['Fold'] = fold
    # Get gene annotation
    annot = annotation(annot_file)
    annot = annot.rename(columns={'seqname': 'chrom', 'names': 'gene_id'})
    annot['ensemblID'] = annot.gene_id.str.replace("\\..*", "", regex=True)
    pred_df = dft.merge(annot, how='left', left_on='Geneid', right_on='gene_id')
    pred_df.to_csv('%s/important_features.txt' % path,
                   sep='\t', mode='a', index=False,
                   header=True if fold == 0 else False)


def rf_run(X_train, X_test, Y_train, Y_test, fold_num, outdir,
           estimator, frac, step_size, feature, tissue):
    # Apply random forest
    features = X_train.columns
    d, pfirst = dRFEtools.rf_rfe(estimator, X_train.to_numpy(),
                                 Y_train.to_numpy().ravel(),
                                 features, fold_num, outdir,
                                 elimination_rate=0.1)
    df_elim = pd.DataFrame([{'fold':fold_num, 'n features':k,
                             'normalized mutual information':d[k][1],
                             'accuracy':d[k][2],
                             'ROC AUC':d[k][3]} for k in d.keys()])
    dRFEtools.plot_nmi(d, fold_num, outdir)
    dRFEtools.plot_acc(d, fold_num, outdir)
    dRFEtools.plot_roc(d, fold_num, outdir)
    n_features_max = max(d, key=lambda x: d[x][1])
    try:
        # Max features based on lowess curve
        n_features,_ = dRFEtools.extract_max_lowess(d, frac=frac)
        n_redundant,_ = dRFEtools.extract_redundant_lowess(d, frac=frac,
                                                           step_size=step_size)
        dRFEtools.plot_with_lowess_vline(d, fold_num, outdir, frac=frac,
                                         step_size=step_size, classify=True)
    except ValueError:
        n_features = n_features_max
        n_redundant = n_features
    # Fit model
    estimator.fit(X_train, Y_train)
    all_fts = estimator.predict(X_test)
    estimator.fit(X_train.values[:, d[n_redundant][4]], Y_train)
    labels_pred_redundant = estimator.predict(X_test.values[:, d[n_redundant][4]])
    estimator.fit(X_train.values[:, d[n_features][4]], Y_train)
    labels_pred = estimator.predict(X_test.values[:, d[n_features][4]])
    # Output test predictions
    pd.DataFrame({'fold': fold_num,
                  'real': Y_test,
                  'predict_all': all_fts,
                  'predict_max': labels_pred,
                  'predict_redundant': labels_pred_redundant},
                 index=X_test.index)\
      .to_csv('%s/test_predictions.txt' % outdir,
              sep='\t', mode='a', index=True,
              header=True if fold_num == 0 else False)
    # Annotate features
    pred_features = sorted(list(zip(estimator.feature_importances_,
                                    X_train.columns[d[n_features][4]])),
                           reverse=True)
    extract_feature_annotation(pred_features, outdir, fold_num, annot_file)
    # Save output data
    output = dict()
    output['n_max'] = n_features_max
    output['n_features'] = n_features
    output['n_redundant'] = n_redundant
    output['n_features_all_features'] = pfirst[0]
    output['train_oob_score_nmi_all_features'] = pfirst[1]
    output['train_oob_score_acc_all_features'] = pfirst[2]
    output['train_oob_score_roc_all_features'] = pfirst[3]
    output['train_oob_score_nmi'] = dRFEtools.oob_score_nmi(estimator, Y_train)
    output['train_oob_score_acc'] = dRFEtools.oob_score_accuracy(estimator,
                                                                 Y_train)
    output['train_oob_score_roc'] = dRFEtools.oob_score_roc(estimator, Y_train)
    output['test_score_acc'] = accuracy_score(Y_test, labels_pred)
    output['test_score_nmi'] = nmi(Y_test, labels_pred,
                                   average_method="arithmetic")
    output['test_score_roc'] = roc_auc_score(Y_test, labels_pred)
    output['test_score_acc_redundant'] = accuracy_score(Y_test,
                                                        labels_pred_redundant)
    output['test_score_roc_redundant'] = roc_auc_score(Y_test,
                                                       labels_pred_redundant)
    output['test_score_nmi_redundant'] = nmi(Y_test, labels_pred_redundant,
                                             average_method="arithmetic")
    return output, df_elim


def dev_clean_names(feature, tissue):
    if feature in ["genes", "transcripts"]:
        new_annot = [re.sub("\\|", "_", x) for x in annot]
        new_features = [re.sub("\\|", "_", x) for x in r['gnames']]
        expr_df = pl.DataFrame(data=r['expr'], schema=r['snames'])\
                    .transpose(include_header=True, header_name="RNum",
                               column_names=np.array(new_features))
    else:
        new_annot = [re.sub("\\-$", "minus",
                            re.sub("\\+$", "plus",
                                   re.sub("\:", "_", x))) for x in annot]
        new_features = [re.sub("\\-$", "minus",
                            re.sub("\\+$", "plus",
                                   re.sub("\:", "_", x))) for x in r['gnames']]
        expr_df = pl.DataFrame(data=r['expr'], schema=r['snames'])\
                    .transpose(include_header=True, header_name="RNum",
                               column_names=np.array(new_features))


def load_data(feature, tissue):
    fn = here("differential_expression", tissue, "_m",
              feature, "voomSVA.RData")
    pandas2ri.activate()
    globalenv['fn'] = str(fn)
    r('''
    load(fn); model <- v$design; expr  <- v$E
    gnames <- rownames(v$E); snames <- colnames(v$E)
    cnames <- gsub("Male", "Sex", colnames(model))
    ''')
    annot = list(annotation(tissue, feature)\
                 .select(pl.col("name")).to_numpy().ravel())
    expr_df = pl.DataFrame(data=r['expr'], schema=r['snames'])\
                .transpose(include_header=True, header_name="RNum",
                           column_names=np.array(r['gnames']))
    gnames = list(set(annot) & set(r['gnames']))
    model_df = pl.concat(
        [
            pl.DataFrame({"RNum": np.array(r['snames'])}),
            pl.DataFrame(data=r['model'], schema=r['cnames']),
        ],
        how="horizontal",
    ).select(["RNum", "Sex"]).with_columns(pl.col("Sex").cast(pl.Int64))
    return expr_df.select(["RNum", pl.col(gnames)]), model_df


def check_matching_columns(x, y):
    new_x = x.select(pl.col("RNum")); new_y = y.select(pl.col("RNum"))
    if (new_x != new_y).sum().to_numpy() != 0:
        ## reorganize the dataframes to match!
        x = x.select([pl.all().sort_by("RNum")])
        y = y.select([pl.all().sort_by("RNum")])
    return x, y


def main_loop(feature, tissue):
    cla = dRFEtools.RandomForestClassifier(n_estimators=100,
                                           oob_score=True,
                                           n_jobs=2)
    X, Y = load_data(feature, tissue)
    X, Y = check_matching_columns(X, Y)
    X = X.select([pl.all().exclude("RNum")])
    Y = Y.select(pl.col("Sex"))
    features = X.columns
    skf = StratifiedKFold(n_splits=10, shuffle=True, random_state=20230222)
    skf.get_n_splits(X.to_numpy(), Y.to_numpy())


    ### working here!!!!, problem with optimization function
    optimize_rf(X.to_numpy(), Y.to_numpy(), features, cla,
                skf, outdir, feature)
    frac = 0.30; step_size = 0.04; fold = 0
    fields = ['n_features_all_features', 'train_oob_score_acc_all_features',
              'train_oob_score_nmi_all_features', 'train_oob_score_roc_all_features',
              'n_max', 'n_features', 'train_oob_score_nmi', 'train_oob_score_acc',
              'train_oob_score_roc', 'test_score_nmi', 'test_score_acc',
              'test_score_roc', 'n_redundant', 'test_score_nmi_redundant',
              'test_score_acc_redundant', 'test_score_roc_redundant']
    df_dict = pd.DataFrame()
    with open("%s/dRFEtools_10folds.txt" % (outdir), "w") as f:
        print("\t".join(["fold"] + fields), file=f, flush=True)
        for train_index, test_index in skf.split(X, Y):
            Y_train, Y_test = Y[train_index], Y[test_index]
            X_train, X_test = residualize(train_index, test_index)
            X_train = X_train.select(pl.col(features))
            X_test = X_test.select(pl.col(features))
            o, df_elim = rf_run(X_train, X_test, Y_train, Y_test, fold, outdir,
                                cla, frac, step_size, feature, tissue)
            df_dict = pd.concat([df_dict, df_elim], axis=0)
            print("\t".join([str(fold)] + [str(o[x]) for x in fields]),
                  flush=True, file=f)
            fold += 1
        df_dict.to_csv("%s/feature_elimination_allFolds_metrics.txt" % outdir,
                       sep="\t", index=False)


def main():
    # Parser inputs
    parser = argparse.ArgumentParser()
    parser.add_argument('--feature', type=str, default="genes")
    parser.add_argument('--tissue', type=str, default="caudate")
    args=parser.parse_args()
    # Set environment
    os.environ['NUMEXPR_MAX_THREADS'] = "5"
    rcParams.update({'figure.max_open_warning': 0})
    # Run analysis
    feature = args.feature; tissue = args.tissue
    outdir = f"{feature}"; mkdir_p(outdir)
    main_loop(feature, tissue)


if __name__ == '__main__':
    main()

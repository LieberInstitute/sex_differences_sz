## This script calculates model statistics across N folds.

suppressPackageStartupMessages({    
    library(ggpubr)
    library(dplyr)
})

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".png", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

extract_rank <- function(tissue, feature){
    new_feature <- paste0(tolower(feature), "s")
    ml_file     <- paste('../../..',tolower(tissue), '_m',
                         new_feature,'rank_features.txt', sep="/")
    ml_df       <- data.table::fread(ml_file) |>
        rename('feature_id'='V1', 'fold'='V2', 'Rank'='V3') |>
        tidyr::pivot_wider(names_from = fold, values_from = Rank) |>
        rowwise() |> mutate(median_all = median(c_across(`0`:`9`))) |>
        arrange(median_all) |> as.data.frame() |>
        mutate(rank_id = rank(median_all)) |>
        select(feature_id, rank_id)
    ##fwrite(ml_df, file = fn, sep='\t')
    return(ml_df)
}

summary_max_data <- function(tissue, feature){
    new_feature <- paste0(tolower(feature), "s")
    ml_file     <- paste('../../..',tolower(tissue),'_m', new_feature,
                         'dRFEtools_10folds.txt', sep='/')
    ml_df       <- data.table::fread(ml_file) |>
        mutate_at("fold", as.character) |>
        select(fold, n_features, train_oob_score_acc, train_oob_score_nmi,
               train_oob_score_roc, test_score_acc, test_score_nmi,
               test_score_roc) |>
        tidyr::pivot_longer(-fold, names_to='test', values_to='score') |>
        group_by(test) |>
        summarise(mean_score = mean(score), median_score = median(score),
                  std_score = sd(score), .groups = "keep") |>
        mutate(region=tissue, feature_type=feature)
    return(ml_df)
}


summary_redundant_data <- function(tissue, feature){
    new_feature <- paste0(tolower(feature), "s")
    ml_file     <- paste('../../..',tolower(tissue),'_m', new_feature,
                         'dRFEtools_10folds.txt', sep='/')
    ml_df       <- data.table::fread(ml_file) |>
        mutate_at("fold", as.character) |>
        select(fold, ends_with("_redundant")) |>
        pivot_longer(-fold, names_to='test', values_to='score') |>
        group_by(test) |>
        summarise(mean_score = mean(score), median_score = median(score),
                  std_score = sd(score), .groups = "keep") |>
        mutate(region=tissue, feature_type=feature)
    return(ml_df)
}

extract_n_save <- function(){
    maxlist1 = list(); pred_featlist1 = list()
    maxlist2 = list(); pred_featlist2 = list()
    for(feature in c("Gene", "Transcript", "Exon", "Junction")){
        for(tissue in c("Caudate", "DLPFC", "Hippocampus")){
            nf = filter(summary_max_data(tissue, feature),
                        test == "n_features")$median_score
            dat0 = extract_rank(tissue, feature) |> filter(rank_id < nf) |>
                mutate(region=tissue, feature_type=feature)
            dat1 <- summary_max_data(tissue, feature)
            pred_featlist1[[tissue]] <- dat0
            maxlist1[[tissue]] <- dat1
        }
        pred_featlist2[[feature]] <- bind_rows(pred_featlist1)
        maxlist2[[feature]] <- bind_rows(maxlist1)
    }
    bind_rows(pred_featlist2) |>
        data.table::fwrite("dRFE_predictive_features.txt.gz", sep='\t')
    bind_rows(maxlist2) |>
        data.table::fwrite("dRFE_minimal_subset_summary.txt", sep='\t')
}

extract_ml_data <- function(tissue, feature){
    new_feature <- paste0(tolower(feature), "s")
    ml_file     <- paste('../../..',tolower(tissue),'_m', new_feature,
                         'dRFEtools_10folds.txt', sep='/')
    ml_df1      <- data.table::fread(ml_file) |>
        mutate_at("fold", as.character) |>
        select(fold, train_oob_score_acc) |>
        tidyr::pivot_longer(-fold, names_to='metric', values_to='score') |>
        mutate(region=tissue, feature_type=feature, dataset_type="train")
    ml_df2      <- data.table::fread(ml_file) |>
        mutate_at("fold", as.character) |>
        select(fold, test_score_acc) |>
        tidyr::pivot_longer(-fold, names_to='metric', values_to='score') |>
        mutate(region=tissue, feature_type=feature, dataset_type="test")
    ml_df       <- bind_rows(ml_df1, ml_df2)
    return(ml_df)
}

plot_n_save <- function(){
    datalist1 = list(); datalist2 = list()
    for(feature in c("Gene", "Transcript", "Exon", "Junction")){
        for(tissue in c("Caudate", "DLPFC", "Hippocampus")){
            datalist1[[tissue]] <- extract_ml_data(tissue, feature)
        }
        datalist2[[feature]] <- bind_rows(datalist1)
    }
    dt <- bind_rows(datalist2) |> mutate_if(is.character, as.factor) |>
        mutate(feature_type=factor(feature_type,
                                   levels=c("Gene", "Transcript",
                                            "Exon", "Junction")),
               dataset_type=factor(dataset_type, levels=c("train", "test")))
    bxp <- ggboxplot(dt, x="region", y="score", color="dataset_type",
                     facet.by="feature_type", palette="npg",
                     ylab="Model Accuracy", xlab="", add="jitter",
                     panel.labs.font=list(face='bold'),
                     legend="bottom", add.params=list(alpha=0.5),
                     ylim=c(0.5, 1.25), outlier.shape=NA,
                     ggtheme=theme_pubr(base_size=20, border=TRUE)) +
        geom_hline(yintercept=0.90, linetype=2) + rotate_x_text(45)
    save_plot(bxp, "model_accuracy_plot", 8, 10)
}

#### Main

extract_n_save()
plot_n_save()

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

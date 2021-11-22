## This script calculates model statistics across N folds.

library(ggpubr)
library(tidyverse)

map_feature <- function(feature){
    return(list("genes"="Gene", "transcripts"="Transcript",
                "exons"="Exon", "junctions"="Junction")[[feature]])
}


map_tissue <- function(tissue){
    return(list("caudate"="Caudate", "dlpfc"="DLPFC",
                "hippocampus"="Hippocampus")[[tissue]])
}


save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".png", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}


extract_rank <- function(tissue, feature){
    ml_file = paste0('../../../',tissue,'/_m/',feature,'/rank_features.txt')
    ##fn = paste0("median_rank_",feature,"_list_", tissue, ".txt")
    ml_df = data.table::fread(ml_file) %>%
        rename('Geneid'='V1', 'Fold'='V2', 'Rank'='V3') %>%
        pivot_wider(names_from = Fold, values_from = Rank) %>%
        mutate(median_all = apply(., 1, median)) %>%
        arrange(median_all) %>% mutate(Rank = rank(median_all)) %>%
        select(Geneid, Rank)
    ##fwrite(ml_df, file = fn, sep='\t')
    return(ml_df)
}


summary_max_data <- function(tissue, feature){
    ml_file = paste0('../../../',tissue,'/_m/',feature,'/dRFEtools_10folds.txt')
    ml_df = data.table::fread(ml_file) %>% mutate_at("fold", as.character) %>%
        select(fold, n_features, train_oob_score_acc, train_oob_score_nmi,
               train_oob_score_roc, test_score_acc, test_score_nmi,
               test_score_roc) %>%
        pivot_longer(-fold, names_to='Test', values_to='Score') %>%
        group_by(Test) %>%
        summarise(Mean = mean(Score), Median = median(Score),
                  Std = sd(Score), .groups = "keep") %>%
        mutate(Tissue=map_tissue(tissue), Feature=map_feature(feature))
    return(ml_df)
}


summary_redundant_data <- function(tissue, feature){
    ml_file = paste0('../../../',tissue,'/_m/',feature,'/dRFEtools_10folds.txt')
    ml_df = data.table::fread(ml_file) %>%
        mutate_at("fold", as.character) %>%
        select(fold, ends_with("_redundant")) %>%
        pivot_longer(-fold, names_to='Test', values_to='Score') %>%
        group_by(Test) %>%
        summarise(Mean = mean(Score), Median = median(Score),
                  Std = sd(Score), .groups = "keep") %>%
        mutate(Tissue=map_tissue(tissue), Feature=map_feature(feature))
    return(ml_df)
}

extract_n_save <- function(){
    maxlist1 = list(); pred_featlist1 = list()
    maxlist2 = list(); pred_featlist2 = list()
    for(feature in c("genes", "transcripts", "exons", "junctions")){
        for(tissue in c("caudate", "dlpfc", "hippocampus")){
            nf = filter(summary_max_data(tissue, feature),
                        Test == "n_features")$Median
            dat0 = extract_rank(tissue, feature) %>% filter(Rank < nf) %>%
                mutate(Tissue=map_tissue(tissue), Feature=map_feature(feature))
            dat1 <- summary_max_data(tissue, feature)
            pred_featlist1[[tissue]] <- dat0
            maxlist1[[tissue]] <- dat1
        }
        pred_featlist2[[feature]] <- bind_rows(pred_featlist1)
        maxlist2[[feature]] <- bind_rows(maxlist1)
    }
    bind_rows(pred_featlist2) %>%
        data.table::fwrite("dRFE_predictive_features.txt.gz", sep='\t')
    bind_rows(maxlist2) %>%
        data.table::fwrite("dRFE_minimal_subset_summary.txt", sep='\t')
}

extract_ml_data <- function(tissue, feature){
    ml_file = paste0('../../../',tissue,'/_m/',feature,'/dRFEtools_10folds.txt')
    ml_df1 = data.table::fread(ml_file) %>% mutate_at("fold", as.character) %>%
        select(fold, train_oob_score_acc) %>%
        pivot_longer(-fold, names_to='Metric', values_to='Score') %>%
        mutate(Tissue=map_tissue(tissue), Feature=map_feature(feature), Type="Train")
    ml_df2 = data.table::fread(ml_file) %>% mutate_at("fold", as.character) %>%
        select(fold, test_score_acc) %>%
        pivot_longer(-fold, names_to='Metric', values_to='Score') %>%
        mutate(Tissue=map_tissue(tissue), Feature=map_feature(feature), Type="Test")
    ml_df = bind_rows(ml_df1, ml_df2)
    return(ml_df)
}


plot_n_save <- function(){
    datalist1 = list(); datalist2 = list()
    for(feature in c("genes", "transcripts", "exons", "junctions")){
        for(tissue in c("caudate", "dlpfc", "hippocampus")){
            datalist1[[tissue]] <- extract_ml_data(tissue, feature)
        }
        datalist2[[feature]] <- bind_rows(datalist1)
    }
    dt = bind_rows(datalist2) %>% mutate_if(is.character, as.factor) %>%
        mutate(Feature=factor(Feature, levels=c("Gene", "Transcript", "Exon", "Junction")),
               Type=factor(Type, levels=c("Train", "Test")))
    bxp = ggboxplot(dt, x="Tissue", y="Score", color="Type", facet.by="Feature",
                    palette="npg", ylab="Model Accuracy", xlab="", add="jitter",
                    outlier.shape=NA, panel.labs.font=list(face='bold'),
                    legend="bottom", add.params=list(alpha=0.5), ylim=c(0.5, 1),
                    ggtheme=theme_pubr(base_size=20, border=TRUE)) +
        geom_hline(yintercept=0.85, linetype=2) + rotate_x_text(45)
    save_plot(bxp, "model_accuracy_plot", 8, 10)
}

extract_n_save()
plot_n_save()

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

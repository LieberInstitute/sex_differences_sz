## This script generates median rank across N folds for each feature.

library(tidyverse)

map_tissue <- function(tissue){
    return(list("caudate"="Caudate", "dlpfc"="DLPFC",
                "hippocampus"="Hippocampus")[[tissue]])
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


rank_n_save <- function(){
    ranking_list = list()
    for(feature in c("genes", "transcripts", "exons", "junctions")){
        for(tissue in c("caudate", "dlpfc", "hippocampus")){
            ranking_list[[tissue]] <- extract_rank(tissue, feature) %>%
                mutate(Tissue=map_tissue(tissue))
        }
        fn = paste0("BrainSeq_sex_prediction_median_rank_", feature, ".txt.gz")
        bind_rows(ranking_list) %>% data.table::fwrite(fn, sep='\t')
    }
}

rank_n_save()

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

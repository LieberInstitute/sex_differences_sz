## This script generates median rank across N folds for each feature.

suppressPackageStartupMessages({    
    library(dplyr)
})

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

rank_n_save <- function(){
    ranking_list = list()
     for(feature in c("Gene", "Transcript", "Exon", "Junction")){
         new_feature <- paste0(tolower(feature), "s")
        for(tissue in c("Caudate", "DLPFC", "Hippocampus")){
           ranking_list[[tissue]] <- extract_rank(tissue, feature) |>
               mutate(region=tissue)
        }
        fn = paste0("BrainSeq_sex_prediction_median_rank_",
                    new_feature, ".txt.gz")
        bind_rows(ranking_list) |>
            data.table::fwrite(fn, sep='\t')
    }
}

## Main
rank_n_save()

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

## This script extracts results from mash model
library(mashr)
library(dplyr)

load("mashr_meta_results.RData")
m2$result$lfsr |> as.data.frame() |>
    tibble::rownames_to_column("effect") |>
    data.table::fwrite("lfsr.strong_signals.tsv", sep='\t')

m2$result$PosteriorMean |> as.data.frame() |>
    tibble::rownames_to_column("effect") |>
    data.table::fwrite("posterior_mean.strong_signals.tsv", sep='\t')

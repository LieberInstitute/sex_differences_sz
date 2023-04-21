## This script runs a specific chunk
suppressPackageStartupMessages({
    library(argparse)
    library(mashr)
    library(dplyr)
})

save_results_chunk <- function(chunk_num, outdir){
    ## Load data
    load("model_variables.RData")
    load(paste0(outdir,"/chunk_",chunk_num,"_bhat_shat.RData"))
    ## Run mash on chunk
    data.chunk = mash_set_data(bhat_chunk, shat_chunk, V=Vhat)
    m.chunk = mash(data.chunk, g=get_fitted_g(m), fixg=TRUE)
    ## Save chunk
                                        # lfsr
    fn1 = paste0(outdir,"/lfsr_",chunk_num,"_ancestry.txt.gz")
    data.frame(m.chunk$result$lfsr) |>
        tibble::rownames_to_column("effect") |>
        tidyr::separate(effect, c("gene_id", "variant_id"),
                        remove=FALSE, sep="_", extra="merge") |>
        data.table::fwrite(fn1, sep='\t')
                                        # posterior
    fn2 = paste0(outdir,"/posterior_mean_",chunk_num,"_ancestry.txt.gz")
    data.frame(m.chunk$result$PosteriorMean) |>
        rename_all(list(~stringr::str_replace_all(., '.mean', ''))) |>
        tibble::rownames_to_column("effect") |>
        tidyr::separate(effect, c("gene_id", "variant_id"),
                        remove=FALSE, sep="_", extra="merge") |>
        data.table::fwrite(fn2, sep='\t')
}

## Create parser object
parser <- ArgumentParser()
parser$add_argument("-c", "--chunk_num", type="integer",
                    help="Chunk num to run")
parser$add_argument("-o", "--output", type="character", default="output",
                    help="Output directory for files [default: %default]")
args <- parser$parse_args()

## Run mashr for specific feature
print(paste("Run chunk:", args$chunk_num))
save_results_chunk(args$chunk_num, args$output)

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()


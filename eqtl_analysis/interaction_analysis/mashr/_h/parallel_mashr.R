suppressPackageStartupMessages({
    library(argparse)
    library(doParallel)
    library(foreach)
    library(mashr)
    library(dplyr)
})


save_results_chunk <- function(feature, chunk_size, threads){
    load(paste0(feature, "/chunk.RData"))
    numCores <- threads #detectCores()
    registerDoParallel(numCores)
    fn = paste0(feature,"/lfsr_allpairs_3tissues.txt.gz")
    chunks = split(1:dim(bhat)[1],
                   cut(1:dim(bhat)[1], chunk_size, labels=FALSE))
    dt <- foreach(i=1:chunk_size, .combine=rbind) %dopar% {
        data.chunk = mash_set_data(bhat[chunks[[i]],],
                                   shat[chunks[[i]],], V=Vhat)
        m.chunk = mash(data.chunk, g=get_fitted_g(m), fixg=TRUE)
        data.frame(m.chunk$result$lfsr) %>%
            tibble::rownames_to_column("effect") %>%
            tidyr::separate(effect, c("gene_id", "variant_id"),
                            remove=FALSE, sep="_")
    }
    dt %>% data.table::fwrite(fn, sep='\t')
}

## Create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
parser$add_argument("-r", "--run_chunk", action="store_true", default=FALSE,
                    help="Run this as a chunk [default]")
parser$add_argument("-c", "--chunk_size", type="integer", default=250,
                    help="Chunk size used for parallel run [default: %default]")
parser$add_argument("-t", "--threads", type="integer", default=10,
                    help="Number of threads to run on [default: %default]")
args <- parser$parse_args()

## Run mashr for specific feature
if(args$run_chunk){
    save_results_chunk(args$feature, args$chunk_size, args$threads)
    ## print(paste("Run chunk:", args$chunk))
}

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

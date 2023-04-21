suppressPackageStartupMessages({
    library(argparse)
})

chunk_data <- function(chunk_size, outdir){
    load("bhat_shat.RData")
    chunks = split(1:dim(bhat)[1],
                   cut(1:dim(bhat)[1], chunk_size, labels=FALSE))
    for(chunk_num in 1:chunk_size){
        bhat_chunk <- bhat[chunks[[chunk_num]],]
        shat_chunk <- shat[chunks[[chunk_num]],]
        save(bhat_chunk, shat_chunk,
             file=paste0(outdir, "/chunk_",chunk_num,
                         "_bhat_shat.RData"))
    }
}


## Create parser object
parser <- ArgumentParser()
parser$add_argument("-c", "--chunk_size", type="integer", default=250,
                    help="Chunk size used for parallel run [default: %default]")
parser$add_argument("-o", "--output", type="character", default="output",
                    help="Output directory for files [default: %default]")
args <- parser$parse_args()

## Run mashr for specific feature
chunk_data(args$chunk_size, args$output)

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()


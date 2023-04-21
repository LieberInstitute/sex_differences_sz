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
    fn = paste0(feature,"/lfsr_allpairs_ancestry.txt.gz")
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

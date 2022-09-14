## Generate supplemental material for colocalization results
suppressPackageStartupMessages({
    library(doParallel)
    library(foreach)
    library(dplyr)
})


tissue_map <- function(tissue){
    return(list("caudate"="Caudate", "dlpfc"="DLPFC",
                "hippocampus"="Hippocampus")[[tissue]])
}


get_coloc_signal <- function(tissue){
    cols = c("cluster_name", "N_SNPS", "eQTL_PIP", "GWAS_PIP", "GWAS_PIP_prior",
             "RCP")
    fn = paste0("../../../../prep_eqtl_analysis/", tissue,
                "/genes/prepare_expression/fastqtl_nominal/",
                "dapg_finemap/fastENLOC/_m/pgc3.enloc.sig.out")
    return(data.table::fread(fn, col.names=cols) %>%
           mutate(Tissue=tissue_map(tissue)))
}


get_coloc_snps <- function(tissue){
    cols = c("cluster_name", "variant_id", "eQTL_PIP", "GWAS_PIP",
             "GWAS_PIP_prior", "SNP_RCP")
    fn = paste0("../../../../prep_eqtl_analysis/", tissue,
                "/genes/prepare_expression/fastqtl_nominal/",
                "dapg_finemap/fastENLOC/_m/pgc3.enloc.snp.out")
    return(data.table::fread(fn, col.names=cols) %>%
           mutate(Tissue=tissue_map(tissue)))
}

numCores <- 5
registerDoParallel(numCores)
tissues <- c("caudate", "dlpfc", "hippocampus")

## Signal level
dt <- foreach(i=seq_along(tissues), .combine=rbind) %dopar% {
    get_coloc_signal(tissues[i])
}
print(table(dt$Tissue))
data.table::fwrite(dt, "BrainSeq_colocalization_3tissues_signal.tsv", sep='\t')

## SNP level
dt <- foreach(i=seq_along(tissues), .combine=rbind) %dopar% {
    get_coloc_snps(tissues[i])
}
print(table(dt$Tissue))
data.table::fwrite(dt, "BrainSeq_colocalization_3tissues_snps.tsv", sep='\t')

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

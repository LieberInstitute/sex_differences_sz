## Generate supplemental material for colocalization results
suppressPackageStartupMessages({
    library(doParallel)
    library(foreach)
    library(dplyr)
    library(here)
})


get_coloc_signal <- function(tissue){
    cols = c("cluster_name", "N_SNPS", "eQTL_PIP", "GWAS_PIP", "GWAS_PIP_prior",
             "RCP")
    fn = here("prep_eqtl_analysis", tolower(tissue),
              "genes/interaction_model/dapg_finemap",
              "fastENLOC/_m/pgc3.enloc.sig.out")
    return(data.table::fread(fn, col.names=cols) |>
           mutate(Tissue=tissue))
}


get_coloc_snps <- function(tissue){
    cols = c("cluster_name", "variant_id", "eQTL_PIP", "GWAS_PIP",
             "GWAS_PIP_prior", "SNP_RCP")
    fn = here("prep_eqtl_analysis", tolower(tissue),
              "genes/interaction_model/dapg_finemap",
              "fastENLOC/_m/pgc3.enloc.snp.out")
    return(data.table::fread(fn, col.names=cols) |>
           mutate(Tissue=tissue))
}

numCores <- 5
registerDoParallel(numCores)
tissues <- c("Caudate", "DLPFC", "Hippocampus")

## Signal level
dt1 <- foreach(i=seq_along(tissues), .combine=rbind) %dopar% {
    get_coloc_signal(tissues[i])
}
print(table(dt1$Tissue))

print("RCP > 0.5")
print(filter(dt1, RCP > 0.5) |> group_by(Tissue) |>
      summarize(n=n()) |> as.data.frame())
print(filter(dt1, RCP > 0.5))

print("RCP > 0.4")
print(filter(dt1, RCP > 0.4) |> group_by(Tissue) |>
      summarize(n=n()) |> as.data.frame())
print(filter(dt1, RCP > 0.4))

## SNP level
dt2 <- foreach(i=seq_along(tissues), .combine=rbind) %dopar% {
    get_coloc_snps(tissues[i])
}
print(table(dt2$Tissue))
print(filter(dt2, SNP_RCP > 0.5))

## Save file
data.table::fwrite(dt1, "BrainSeq_colocalization_3tissues_signal.tsv", sep='\t')
data.table::fwrite(dt2, "BrainSeq_colocalization_3tissues_snps.tsv", sep='\t')

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

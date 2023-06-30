#### Functional enrichment analysis with g:Profiler

library(dplyr)
library(gprofiler2)

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png', '.svg')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

load_eqtl <- function(){
    fn <- "../../_m/BrainSeq_sexGenotypes_4features_3regions.txt.gz"
    return(data.table::fread(fn) |>
           filter(feature_type == "Gene") |>
           mutate(ensembl_id = gsub("\\..*", "", gene_id)))
}

extract_gostres <- function(df, label){
    gostres <- gost(query=df$ensembl_id, organism="hsapiens")
    data.table::fwrite(gostres$result,
                       file=paste(label, "functional_enrichment.txt", sep="."),
                       sep="\t")
}

#### Main
                                        # Load eQTL results
eqtl_df <- load_eqtl()
caudate <- filter(eqtl_df, region == "Caudate") |>
    distinct(gene_id, .keep_all=TRUE)
dlpfc   <- filter(eqtl_df, region == "DLPFC") |>
    distinct(gene_id, .keep_all=TRUE)
hippo   <- filter(eqtl_df, region == "Hippocampus") |>
    distinct(gene_id, .keep_all=TRUE)
shared  <- intersect(intersect(caudate$ensembl_id, dlpfc$ensembl_id),
                     hippo$ensembl_id)

                                        # Calculated enrichment
extract_gostres(caudate, "caudate")
extract_gostres(dlpfc, "dlpfc")
extract_gostres(hippo, "hippocampus")

                                        # Shared si-eQTL genes
gostres <- gost(query=shared, organism="hsapiens")
data.table::fwrite(gostres$result,
                   file="shared.functional_enrichment.txt",
                   sep="\t")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

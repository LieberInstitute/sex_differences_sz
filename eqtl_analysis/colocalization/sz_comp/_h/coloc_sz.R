#### Extract colocaliztion results
library(here)
library(dplyr)

load_enloc_results <- function(tissue){
    fn <- here("prep_eqtl_analysis", tolower(tissue),
               "genes/interaction_model/dapg_finemap",
               "fastENLOC/_m/pgc3.enloc.sig.out")
    return(data.table::fread(fn) |> rename("RCP"=V6) |>
           mutate(ensembl_id = gsub("\\:.*", "", V1),
                  region = tissue) |> arrange(desc(RCP)) |>
           select(region, ensembl_id, RCP))
}

load_sex_specific <- function(){
    fn <- here("interaction_sex_sz/by_sex_sz/tissue_comparison/summary_table/_m",
               "differential_expression_schizophrenia_by_sex_4features.txt.gz")
    return(data.table::fread(fn) |> filter(Type == "Gene"))
}

#### MAIN
                                        # Load data
dt <- purrr::map(c("Caudate", "DLPFC", "Hippocampus"), load_enloc_results) |>
    bind_rows()
data.table::fwrite(dt, file="brainseq_siQTL_colocalization.enloc.tsv", sep="\t")

                                        # Merge data
dx <- inner_join(dt, load_sex_specific(), by=c("ensembl_id"="ensemblID"),
                 relationship="many-to-many") |> filter(RCP > 0.4) |>
    select(region, Tissue, Feature, ensembl_id, Symbol, RCP, logFC, `adj.P.Val`, Sex)
data.table::fwrite(dx, file="brainseq_siQTL_colocalization.enloc.degs.tsv", sep="\t")

print(filter(dx, RCP > 0.5) |> select(-region) |>  distinct())

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

## This script estimates pi1 across brain regions.

suppressPackageStartupMessages({
    library(dplyr)
    library(qvalue)
})

load_eqtls <- function(){
    fn <- "../../_m/lfsr.sex_interaction.txt.gz"
    return(data.table::fread(fn))
}

load_nominal_eqtls <- function(tissue, egenes){
    fn <- here::here("prep_eqtl_analysis", tolower(tissue),
                     "genes/interaction_model/_m",
                     "BrainSEQ_TOPMed.interaction.txt.gz")
    return(data.table::fread(fn) |>
           filter(phenotype_id %in% egenes$gene_id,
                  variant_id %in% egenes$variant_id) |>
           mutate(effect=paste(phenotype_id, variant_id, sep="_")) |>
           filter(effect %in% egenes$effect))
}

estimate_pi1 <- function(pvals){
    y_max <- round(max(pvals), 2)
    y_max <- min(c(y_max, 0.9))
    return(1 - pi0est(pvals, lambda=seq(0.05, y_max, 0.05))$pi0)
}

plot_dist <- function(pvals, label, title){
    pdf(paste("pi1.histogram", label, "pdf", sep="."))
    hist(pvals,xlim=c(0, 1),xlab="Nominal P-values",
         main=title, col="gray")
    dev.off()
}

pi1_analysis <- function(tissue1, tissue2){
    pvals <- pull(load_nominal_eqtls(tissue1, config[[tissue2]]), "pval_gi")
    pi1   <- estimate_pi1(pvals)
    plot_dist(pvals, tolower(paste(tissue1, tissue2, sep="_")),
              paste("pi1 =", round(pi1, 3)))
}

#### Main
                                        # Load significant results
eqtl_df <- load_eqtls()

                                        # Filter for brain region
c_egene <- filter(eqtl_df, Caudate < 0.05)
d_egene <- filter(eqtl_df, DLPFC < 0.05)
h_egene <- filter(eqtl_df, Hippocampus < 0.05)

                                        # pi1 analysis
config <- list("Caudate"=c_egene, "DLPFC"=d_egene,
               "Hippocampus"=h_egene)

pi1_analysis("Caudate", "DLPFC"); pi1_analysis("Caudate", "Hippocampus")
pi1_analysis("DLPFC", "Hippocampus"); pi1_analysis("DLPFC", "Caudate")
pi1_analysis("Hippocampus", "Caudate"); pi1_analysis("Hippocampus", "DLPFC");

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

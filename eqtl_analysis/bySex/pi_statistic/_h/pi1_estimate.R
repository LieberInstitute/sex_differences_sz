## This script estimates pi1 across brain regions.

suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(qvalue)
})

load_eqtls <- function(tissue, sex){
    fn <- here("prep_eqtl_analysis/by_sex", tolower(tissue), tolower(sex),
               "cis_model/_m/BrainSEQ_TOPMed.signif_variants.txt.gz")
    return(data.table::fread(fn) |>
           mutate(effect=paste(phenotype_id, variant_id, sep="_")))
}

load_nominal_eqtls <- function(tissue, sex, egenes){
    fn <- here("prep_eqtl_analysis/by_sex", tolower(tissue), tolower(sex),
               "cis_model/_m/BrainSEQ_TOPMed.allpairs.txt.gz")
    return(data.table::fread(fn) |>
           filter(phenotype_id %in% egenes$phenotype_id,
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

pi1_analysis <- function(tissue, sex1, sex2){
    egenes  <- load_eqtls(tissue, sex1)
    pvals   <- pull(load_nominal_eqtls(tissue, sex2, egenes), "pval_nominal")
    pi1     <- estimate_pi1(pvals)
    sex_lab <- paste(sex1, sex2, sep="_")
    plot_dist(pvals, tolower(paste(tissue, sex_lab, sep=".")),
              paste("pi1 =", round(pi1, 3)))
}

pi1_analysis_across <- function(tissue1, tissue2, sex){
    egenes <- load_eqtls(tissue1, sex)
    pvals  <- pull(load_nominal_eqtls(tissue2, sex, egenes), "pval_nominal")
    pi1    <- estimate_pi1(pvals)
    lab    <- paste(tissue1, tissue2, sep="_")
    plot_dist(pvals, tolower(paste(sex, lab, sep=".")),
              paste("pi1 =", round(pi1, 3)))
}

#### Main
                                        # pi1 analysis within region
for(tissue in c("Caudate", "DLPFC", "Hippocampus")){
    pi1_analysis(tissue, "Female", "Male")
    pi1_analysis(tissue, "Male", "Female")
}

                                        # pi analysis between region
for(sex in c("Female", "Male")){
    pi1_analysis_across("Caudate", "DLPFC", sex)
    pi1_analysis_across("Caudate", "Hippocampus", sex)
    pi1_analysis_across("DLPFC", "Caudate", sex)
    pi1_analysis_across("DLPFC", "Hippocampus", sex)
    pi1_analysis_across("Hippocampus", "Caudate", sex)
    pi1_analysis_across("Hippocampus", "DLPFC", sex)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

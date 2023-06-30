#############################################
##
## mash modeling in R
## Author: Kynon J Benjamin
##
#############################################

suppressPackageStartupMessages({
    library(mashr)
    library(dplyr)
    library(ggpubr)
})

save_img <- function(image, fn, w, h){
    for(ext in c(".svg", ".pdf")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

get_bhat <- function(){
    fn = "bhat_interaction_3regions.txt.gz"
    cClasses = c("character", "character", "numeric", "numeric", "numeric")
    return(data.table::fread(fn, header=TRUE, sep="\t", colClasses=cClasses) |>
           mutate(effect=paste(phenotype_id, variant_id, sep="_")) |>
	   distinct(effect, .keep_all=TRUE))
}

get_shat <- function(){
    fn = "shat_interaction_3regions.txt.gz"
    cClasses = c("character", "character", "numeric", "numeric", "numeric")
    return(data.table::fread(fn, header=TRUE, sep='\t', colClasses=cClasses) |>
           mutate(effect=paste(phenotype_id, variant_id, sep="_")) |>
	   distinct(effect, .keep_all=TRUE))
}

plot_mixture_prop <- function(m){
    fn = "barplot_estimated_pi"
    df = get_estimated_pi(m) |> as.data.frame() |>
        tibble::rownames_to_column("Model")
    colnames(df)[2] <- "Estimated pi"
    brp = ggbarplot(df, x="Model", y="Estimated pi", fill="gray",
                    ggtheme=theme_pubr(base_size=20), xlab="",
                    label=TRUE, label.pos="out", lab.nb.digits=2) +
        font("y.title", face="bold") + rotate_x_text(45)
    save_img(brp, fn, 7, 7)
}

run_mashr <- function(SEED, percentage){
    set.seed(SEED)
    ## Load prepared data
    print("Load prepared data!")
    bhat <- get_bhat() |> tibble::column_to_rownames("effect") |>
        select(-phenotype_id, -variant_id) |> as.matrix()
    shat <- get_shat() |> tibble::column_to_rownames("effect") |>
        select(-phenotype_id, -variant_id) |> as.matrix()
    save(bhat, shat, file="bhat_shat.RData")
    ## Prepare for mashr
    print("Prepare random and strong subsets!")
                                        # Random subset
    rand_n <- round(dim(bhat)[1] * percentage)
    random.subset <- sample(1:nrow(bhat), rand_n)
                                        # Strong subset
    m.1by1 = mash_1by1(mash_set_data(bhat, shat))
    strong.subset = get_significant_results(m.1by1, 0.20)
    rm(m.1by1); gc(verbose=TRUE)
    ## Correlation structure
    print("Correlate structures!")
    data.init = mash_set_data(bhat[random.subset,], shat[random.subset,])
    Vhat = estimate_null_correlation_simple(data.init)
    rm(data.init); gc(verbose=TRUE)
                                        # Apply correlation structure
    print("Apply to random and strong subset!")
    data.random = mash_set_data(bhat[random.subset,],shat[random.subset,],V=Vhat)
    data.strong = mash_set_data(bhat[strong.subset,],shat[strong.subset,],V=Vhat)
    rm(bhat, shat); gc(verbose=TRUE)
    ## Calculate data driven covariances
    print("Calculated covariances!")
    U.pca = cov_pca(data.strong, 3)
    U.ed = cov_ed(data.strong, U.pca)
    rm(U.pca); gc(verbose=TRUE)
    ## Fit mash model
    print("Fit model!")
    U.c = cov_canonical(data.random)
    m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel=1)
    rm(U.ed, U.c, data.random); gc(verbose=TRUE)
    ## Compute posterior summaries
    print("Compute posterior summaries!")
    m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
                                        # Save model variables
    save(m, Vhat, file=paste0("model_variables.RData"))
    rm(data.strong, m, Vhat); gc(verbose=TRUE)
    ## Examine pairwise sharing
    print("Examine pairwise sharing!")
    print(get_pairwise_sharing(m2))
    ## Save significant results
    print("Save results!")
    print(get_significant_results(m2) |> length())
    save(m2, file=paste0("mashr_meta_results.RData"))
    ## Estimate mixture proportions
    print("Estimate mixture of proportions!")
    print(get_estimated_pi(m2))
    plot_mixture_prop(m2)
    
}

## Run mashr for specific feature
SEED = 20220422; percentage = 0.05
run_mashr(SEED, percentage)

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

library(argparse)
library(mashr)
library(dplyr)
library(ggpubr)

save_img <- function(image, fn, w, h){
    for(ext in c(".svg", ".pdf", ".png")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

get_pvals <- function(feature){
    fn = paste0("../../_m/", feature, "/pvalue_fastqtl_3tissues.tsv")
    pval <- data.table::fread(fn) %>%
        mutate(effect=paste(gene_id, variant_id, sep="_")) %>%
        group_by(gene_id) %>%
        arrange(Caudate, DLPFC, Hippocampus, .by_group = TRUE) %>% slice(1)
    return(pval)
}

get_bhat <- function(feature){
    fn = paste0("../../_m/", feature, "/bhat_fastqtl_3tissues.tsv")
    return(data.table::fread(fn) %>%
           mutate(effect=paste(gene_id, variant_id, sep="_")) %>%
           distinct(effect, .keep_all=TRUE))
}

get_shat <- function(feature){
    fn = paste0("../../_m/", feature, "/shat_fastqtl_3tissues.tsv")
    return(data.table::fread(fn) %>%
           mutate(effect=paste(gene_id, variant_id, sep="_")) %>%
           distinct(effect, .keep_all=TRUE))
}

get_rand_set <- function(bhat0, percentage, SEED){
    set.seed(SEED)
    rand_n = round(dim(bhat0)[1] * percentage)
    return(sample(bhat0$effect, rand_n))
}

save_results <- function(pval, m2, feature){
    fn = paste0(feature,"/significant_geneSNP_pairs_3tissues.tsv")
    sig_eqtls <- get_significant_results(m2)
    df = get_n_significant_conditions(m2, thresh = 0.05,
                                      conditions = NULL, sig_fn = get_lfsr)
    cc = get_n_significant_conditions(m2, thresh = 0.05,
                                      conditions = 1, sig_fn = get_lfsr)
    dd = get_n_significant_conditions(m2, thresh = 0.05,
                                      conditions = 2, sig_fn = get_lfsr)
    hh = get_n_significant_conditions(m2, thresh = 0.05,
                                      conditions = 3, sig_fn = get_lfsr)
    dt = data.frame(N_Regions_Shared=df, Caudate=cc, DLPFC=dd, Hippocampus=hh) %>%
        tibble::rownames_to_column("effect")
    ## Total number of tissue specific eQTL
    print(paste("There are",sum(df == 1), "genes with specific expression!"))
    pval %>% select(gene_id, variant_id, effect) %>%
        filter(effect %in% names(sig_eqtls)) %>%
        inner_join(dt, by="effect") %>% select(-effect) %>%
        data.table::fwrite(fn, sep='\t')
}

plot_mixture_prop <- function(m2, feature){
    fn = paste0(feature,"/barplot_estimated_pi")
    df = get_estimated_pi(m2) %>% as.data.frame %>%
        tibble::rownames_to_column("Model")
    colnames(df)[2] <- "Estimated pi"
    brp = ggbarplot(df, x="Model", y="Estimated pi", fill="gray",
                    ggtheme=theme_pubr(base_size=20), xlab="",
                    label=TRUE, label.pos="out", lab.nb.digits=2) +
        font("y.title", face="bold") + rotate_x_text(45)
    save_img(brp, fn, 7, 7)
}

run_mashr <- function(feature, percentage, SEED){
    dir.create(feature)
    ## Load prepared data
    bhat0 <- get_bhat(feature)
    shat0 <- get_shat(feature)
    pval <- get_pvals(feature)
    ## Get random and strong set
    strong_set <- pval$effect
    random_set <- get_rand_set(bhat0, percentage, SEED)
    ## Prepared data for mashr
    bhat <- bhat0 %>% tibble::column_to_rownames("effect") %>%
        select(-gene_id, -variant_id) %>% as.matrix
    shat <- shat0 %>% tibble::column_to_rownames("effect") %>%
        select(-gene_id, -variant_id) %>% as.matrix
    ## Calculate correlation matrix
    data.temp = mash_set_data(bhat[random_set,],shat[random_set,])
    Vhat = estimate_null_correlation_simple(data.temp)
    rm(data.temp)
    ## Initial mash
    data.random = mash_set_data(bhat[random_set,],shat[random_set,], V=Vhat)
    data.strong = mash_set_data(bhat[strong_set,],shat[strong_set,], V=Vhat)
    ## Calculate data driven covariances
    U.pca = cov_pca(data.strong,3)
    U.ed = cov_ed(data.strong, U.pca)
    U.c = cov_canonical(data.random)
    ## Fit mash model
    m = mash(data.random, Ulist = c(U.ed,U.c), outputlevel = 1)
    ## Compute posterior summaries
    m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)
    ## Examine pairwise sharing
    print(get_pairwise_sharing(m2))
    ## Save significant results
    print(get_significant_results(m2) %>% length)
    save_results(pval, m2, feature)
    save(m2, file=paste0(feature, "/mashr_meta_results.RData"))
    ## Estimate mixture proportions
    print(get_estimated_pi(m2))
    plot_mixture_prop(m2, feature)
}

## Create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to be analyzed [default: %default]")
parser$add_argument("-s", "--seed", type="integer", default=13131313,
                    help="Seed for randomization [default: %default]")
args <- parser$parse_args()

if(args$feature == "genes"){
    percentage = 0.05
} else {
    percentage = 0.01
}

## Run mashr for specific feature
run_mashr(args$feature, percentage, args$seed)

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

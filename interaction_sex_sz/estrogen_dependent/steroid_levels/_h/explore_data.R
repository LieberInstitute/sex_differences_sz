#### This script explores sex steroid hormone levels
#### in the brain

suppressPackageStartupMessages({
    library(dplyr)
})

get_expression <- function(region){
    fn <- here::here("differential_expression",
                     tolower(region), "_m/genes/voomSVA.RData")
    load(fn)
    return(v)
}
memVOOM <- memoise::memoise(get_expression)

get_covs <- function(region){
    v      <- memVOOM(region)
    sample <- select(as.data.frame(v$targets), "Sex", "Dx") |>
        tibble::rownames_to_column("RNum")
    pheno  <- as.data.frame(v$design) |> select(-Intercept)
    return(pheno)
}

extract_sex_hormones <- function(region){
    v      <- memVOOM(region)
    return(v$genes |> as.data.frame() |>
           filter(Symbol %in% c("AR", "ESR1", "ESR2", "PGR")) |>
           select(gencodeID, Symbol) |>
           tibble::rownames_to_column("feature_id"))
}

get_norm_expr <- function(region){
    v     <- memVOOM(region)
    genes <- extract_sex_hormones(region)
    return(as.data.frame(t(v$E[genes$feature_id, ])))
}

get_res_expr <- function(region, protect_vars){
    ##protect_vars = c("Male", "SCZD", "Age")
    genes <- extract_sex_hormones(region)
    covs  <- get_covs(region) |> select(-any_of(protect_vars))
    norm  <- get_norm_expr(region)
    res   <- norm
    for(ii in 1:length(genes$feature_id)){
        d       <- as.data.frame(cbind(y=norm[,ii],covs))
	model   <- lm(y ~ .,data=d)
	res[,ii] <- resid(model)
    }
    return(res)
}

merge_data <- function(region, protect_vars){
    genes <- extract_sex_hormones(region)
    res   <- get_res_expr(region, protect_vars) |>
        tibble::rownames_to_column("RNum") |>
        tidyr::pivot_longer(-RNum, names_to="feature_id",
                            values_to="expression")
    covs  <- get_covs(region) |>
        select(any_of(protect_vars)) |>
        tibble::rownames_to_column("RNum")
    return(inner_join(res, covs, by="RNum") |>
           inner_join(genes, by="feature_id"))
}

fit_model <- function(region, protect_var){
    model = paste("expression ~", protect_var)
    est_fit <- merge_data(region, protect_var) |>
        group_by(Symbol) |>
        reframe(fitEST = broom::tidy(lm(as.formula(model)))) |>
        tidyr::unnest(fitEST) |>
        filter(term != "(Intercept)") |>
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.fdr  = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf+10**(-100)))
    return(as.data.frame(est_fit))
}

fit_model_interaction <- function(region, var1, var2){
    protect_vars <- c(var1, var2)
    model = paste("expression ~", paste(var1, var2, sep="*"))
    est_fit <- merge_data(region, protect_vars) |>
        group_by(Symbol) |>
        reframe(fitEST = broom::tidy(lm(as.formula(model)))) |>
        tidyr::unnest(fitEST) |>
        filter(term != "(Intercept)") |>
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.fdr  = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf+10**(-100)))
    return(as.data.frame(est_fit))
}

loop_single_vars <- function(region){
    data1 <- list()
    for(var in c("Age", "Male", "SCZD")){
        data1[[var]] <- fit_model(region, var)
    }
    return(bind_rows(data1) |> mutate(tissue = region))
}

loop_interaction <- function(region){
    data2 <- list()
    for(var in c("Male", "SCZD")){
        data2[[var]] <- fit_model_interaction(region, "Age", var) |>
            filter(term == paste("Age", var, sep=":"))
    }
    data2[["Age"]] <- fit_model_interaction(region, "Male", "SCZD") |>
            filter(term == "Male:SCZD")
    return(bind_rows(data2) |> mutate(tissue = region))
}

#### Main
                                        # Linear model
df1 <- purrr::map(c("Caudate", "DLPFC", "Hippocampus"),
                  loop_single_vars) |>
    bind_rows()

                                        # Interaction model
df2 <- purrr::map(c("Caudate", "DLPFC", "Hippocampus"),
                  loop_interaction) |>
    bind_rows()

                                        # Combine models
bind_rows(df1, df2) |>
    data.table::fwrite("linear_model.sex_hormones.tsv",
                       sep='\t')

#### Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

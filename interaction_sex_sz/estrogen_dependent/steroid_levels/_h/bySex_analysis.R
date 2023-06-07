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

get_res_expr <- function(region, protect_vars, FEMALE=TRUE){
    sex   <- ifelse(FEMALE, 0, 1)
    genes <- extract_sex_hormones(region)
    covs  <- get_covs(region) |> select(-any_of(protect_vars))
    norm  <- get_norm_expr(region)
    nsize <- dim(select(filter(covs, Male == sex), -Male))[1]
    res   <- norm[1:nsize, ]
    for(ii in 1:length(genes$feature_id)){
        d       <- as.data.frame(cbind(y=norm[,ii],covs)) |>
            filter(Male == sex) |> select(-Male) ## Select female
	model   <- lm(y ~ .,data=d)
	res[,ii] <- resid(model)
    }
    return(res)
}

merge_data <- function(region, protect_vars, FEMALE=TRUE){
    genes <- extract_sex_hormones(region)
    res   <- get_res_expr(region, protect_vars, FEMALE) |>
        tibble::rownames_to_column("RNum") |>
        tidyr::pivot_longer(-RNum, names_to="feature_id",
                            values_to="expression")
    covs  <- get_covs(region) |>
        select(any_of(protect_vars)) |>
        tibble::rownames_to_column("RNum")
    return(inner_join(res, covs, by="RNum") |>
           inner_join(genes, by="feature_id"))
}

fit_model <- function(region, protect_var, FEMALE=TRUE){
    model = paste("expression ~", protect_var)
    est_fit <- merge_data(region, protect_var, FEMALE) |>
        group_by(Symbol) |>
        reframe(fitEST = broom::tidy(lm(as.formula(model)))) |>
        tidyr::unnest(fitEST) |>
        filter(term != "(Intercept)") |>
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.fdr  = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf+10**(-100)))
    return(as.data.frame(est_fit))
}

fit_model_interaction <- function(region, var1, var2, FEMALE=TRUE){
    protect_vars <- c(var1, var2)
    model = paste("expression ~", paste(var1, var2, sep="*"))
    est_fit <- merge_data(region, protect_vars, FEMALE) |>
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
    data1 <- list(); data2 <- list()
    for(var in c("Age", "SCZD")){
        data1[[var]] <- fit_model(region, var, TRUE)
        data2[[var]] <- fit_model(region, var, FALSE)
    }
    return(bind_rows(mutate(bind_rows(data1), sex="Female"),
                     mutate(bind_rows(data2), sex="Male")) |>
           mutate(tissue = region))
}

loop_interaction <- function(region){
    data1 <- fit_model_interaction(region, "Age", "SCZD", TRUE) |>
        filter(term == "Age:SCZD") |> mutate(sex="Female")
    data2 <- fit_model_interaction(region, "Age", "SCZD", FALSE) |>
        filter(term == "Age:SCZD") |> mutate(sex="Male")
    return( bind_rows(data1, data2) |> mutate(tissue = region) )
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
    data.table::fwrite("linear_model.sex_hormones.bySex.tsv",
                       sep='\t')

#### Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

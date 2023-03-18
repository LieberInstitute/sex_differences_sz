## This script plots correlation between potential confounders.
suppressPackageStartupMessages({
    library(here)
    library(tidyverse)
})

## Functions
save_img <- function(image, fn, w=14, h=7){
    for(ext in c(".svg", ".pdf")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}

pca_norm_data <- function(region){
    ## if(region == "Hippocampus"){
    ##     new_region = "HIPPO"
    ## } else {
    ##     new_region = tolower(region)
    ## }
    new_region = tolower(region)
    ## Load voom normalized data
    load(here("differential_expression", new_region,
              "_m/genes/voomSVA.RData"))
    ## Transpose expression
    norm_df = v$E |> t()
    ## Calculate PCA
    pca_df = prcomp(norm_df, center=TRUE)$x[, 1:20]
    ## Convert to data frame
    norm_dt = pca_df |> as.data.frame() |> rownames_to_column("RNum") |>
        pivot_longer(-RNum, names_to="PC", values_to="PC_values") |>
        mutate_if(is.character, as.factor)
    return(norm_dt)
}
memNORM <- memoise::memoise(pca_norm_data)

pca_res_data <- function(region){
    ## if(region == "Hippocampus"){
    ##     new_region = "HIPPO"
    ## } else {
    ##     new_region = tolower(region)
    ## }
    new_region = tolower(region)
    ## Read in residualized data
    fname = here("differential_expression", new_region,
                 "_m/genes/residualized_expression.tsv")
    res_df = data.table::fread(fname) |>
        column_to_rownames("feature_id") |> t()
    ## Calculate PCA
    pca_df = prcomp(res_df, center=TRUE)$x[, 1:20]
    res_dt = pca_df |> as.data.frame() |> rownames_to_column("RNum") |>
        pivot_longer(-RNum, names_to="PC", values_to="PC_values") |>
        mutate_if(is.character, as.factor)
    return(res_dt)
}
memRES <- memoise::memoise(pca_res_data)

get_pheno <- function(region){
    ## if(region == "Hippocampus"){
    ##     new_region = "HIPPO"
    ## } else {
    ##     new_region = tolower(region)
    ## }
    new_region = tolower(region)
                                        # Load shared variables
    vars <- data.table::fread("../../_m/shared_variables.tsv")
    ## Load voom normalized data
    load(here("differential_expression", new_region, "_m/genes/voomSVA.RData"))
    dfx = v$design |> as.data.frame() |> select(starts_with("sv")) |>
        rownames_to_column("RNum") |>
        pivot_longer(!RNum, names_to="Covariate", values_to="Variable")
    dfx$Covariate <- factor(dfx$Covariate, levels = paste0('sv', 1:20))
    dfy = v$targets |> as.data.frame() |> distinct(RNum, .keep_all = TRUE) |>
        select(all_of(vars$Variables)) |> 
        select(-c("SAMPLE_ID", "RNum", "BrNum", "Region", "Dataset", "IID", 
                  "SOL", "Read_Length", "Protocol")) |>
        mutate(across(where(is.character), as.factor)) |>
        mutate(across(where(is.factor), as.numeric)) |>
        mutate(across(where(is.logical), as.numeric))
    dfy = dfy |> rownames_to_column("RNum") |>
        pivot_longer(!RNum, names_to="Covariate", values_to="Variable")
    return(bind_rows(dfy, dfx) |> mutate_if(is.character, as.factor))
}
memPHENO <- memoise::memoise(get_pheno)

merge_expr <- function(region, fnc){
    df <- inner_join(memPHENO(region), fnc(region), by="RNum", multiple = "all")
    df$PC <- factor(df$PC, levels = paste0('PC', 1:50))
    dfx = df |> filter(!grepl("sv", Covariate)) |>
        mutate(Covariate=droplevels(Covariate))
    dfy = df |> filter(grepl("sv", Covariate)) |>
        mutate(Covariate=droplevels(Covariate))
    dfy$Covariate = factor(dfy$Covariate, levels = paste0("sv", 1:10))
    levels(dfy$Covariate) <- paste0("sv", 1:10)
    return(bind_rows(dfx,dfy) |> mutate_if(is.character, as.factor))
}
memEXPR <- memoise::memoise(merge_expr)

merge_covars <- function(region){
    df0 <- memPHENO(region)
    dfx = df0 |> filter(!grepl("sv", Covariate)) |>
        mutate(Covariate=droplevels(Covariate))
    dfy = df0 |> filter(grepl("sv", Covariate)) |>
        mutate(Covariate=droplevels(Covariate))
    dfy$Covariate = factor(dfy$Covariate, levels = paste0("sv", 1:10))
    levels(dfy$Covariate) <- paste0("sv", 1:10)
    df = bind_rows(dfx,dfy)
    return(inner_join(df, df, by="RNum", multiple = "all") |> 
           mutate_if(is.character, as.factor))
}
memCOVARS <- memoise::memoise(merge_covars)

fit_model <- function(region, norm, identity){
    if(identity){
        est_fit0 <- memCOVARS(region) |>
            group_by(Covariate.x, Covariate.y) |>
            reframe(fitEST = broom::tidy(lm(Variable.x ~ Variable.y)))
    } else {
        if(norm){
            est_fit0 <- memEXPR(region, memNORM) |>
                group_by(Covariate, PC) |>
                reframe(fitEST = broom::tidy(lm(Variable ~ PC_values)))
        } else {
            est_fit0 <- memEXPR(region, memRES) |>
                group_by(Covariate, PC) |>
                reframe(fitEST = broom::tidy(lm(Variable ~ PC_values)))
        }
    }
    ## Calculate p-values
    est_fit <- est_fit0 |> tidyr::unnest(fitEST) |>
        filter(term != "(Intercept)") |>
        mutate(p.bonf = p.adjust(p.value, "bonf"),
               p.bonf.sig = p.bonf < 0.05,
               p.bonf.cat = cut(p.bonf, breaks = c(1,0.05, 0.01, 0.005, 0),
                                labels = c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                                include.lowest = TRUE),
               p.fdr = p.adjust(p.value, "fdr"),
               log.p.bonf = -log10(p.bonf+10**(-100)))
    print(est_fit |> count(p.bonf.cat))
    return(est_fit)
}
memFIT <- memoise::memoise(fit_model)

tile_plot <- function(region, norm=TRUE, identity=TRUE){
    ## Tile plot (heatmap)
    my_breaks <- c(0.05, 0.01, 0.005, 0)
    xlabel = "Covariate"
    if(identity){
        ylabel = "Covariate"; out = "covars"
        tile_plot0 <- memFIT(region, norm, identity) |>
            mutate_if(is.character, as.factor) |> rowwise() |>
            mutate(pair=sort(c(Covariate.x,Covariate.y)) |>
                       paste(collapse=",")) |>
            group_by(pair) |> distinct(pair, .keep_all=TRUE) |>
            ggplot(aes(x = Covariate.x, y = Covariate.y, fill = log.p.bonf,
                       label=ifelse(p.bonf.sig,format(round(log.p.bonf,1),
                                                      nsmall=1), "")))
        h = 15; w = 15; limits = c(0, 100)
    } else {
        tile_plot0 <- memFIT(region, norm, identity) |>
            ggplot(aes(x = Covariate, y = PC, fill = log.p.bonf,
                       label=ifelse(p.bonf.sig,format(round(log.p.bonf,1),
                                                      nsmall=1), "")))
        h = 10; w = 12; limits = c(0, 100)
        if(norm){
            ylabel = "Normalized Expression"; out = "norm"
        } else {
            ylabel = "Residualized Expression"; out = "res"
        }
    }
    tile_plot <- tile_plot0 + geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE, aes(fontface="bold")) +
        viridis::scale_color_viridis(option = "magma") +
        viridis::scale_fill_viridis(name="-log10(p-value Bonf)",
                                    option="magma",
                                    direction=-1, limits=limits) +
        labs(x=xlabel, y=ylabel) + ggpubr::theme_pubr(base_size = 15) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              panel.grid = element_blank(),
              axis.title = element_text(size=18, face="bold"))
    fn = paste0(region, "_tilePlot_",out,"_covariates")
    save_img(tolower(tile_plot), fn, w, h)
    return(tile_plot)
}

#### Correlation with expression PCs ####
for(region in c("Caudate", "DLPFC", "Hippocampus")){
    ## Plotting
    tile_plot(region, FALSE, TRUE)  # Identity
    tile_plot(region, TRUE, FALSE)  # Normalized
    tile_plot(region, FALSE, FALSE) # Residualized
}

#### Reproducibility information ####
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

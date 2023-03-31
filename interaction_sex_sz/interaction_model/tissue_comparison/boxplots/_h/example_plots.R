## This script plots example boxplots

suppressPackageStartupMessages({
    library(ggpubr)
    library(dplyr)
})

ggplot_save <- function(p, fn, w=7, h=7){
    for(ext in c('.png', '.pdf')){
        ggsave(p, filename=paste0(fn, ext), width=w, height=h)
    }
}

get_de <- function(){
    fn <- paste0("../../summary_table/_m/",
                 "differential_expression_interaction_jxn.txt")
    return(data.table::fread(fn) |> arrange(Tissue, adj.P.Val) |>
           select(Tissue, Feature, Symbol, logFC,
                  adj.P.Val, Direction) |>
           group_by(Tissue, Direction) |> slice(1) |>
           as.data.frame())
}

get_pheno <- function(){
    fn <- here::here("input/phenotypes/_m/phenotypes.csv")
    pheno_df <- data.table::fread(fn) |>
        select(Region, RNum, Sex, Dx, Race, Age) |>
        distinct() |> mutate_if(is.character, as.factor)
    levels(pheno_df$Sex) <- c("Female", "Male")
    levels(pheno_df$Region)[3] <- "Hippocampus"
    return(pheno_df)
}

get_resid <- function(tissue){
    fn <- here::here("interaction_sex_sz/interaction_model",
                     tolower(tissue),
                     "_m/junctions/residualized_expression.tsv")
    return(data.table::fread(fn) |>
           filter(feature_id %in% get_de()$Feature) |>
           tibble::column_to_rownames("feature_id") |>
           t() |> as.data.frame() |>
           tibble::rownames_to_column("RNum"))
}

merge_data <- function(){
    dt <- list()
    for(region in c("Caudate", "DLPFC", "Hippocampus")){
      dt[[region]] <- get_resid(region)  
    }
    return(bind_rows(dt) |>
           tidyr::pivot_longer(-RNum, names_to="feature_id",
                               values_to="resid") |>
           inner_join(get_pheno(), by="RNum") |>
           inner_join(get_de(),
                      by=c("feature_id"="Feature",
                           "Region"="Tissue")) |>
           as.data.frame())
}


plot_results <- function(df, tissue){
    y0  <- round(min(df$resid)) - 0.5;
    y1  <- round(max(df$resid)) + 0.25;
    tmp <- get_de() |> filter(Tissue == tissue) |>
        mutate(group1="Control", group2="SCZD", y_pos=y1+0.25) |>
        mutate_if(is.character, as.factor)
    bxp <- df |> filter(Region == tissue) |>
        ggboxplot(x="Dx", y="resid", color="Sex", add="jitter",
                  panel.labs.font=list(face='bold', size=14),
                  facet.by=c("feature_id"), palette="npg",
                  xlab="", ylab="Residualized Expression",
                  ylim=c(y0, y1+1), add.params=list(alpha=0.5),
                  legend="bottom", outlier.shape=NA,
                  ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        geom_signif(data=tmp, aes(xmin=group1, xmax=group2,
                                  annotations=format.pval(adj.P.Val, digits=2),
                                  y_position=y_pos),
                    manual=TRUE) +
        font("xy.title", face="bold", size=16) +
        font("legend.title", face="bold", size=16)
    return(bxp)
}

#### MAIN
df  <- merge_data()
for(region in c("Caudate", "DLPFC", "Hippocampus")){
    bxp <- plot_results(df, region)
    ggplot_save(bxp, paste0(tolower(region), '.jxn_interaction'), 6, 4.2)
}

#### Reproducility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

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

get_sex_de <- function(){
    fn <- here::here("differential_expression/tissue_comparison",
                     "summary_table/_m",
                     "differential_expression_analysis_4features_sex.txt.gz")
    return(data.table::fread(fn) |> filter(Type == "Gene") |>
           group_by(Feature) |> summarise(n = n()) |>
           filter(n == 3) |> as.data.frame())
}
memSEX <- memoise::memoise(get_sex_de)

get_de <- function(){
    fn <- paste0("../../metrics_summary/_m/",
                 "differential_expression_region_interaction_sex_4features.txt")
    return(data.table::fread(fn) |>
           filter(Type == "Gene", Feature %in% memSEX()$Feature) |>
           arrange(Comp, adj.P.Val) |> mutate(Direction = sign(t)) |>
           select(Comp, Feature, Symbol, logFC,
                  adj.P.Val, Direction) |>
           group_by(Comp, Direction) |> slice(1) |>
           as.data.frame() |> select(Feature, Symbol) |>
           distinct())
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

get_resid <- function(){
    fn <- "../../_m/genes/residualized_expression.tsv"
    return(data.table::fread(fn) |>
           filter(feature_id %in% get_de()$Feature) |>
           tibble::column_to_rownames("feature_id") |>
           t() |> as.data.frame() |>
           tibble::rownames_to_column("RNum"))
}

merge_data <- function(){
    return(get_resid() |>
           tidyr::pivot_longer(-RNum, names_to="feature_id",
                               values_to="resid") |>
           inner_join(get_pheno(), by="RNum") |>
           inner_join(get_de(), by=c("feature_id"="Feature")) |>
           as.data.frame())
}


plot_results <- function(){
    df  <- merge_data()
    y0  <- round(min(df$resid)) - 0.5;
    y1  <- round(max(df$resid)) + 0.5;
    bxp <- ggboxplot(df, x="Region", y="resid", color="Sex",
                     add="jitter", facet.by=c("Symbol"), palette="npg",
                     panel.labs.font=list(face='bold', size=14),
                     ylim=c(y0, y1), xlab="", ylab="Residualized Expression",
                     add.params=list(alpha=0.5), #scales="free_y",
                     legend="bottom", outlier.shape=NA,
                     ggtheme=theme_pubr(base_size=15, border=TRUE)) +
        font("xy.title", face="bold", size=16) +
        font("legend.title", face="bold", size=16)
    return(bxp)
}

#### MAIN
bxp     <- plot_results()
outfile <- 'region_sex_interaction'
ggplot_save(bxp, outfile, 11, 5)


#### Reproducility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

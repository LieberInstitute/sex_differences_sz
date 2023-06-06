#### This script compares pi1 statistics with CMC

library(ggpubr)
library(dplyr)

save_plot <- function(p, fn, w, h){
    for(ext in c('.png', '.pdf')){
        ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_cmc <- function(){
    fn <- here::here("input/public_results/cmc_hoffman/_m",
                     "joint_analysis/df_DE_results.inter_sex.disease.RDS")
    cmc_df <- readRDS(fn)
    return(list("MSSM"=cmc_df$`MSSM-Penn-Pitt`$`DxSCZ:Reported_GenderMale`,
                "NIMH"=cmc_df$`NIMH-HBCC`$`DxSCZ:Reported_GenderMale`))
}

get_deg <- function(region){
    fn <- paste0("../../../", tolower(region),
                 "/_m/genes/diffExpr_interaction_full.txt")
    return(data.table::fread(fn) |> filter(P.Value < 0.05) |>
           select(gencodeID, Symbol, logFC, P.Value) |>
           tibble::column_to_rownames("gencodeID") |>
           tidyr::drop_na(Symbol))
}

merge_data <- function(region, cmc){
    cmc_df <- load_cmc()[[cmc]] |> #filter(P.Value < 0.05) |>
        select(Symbol, logFC, P.Value) |> tidyr::drop_na(Symbol)
    return(inner_join(cmc_df, get_deg(region), by="Symbol",
                      suffix=c("_CMC", "_BS"),
                      relationship="many-to-many"))
}

pi1_stat <- function(region, cmc){
    df  <- merge_data(region, cmc)
    pi1 <- 1 - qvalue::qvalue(df$P.Value_CMC)$pi0
    return(pi1)
}

plot_corr <- function(region, cmc){
    df <- merge_data(region, cmc)
    xlab = paste0("Effect Size (", region, ")")
    ylab = paste0("Effect Size (CMC DLPFC)")
    fn = tolower(paste("effectsize_scatter", region, cmc, sep="."))
    pp = ggscatter(df, x="logFC_BS", y="logFC_CMC", add="reg.line",
                   xlab=xlab, ylab=ylab, size=1,
                   add.params=list(color="blue", fill="lightgray"),
                   conf.int=TRUE, cor.coef=TRUE, cor.coef.size=3,
                   cor.coeff.args=list(method="spearman", label.sep="\n"),
                   ggtheme=ggpubr::theme_pubr(base_size=15)) +
         ggpubr::font("xy.title", face="bold", size=16)
    save_plot(pp, fn, 6, 6)
}

#### Main
for(cmc in c("MSSM", "NIMH")){
    for(region in c("Caudate", "DLPFC", "Hippocampus")){
        print(paste(paste(cmc, region, sep=" vs "),
                    pi1_stat(region, cmc), sep=": "))
        plot_corr(region, cmc)
    }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

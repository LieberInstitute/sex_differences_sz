## This script generates plots for GSEA analysis.
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

save_plot <- function(p, fn, w, h){
    for(ext in c('.svg', '.pdf')){
        ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}


get_top_GO <- function(comp){
    new_comp  <- tolower(comp)
    filenames <- list.files(path="./", full.names=TRUE,
                            pattern=paste0(new_comp, "_gsea*"))
    filenames <- filenames[!grepl("CC", filenames)]
    df_list   <- list()
    for(i in seq_along(filenames)){
        df_list[[i]] <- data.table::fread(filenames[i])
    }
    return(bind_rows(df_list) |> arrange(pvalue) |> head(10) |>
           mutate(`Log10`=-log10(qvalue), Comparison=comp))
}

generate_dataframe <- function(){
    df_list <- list()
    comps <- c("CvD", "CvH", "DvH")
    for(jj in seq_along(comps)){
        df_list[[jj]] <- get_top_GO(comps[jj])
    }
    return( bind_rows(df_list) )
}

plot_GO <- function(){
    dt <- generate_dataframe()
    cbPalette <- ggpubr::get_palette(palette = "npg", 4)
    gg1 = ggplot(dt, aes(x=NES, y=Description, color=Comparison,
                         size=Log10)) +
        geom_point(shape=18, alpha=0.8) + xlim(-3, 3) +
        labs(y='', size='-log10 (qvalue)') +
        scale_colour_manual(name="Brain Regions", values=cbPalette,
                            labels=c("Caudate vs DLPFC",
                                     "Caudate vs Hippocampus",
                                     "DLPFC vs Hippocampus")) +
        scale_size_continuous(range = c(2, 10)) +
        theme_bw(base_size=15) +
        theme(axis.title=element_text(face='bold'),
              strip.text=element_text(face='bold'))
    return(gg1)
}

plot_GO_sig <- function(){
    dt <- generate_dataframe() |> filter(Comparison == "CvD")
    cbPalette <- ggpubr::get_palette(palette = "npg", 4)
    gg <- ggplot(dt, aes(x=NES, y=Description, size=Log10)) +
        geom_point(shape=18, alpha=0.8) + xlim(-3, 3) +
        labs(y='', size='-log10 (qvalue)') +
        scale_size_continuous(range = c(2, 10)) +
        theme_bw(base_size=15) +
        theme(axis.title=element_text(face='bold'),
              strip.text=element_text(face='bold'))
    return(gg)
}

#### MAIN
gg1 <- plot_GO()
gg2 <- plot_GO_sig()
save_plot(gg1, "sex_region_GSEA_top10_stacked", 10, 6)
save_plot(gg2, "sex_region_GSEA_top10_cvd", 9, 4.5)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

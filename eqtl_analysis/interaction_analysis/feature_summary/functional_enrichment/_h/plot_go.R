#### Visualization of GO plots
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

save_plot <- function(p, fn, w, h){
    for(ext in c('.svg', '.pdf')){
        ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}


get_top_GO <- function(tissue){
    fname <- paste0(tolower(tissue), ".functional_enrichment.txt")
    return(data.table::fread(fname) |>
           filter(source %in% c("KEGG", "REAC", "GO:BP")) |>
           arrange(p_value) |> head(10) |>
           mutate(`-Log10`=-log10(p_value), Tissue=tissue))
}

generate_dataframe <- function(){
    df_list <- list()
    tissues <- c("Caudate", "DLPFC", "Hippocampus")
    for(jj in seq_along(tissues)){
        df_list[[jj]] <- get_top_GO(tissues[jj])
    }
    return( bind_rows(df_list) )
}

plot_GO <- function(){
    dt <- generate_dataframe()
    cbPalette <- ggpubr::get_palette(palette = "jco", 3)
    gg1 = ggplot(dt, aes(x=`-Log10`, y=term_name, color=Tissue)) +
        geom_point(shape=18, alpha=0.8, size=5) +
        labs(y='', size='-Log10 (p-value)') +
        theme_bw(base_size=15) +
        scale_colour_manual(name="Brain Region", values=cbPalette,
                            labels=c("Caudate","DLPFC","Hippocampus")) +
        geom_vline(xintercept = -log10(0.05), linetype = "dotted") +
        theme(axis.title=element_text(face='bold', size=18),
              strip.text=element_text(face='bold', size=18))
    return(gg1)
}

#### MAIN
gg <- plot_GO()
save_plot(gg, "GOBP_top10_stacked", 7, 6)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

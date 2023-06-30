library(ggplot2)
library(tidyverse)

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

load_DEG_enrichment <- function(){
    return(data.table::fread("clincial_phenotypes_enrichment_analysis_DEGs.tsv"))
}
memDEG <- memoise::memoise(load_DEG_enrichment)

load_TWAS_enrichment <- function(){
    return(data.table::fread("clincial_phenotypes_enrichment_analysis_TWAS.tsv") |>
           mutate_at("OR", ~replace(., is.na(.), 1)))
}
memTWAS <- memoise::memoise(load_TWAS_enrichment)

gen_data <- function(fnc){
    err = 0.0000001
    dt <- fnc() |> mutate_if(is.character, as.factor) |>
        mutate(Tissue=fct_relevel(Tissue, rev), `-log10(FDR)`= -log10(FDR),
               `OR Percentile`= OR / (1+OR), p.fdr.sig=FDR < 0.05,
               `log2(OR)` = log2(OR+err),
               p.fdr.cat=cut(FDR, breaks=c(1,0.05,0.01,0.005,0),
                             labels=c("<= 0.005","<= 0.01","<= 0.05","> 0.05"),
                             include.lowest=TRUE)) |>
        filter(Direction == "All")
    return(dt)
}
memDF <- memoise::memoise(gen_data)

plot_tile <- function(fnc, label, w, h){
    y0 <- min(memDF(fnc)$`log2(OR)`)-0.1
    y1 <- max(memDF(fnc)$`log2(OR)`)+0.1
    tile_plot <- memDF(fnc) |>
        ggplot(aes(x = Dataset, y = Tissue, fill = `log2(OR)`,
                   label = ifelse(p.fdr.sig,
                                  format(round(`-log10(FDR)`,1), nsmall=1), ""))) +
        ylab('si-eQTL (eGenes)') + xlab(label) + geom_tile(color = "grey") +
        ggfittext::geom_fit_text(contrast = TRUE) +
        scale_fill_gradientn(colors=c("blue", "white", "red"),
                             values=scales::rescale(c(y0, 0, y1)),
                             limits=c(y0,y1)) +
        ggpubr::theme_pubr(base_size = 20, border=TRUE) +
        theme(axis.text.x = element_text(angle = 45, hjust=1),
              axis.title  = element_text(face="bold"),
              axis.text.y = element_text(face="bold"),
              strip.text  = element_text(face="bold"))
    save_plot(tile_plot, paste0("tileplot_enrichment_",tolower(label)), w, h)
}

## Run script
plot_tile(memDEG, "DEGs", 10, 6)
plot_tile(memTWAS, "TWAS", 7, 6)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

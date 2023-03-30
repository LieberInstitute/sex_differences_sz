## This script plots the Rsq for N DEGs

library(ggpubr)

save_plot <- function(p, fn, w, h){
    for(ext in c('.png', '.pdf')){
        ggplot2::ggsave(file=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

get_variance <- function(tissue){
    new_tissue <- tolower(tissue)
    fn <- paste0(new_tissue,
                 "/variance_explained_rsq_recursive_autosomes.csv")
    return(data.table::fread(fn) |>
           dplyr::mutate("Region"=tissue))
}

merge_data <- function(){
    dt <- list()
    for(region in c("Caudate", "DLPFC", "Hippocampus")){
        dt[[region]] <- get_variance(region)
    }
    return( dplyr::bind_rows(dt) )
}

#### Main

sca <- merge_data() |>
    ggscatter(x="DEGs", y="Rsq", facet.by="Region", scales="free_x",
              size=2, panel.labs.font=list(face="bold", size=16),
              alpha=0.75, xlab="# of Autosomal DEGs",
              ylab="R2 of PC1\n(Residualized Expression)",
              ggtheme=theme_pubr(base_size=15, border=TRUE),
              ncol=1) +
    font("xy.title", face="bold")

save_plot(sca, "correlation_sex_autosomal_DEGs", 4, 10)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

## This script plots the dRFE curve

suppressPackageStartupMessages({    
    library(dplyr)
    library(ggpubr)
})

read_data <- function(tissue){
    fn <- paste0("../../../", tolower(tissue), "/_m/genes/",
                 "feature_elimination_allFolds_metrics.txt")
    return(data.table::fread(fn) |> group_by(`n features`) |>
           summarise(nmi=mean(`normalized mutual information`)) |>
           as.data.frame() |> rename("n_features"="n features") |>
           mutate(region=tissue))
}

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

#### Main
tissues <- c("Caudate", "DLPFC", "Hippocampus")
dt <- purrr::map(tissues, read_data) |> bind_rows()

sca <- ggscatter(dt, x="n_features", y="nmi", facet.by="region",
                 palette="npg", ylab="Mean NMI", xlab="N Features",
                 legend="bottom", ylim=c(0,1), xscale="log10",
                 panel.labs.font=list(face='bold'), format.scale=TRUE,
                 ggtheme=theme_pubr(base_size=20, border=TRUE)) +
    font("xy.title", face="bold")
save_plot(sca, "drfe_curve", 12, 5)

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

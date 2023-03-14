## This script generates a clustered dotplot for
## enrichment analysis (XCI)

library(dplyr)
library(ggplot2)

save_plot <- function(p, fn, w, h){
    for(ext in c(".pdf", ".png", ".svg")){
        ggsave(filename=paste0(fn,ext), plot=p, width=w, height=h)
    }
}

#### Main
xci <- data.table::fread("xci_enrichment_analysis.txt")
print(xci |> filter(Bonferroni < 0.05))

err <- 0.0000001
df  <- xci |> 
    mutate(`-log10(Bonferroni)` = -log10(Bonferroni), 
           `OR Percentile` = OR / (1+OR), 
           `log10(OR)` = log10(OR+err),
           Direction=as.factor(Direction)) |>
    tidyr::drop_na()
levels(df$Direction) <- c("All", "Female Bias", "Male Bias")

y1 <- max(df$`log10(OR)`, na.rm=TRUE)+0.1
y0 <- min(df$`log10(OR)`, na.rm=TRUE)-0.1

dotplot <- df |> 
    ggplot(aes(x=`XCI_status`, y=Direction, color=`log10(OR)`, 
               size=`-log10(Bonferroni)`)) + 
    geom_point() + ylab('') + xlab('') + facet_grid(~Tissue) +
    scale_color_gradientn(colors=c("blue", "grey", "red"), 
                          values=scales::rescale(c(y0, 0, y1)), 
                          limits=c(y0, y1)) +
    theme_bw() + 
    theme(axis.line  = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust=1), 
          axis.text = element_text(size=14), 
          axis.ticks = element_blank(), 
          legend.position="right", 
          panel.grid = element_blank(), 
          strip.text=element_text(size=14, face="bold"))

save_plot(dotplot, "dotplot_enrichment_xci", 9, 4)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

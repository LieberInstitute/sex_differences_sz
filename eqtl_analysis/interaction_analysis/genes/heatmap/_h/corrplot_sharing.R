## Example sharing using posterior sampling.
library(mashr)
library(corrplot)

plot_heatmap <- function(feature, factor){
    config <- list("0"="signOnly", "0.5"="general", "0.99"="strict")
    fn     <- paste("../../_m/mashr_meta_results.RData", sep="/")
    load(fn)
    pdf(file = paste0("sharing_heatmap_",config[[as.character(factor)]],
                      "_all_",feature,".pdf"), width = 6, height = 6)
    x = get_pairwise_sharing(m2, factor=factor, lfsr_thresh = 1)
    corrplot(x,method='color',is.corr=FALSE,col.lim=c(0,1),type='upper',
             addCoef.col="white",tl.col="black", tl.srt=45,
             title='Pairwise Sharing by Magnitude',
             mar=c(4,0,4,0), col=viridisLite::mako(25, direction=-1))
    dev.off()

    pdf(file = paste0("sharing_heatmap_",config[[as.character(factor)]],
                      "_sig_",feature,".pdf"), width = 6, height = 6)
    x = get_pairwise_sharing(m2, factor=factor, lfsr_thresh = 0.05)
    corrplot(x,method='color',is.corr=FALSE,col.lim=c(0,1),type='upper',
             addCoef.col="white",tl.col="black",tl.srt=45,
             title='Pairwise Sharing by Magnitude',
             mar=c(4,0,4,0), col=viridisLite::mako(25, direction=-1))
    dev.off()
}

#### MAIN
for(feature in c("genes")){
    for(factor in c(0, 0.5, 0.99)){
        plot_heatmap(feature, factor)
    }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

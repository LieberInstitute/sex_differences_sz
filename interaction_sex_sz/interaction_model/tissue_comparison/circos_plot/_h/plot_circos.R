## Generate circlized plot
suppressPackageStartupMessages({
    library(dplyr)
    library(circlize)
    library(ComplexHeatmap)
})

extract_bed <- function(tissue){
    fn  <- paste0("../../summary_table/_m/",
                  "differential_expression_interaction_jxn.txt")
    bed <- data.table::fread(fn) %>%
        filter(Type == "Junction", Tissue == tissue) %>%
        dplyr::select(seqnames,start,end,logFC,`adj.P.Val`)
    bed_up   <- bed %>% filter(logFC > 0, adj.P.Val < 0.05)
    bed_down <- bed %>% filter(logFC < 0, adj.P.Val < 0.05)
    return(list("up"=bed_up, "down"=bed_down))
}

plot_circos <- function(caudate, dlpfc, hippo){
    lgd_points = Legend(at=c("Downregulated", "Upregulated"),
                        type="points", legend_gp=gpar(col=c("blue", "red")),
                        title_position="topleft", title="Interaction",
                        background="#FFFFFF")
    circos.clear() # clear plot if there is any
    circos.par("start.degree" = 0, "cell.padding" = c(0, 0, 0, 0),
               "track.height" = 0.15) # rotate 90 degrees
    # initialize with ideogram
    # use hg38, default is hg19
    circos.initializeWithIdeogram(species="hg38")
    circos.genomicTrack(caudate, bg.border="#000080",
                        bg.col=add_transparency("#000080", transparency=0.6),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.7, col = c("blue", "red")[i], ...)
    })
    circos.genomicTrack(dlpfc, bg.border="#8B0000",
                        bg.col=add_transparency("#8B0000", transparency=0.6),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.7, col = c("blue", "red")[i], ...)
    })
    circos.genomicTrack(hippo, bg.border="#006400",
                        bg.col=add_transparency("#006400", transparency=0.6),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.7, col = c("blue", "red")[i], ...)
    })
    draw(lgd_points, x=unit(6, "mm"), y=unit(6, "mm"), just=c("left", "bottom"))
}

####### MAIN
caudate <- extract_bed("Caudate")
dlpfc   <- extract_bed("DLPFC")
hippo   <- extract_bed("Hippocampus")
                                        # plot
pdf(file = paste0("significant_circos_plot.pdf"),
    width = 6, height = 6)
plot_circos(caudate, dlpfc, hippo)
dev.off()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

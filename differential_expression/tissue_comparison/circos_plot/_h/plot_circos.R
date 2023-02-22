## Generate circlized plot
suppressPackageStartupMessages({
    library(dplyr)
    library(circlize)
    library(ComplexHeatmap)
})

extract_bed <- function(tissue){
    fn  <- paste0("../../summary_table/_m/",
                  "differential_expression_analysis_4features_sex.txt.gz")
    bed <- data.table::fread(fn) %>%
        filter(Type == "Gene", Tissue == tissue) %>%
        dplyr::select(seqnames,start,end,logFC,`adj.P.Val`)
    bed_male   <- bed %>% filter(logFC > 0, adj.P.Val < 0.05)
    bed_female <- bed %>% filter(logFC < 0, adj.P.Val < 0.05)
    return(list("male"=bed_male, "female"=bed_female))
}

plot_circos <- function(caudate, dlpfc, hippo){
    lgd_points = Legend(at=c("Female Bias", "Male Bias"),
                        type="points", legend_gp=gpar(col = c("red", "blue")),
                        title_position="topleft", title="Sex Bias",
                        background="#FFFFFF")
    circos.clear() # clear plot if there is any
    circos.par("start.degree" = 0, "cell.padding" = c(0, 0, 0, 0),
               "track.height" = 0.15) # rotate 90 degrees
    # initialize with ideogram
    # use hg38, default is hg19
    circos.initializeWithIdeogram(species="hg38")
    circos.genomicTrack(caudate, bg.border="#000080",
                        bg.col=add_transparency("#000080", transparency=0.8),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.6, col = c("blue", "red")[i], ...)
    })
    circos.genomicTrack(dlpfc, bg.border="#8B0000",
                        bg.col=add_transparency("#8B0000", transparency=0.8),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.6, col = c("blue", "red")[i], ...)
    })
    circos.genomicTrack(hippo, bg.border="#006400",
                        bg.col=add_transparency("#006400", transparency=0.8),
                        panel.fun = function(region, value, ...) {
                            i = getI(...)
                            circos.genomicPoints(region, value, pch = 16,
                                                 cex = 0.6, col = c("blue", "red")[i], ...)
    })
    draw(lgd_points, x=unit(5, "mm"), y=unit(5, "mm"), just=c("left", "bottom"))
}

####### MAIN
caudate <- extract_bed("Caudate")
dlpfc   <- extract_bed("DLPFC")
hippo   <- extract_bed("Hippocampus")
                                        # plot
pdf(file = paste0("significant_circos_plot.pdf"),
    width = 10, height = 10)
plot_circos_4tissue(caudate, dlpfc, hippo)
dev.off()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

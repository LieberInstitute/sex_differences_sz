## This script generates the upset plot for mash sig results.

library(magrittr)
library(ComplexHeatmap)

get_sig_de <- function(){
    fn  <- paste0("../../_m/",
                  "differential_expression_analysis_4features_sex.txt.gz")
    return(data.table::fread(fn) %>% dplyr::filter(adj.P.Val < 0.05))
}
memDE <- memoise::memoise(get_sig_de)

subset_region <- function(region){
    return(dplyr::filter(memDE(), Tissue == region))
}

prep_heatmap <- function(region){
    df   <- subset_region(region)
    gene <- dplyr::filter(df, Type == "Gene")
    tx   <- dplyr::filter(df, Type == "Transcript")
    exon <- dplyr::filter(df, Type == "Exon")
    jxn  <- dplyr::filter(df, Type == "Junction", gencodeID != "")
    lt = list(Gene       = gene$gencodeID,
              Transcript = tx$gencodeID,
              Exon       = exon$gencodeID,
              Junction   = jxn$gencodeID)
    m = make_comb_mat(lt)
    return(m)
}

plot_upset <- function(region, atop, aright){
    m <- prep_heatmap(region)
    cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                        # Right annotation
    right_annot = upset_right_annotation(
        m, ylim = c(0, aright),
        gp = gpar(fill = "black"),
        annotation_name_side = "top",
        axis_param = list(side = "top"))
                                        # Top annotation
    top_annot = upset_top_annotation(
        m, height=unit(7, "cm"),
        ylim = c(0, atop),
        gp=gpar(fill=cbb_palette[comb_degree(m)]),
        annotation_name_rot = 90)
                                        # Plotting
    outfile <- tolower(paste0(region,'.upsetplot.sex.feature_comparison.pdf'))
    pdf(outfile, width=8, height=4)
    ht = draw(UpSet(m, pt_size=unit(4, "mm"), lwd=3,
                    comb_col=cbb_palette[comb_degree(m)],
                    set_order = c("Gene", "Transcript", "Exon", "Junction"),
                    comb_order = order(-comb_size(m)),
                    row_names_gp = gpar(fontsize = 14, fontface='bold'),
                    right_annotation = right_annot, top_annotation = top_annot))
    od = column_order(ht); cs = comb_size(m)
    decorate_annotation("intersection_size", {
        grid.text(cs[od], x = seq_along(cs),
                  y = unit(cs[od], "native") + unit(6, "pt"),
                  default.units="native", just="bottom", gp=gpar(fontsize=11))
    })
    dev.off()
}

###### MAIN
plot_upset("Caudate", 500, 800)
plot_upset("DLPFC", 400, 600)
plot_upset("Hippocampus", 400, 600)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

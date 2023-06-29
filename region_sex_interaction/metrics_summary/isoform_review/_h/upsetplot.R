## This script generates the upset plot for mash sig results.

library(magrittr)
library(ComplexHeatmap)

get_sig_de <- function(){
    fn  <- paste0("../../_m/",
                  "differential_expression_region_interaction_sex_4features.txt")
    return(data.table::fread(fn) %>% dplyr::filter(adj.P.Val < 0.05))
}
memDE <- memoise::memoise(get_sig_de)

subset_comp <- function(comp){
    return(dplyr::filter(memDE(), Comp == comp))
}

prep_heatmap <- function(comp){
    df   <- subset_comp(comp)
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

plot_upset <- function(comp, atop, aright){
    m <- prep_heatmap(comp)
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
    outfile <- tolower(paste0(comp,'.upsetplot.sex.feature_comparison.pdf'))
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
plot_upset("CvD", 500, 800)
plot_upset("CvH", 100, 150)
plot_upset("DvH", 40, 50)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

## This script generates the upset plot for mash sig results.

library(magrittr)
library(ComplexHeatmap)

get_sig_genes <- function(){
    fn  <- paste0("../../summary_table/_m/",
                  "differential_expression_analysis_4features_sex.txt.gz")
    return(data.table::fread(fn) %>%
           dplyr::filter(Type == "Gene", adj.P.Val < 0.05))
}
memGENES <- memoise::memoise(get_sig_genes)

prep_heatmap <- function(){
    caudate <- memGENES() %>% dplyr::filter(Tissue == "Caudate")
    dlpfc <- memGENES() %>% dplyr::filter(Tissue == "DLPFC")
    hippo <- memGENES() %>% dplyr::filter(Tissue == "Hippocampus")
    lt = list(Caudate = caudate$Feature,
              DLPFC=dlpfc$Feature,
              Hippocampus=hippo$Feature)
    m = make_comb_mat(lt)
    return(m)
}

plot_upset <- function(){
    m <- prep_heatmap()
    cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                        # Right annotation
    right_annot = upset_right_annotation(
        m, ylim = c(0, 750),
        gp = gpar(fill = "black"),
        annotation_name_side = "top",
        axis_param = list(side = "top"))
                                        # Top annotation
    top_annot = upset_top_annotation(
        m, height=unit(7, "cm"),
        ylim = c(0, 750),
        gp=gpar(fill=cbb_palette[comb_degree(m)]),
        annotation_name_rot = 90)
                                        # Plotting
    pdf('BrainSeq_sex_upsetR.pdf', width=8, height=4)
    ht = draw(UpSet(m, pt_size=unit(4, "mm"), lwd=3,
                    comb_col=cbb_palette[comb_degree(m)],
                    set_order = c("Caudate", "DLPFC", "Hippocampus"),
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

plot_upset_transposed <- function(){
    m <- prep_heatmap()
    cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                        # Right annotation
    right_ha = rowAnnotation(
        "Intersection\nsize" = anno_barplot(comb_size(m), border=F,
                                            ylim = c(0, 750),
                                            gp=gpar(fill=cbb_palette[comb_degree(m)]),
                                            width = unit(7, "cm")))
                                        # Top annotation
    top_ha = HeatmapAnnotation(
        "Set size" = anno_barplot(set_size(m), border=F,
                                  ylim = c(0, 750),
                                  gp = gpar(fill = "black"),
                                  height = unit(2, "cm")),
        gap = unit(2, "mm"), annotation_name_side = "left",
        annotation_name_rot = 90)
                                        # Plotting
    pdf('BrainSeq_sex_upsetR_tranposed.pdf', width=5, height=10)
    ht = draw(UpSet(t(m), pt_size=unit(5, "mm"), lwd=3,
                    comb_order = order(-comb_size(m)),
                    comb_col=cbb_palette[comb_degree(m)],
                    set_order = c("Caudate", "DLPFC", "Hippocampus"),
                    row_names_gp = gpar(fontsize = 16, fontface='bold'),
                    right_annotation = right_ha, top_annotation = top_ha))
    od = rev(row_order(ht)); cs = comb_size(m)
    decorate_annotation("Intersection\nsize", {
        grid.text(cs[od], y = seq_along(cs),
                  x = unit(cs[od], "native") + unit(6, "pt"),
                  default.units="native", just="left", gp=gpar(fontsize=11))
    })
    dev.off()
}

###### MAIN
plot_upset()
plot_upset_transposed()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

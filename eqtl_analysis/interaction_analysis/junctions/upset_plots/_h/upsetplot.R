## This script generates the upset plot for mash sig results.

library(ComplexHeatmap)

get_eqtls <- function(){
    fn = "../../_m/lfsr.sex_interaction.txt.gz"
    return(data.table::fread(fn))
}
memEQTLS <- memoise::memoise(get_eqtls)

prep_heatmap <- function(df){
    caudate <- df |> dplyr::filter(Caudate < 0.05)
    dlpfc <- df |> dplyr::filter(DLPFC < 0.05)
    hippo <- df |> dplyr::filter(Hippocampus < 0.05)
    lt = list(Caudate = unique(caudate$gene_id),
              DLPFC=unique(dlpfc$gene_id),
              Hippocampus=unique(hippo$gene_id))
    m = make_comb_mat(lt)
    return(m)
}

plot_upset <- function(df){
    m <- prep_heatmap(df)
    cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                        # Right annotation
    right_annot = upset_right_annotation(
        m, ylim = c(0, 400),
        gp = gpar(fill = "black"),
        annotation_name_side = "top",
        axis_param = list(side = "top"))
                                        # Top annotation
    top_annot = upset_top_annotation(
        m, height=unit(7, "cm"),
        ylim = c(0, 400),
        gp=gpar(fill=cbb_palette[comb_degree(m)]),
        annotation_name_rot = 90)
                                        # Plotting
    pdf('BrainSeq_eQTL_upsetR_mash.pdf', width=4, height=4)
    ht = draw(UpSet(m, pt_size=unit(4, "mm"), lwd=3,
                    comb_col=cbb_palette[comb_degree(m)],
                    set_order = c("Caudate","DLPFC", "Hippocampus"),
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

plot_upset_transposed <- function(df){
    m <- prep_heatmap(df)
    cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                        # Right annotation
    right_ha = rowAnnotation(
        "Intersection\nsize" = anno_barplot(comb_size(m), border=F,
                                            ylim = c(0, 400),
                                            gp=gpar(fill=cbb_palette[comb_degree(m)]),
                                            width = unit(7, "cm")))
                                        # Top annotation
    top_ha = HeatmapAnnotation(
        "Set size" = anno_barplot(set_size(m), border=F,
                                  ylim = c(0, 400),
                                  gp = gpar(fill = "black"),
                                  height = unit(2, "cm")),
        gap = unit(2, "mm"), annotation_name_side = "left",
        annotation_name_rot = 90)
                                        # Plotting
    pdf('BrainSeq_eQTL_upsetR_mash_tranposed.pdf', width=5, height=10)
    ht = draw(UpSet(t(m), pt_size=unit(5, "mm"), lwd=3,
                    comb_order = order(-comb_size(m)),
                    comb_col=cbb_palette[comb_degree(m)],
                    set_order = c("Caudate","DLPFC","Hippocampus"),
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
df = memEQTLS()
plot_upset(df)
plot_upset_transposed(df)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

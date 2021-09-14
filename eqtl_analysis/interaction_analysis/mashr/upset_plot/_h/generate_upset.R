## Generate an Upset Plot using the output from mashr
suppressPackageStartupMessages({
    library(dplyr)
    library(optparse)
    library(ComplexHeatmap)
})
                                        # Create parser object
option_list <- list(
    make_option(c("-f", "--feature"), type="character", default="genes",
                help="Feature to extract covariates [default: %(default)]",
                metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

right_size <- function(feature){
    size = 1000
    if(feature == "genes"){
        r_size = size
    } else if(feature == "transcripts"){
        r_size = size + 400
    } else if(feature == "exons"){
        r_size = size + 1000
    } else {
        r_size = size - 400
    }
    return(r_size)
}

top_size <- function(feature){
    size = 1000
    if(feature == "genes"){
        t_size = size
    } else if(feature == "transcripts"){
        t_size = size + 400
    } else if(feature == "exons"){
        t_size = size + 1000
    } else {
        t_size = size - 400
    }
    return(t_size)
}

generate_upset <- function(feature){
    fn = paste0("../../_m/",feature,"/lfsr_allpairs_3tissues.txt.gz")
    outfile = paste0("brainseq_e",toupper(feature),"_upsetplot.pdf")
    df = data.table::fread(fn)
                                        # Extract eFeatures by brain region
    caudate = filter(df, Caudate < 0.05)$`gene_id` %>% unique
    dlpfc = filter(df, DLPFC < 0.05)$`gene_id` %>% unique
    hippo = filter(df, Hippocampus < 0.05)$`gene_id` %>% unique
                                        # Generate list for ComplexHeatmap
    lt = list(Caudate = caudate, DLPFC = dlpfc, Hippocampus = hippo)
    m = make_comb_mat(lt)
    cbb_palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                     "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
                                        # Plot annotation
    r_size = right_size(feature)
    t_size = top_size(feature)
    right_annot = upset_right_annotation(
        m, ylim = c(0, r_size),
        gp = gpar(fill = "black"),
        annotation_name_side = "top",
        axis_param = list(side = "top"))
    top_annot = upset_top_annotation(
        m, height=unit(8, "cm"),
        ylim = c(0, t_size),
        gp=gpar(fill=cbb_palette[comb_degree(m)]),
        annotation_name_rot = 90)
                                        # Export plot
    pdf(outfile, width=8, height=4)
    ht = draw(UpSet(m, pt_size=unit(4, "mm"), lwd=3,
                    comb_col=cbb_palette[comb_degree(m)],
                    set_order = c("Caudate", "DLPFC", "Hippocampus"),
                    comb_order = order(-comb_size(m)),
                    row_names_gp = gpar(fontsize = 15, fontface='bold'),
                    right_annotation = right_annot,
                    top_annotation = top_annot))
    od = column_order(ht)
    cs = comb_size(m)
    decorate_annotation("intersection_size", {
        grid.text(cs[od], x=seq_along(cs), y=unit(cs[od], "native") +
                                                 unit(6, "pt"),
                  default.units="native", just="bottom", gp=gpar(fontsize=12))
    })
    dev.off()
}

generate_upset(opt$feature)

                                        # Reproducibility
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

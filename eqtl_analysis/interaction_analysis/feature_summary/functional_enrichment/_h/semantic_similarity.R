#### Gene ontology semantic similarity analysis
suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
    library(GOSemSim)
    library(org.Hs.eg.db)
})

get_semData <- function(ont){
    return(godata('org.Hs.eg.db', ont=ont))
}

save_img <- function(image, fn, w=7, h=7){
    for(ext in c(".pdf", ".png")){
        ggsave(file=paste0(fn, ext), plot=image, width=w, height=h)
    }
}   

GO_semantic_similarity <- function(ont){
    hsGO = get_semData(ont)
    tissues = c("Caudate", "DLPFC", "Hippocampus")
    t1 = c(); t2 = c(); ss = c()
    for(tissue1 in c("Caudate", "DLPFC")){
        fn1 = paste0(tolower(tissue1), ".functional_enrichment.txt")
        for(tissue2 in c("DLPFC", "Hippocampus")){
            fn2 = paste0(tolower(tissue2), ".functional_enrichment.txt")
            if(tissue1 != tissue2 & file.exists(fn1) & file.exists(fn2)){
                df1 = data.table::fread(fn1) |> filter(source == paste0("GO:", ont))
                df2 = data.table::fread(fn2) |> filter(source == paste0("GO:", ont))
                sim = mgoSim(df1$term_id, df2$term_id, semData=hsGO, 
                             measure="Wang", combine="BMA")
                t1 = c(t1, tissue1); t2 = c(t2, tissue2);ss = c(ss, sim)
            }
        }
    }
    return(data.frame("Tissue_1"=t1, "Tissue_2"=t2, "Semantic_Similarity"=ss, "Ont"=ont))
}

#### Main
                                        # calculate similarity
datalist <- list()
for(ONT in c("MF", "BP", "CC")){
    datalist[[ONT]] <- GO_semantic_similarity(ONT)
}
dt <- bind_rows(datalist) |> mutate_if(is.character, as.factor)
data.table::fwrite(dt, file="go_semantic_similarity.tsv", sep="\t")

                                        # plot similarity
tile_plot <- dt |> #tidyr::drop_na() |> 
    ggplot(aes(x=Tissue_2, y=Tissue_1, fill=Semantic_Similarity, 
               label=format(round(Semantic_Similarity, 2)))) + 
    geom_tile(color="grey") + ggfittext::geom_fit_text(contrast=TRUE) + 
    viridis::scale_color_viridis(option="magma") + facet_wrap("~Ont") +
    viridis::scale_fill_viridis(name="Semantic Similarity", limits=c(0.5,1),
                                direction=-1, option="magma") +
    labs(x="", y="") + ggpubr::theme_pubr(base_size=15, border=TRUE) +
    scale_y_discrete(limits = rev(levels(dt$Tissue_1))) + 
    theme(axis.text.x=element_text(angle = 45, hjust=1), 
          strip.text=element_text(face="bold"), 
          legend.key.width=unit(2, 'cm'))

save_img(tile_plot, "GO_semantic_similarity", w=8, h=5)

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()



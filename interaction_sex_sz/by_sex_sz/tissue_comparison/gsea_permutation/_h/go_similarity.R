#### This script measures sematic similarity between GO annotations
suppressPackageStartupMessages({
    library(dplyr)
    library(GOSemSim)
    library(corrplot)
})

load_go_data <- function(tissue, sex, ont){
    fname <- paste0(tolower(sex), "/", tolower(tissue),
                    "_gsea_", ont, ".tsv")
    return(data.table::fread(fname) |> pull(ID))
}

get_go_data <- function(ONT){
    return(godata("org.Hs.eg.db", ont=ONT))
}
memGO <- memoise::memoise(get_go_data)

cal_sem_sim <- function(ont, tissue){
    go1  <- load_go_data(tissue, "Female", ont)
    go2  <- load_go_data(tissue, "Male", ont)
    simM <- mgoSim(go1, go2, semData=memGO(ont), measure="Wang")
    return(data.frame("Tissue"=tissue, "ONT"=ont, SemSim=simM))
}

#### Main

dt_list <- list()
for(ont in c("BP", "CC", "MF")){
    df_list <- list(); 
    for(tissue in c("Caudate", "DLPFC", "Hippocampus")){
        df_list[[tissue]] <- cal_sem_sim(ont, tissue)
    }
    dt_list[[ont]] <- bind_rows(df_list)
}

df  <- bind_rows(dt_list)
data.table::fwrite(df, file="semantic_similarity.tsv", sep='\t')


mat <- tidyr::pivot_wider(df, names_from=ONT, values_from=SemSim) |>
    tibble::column_to_rownames("Tissue") |> as.matrix()
print(mat)

pdf(file = "semantic_similarity.heatmap.pdf", width = 6, height = 6)
corrplot(mat, method='color', is.corr=FALSE, col.lim=c(0.20,0.90),
         addCoef.col="white", tl.col="black", tl.srt=0,
         title='GO Semantic Similarity (Female vs Male)',
         mar=c(4,0,4,0), col=viridisLite::magma(25, direction=-1))
dev.off()

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

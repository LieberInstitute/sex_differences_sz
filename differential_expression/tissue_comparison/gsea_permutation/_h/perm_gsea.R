## This script analyzsis gene set enrichment with permutation
suppressPackageStartupMessages({
    library(DOSE)
    library(here)
    library(dplyr)
    library(biomaRt)
    library(org.Hs.eg.db)
    library(clusterProfiler)
})

get_biomart <- function(){
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mapping <- getBM(attributes = c('ensembl_gene_id', 'entrezgene_id'),
                     mart = ensembl)
    return(mapping)
}
memBIOMART <- memoise::memoise(get_biomart)

get_features <- function(feature, tissue){
    fn <- paste0("../../summary_table/_m/",
                 "differential_expression_analysis_4features_sex.txt.gz")
    return(data.table::fread(fn) |>
           filter(Type == feature, Tissue == tissue) |>
           mutate(DE=ifelse(`adj.P.Val` < 0.05, 1, 0)))
}
memFEAT <- memoise::memoise(get_features)

get_annotation <- function(feature, tissue){
    dt <- inner_join(memFEAT(feature, tissue), memBIOMART(),
                     by=c("ensemblID"="ensembl_gene_id"),
                     multiple="all") |>
        distinct(Feature, .keep_all = TRUE)
    return(dt)
}
memDF <- memoise::memoise(get_annotation)

get_genes <- function(feature, tissue){
    return(memDF(feature, tissue) |> dplyr::filter(DE == 1) |>
           dplyr::select(entrezgene_id, logFC) |>
           mutate(abs_beta=abs(logFC)) |> group_by(entrezgene_id) |>
           arrange(desc(abs_beta), .by_group=TRUE) |> slice(1) |>
           tidyr::drop_na() |> pull("entrezgene_id") |> as.character())
}

get_geneList <- function(feature, tissue){
    tmp <- memDF(feature, tissue) |> na.exclude() |>
        dplyr::select(entrezgene_id, logFC) |>
        mutate(abs_beta=abs(logFC)) |> group_by(entrezgene_id) |>
        arrange(desc(abs_beta), .by_group=TRUE) |> slice(1)
    geneList <- tmp |> pull("logFC")
    names(geneList) <- tmp$entrezgene_id
    return(geneList |> sort(decreasing=TRUE))
}

run_gseGO <- function(feature, tissue, ont, label){
                                        # GSEA with GO terms
    ego <- gseGO(geneList     = get_geneList(feature, tissue),
                 OrgDb        = org.Hs.eg.db,
                 ont          = ont,
                 keyType      = "ENTREZID",
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE,
                 seed         = TRUE)
                                        # Save results
    if(dim(ego)[1] != 0){
        fn = paste0(label, "_", ont, ".tsv")
        ego <- setReadable(ego, "org.Hs.eg.db")
        ego |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

run_enichmentDGN <- function(feature, tissue, label){
    dgn <- enrichDGN(gene      = get_genes(feature, tissue),
                     universe  = names(get_geneList(feature, tissue)),
                     minGSSize = 5)
      # Save results
    if(dim(dgn)[1] != 0){
        fn = paste0(label, "_DGN.tsv")
        dgn <- setReadable(dgn, "org.Hs.eg.db")
        dgn |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

run_gseDGN <- function(feature, tissue, ont, label){
    dgn <- gseDGN(geneList     = get_geneList(feature, tissue),
                  minGSSize    = 5,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE,
                  seed         = TRUE)
      # Save results
    if(dim(dgn)[1] != 0){
        fn = paste0(label, "_", ont, ".tsv")
        dgn <- setReadable(dgn, "org.Hs.eg.db")
        dgn |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

##### MAIN
                                        # Set seed
set.seed(13131313)
for(feature in c("Gene", "Transcript")){
    new_feature <- paste0(tolower(feature), "s")
    if(!dir.exists(new_feature)){dir.create(new_feature)}
    for(tissue in c("Caudate", "DLPFC", "Hippocampus")){
        ##print(tissue)
        label1 = paste0(new_feature,"/",tolower(tissue),"_gsea")
        label2 = paste0(new_feature,"/",tolower(tissue),"_enrich")
        for(ONT in c("BP", "MF", "CC")){
            run_gseGO(feature, tissue, ONT, label1)
        }
        run_gseDGN(feature, tissue, "DGN", label1)
        run_enichmentDGN(feature, tissue, label2)
    }
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

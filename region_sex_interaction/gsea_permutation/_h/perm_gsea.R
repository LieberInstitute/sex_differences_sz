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

get_features <- function(comp){
    fn <- paste0("../../metrics_summary/_m/",
                 "differential_expression_region_",
                 "interaction_sex_4features_full.txt.gz")
    return(data.table::fread(fn) |>
           filter(Type == "Gene", Comp == comp) |>
           mutate(DE=ifelse(`adj.P.Val` < 0.05, 1, 0)))
}
memFEAT <- memoise::memoise(get_features)

get_annotation <- function(comp){
    dt <- inner_join(memFEAT(comp), memBIOMART(),
                     by=c("ensemblID"="ensembl_gene_id"),
                     multiple="all") |>
        distinct(Feature, .keep_all = TRUE)
    return(dt)
}
memDF <- memoise::memoise(get_annotation)

get_genes <- function(comp){
    return(memDF(comp) |> dplyr::filter(DE == 1) |>
           dplyr::select(entrezgene_id, logFC) |>
           mutate(abs_beta=abs(logFC)) |> group_by(entrezgene_id) |>
           arrange(desc(abs_beta), .by_group=TRUE) |> slice(1) |>
           tidyr::drop_na() |> pull("entrezgene_id") |> as.character())
}

get_geneList <- function(comp){
    tmp <- memDF(comp) |> na.exclude() |>
        dplyr::select(entrezgene_id, logFC) |>
        mutate(abs_beta=abs(logFC)) |> group_by(entrezgene_id) |>
        arrange(desc(abs_beta), .by_group=TRUE) |> slice(1)
    geneList <- tmp |> pull("logFC")
    names(geneList) <- tmp$entrezgene_id
    return(geneList |> sort(decreasing=TRUE))
}

run_gseGO <- function(comp, ont, label){
                                        # GSEA with GO terms
    ego <- gseGO(geneList     = get_geneList(comp),
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

run_enichmentDGN <- function(comp, label){
    dgn <- enrichDGN(gene      = get_genes(comp),
                     universe  = names(get_geneList(comp)),
                     minGSSize = 5)
      # Save results
    if(dim(dgn)[1] != 0){
        fn = paste0(label, "_DGN.tsv")
        dgn <- setReadable(dgn, "org.Hs.eg.db")
        dgn |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

run_gseDGN <- function(comp, ont, label){
    dgn <- gseDGN(geneList     = get_geneList(comp),
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

for(comp in c("CvD", "CvH", "DvH")){
    ##print(comp)
    label1 = paste0(tolower(comp),"_gsea")
    label2 = paste0(tolower(comp),"_enrich")
    for(ONT in c("BP", "MF", "CC")){
        run_gseGO(comp, ONT, label1)
    }
    run_gseDGN(comp, "DGN", label1)
    run_enichmentDGN(comp, label2)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

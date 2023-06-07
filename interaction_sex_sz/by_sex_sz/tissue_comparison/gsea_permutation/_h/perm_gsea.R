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

get_significant <- function(tissue, sex){
    fn <- paste0("../../summary_table/_m/",
                 "differential_expression_schizophrenia_",
                 "by_sex_4features.sig.txt.gz")
    return(data.table::fread(fn) |>
           filter(Type == "Gene", Tissue == tissue, Sex == sex))
}
memSIG <- memoise::memoise(get_significant)

get_features <- function(tissue, sex){
    fn <- paste0("../../summary_table/_m/",
                 "differential_expression_schizophrenia_",
                 "by_sex_4features.txt.gz")
    return(data.table::fread(fn) |>
           filter(Type == "Gene", Tissue == tissue, Sex == sex) |>
           mutate(DE=ifelse(Feature %in% memSIG(tissue, sex)$Feature, 1, 0)))
}
memFEAT <- memoise::memoise(get_features)

get_annotation <- function(tissue, sex){
    dt <- inner_join(memFEAT(tissue, sex), memBIOMART(),
                     by=c("ensemblID"="ensembl_gene_id"),
                     multiple="all") |>
        distinct(Feature, .keep_all = TRUE)
    return(dt)
}
memDF <- memoise::memoise(get_annotation)

get_genes <- function(tissue, sex){
    return(memDF(tissue, sex) |> dplyr::filter(DE == 1) |>
           dplyr::select(entrezgene_id, logFC) |>
           mutate(abs_beta=abs(logFC)) |> group_by(entrezgene_id) |>
           arrange(desc(abs_beta), .by_group=TRUE) |> slice(1) |>
           tidyr::drop_na() |> pull("entrezgene_id") |> as.character())
}

get_geneList <- function(tissue, sex){
    tmp <- memDF(tissue, sex) |> na.exclude() |>
        dplyr::select(entrezgene_id, logFC) |>
        mutate(abs_beta=abs(logFC)) |> group_by(entrezgene_id) |>
        arrange(desc(abs_beta), .by_group=TRUE) |> slice(1)
    geneList <- tmp |> pull("logFC")
    names(geneList) <- tmp$entrezgene_id
    return(geneList |> sort(decreasing=TRUE))
}

run_gseGO <- function(tissue, sex, ont, label){
                                        # GSEA with GO terms
    ego <- gseGO(geneList     = get_geneList(tissue, sex),
                 OrgDb        = org.Hs.eg.db,
                 ont          = ont,
                 keyType      = "ENTREZID",
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE,
                 seed         = TRUE,
                 eps          = 0)
                                        # Save results
    if(dim(ego)[1] != 0){
        fn = paste0(label, "_", ont, ".tsv")
        ego <- setReadable(ego, "org.Hs.eg.db")
        ego |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

run_enrichGO <- function(tissue, sex, ont, label){
    ego <- enrichGO(gene      = get_genes(tissue, sex),
                    OrgDb        = org.Hs.eg.db,
                    ont          = ont,
                    keyType      = "ENTREZID",
                    universe  = names(get_geneList(tissue, sex)),
                    minGSSize    = 10,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.05)

      # Save results
    if(dim(ego)[1] != 0){
        fn = paste0(label, "_", ont, ".tsv")
        ego <- setReadable(ego, "org.Hs.eg.db")
        ego |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

run_enrichDGN <- function(tissue, sex, label){
    dgn <- enrichDGN(gene      = get_genes(tissue, sex),
                     universe  = names(get_geneList(tissue, sex)),
                     minGSSize = 5)
      # Save results
    if(dim(dgn)[1] != 0){
        fn = paste0(label, "_DGN.tsv")
        dgn <- setReadable(dgn, "org.Hs.eg.db")
        dgn |> as.data.frame() |> data.table::fwrite(fn, sep='\t')
    }
}

run_gseDGN <- function(tissue, sex, ont, label){
    dgn <- gseDGN(geneList     = get_geneList(tissue, sex),
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

                                        # Female analysis
sex = "Female"
if(!dir.exists(tolower(sex))){dir.create(tolower(sex))}
for(tissue in c("Caudate", "DLPFC", "Hippocampus")){
    ##print(tissue)
    label1 = paste0(tolower(sex),"/",tolower(tissue),"_gsea")
    label2 = paste0(tolower(sex),"/",tolower(tissue),"_enrich")
    for(ONT in c("BP", "MF", "CC")){
        run_gseGO(tissue, sex, ONT, label1)
        if(tissue == "Caudate"){
            run_enrichGO(tissue, sex, ONT, label2)
        }
    }
    run_gseDGN(tissue, sex, "DGN", label1)
    if(tissue == "Caudate"){
        run_enrichDGN(tissue, sex, label2)
    }
}

                                        # Male analysis
sex = "Male"
if(!dir.exists(tolower(sex))){dir.create(tolower(sex))}

for(tissue in c("Caudate", "DLPFC", "Hippocampus")){
    ##print(tissue)
    label1 = paste0(tolower(sex),"/",tolower(tissue),"_gsea")
    label2 = paste0(tolower(sex),"/",tolower(tissue),"_enrich")
    for(ONT in c("BP", "MF", "CC")){
        run_gseGO(tissue, sex, ONT, label1)
        run_enrichGO(tissue, sex, ONT, label2)
    }
    run_gseDGN(tissue, sex, "DGN", label1)
    run_enrichDGN(tissue, sex, label2)
}

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

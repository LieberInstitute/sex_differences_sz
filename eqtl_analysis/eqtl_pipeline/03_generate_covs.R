## Prepare covariates for FastQTL / tensorQTL analysis
suppressPackageStartupMessages({
    library(here)
    library(argparse)
    library(tidyverse)
    library(SummarizedExperiment)
})

#### Functions
getRPKM <- function(rse, length_var = "bp_length", mapped_var = NULL) {
    mapped <- if (!is.null(mapped_var)) colData(rse)[, mapped_var] else colSums(assays(rse)$counts)
    bg <- matrix(mapped, ncol = ncol(rse), nrow = nrow(rse), byrow = TRUE)
    len <- if (!is.null(length_var)) rowData(rse)[, length_var] else width(rowRanges(rse))
    wid <- matrix(len, nrow = nrow(rse), ncol = ncol(rse), byrow = FALSE)
    return(assays(rse)$counts / (wid / 1000) / (bg / 1e6))
}

load_counts <- function(feature){
    counts_lt = list("genes"=here("input/counts/_m",
                                  "rse_gene.BrainSeq_Phase3.Gencode41.n464.Rdata"), 
                     "transcripts"=here("input/counts/_m",
                                        "rse_tx.BrainSeq_Phase3.Gencode41.n464.Rdata"),
                     "exons"=here("input/counts/_m",
                                  "rse_exon.BrainSeq_Phase3.Gencode41.n464.Rdata"),
                     "junctions"=here("input/counts/_m",
                                      "rse_jxn.BrainSeq_Phase3.Gencode41.n464.Rdata"))
    load(counts_lt[[feature]])
    if(exists("rse_gene")){
        rse_df = rse_gene
    } else if (exists("rse_tx")){
        rse_df = rse_tx
    } else if (exists("rse_exon")){
        rse_df = rse_exon
    } else {
        rse_df = rse_jxn
    }
    return(rse_df)
}

filter_data <- function(feature){
                                        # Load counts
    rse_df  <- load_counts(feature)
                                        # Sample selection
    keepInd <- which((colData(rse_df)$Age > 13) &
                      (colData(rse_df)$Race %in% c("AA", "CAUC")) &
                      (colData(rse_df)$Dx %in% c("Control", "SCZD")))
    rse_df  <- rse_df[,keepInd]
    if(feature == "junction"){
        rowRanges(rse_df)$Length <- 100
        rpkm_df   <- getRPKM(rse_df, 'Length')
        rse_df  <- rse_df[rowMeans(rpkm_df) > 0.4,]; rm(rpkm_df)
    } else if (feature == "transcript"){
        rse_df   <- rse_df[rowMeans(assays(rse_df)$tpm) > 0.2,]
    } else {
        rpkm_df <- getRPKM(rse_df, "Length")
        rse_df  <- rse_df[rowMeans(rpkm_df) > 0.2,]; rm(rpkm_df)
    }
    return(rse_df)
}
memDF <- memoise::memoise(filter_data)

normalize_data <- function(feature){
    rse_df <- memDF(feature)
    if(feature == "transcript"){
        norm_df <- assays(rse_df)$tpm
    } else {
        norm_df <- getRPKM(rse_df, 'Length')
    }
    colnames(norm_df) <- gsub("\\_.*", "", colnames(norm_df))
    return(norm_df)
}

load_mds <- function(){
    gfile  <- here::here("input/genotypes/mds/_m/LIBD_Brain_TopMed.mds")
    mds    <- data.table::fread(gfile) |>
        rename_at(.vars = vars(starts_with("C")),
                  function(x){sub("C", "snpPC", x)}) |>
        mutate_if(is.character, as.factor)
    return(mds)
}

get_pheno <- function(feature){    
    return(colData(memDF(feature)) |> as.data.frame() |>
           inner_join(load_mds(), by=c("BrNum"="FID"), multiple="all") |>
           distinct(RNum, .keep_all = TRUE) |> mutate(ids=RNum) |>
           column_to_rownames("ids"))
}

cal_pca <- function(norm_df, new_pd, mod){
    pca_df <- prcomp(t(log2(norm_df[, new_pd$RNum]+1)))
    if(dim(norm_df)[1] > 50000){
        k  <- sva::num.sv(log2(norm_df[, new_pd$RNum]+1),
                              mod, vfilter=50000)
    } else {
        k  <- sva::num.sv(log2(norm_df[, new_pd$RNum]+1), mod)
    }
    pcs    <- pca_df$x[,1:k]
    return(pcs)
}

#### Main
                                        # Create parser object
parser <- ArgumentParser()
parser$add_argument("-f", "--feature", type="character", default="genes",
                    help="Feature to extract [default: %default]")
parser$add_argument("-r", "--region", type="character", default="caudate",
                    help="Brain region to extract [default: %default]")
args   <- parser$parse_args()
region <- args$region; feature <- args$feature

                                        # Load data
norm_df <- normalize_data(feature)
if(feature == "junction"){ norm_df <- as.matrix(norm_df)}

                                        # Full model
new_pd  <- get_pheno(feature)
mod     <- model.matrix(~Sex + Dx + Age + snpPC1 + snpPC2 + snpPC3,
                        data = new_pd)

                                        # PCA
pcs    <- cal_pca(norm_df, new_pd, mod)

                                        # Subset samples
sample_df <- data.table::fread("../../_m/sample_id_to_brnum.tsv")

                                        # Extract covariates
covs <- cbind(mod[,c(-1, -2)], pcs) |> as.data.frame() |>
    rownames_to_column("RNum") |>
    inner_join(sample_df, by=c("RNum")) |> select(-"RNum") |>
    rename("BrNum"="ID") |> column_to_rownames("ID") |> t() |>
    as.data.frame() |> rownames_to_column("ID")
covs <- covs[, c("ID", sample_df$BrNum)]

                                        # Save file
data.table::fwrite(covs, paste0(feature, ".combined_covariates.txt"),
                   sep='\t')

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

#### Residualize data
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(tidyverse)
    library(here)
})

## Helper functions for loading and combining data
load_rse <- function(filename){
    load(filename)
    rse_df = rse_gene
    colData(rse_df) = colData(rse_df)[,c('BrNum', 'RNum', 'Region', 'Dx',
                                         'Age', 'Sex', 'Race')]
    rowData(rse_df)$meanExprs = NULL
    return(rse_df)
}

combine_expression <- function(){
    caudate <- load_rse(here("input/counts/_m",
                             "rse_gene.BrainSeq_Phase3.Gencode41.n464.Rdata"))
    dlpfc   <- load_rse(here("input/counts/_m",
                             "rse_gene.BrainSeq_Phase2_DLPFC.Gencode41.n453.Rdata"))
    hippo   <- load_rse(here("input/counts/_m",
                             "rse_gene.BrainSeq_Phase2_HIPPO.Gencode41.n447.Rdata"))
    rse_df  <- cbind(caudate, dlpfc, hippo)
    colnames(rse_df) <- gsub("\\_.*", "", colnames(rse_df))
    return(list("gene"=rse_df))
}

filter_data <- function(){
                                        # Load combined data
    rse_df  <- combine_expression()[["gene"]]
                                        # Sample selection
    keepInd <- which((colData(rse_df)$Age > 13) &
                      (colData(rse_df)$Race %in% c("AA", "CAUC")) &
                      (colData(rse_df)$Dx %in% c("Control", "SCZD")))
    rse_df  <- rse_df[,keepInd]
                                        # Filter based on gene expression
    rpkm_df <- recount::getRPKM(rse_df, "Length")
    rse_df  <- rse_df[rowMeans(rpkm_df) > 0.2,]; rm(rpkm_df)
    return(rse_df)
}
memDF <- memoise::memoise(filter_data)

normalize_data <- function(){
    rse_df  <- memDF()
    norm_df <- recount::getRPKM(rse_df, 'Length')
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

get_pheno <- function(){    
    return(colData(memDF()) |> as.data.frame() |>
           inner_join(load_mds(), by=c("BrNum"="FID"),
                      multiple="all", relationship="many-to-many") |>
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
                                        # Load data
norm_df <- normalize_data()
expr_df <- log2(norm_df+1)
                                        # Full model
new_pd  <- get_pheno()
mod     <- model.matrix(~0+Dx + Sex + Age + Region + snpPC1 +
                            snpPC2 + snpPC3, data = new_pd)

                                        # PCA
pcs     <- cal_pca(norm_df, new_pd, mod)

                                        # Covariates
covs0   <- t(cbind(mod[,-1], pcs))

covs    <- covs0 |> as.data.frame() |>
    rownames_to_column("Vars") |>
    filter(!str_detect(Vars, "Region"), !str_detect(Vars, "Sex")) |>
    column_to_rownames("Vars") |> as.matrix()

                                        # Write residualized expression
null_mod <- limma::lmFit(expr_df[, colnames(covs)],
                         t(rbind(intercept=1, covs)))
res_df   <- residuals(null_mod, expr_df[, colnames(covs)]) |>
    as.data.frame() |> rownames_to_column("gene_id")

                                        # Save file
data.table::fwrite(res_df, "genes_residualized_expression.csv", sep=',')

#### Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

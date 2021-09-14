suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(tidyverse)
    library(jaffelab)
    library(sva)
})

## Helper functions for loading and combining data
load_rse <- function(filename, feature){
    load(filename)
    if(feature == "gene"){
        rse_df = rse_gene
    } else if(feature == "tx"){
        rse_df = rse_tx
    } else if(feature == "exon"){
        rse_df = rse_exon
    } else {
        rse_df = rse_jxn
    }
    colData(rse_df) = colData(rse_df)[,c('BrNum', 'RNum', 'Region', 'Dx',
                                         'Age', 'Sex', 'Race')]
    rowData(rse_df)$meanExprs = NULL
    return(rse_df)
}

common_rse <- function(x, y, z){
    common_rows = intersect(intersect(rownames(x), rownames(y)),
                            rownames(z))
    new_row_ranges = rowRanges(x)[common_rows,]
    new_colData = rbind(colData(x), colData(y), colData(z))
    common_assays = intersect(intersect(names(assays(x)), names(assays(y))),
                              names(assays(z)))
    new_assays = SimpleList()
    for( nn in common_assays ){
        new_assay = SimpleList(cbind(assays(x)[[nn]][common_rows,],
                                     assays(y)[[nn]][common_rows,],
                                     assays(z)[[nn]][common_rows,]))
        new_assays = c(new_assays, new_assay)
    }
    names(new_assays) = common_assays
    SummarizedExperiment(rowRanges = new_row_ranges,
                         colData = new_colData,
                         assays = new_assays)
}

combine_expression <- function(){
    counts_dir = "/ceph/projects/v4_phase3_paper/inputs/counts/_m/"
    config_caudate <- list(
        "gene"=paste0(counts_dir,
                      "caudate_brainseq_phase3_hg38_rseGene_merged_n464.rda"),
        "tx"=paste0(counts_dir,
                    "caudate_brainseq_phase3_hg38_rseTx_merged_n464.rda"),
        "exon"=paste0(counts_dir,
                      "caudate_brainseq_phase3_hg38_rseExon_merged_n464.rda"),
        "jxn"=paste0(counts_dir,
                     "caudate_brainseq_phase3_hg38_rseJxn_merged_n464.rda")
    )
    config_dlpfc = list(
        "gene"=paste0(counts_dir,
                      "dlpfc_ribozero_brainseq_phase2_hg38_rseGene_merged_n453.rda"),
        "tx"=paste0(counts_dir,
                    "dlpfc_ribozero_brainseq_phase2_hg38_rseTx_merged_n453.rda"),
        "exon"=paste0(counts_dir,
                      "dlpfc_ribozero_brainseq_phase2_hg38_rseExon_merged_n453.rda"),
        "jxn"=paste0(counts_dir,
                     "dlpfc_ribozero_brainseq_phase2_hg38_rseJxn_merged_n453.rda")
        )
    config_hippo = list(
        "gene"=paste0(counts_dir,
                      "hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda"),
        "tx"=paste0(counts_dir,
                    "hippo_brainseq_phase2_hg38_rseTx_merged_n447.rda"),
        "exon"=paste0(counts_dir,
                      "hippo_brainseq_phase2_hg38_rseExon_merged_n447.rda"),
        "jxn"=paste0(counts_dir,
                     "hippo_brainseq_phase2_hg38_rseJxn_merged_n447.rda")
    )
    rse_gene = common_rse(load_rse(config_caudate[["gene"]], "gene"),
                          load_rse(config_dlpfc[["gene"]], "gene"),
                          load_rse(config_hippo[["gene"]], "gene"))
    rse_tx = common_rse(load_rse(config_caudate[["tx"]], "tx"),
                        load_rse(config_dlpfc[["tx"]], "tx"),
                        load_rse(config_hippo[["tx"]], "tx"))
    rse_exon = common_rse(load_rse(config_caudate[["exon"]], "exon"),
                          load_rse(config_dlpfc[["exon"]], "exon"),
                          load_rse(config_hippo[["exon"]], "exon"))
    rse_jxn = common_rse(load_rse(config_caudate[["jxn"]], "jxn"),
                         load_rse(config_dlpfc[["jxn"]], "jxn"),
                         load_rse(config_hippo[["jxn"]], "jxn"))
    return(list("gene"=rse_gene, "tx"=rse_tx, "exon"=rse_exon, "jxn"=rse_jxn))
}


filter_convert <- function(){
    ## Load combined data
    rse_lt <- combine_expression()
    ## Filter based on gene expression
    rse_gene <- rse_lt[["gene"]]
    tmp_gene_rpkm = log2(recount::getRPKM(rse_gene, "Length") + 1)
    rse_gene = rse_gene[rowMeans(tmp_gene_rpkm) > 0.2,]
    rm(tmp_gene_rpkm)
    rse_exon <- rse_lt[["exon"]]
    tmp_exon_rpkm = log2(recount::getRPKM(rse_exon, "Length") + 1)
    rse_exon = rse_exon[rowMeans(tmp_exon_rpkm) > 0.2,]
    rm(tmp_exon_rpkm)
    rse_jxn = rse_lt[["jxn"]]
    rowRanges(rse_jxn)$Length <- 100
    jRp10m = recount::getRPKM(rse_jxn, 'Length')
    rse_jxn = rse_jxn[rowMeans(jRp10m) > 0.4,]
    rse_tx <- rse_lt[["tx"]]
    rse_tx = rse_tx[rowMeans(assays(rse_tx)$tpm) > 0.4,]
    ## Sample selection
    keepInd = which((colData(rse_gene)$Age > 13) &
                    (colData(rse_gene)$Dx %in% c('Schizo', 'Control')))
    rse_gene = rse_gene[,keepInd]
    rse_exon = rse_exon[,keepInd]
    rse_jxn = rse_jxn[,keepInd]
    rse_tx = rse_tx[,keepInd]
    ## extract pd and rpkms
    pd = colData(rse_gene)
    geneRpkm = recount::getRPKM(rse_gene, "Length")
    exonRpkm = recount::getRPKM(rse_exon, "Length")
    jxnRp10m = recount::getRPKM(rse_jxn, 'Length')
    txTpm = assays(rse_tx)$tpm
                                        # save expression
    geneExpression = log2(geneRpkm+1)
    exonExpression = log2(exonRpkm+1)
    jxnExpression = log2(jxnRp10m+1)
    txExpression = log2(txTpm+1)
    save(geneExpression, exonExpression, jxnExpression, txExpression,
         file="expression.rda")
    return(pd)
}

generate_model <- function(){
    ## Statistical model
                                        # Add MDS to phenotyeps
    pd <- filter_convert()
    mds_file = paste0("/ceph/projects/v4_phase3_paper/inputs/genotypes/mds/",
                      "_m/LIBD_Brain_TopMed.mds")
    mds = data.table::fread(mds_file) %>%
        rename_at(.vars = vars(starts_with("C")),
                  function(x){sub("C", "snpPC", x)}) %>%
        mutate_if(is.character, as.factor)
    new_pd = pd %>% as.data.frame %>%
        inner_join(mds, by=c("BrNum"="FID")) %>%
        distinct(RNum, .keep_all = TRUE) %>% mutate(ids=RNum) %>%
        column_to_rownames("ids")
                                        # Model
    new_pd$Dx = factor(new_pd$Dx,levels = c("Control", "Schizo"))
    mod = model.matrix(~0+Dx + Sex + Region + snpPC1 + snpPC2 + snpPC3 +
                           snpPC4 + snpPC5, data = new_pd)
    return(list("model"=mod, "phenotypes"=new_pd))
}

calculate_pca <- function(){
    ## Load data
    model_lt = generate_model()
    load("expression.rda")
    new_pd = model_lt[["phenotypes"]]
    mod = model_lt[["model"]]
    ## PCA
                                        # Gene
    pcaGene = prcomp(t(geneExpression[, new_pd$RNum]))
    kGene = num.sv(geneExpression[, new_pd$RNum], mod)
    genePCs = pcaGene$x[,1:kGene]
                                        # Exon
    pcaExon = prcomp(t(exonExpression[, new_pd$RNum]))
    kExon = num.sv(exonExpression[, new_pd$RNum], mod, vfilter=50000)
    exonPCs = pcaExon$x[,1:kExon]
                                        # Junction
    pcaJxn = prcomp(t(jxnExpression[, new_pd$RNum]))
    kJxn = num.sv(jxnExpression[, new_pd$RNum], mod, vfilter=50000)
    jxnPCs = pcaJxn$x[,1:kJxn]
                                        # Transcript
    pcaTx = prcomp(t(txExpression[, new_pd$RNum]))
    kTx = num.sv(txExpression[, new_pd$RNum], mod, vfilter=50000)
    txPCs = pcaTx$x[,1:kTx]
                                        # Save PCs
    save(mod, genePCs, exonPCs, jxnPCs, txPCs,
         file="pcs_caudate_4features_filtered_over13.rda")
}

calculate_covs <- function(){
                                        # Calculate PCs
    calculate_pca()
    load("pcs_caudate_4features_filtered_over13.rda")
                                        # Combind Covariates
    covsGene0 = t(cbind(mod[,-1],genePCs))
    covsExon0 = t(cbind(mod[,-1],exonPCs))
    covsJxn0 = t(cbind(mod[,-1],jxnPCs))
    covsTx0 = t(cbind(mod[,-1],txPCs))
                                        # Save covariates
    save(covsGene0, covsExon0, covsJxn0, covsTx0, file="covariates.rda")
}


write_residualized <- function(){
                                        # Run covariates prep
    calculate_covs()
                                        # Load data
    load("covariates.rda")
    load("expression.rda")
                                        # Genes
    covsGene = covsGene0 %>% as.data.frame %>%
        rownames_to_column("Vars") %>% filter(!str_detect(Vars, "Region"),
                                              !str_detect(Vars, "Sex")) %>%
        column_to_rownames("Vars") %>% as.matrix
    mod = limma::lmFit(geneExpression[, colnames(covsGene)],
                       t(rbind(intercept=1, covsGene)))
    residuals(mod, geneExpression[, colnames(covsGene)]) %>%
        as.data.frame %>% rownames_to_column("gene_id") %>%
        data.table::fwrite("genes_residualized_expression.csv", sep=',')
                                        # Transcripts
    covsTx = covsTx0 %>% as.data.frame %>%
        rownames_to_column("Vars") %>% filter(!str_detect(Vars, "Region"),
                                              !str_detect(Vars, "Sex")) %>%
        column_to_rownames("Vars") %>% as.matrix
    mod = limma::lmFit(txExpression[, colnames(covsTx)],
                       t(rbind(intercept=1, covsTx)))
    residuals(mod, txExpression[, colnames(covsTx)]) %>%
        as.data.frame %>% rownames_to_column("gene_id") %>%
        data.table::fwrite("transcripts_residualized_expression.csv", sep=',')
                                        # Exons
    covsExon = covsExon0 %>% as.data.frame %>%
        rownames_to_column("Vars") %>% filter(!str_detect(Vars, "Region"),
                                              !str_detect(Vars, "Sex")) %>%
        column_to_rownames("Vars") %>% as.matrix
    mod = limma::lmFit(exonExpression[, colnames(covsExon)],
                       t(rbind(intercept=1, covsExon)))
    residuals(mod, exonExpression[, colnames(covsExon)]) %>%
        as.data.frame %>% rownames_to_column("gene_id") %>%
        data.table::fwrite("exons_residualized_expression.csv", sep=',')
                                        # Junctions
    covsJxn = covsJxn0 %>% as.data.frame %>%
        rownames_to_column("Vars") %>% filter(!str_detect(Vars, "Region")) %>%
        column_to_rownames("Vars") %>% as.matrix
    mod = limma::lmFit(jxnExpression[, colnames(covsJxn)],
                       t(rbind(intercept=1, covsJxn)))
    residuals(mod, jxnExpression[, colnames(covsJxn)]) %>%
        as.data.frame %>% rownames_to_column("gene_id") %>%
        data.table::fwrite("junctions_residualized_expression.csv", sep=',')
}

write_residualized()

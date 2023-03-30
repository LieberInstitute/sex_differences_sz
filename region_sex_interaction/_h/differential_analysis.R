## Differential expression analysis with limma-voom
## Sex-specific analysis for specific brain region
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(variancePartition)
    library(BiocParallel)
    library(optparse)
    library(dplyr)
    library(limma)
    library(edgeR)
    library(here)
    library(glue)
    library(sva)
})

load_counts <- function(region, feature){
    if(region == "caudate"){
        counts_lt = list("genes"=here("input/counts/_m",
                                      paste0("rse_gene.BrainSeq_Phase3",
                                             ".Gencode41.n464.Rdata")), 
                         "transcripts"=here("input/counts/_m",
                                            paste0("rse_tx.BrainSeq_Phase3",
                                                   ".Gencode41.n464.Rdata")),
                         "exons"=here("input/counts/_m",
                                      paste0("rse_exon.BrainSeq_Phase3.",
                                             "Gencode41.n464.Rdata")),
                         "junctions"=here("input/counts/_m",
                                          paste0("rse_jxn.BrainSeq_Phase3.",
                                                 "Gencode41.n464.Rdata")))
    } else if(region == "dlpfc"){
        counts_lt = list("genes"=here("input/counts/_m",
                                      paste0("rse_gene.BrainSeq_Phase2_DLPFC",
                                             ".Gencode41.n453.Rdata")), 
                         "transcripts"=here("input/counts/_m",
                                            paste0("rse_tx.BrainSeq_Phase2_",
                                                   "DLPFC.Gencode41.n453.Rdata")),
                         "exons"=here("input/counts/_m",
                                      paste0("rse_exon.BrainSeq_Phase2_DLPFC",
                                             ".Gencode41.n453.Rdata")),
                         "junctions"=here("input/counts/_m",
                                          paste0("rse_jxn.BrainSeq_Phase2_",
                                                 "DLPFC.Gencode41.n453.Rdata")))
    } else {
        counts_lt = list("genes"=here("input/counts/_m",
                                      paste0("rse_gene.BrainSeq_Phase2_HIPPO",
                                             ".Gencode41.n447.Rdata")), 
                         "transcripts"=here("input/counts/_m",
                                            paste0("rse_tx.BrainSeq_Phase2_HIPPO",
                                                   ".Gencode41.n447.Rdata")),
                         "exons"=here("input/counts/_m",
                                      paste0("rse_exon.BrainSeq_Phase2_HIPPO",
                                             ".Gencode41.n447.Rdata")),
                         "junctions"=here("input/counts/_m",
                                          paste0("rse_jxn.BrainSeq_Phase2_HIPPO",
                                                 ".Gencode41.n447.Rdata")))
    }
    return(list("count_file"=counts_lt[[feature]]))
}

get_juncs <- function(region){
    lt     <- load_counts(region, "junctions")
    load(lt[['count_file']])
    rse_df <- rse_jxn
    genes  <- rowData(rse_df) |> as.data.frame()
    return(rownames(genes))
}

get_model <- function(){
    covs = c('Dx', 'Age', 'mitoRate', 'rRNA_rate', 'totalAssignedGene',
             'RIN',  'Mapping_Rate', 'Mean3Bias', 'snpPC1', 
             'snpPC2', 'snpPC3', 'Sex', 'Region', 'BrNum')
    adjust_covs = covs[1:11]
    random_effect = covs[14]
    formula <- glue("~ Sex * Region + ", 
                    glue_collapse(adjust_covs, sep = " + ")) %>% 
        glue(" + (1|{random_effect})")
    return(formula)
}

get_null <- function(){
    covs = c('Dx', 'Age', 'mitoRate', 'rRNA_rate', 'totalAssignedGene',
             'RIN',  'Mapping_Rate', 'Mean3Bias', 'snpPC1', 
             'snpPC2', 'snpPC3', 'Sex', 'Region', 'BrNum')
    adjust_covs = covs[1:11]
    random_effect = covs[14]
    null_form <- glue("~ ", glue_collapse(adjust_covs[2:11], sep = " + ")) %>% 
      glue(" + (1|{random_effect}) + (1|Dx)")
    return(null_form)
}

save_volcanoPlot <- function(top, label, feature){
    pdf(file=paste0(feature, "/volcanoPlot_", label, ".pdf"), 8, 6)
    with(top, plot(logFC, -log10(P.Value), pch=20, cex=0.6))
    with(subset(top, adj.P.Val<=0.05), points(logFC, -log10(P.Value),
                                              pch=20, col='red', cex=0.6))
    with(subset(top, abs(logFC)>0.50), points(logFC, -log10(P.Value),
                                              pch=20, col='orange', cex=0.6))
    with(subset(top, adj.P.Val<=0.05 & abs(logFC)>0.50),
         points(logFC, -log10(P.Value), pch=20, col='green', cex=0.6))
    dev.off()
}

save_MAplot <- function(top, label, feature){
    pdf(file=paste0(feature, "/MAplot_", label, ".pdf"), 8, 6)
    with(top, plot(AveExpr, logFC, pch=20, cex=0.5))
    with(subset(top, adj.P.Val<0.05),
         points(AveExpr, logFC, col="red", pch=20, cex=0.5))
    dev.off()
}

extract_de <- function(contrast, label, efit, feature){
    top     <- topTable(efit, coef=contrast, number=Inf, sort.by="P")
    if(feature == "junctions"){
        top <- top |> select(!contains("Exon"))
    }
    top$SE  <- efit$sigma * pull(as.data.frame(efit$stdev.unscaled), contrast)
    top     <- top[order(top$P.Value), ]
    top.fdr <- top |> filter(adj.P.Val<=0.05)
    print(paste("Comparison for:", label))
    print(paste('There are:', dim(top.fdr)[1], 'DE features!'))
    data.table::fwrite(top, file=paste0(feature, "/diffExpr_", label, "_full.txt"), 
                       sep='\t', row.names=TRUE)
    data.table::fwrite(top.fdr, file=paste0(feature, "/diffExpr_", label,
                                            "_FDR05.txt"), 
                       sep='\t', row.names=TRUE)
    save_volcanoPlot(top, label, feature)
    save_MAplot(top, label, feature)
}

get_mds <- function(){
    mds = data.table::fread(here("input/genotypes/mds/_m/LIBD_Brain_TopMed.mds")) |>
        rename_at(.vars = vars(starts_with("C")),
                  function(x){sub("C", "snpPC", x)}) |>
        mutate_if(is.character, as.factor)
    return(mds)
}
memMDS <- memoise::memoise(get_mds)

get_antipsychotics <- function(){
    anti_file <- here("input/phenotypes/_h/antipsychotics_sample_annotation.csv")
    anti_df  <- data.table::fread(anti_file) |>
        select(c("BrNum", "New_Dx", "antipsychotics", "lifetime_antipsych"))
    return(anti_df)
}
memANTI <- memoise::memoise(get_antipsychotics)

get_overlapping_juncs <- function(){
    common_genes <- intersect(intersect(get_juncs("caudate"),
                                        get_juncs("dlpfc")),
                              get_juncs("hippocampus"))
    return(common_genes)
}
memJUNCS <- memoise::memoise(get_overlapping_juncs)

prep_data <- function(feature){
    datalist1 <- list(); datalist2 <- list()
    for(region in c("dlpfc", "caudate", "hippocampus")){
        print(region)
        flush.console()
        lt <- load_counts(region, feature)
        load(lt[["count_file"]])
        if(exists("rse_gene")){
            rse_df = rse_gene
            rm(rse_gene)
        } else if (exists("rse_tx")){
            rse_df = rse_tx
            rm(rse_tx)
        } else if (exists("rse_exon")){
            rse_df = rse_exon
            rm(rse_exon)
        } else {
            rse_df = rse_jxn
            rm(rse_jxn)
        }
        keepIndex = which((rse_df$Dx %in% c("Control", "SCZD")) & 
                          rse_df$Age > 17 & rse_df$Race %in% c("AA", "CAUC"))
        rse_df = rse_df[, keepIndex]
        rse_df$Dx = factor(rse_df$Dx, levels = c("Control", "SCZD"))
        rse_df$Sex <- factor(rse_df$Sex)
        rownames(colData(rse_df)) <- sapply(strsplit(rownames(colData(rse_df)), "_"),
                                            "[", 1)
        cols <- c('Dx', 'Age', 'mitoRate', 'rRNA_rate', 'totalAssignedGene',
                  'RIN',  'Mapping_Rate', 'Mean3Bias', 'snpPC1', 
                  'snpPC2', 'snpPC3', 'Sex', 'Region', 'BrNum', "RNum")
        pheno <- colData(rse_df) |> as.data.frame() |> 
            inner_join(memMDS(), by=c("BrNum"="FID"), multiple="all") |> 
            distinct(RNum, .keep_all = TRUE) |>
            left_join(memANTI(), by="BrNum") |>
            select(all_of(cols))
        if(feature == "junctions"){
            datalist1[[region]] <-  assays(rse_df)$counts[memJUNCS(),
                                                          pheno$RNum]
        } else {
            datalist1[[region]] <- assays(rse_df)$counts[, pheno$RNum]
        }
        datalist2[[region]] <- pheno
    }
                                        # Generate DGE list
    samples  <- bind_rows(datalist2) |> as.data.frame() |>
        mutate_if(is.numeric, scales::rescale) |>
        tibble::column_to_rownames("RNum")
    subjects <- samples$BrNum
    if(feature == "junctions"){
        genes <- rowData(rse_df)[memJUNCS(), ]
    } else {
        genes <- rowData(rse_df)
    }
    counts   <- do.call(cbind, datalist1)
    rm(datalist1, datalist2, rse_df, pheno)
    x <- DGEList(counts=counts, genes=genes, samples=samples)
    # Filter by expression
    design0 <- model.matrix(~Sex*Region, data=x$samples)
    keep.x <- filterByExpr(x, design=design0)
    print(paste('There are:', sum(keep.x), 'features left!', sep=' '))
    x <- x[keep.x, , keep.lib.sizes=FALSE]
    # Normalize library size
    x <- calcNormFactors(x, method="TMM")
    return(x)
}
memPREP <- memoise::memoise(prep_data)

plot_corr_matrix <- function(feature){
    x    <- memPREP(feature)
    covs <- c('Dx', 'Age', 'mitoRate', 'rRNA_rate', 'totalAssignedGene',
             'RIN',  'Mapping_Rate', 'Mean3Bias', 'snpPC1', 
             'snpPC2', 'snpPC3', 'Sex', 'Region', 'BrNum')
    form <- paste("~", glue_collapse(covs, sep=" + "))
    C    <- canCorPairs(form, x$samples)
    pdf(file = paste0(feature, "/correlation_matrix.pdf"),
        width = 6, height = 6)
    plotCorrMatrix(C)
    dev.off()
}
memCORR <- memoise::memoise(plot_corr_matrix)

get_voom <- function(feature){
### Preform voom
                                        # Plot correlation of covariates
    memCORR(feature)
                                        # Load features
    x    <- memPREP(feature)
    vobj <- voomWithDreamWeights(x, get_model(),
                                 x$samples, plot=FALSE)
    return(vobj)
}
memVOOM <- memoise::memoise(get_voom)

fit_dream <- function(feature){
    x    <- memPREP(feature)
    vobj <- memVOOM(feature)
                                        # define contrasts
    L1   <- getContrast(vobj, get_model(), x$samples,
                        c("RegionDLPFC", "RegionHIPPO"))
    L2   <- getContrast(vobj, get_model(), x$samples,
                        c("SexM:RegionDLPFC", "SexM:RegionHIPPO"))
                                        # combine contrasts
    L    <- cbind(L1, L2)
                                        # visualize
    pdf(file = paste0(feature, "/plot_contrasts.pdf"),
        width = 6, height = 6)
    plotContrasts(L)
    dev.off()
                                        # fit model
    fit  <- dream(vobj, get_model(), x$samples, L=L)
    return(fit)
}
memDREAM <- memoise::memoise(fit_dream)

cal_res <- function(feature){
    ### Calculate residuals
    vobj      <- memVOOM(feature)
    residList <- fitVarPartModel(vobj, get_null(), vobj$targets, fxn=residuals)
    do.call(rbind, residList) |> as.data.frame() |> 
        tibble::rownames_to_column("feature_id") |>
        write.table(file=paste0(feature, '/residualized_expression.tsv'), 
                    sep="\t", quote=FALSE, row.names=FALSE)
}
memRES <- memoise::memoise(cal_res)

#### MAIN
                                        # Input parser
option_list <- list(
    make_option(c("--feature"), type="character", default="genes",
                help="feature to analyze [default=%default]",
                metavar="character"),
    make_option(c("--threads"), type="integer", default=16,
                help="threads to use [default=%default]",
                metavar="integer")
)
opt_parser  <- OptionParser(option_list=option_list)
opt         <- parse_args(opt_parser)
                                        # Set parallel
param <- SnowParam(opt$threads, "SOCK", progressbar=FALSE)
register(param)
                                        # Differential expression by feature
feature <- opt$feature
if(!dir.exists(feature)){dir.create(feature)}
                                        # Preform voom
print("Prepare data!")
x   <- memPREP(feature)
print("Preform voom!")
v   <- memVOOM(feature)
save(v, file=paste0(feature,'/voom_dream.RData'))
                                        # Fit model and apply eBayes
print("Fit model!")
fit <- memDREAM(feature)
                                        # Save differential expression
extract_de("SexM", "sex", fit, feature)
extract_de("SexM:RegionDLPFC", "CvD_sex", fit, feature)
extract_de("SexM:RegionHIPPO", "CvH_sex", fit, feature)
extract_de("L2", "DvH_sex", fit, feature)
                                        # Calculate residuals
print("Residualized expression!")
memRES(feature)

#### Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

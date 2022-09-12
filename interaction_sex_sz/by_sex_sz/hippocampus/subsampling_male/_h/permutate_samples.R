## Differential expression with limma (permutation, subsampling)

suppressMessages({
    library(SummarizedExperiment)
    library(tidyverse)
    library(argparse)
    library(limma)
    library(edgeR)
    library(sva)
})

## Functions
# Function from jaffelab github
merge_rse_metrics <- function(rse) {
    stopifnot(is(rse, 'RangedSummarizedExperiment'))
    stopifnot(
        c('concordMapRate', 'overallMapRate', 'mitoRate', 'rRNA_rate',
          'totalAssignedGene', 'numMapped', 'numReads', 'numUnmapped',
          'mitoMapped', 'totalMapped') %in%
            colnames(SummarizedExperiment::colData(rse))
    )

    stopifnot(all(sapply(c(
        'concordMapRate', 'overallMapRate', 'mitoRate', 'rRNA_rate',
        'totalAssignedGene', 'numMapped', 'numReads', 'numUnmapped',
        'mitoMapped', 'totalMapped'), function(var) {
            is(colData(rse)[, var], 'List')
        })
    ))

    rse$concordMapRate = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$concordMapRate, rse$numReads)
    rse$overallMapRate = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$overallMapRate, rse$numReads)
    rse$mitoRate = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$mitoRate, rse$numMapped)
    rse$rRNA_rate = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$rRNA_rate, rse$numMapped)
    rse$totalAssignedGene = mapply(function(r, n) {
        sum(r*n)/sum(n)
    }, rse$totalAssignedGene, rse$numMapped)

    rse$numMapped = sapply(rse$numMapped, sum)
    rse$numReads = sapply(rse$numReads, sum)
    rse$numUnmapped = sapply(rse$numUnmapped, sum)
    rse$mitoMapped = sapply(rse$mitoMapped, sum)
    rse$totalMapped = sapply(rse$totalMapped, sum)
    return(rse)
}

extract_de <- function(contrast, label, efit, seed_int){
    perm_name <- paste0("permutation_", stringr::str_pad(seed_int, 4, pad = "0"))
    se <- as.data.frame(sqrt(efit$s2.post) * efit$stdev.unscaled) %>%
        rename("SE"="CtrlvsSZ") %>% rownames_to_column()
    top <- topTable(efit, coef=1, number=Inf, sort.by="P") %>%
        rownames_to_column() %>% inner_join(se, by="rowname") %>%
        column_to_rownames("rowname")
    top.fdr <- top %>% filter(adj.P.Val<=0.05)
    print(paste("Comparison for:", label))
    print(paste('There are:', dim(top.fdr)[1], 'DE features!'))
    data.table::fwrite(top, 
                       file=paste0(perm_name, "_diffExpr_", label, "_full.txt"), 
                       sep='\t', row.names=TRUE)
}

load_counts <- function(){
    counts_file = paste0("/dcs04/lieber/statsgen/jbenjami/projects/",
                         "sex_differences_sz/input/counts/_m/", 
                         "hippo_brainseq_phase2_hg38_rseGene_merged_n447.rda")
    load(counts_file)
    rse_df = rse_gene
    return(rse_df)
}
memCounts <- memoise::memoise(load_counts)

get_random_samples <- function(seed_int, new_dir=TRUE){
    set.seed(seed_int + 113) # seed for reproducibility
    rse_df <- memCounts()
    keepIndex = which((rse_df$Dx %in% c("Control", "Schizo")) & 
                      rse_df$Age > 17 & rse_df$Sex == "M" & 
                      rse_df$Race %in% c("AA", "CAUC"))
    snames = sample(keepIndex, 121, replace=FALSE) # subsampling to Female N (sample size)
    return(snames)
}
memSamples <- memoise::memoise(get_random_samples)

get_mds <- function(){
    mds_file = paste0("/dcs04/lieber/statsgen/jbenjami/projects/sex_differences_sz/",
                      "input/genotypes/mds/_m/LIBD_Brain_TopMed.mds")
    mds = data.table::fread(mds_file) %>%
        rename_at(.vars = vars(starts_with("C")),
                  function(x){sub("C", "snpPC", x)}) %>%
        mutate_if(is.character, as.factor)
    return(mds)
}
memMDS <- memoise::memoise(get_mds)

prep_data <- function(seed_int){
    rse_df <- memCounts()
    keepIndex <- memSamples(seed_int)
    rse_df = rse_df[, keepIndex]
    rse_df$Dx = factor(rse_df$Dx, levels = c("Control", "Schizo"))
    rse_df$Sex <- factor(rse_df$Sex)
    rse_df <- merge_rse_metrics(rse_df)
    rse_df$ERCCsumLogErr <- mapply(function(r, n) {
        sum(r * n)/sum(n)
    }, rse_df$ERCCsumLogErr, rse_df$numReads)
    colData(rse_df)$RIN = sapply(colData(rse_df)$RIN,"[",1)
    rownames(colData(rse_df)) <- sapply(strsplit(rownames(colData(rse_df)), "_"), "[", 1)
    pheno = colData(rse_df) %>% as.data.frame %>% 
        inner_join(memMDS(), by=c("BrNum"="FID")) %>% 
        distinct(RNum, .keep_all = TRUE) 
    # Generate DGE list
    x <- DGEList(counts=assays(rse_df)$counts[, pheno$RNum], 
                 genes=rowData(rse_df), samples=pheno)
    # Filter by expression
    design0 <- model.matrix(~Dx, data=x$samples)
    keep.x <- filterByExpr(x, design=design0)
    x <- x[keep.x, , keep.lib.sizes=FALSE]
    print(paste('There are:', sum(keep.x), 'features left!', sep=' '))
    # Normalize library size
    x <- calcNormFactors(x, method="TMM")
    return(x)
}
memRPEP <- memoise::memoise(prep_data)

SVA_model <- function(seed_int){
    x <- memRPEP(seed_int)
    # Design matrix
    mod = model.matrix(~Dx + Age + mitoRate + rRNA_rate + RIN + 
                       totalAssignedGene + overallMapRate + ERCCsumLogErr + 
                       snpPC1 + snpPC2 + snpPC3, data=x$samples)
    colnames(mod) <- gsub("Dx", "", colnames(mod))
    colnames(mod) <- gsub("\\(Intercept\\)", "Intercept", colnames(mod))
    # Calculated SVs
    null.model = mod %>% as.data.frame %>% select(-c("Schizo")) %>% as.matrix
    n.sv <- num.sv(x$counts, mod, method="be")
    svobj <- svaseq(x$counts, mod, null.model, n.sv=n.sv)
    if(svobj$sv == 0){
        modQsva <- mod
    } else {
        modQsva <- cbind(mod, svobj$sv)
        len.d <- length(colnames(modQsva))
        colnames(modQsva)[((len.d - n.sv)+1):len.d] <- make.names(paste0("sv",1:n.sv))
    }
    return(modQsva)
}
memSVA <- memoise::memoise(SVA_model)

get_voom <- function(seed_int){
    ### Preform voom
    x <- memRPEP(seed_int)
    modQsva <- memSVA(seed_int)
    v <- voom(x[, rownames(modQsva)], modQsva)
    return(v)
}
memVOOM <- memoise::memoise(get_voom)

fit_voom <- function(seed_int){
    v <- memVOOM(seed_int)
    modQsva <- memSVA(seed_int)
    fit0 <- lmFit(v, modQsva)
    contr.matrix <- makeContrasts(CtrlvsSZ = Schizo,  
                                  levels=colnames(modQsva))
    fit <- contrasts.fit(fit0, contrasts=contr.matrix)
    esv <- eBayes(fit)
    return(esv)
}
memFIT <- memoise::memoise(fit_voom)

## MAIN
parser <- ArgumentParser()
parser$add_argument("-p", "--perm_num", type="integer",
                    help="Permutation number")
args <- parser$parse_args()

                                        # Preform voom
v <- memVOOM(args$perm_num)
                                        # Fit model and apply eBayes
efit = memFIT(args$perm_num)
                                        # Save differential expression
extract_de(1, "szVctl", efit, args$perm_num)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

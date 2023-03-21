## Differential expression with limma (permutation, subsampling)

suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(optparse)
    library(dplyr)
    library(limma)
    library(edgeR)
    library(here)
    library(sva)
})

## Functions

extract_de <- function(contrast, label, efit, seed_int){
    nperm   <- paste0("permutation_",
                      stringr::str_pad(seed_int, 4, pad = "0"))
    top     <- topTable(efit, coef=contrast, number=Inf, sort.by="P")
    top$SE  <- sqrt(efit$s2.post) * efit$stdev.unscaled
    top     <- top[order(top$P.Value), ]
    top.fdr <- top |> filter(adj.P.Val<=0.05)
    print(paste("Comparison for:", label))
    print(paste('There are:', dim(top.fdr)[1], 'DE features!'))
    data.table::fwrite(top, 
                       file=paste0(nperm, "_diffExpr_", label, "_full.txt"), 
                       sep='\t', row.names=TRUE)
}

load_counts <- function(){
    counts_file <- here("input/counts/_m",
                        "rse_gene.BrainSeq_Phase2_HIPPO.Gencode41.n447.Rdata")
    load(counts_file)
    return(rse_gene)
}
memCounts <- memoise::memoise(load_counts)

get_random_samples <- function(seed_int){
                                        # seed for reproducibility
    set.seed(seed_int + 113)
                                        # load count
    rse_df    <- memCounts()
                                        # subset samples
    keepIndex <- which((rse_df$Dx %in% c("Control", "SCZD")) & 
                       (rse_df$Age > 17) & (rse_df$Sex == "M") &
                       (rse_df$Race %in% c("AA", "CAUC")))
                                        # subsampling to female N (sample size)
    snames    <- sample(keepIndex, 121, replace=FALSE) 
    return(snames)
}
memSamples <- memoise::memoise(get_random_samples)

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

prep_data <- function(seed_int){
    rse_df     <- memCounts()
    keepIndex  <- memSamples(seed_int)
    rse_df     <- rse_df[, keepIndex]
    rse_df$Dx  <- factor(rse_df$Dx, levels = c("Control", "SCZD"))
    rse_df$Sex <- factor(rse_df$Sex)
    rownames(colData(rse_df)) <- sapply(strsplit(rownames(colData(rse_df)), "_"),
                                        "[", 1)
    pheno <- colData(rse_df) |> as.data.frame() |> 
        inner_join(memMDS(), by=c("BrNum"="FID"), multiple="all") |> 
        distinct(RNum, .keep_all = TRUE) |>
        left_join(memANTI(), by="BrNum")
    # Generate DGE list
    x <- DGEList(counts=assays(rse_df)$counts[, pheno$RNum], 
                 genes=rowData(rse_df), samples=pheno)
    # Filter by expression
    design0 <- model.matrix(~Dx, data=x$samples)
    keep.x  <- filterByExpr(x, design=design0)
    print(paste('There are:', sum(keep.x), 'features left!', sep=' '))
    x <- x[keep.x, , keep.lib.sizes=FALSE]
    # Normalize library size
    x <- calcNormFactors(x, method="TMM")
    return(x)
}
memRPEP <- memoise::memoise(prep_data)

SVA_model <- function(seed_int){
    x   <- memRPEP(seed_int)
    # Design matrix
    mod <- model.matrix(~Dx + Age + mitoRate + rRNA_rate + 
                            totalAssignedGene + RIN + Mapping_Rate + 
                            Mean3Bias + snpPC1 + snpPC2 + snpPC3,
                        data = x$samples)
    colnames(mod) <- gsub("Dx", "", colnames(mod))
    colnames(mod) <- gsub("\\(Intercept\\)", "Intercept", colnames(mod))
    # Calculated SVs
    null.model <- mod |> as.data.frame() |>
        select(-"SCZD") |> as.matrix()
    n.sv       <- num.sv(x$counts, mod, method="be")
    svobj <- svaseq(x$counts, mod, null.model, n.sv=n.sv)
    if(n.sv == 0){
        modQsva <- mod
    } else {
        modQsva <- cbind(mod, svobj$sv)
        len.d   <- length(colnames(modQsva))
        colnames(modQsva)[((len.d - n.sv)+1):len.d] <- make.names(paste0("sv",1:n.sv))
    }
    return(modQsva)
}
memSVA <- memoise::memoise(SVA_model)

get_voom <- function(seed_int){
    ### Preform voom
    x       <- memRPEP(seed_int)
    modQsva <- memSVA(seed_int)
    v       <- voom(x[, rownames(modQsva)], modQsva, plot=FALSE)
    return(v)
}
memVOOM <- memoise::memoise(get_voom)

fit_voom <- function(seed_int){
    v       <- memVOOM(seed_int)
    modQsva <- memSVA(seed_int)
    fit0    <- lmFit(v, modQsva)
    contr.matrix <- makeContrasts(CtrlvsSZ = SCZD,  
                                  levels=colnames(modQsva))
    fit <- contrasts.fit(fit0, contrasts=contr.matrix)
    esv <- eBayes(fit)
    return(esv)
}
memFIT <- memoise::memoise(fit_voom)

## MAIN
                                        # Input parser
option_list <- list(
    make_option(c("--perm_num"), type="integer",
                help="Permutation number",
                metavar="integer")
)
opt_parser  <- OptionParser(option_list=option_list)
opt         <- parse_args(opt_parser)

                                        # Preform voom
perm_num <- opt$perm_num
v        <- memVOOM(perm_num)

                                        # Fit model and apply eBayes
efit     <- memFIT(perm_num)

                                        # Save differential expression
extract_de(1, "szVctl", efit, perm_num)

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

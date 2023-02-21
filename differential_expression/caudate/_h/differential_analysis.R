## Differential expression analysis with limma-voom
## Sex-specific analysis for specific brain region
suppressPackageStartupMessages({
    library(SummarizedExperiment)
    library(optparse)
    library(dplyr)
    library(limma)
    library(edgeR)
    library(here)
    library(sva)
})

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
    top$SE  <- sqrt(efit$s2.post) * efit$stdev.unscaled
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

prep_data <- function(feature){
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
    keepIndex = which((rse_df$Dx %in% c("Control", "SCZD")) & 
                      rse_df$Age > 17 & rse_df$Race %in% c("AA", "CAUC"))
    rse_df = rse_df[, keepIndex]
    rse_df$Dx = factor(rse_df$Dx, levels = c("Control", "SCZD"))
    rse_df$Sex <- factor(rse_df$Sex)
    rownames(colData(rse_df)) <- sapply(strsplit(rownames(colData(rse_df)), "_"), "[", 1)
    pheno = colData(rse_df) |> as.data.frame() |> 
        inner_join(memMDS(), by=c("BrNum"="FID"), multiple="all") |> 
        distinct(RNum, .keep_all = TRUE) |>
        left_join(memANTI(), by="BrNum")
    # Generate DGE list
    x <- DGEList(counts=assays(rse_df)$counts[, pheno$RNum], 
                 genes=rowData(rse_df), samples=pheno)
    # Filter by expression
    design0 <- model.matrix(~Sex, data=x$samples)
    keep.x <- filterByExpr(x, design=design0)
    print(paste('There are:', sum(keep.x), 'features left!', sep=' '))
    x <- x[keep.x, , keep.lib.sizes=FALSE]
    # Normalize library size
    x <- calcNormFactors(x, method="TMM")
    return(x)
}
memPREP <- memoise::memoise(prep_data)

SVA_model <- function(feature){
    x <- memPREP(feature)
    # Design matrix
    mod = model.matrix(~Sex + Dx + Age + mitoRate + rRNA_rate + 
                       totalAssignedGene + RIN + Mapping_Rate + 
                       Mean3Bias + snpPC1 + snpPC2 + snpPC3, data = x$samples)
    colnames(mod) <- gsub("Dx", "", colnames(mod))
    colnames(mod) <- gsub("SexM", "Male", colnames(mod))
    colnames(mod) <- gsub("\\(Intercept\\)", "Intercept", colnames(mod))
    # Null model
    null.model = mod |> as.data.frame() |> select(-"Male") |> as.matrix()
    n.sv <- num.sv(x$counts, mod, method="be")
    ## Fit SVA
    svobj <- svaseq(x$counts, mod, null.model, n.sv=n.sv)
    ## Add to model
    print("")
    print(paste('Adding SV to design matrix ...', Sys.time()))
    modQsva <- cbind(mod, svobj$sv)
    len.d <- length(colnames(modQsva))
    colnames(modQsva)[((len.d - n.sv)+1):len.d] <- make.names(paste0("sv",1:n.sv))
    return(modQsva)
}
memSVA <- memoise::memoise(SVA_model)

get_voom <- function(feature){
    ### Preform voom
    x       <- memPREP(feature)
    modQsva <- memSVA(feature)
    v       <- voom(x[, rownames(modQsva)], modQsva, plot=FALSE)
    return(v)
}
memVOOM <- memoise::memoise(get_voom)

cal_res <- function(feature){
    ### Calculate residuals
    v <- memVOOM(feature)
    null_model <- v$design |> as.data.frame() |> select(-c("Male")) |> as.matrix()
    fit_res <- lmFit(v, design=null_model)
    res = v$E - ( fit_res$coefficients %*% t(null_model) )
    res_sd = apply(res, 1, sd)
    res_mean = apply(res, 1, mean)
    res_norm = (res - res_mean) / res_sd
    res_norm |> as.data.frame() |> tibble::rownames_to_column("feature_id") |>
    write.table(file=paste0(feature, '/residualized_expression.tsv'), 
                sep="\t", quote=FALSE, row.names=FALSE)
}
memRES <- memoise::memoise(cal_res)

fit_voom <- function(feature){
    v       <- memVOOM(feature)
    modQsva <- memSVA(feature)
    fit0    <- lmFit(v, modQsva)
    contr.matrix <- makeContrasts(MvsF = Male,  
                                  levels=colnames(modQsva))
    fit <- contrasts.fit(fit0, contrasts=contr.matrix)
    esv <- eBayes(fit)
    return(esv)
}
memFIT <- memoise::memoise(fit_voom)

#### MAIN
                                        # Input parser
option_list <- list(
    make_option(c("--feature"), type="character", default="genes",
                help="feature to analyze [default=%default]",
                metavar="character")
)
opt_parser  <- OptionParser(option_list=option_list)
opt         <- parse_args(opt_parser)
                                        # Differential expression by feature
feature <- opt$feature
dir.create(feature)
                                        # Preform voom
v <- memVOOM(feature)
save(v, file=paste0(feature,'/voomSVA.RData'))
                                        # Fit model and apply eBayes
efit = memFIT(feature)
                                        # Save differential expression
extract_de(1, "maleVfemale", efit, feature)
                                        # Calculate residuals
memRES(feature)

#### Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

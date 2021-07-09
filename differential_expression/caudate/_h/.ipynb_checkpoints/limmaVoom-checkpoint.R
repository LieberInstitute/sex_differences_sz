#!/usr/bin/env Rscript
                                        #@author = 'Kynon J Benjamin'
library(SummarizedExperiment)
library(data.table)
library(calibrate)
library(optparse)
library(limma)
library(edgeR)
library(sva)

arguments <- parse_args(
    OptionParser(
        usage = "%prog [options] counts_file directory_location",
        description="Run differential expression with limma voom for genes. Required inputs:\n <counts_file>: Features counts file PATH.\n <directory_loction>: PATH of output directory.",
        option_list=list(
            make_option(c("-f","--feature"),
                        default = "genes",
                        help = "Feature selection, needed for annotation subsetting [default %default]")
        )
    ),
    positional_arguments = 2
)

prepare_data <- function(counts.files, dir.loc, feature)
{
                                        # Load phenotype data
    print(paste('Starting processing phenotypes ...', Sys.time(), sep=' '))
    pheno.file <- '../../_m/phenotypes/sample_table.csv'
    sample.table <- read.csv(pheno.file, row.names=1)
    pheno0 <- subset(sample.table, Region=='Caudate',
                     select= c('BrNum', 'primarydx', 'sex', 'race',
                               'Protocol', 'RIN', 'Region', 'agedeath'))
    pheno0['RNum'] <- row.names(pheno0)
    print(paste('Finished with phenotype selection ...', Sys.time(),
                sep=' '))
                                        # Load PCA data
    print(paste('Load PCA from plink ..', Sys.time(), sep=' '))
    pca.file <- '../../_m/genotypes/hg38/pca/plink.eigenvec'
    pca0 <- read.delim(pca.file)
    pca <- merge(pca0[, 1:7], pheno0, by.x='FID', by.y='BrNum')
    rownames(pca) <- pca$RNum
    pheno <- pca[order(pca$RNum), ]
                                        # Load Mapping
    mito.file <- '../../_m/phenotypes/mitoMapping_rate.txt'
    mito <- read.delim(mito.file)
    map.file <- '../../_m/phenotypes/geneAlignment_rate.tsv'
    mapping <- read.delim(map.file)
    newMap0 <- merge(mito, mapping, by='SampleID')
    newMap <- newMap0[, c('SampleID', 'PercentRate', 'MitoRate')]
                                        # Genome Mapping
    tot.file <- '../../_m/phenotypes/total_reads.tsv'
    unmap.file <- '../../_m/phenotypes/unmapped_reads.tsv'
    tot.map <- read.delim(tot.file)
    unmap <- read.delim(unmap.file)
    mm <- merge(tot.map, unmap, by='RNum')
    mm['genomeMapping'] <- mm['unmapped_reads'] / mm['total_reads']
    allMap0 <- merge(newMap, mm, by.y='RNum', by.x='SampleID')
    allMap <- subset(allMap0, select=c(SampleID, PercentRate, MitoRate,
                                       genomeMapping))
                                        # Merge final
    new.pheno <- merge(pheno, allMap, by.x=0, by.y='SampleID')
    rownames(new.pheno) <- new.pheno$RNum
    new.pheno <- new.pheno[, -1]
    print(head(new.pheno, 5))
    print(paste('Finished with phenotype selection ...', Sys.time(),
                sep=' '))
                                        # Load counts data
    print(paste('Starting preparing data ...', Sys.time(), sep=' '))
    df.raw0 <- fread(counts.files, header=TRUE, verbose=TRUE)
    print(paste('Dimensions of counts object:', dim(df.raw0)[1],
                dim(df.raw0)[2], sep=' '))
    print(paste('Finished loading count data ...', Sys.time(), sep=' '))
                                        # Get Annotation
    if(feature == "genes"){
        annot <- df.raw0[, 1:6]
        df.raw <- df.raw0[, c(-1:-6)]
    } else if (feature == "exons"){
        annot <- df.raw0[, 1:7]
        df.raw <- df.raw0[, c(-1:-7)]
    } else if (feature == "transcripts"){
        annot <- df.raw0[, 1]
        df.raw <- df.raw0[, -1]
    } else if (feature == "junctions"){
        annot <- df.raw0[, 1:14]
        df.raw <- df.raw0[, c(-1:-14)]
    } else {
        stop("Please enter a correct feature!")
    }
    
    samp <- intersect(colnames(df.raw), rownames(new.pheno))
    new.pheno <- new.pheno[samp, ]
    df  <- df.raw[, ..samp]
    df[is.na(df)] <- 0
    print(dim(df))
                                        # Convert to DGEList
    print(paste('Converting to DGEList ...', Sys.time(), sep=' '))
    x <- DGEList(counts=df, genes=annot)
                                        # Add variables for group design
    x$samples$disorder <- factor(new.pheno$primarydx)
    x$samples$sex <- factor(new.pheno$sex)
    x$samples$race <- factor(new.pheno$race)
    x$samples$region <- factor(new.pheno$Region)
    x$samples$age <- new.pheno$age
    x$samples$brnum <- factor(new.pheno$FID)
    x$samples$rin <- new.pheno$RIN
    x$samples$rz <- factor(new.pheno$Protocol)
    x$samples$pc1 <- new.pheno$PC1
    x$samples$pc2 <- new.pheno$PC2
    x$samples$pc3 <- new.pheno$PC3
    x$samples$pc4 <- new.pheno$PC4
    x$samples$pc5 <- new.pheno$PC5
    x$samples$mito <- new.pheno$MitoRate
    x$samples$rate <- new.pheno$PercentRate
    x$samples$genome <- new.pheno$genomeMapping
    print(paste('Save raw counts ...', Sys.time(), sep=' '))
    save(x, file = paste0(dir.loc, 'rawCounts.RData'))
                                        # Filter low counts
    print(paste('Filtering low counts ...', Sys.time(), sep=' '))
    design0 <- model.matrix(~disorder, data=x$samples)
    keep.x <- filterByExpr(x, design=design0)
    x <- x[keep.x, , keep.lib.sizes=FALSE]
    print(paste('There are:', sum(keep.x), 'features left!', sep=' '))
                                        # Normalize for composition bias
    print(paste('Normalizing with TMM ...', Sys.time(), sep=' '))
    x <- calcNormFactors(x, method="TMM")
    print(paste('Saving normalized counts ...', Sys.time(), sep=' '))
    save(x, file = paste0(dir.loc, 'normCounts.RData'))
}

diff_expr <- function(dir.loc, feature)
{
    print(paste('Loading normalized counts ...', Sys.time(), sep=' '))
    load(paste0(dir.loc, 'normCounts.RData'))
    print(paste('Setting up model matrices ...', Sys.time(), sep=' '))
                                        # Full model
    design <- model.matrix(~disorder + sex + race + age + rin + mito + rate +
                               genome + pc1 + pc2 + pc3, data=x$samples)
                                        # Clean up column names
    colnames(design) <- gsub("disorder", "", colnames(design))
    colnames(design) <- gsub("sexM", "Male", colnames(design))
    colnames(design) <- gsub("race", "", colnames(design))
    colnames(design) <- gsub("\\(Intercept\\)", "Intercept",
                             colnames(design))
                                        # Calculate qSV
    print(paste('Starting qSVA ...', Sys.time(), sep=' '))
    print(paste('Loading degradation matrix ...', Sys.time(), sep=' '))
    load("../_h/degradation_rse_phase3_caudate.rda")

    print(paste('Calculating qSV components ...', Sys.time(), sep=' '))
    dm    <- assays(cov_rse_caudate)$counts
    qSV   <- qsva(dm)

    if("TRUE" %in% grepl("_", rownames(qSV))){# Remove underscore if needed
        rownames(qSV) <- sapply(strsplit(rownames(qSV), "_"), "[", 1)
    }
                                        # Add qSV to design matrix
    designQ  <- merge(design, qSV, by=0, all=FALSE)
    rownames(designQ) <- designQ$Row.names
    designQ  <- designQ[,-1]

    samp2 <- intersect(rownames(designQ), rownames(x$samples))
    x$samples <- x$samples[samp2, ]
    x$counts <- x$counts[ ,samp2]
                                        # Preform Initial Voom with SV model
    print(paste('Preforming inital voom ...', Sys.time(), sep=' '))
    v <- voom(x, designQ)
    print(paste('Finished voom ...', Sys.time(), sep=' '))
    save(v, file=paste0(dir.loc, 'voomSVA.RData'))
                                        # Calculate residuals
    null_model = v$design[, !(names(v$design) %in% c("Schizo"))]
    if(feature == 'junctions')
    {
        v$genes <- cbind(data.frame(ID=paste0('Junction_', seq(1:dim(v$genes)[1]))),
                         v$genes)
    }
    rownames(v$E) = v$genes[, 1]

    fit_res <- lmFit(v, design=null_model)
    res = v$E - ( fit_res$coefficients %*% t(null_model) )
    res_sd = apply(res, 1, sd)
    res_mean = apply(res, 1, mean)
					# Normalize residuals
    res_norm = (res - res_mean) / res_sd

    write.table(res_norm, file=paste0(dir.loc, 'residualized_expression.tsv'),
		sep="\t", quote=FALSE)
                                        # Fitting model with limma
    print(paste('Fitting model with limma ...', Sys.time(), sep=' '))
    fit0 <- lmFit(v, designQ)
                                        # Computing contrasts with qSV
    print(paste('Computing contrasts ...', Sys.time(), sep=' '))
    contr.matrix <- makeContrasts(
	CtrlvsSZ  = Schizo,
	levels    = colnames(designQ))
                                        # Fit contrast
    fit <- contrasts.fit(fit0, contrasts=contr.matrix)
    save(fit, file=paste0(dir.loc, 'fitSVA.RData'))
    save(designQ, contr.matrix, file=paste0(dir.loc, 'modelVariables.RData'))
                                        # Calculate diff genes with eBayes
    print(paste('Fitting statistical model eBayes ...', Sys.time(), sep=' '))
    esv <- eBayes(fit)

    options(width=200)
    print(paste('Getting statistics with topTable ...', Sys.time(), sep=' '))
    top0 <- topTable(esv, coef=1, number=Inf, sort.by="P")# SZ vs CTL
    sigTest <- decideTests(esv)
    top <- merge(top0, sigTest, by=0)
    rownames(top) <- top$Row.names
    top <- top[,-1]
    top <- top[order(top$P.Value), ]
                                        # Save differential expression data
    save(top, file=paste0(dir.loc, 'topTable_allData.RData'))
                                        # Write full differential expression
    print(paste('Writing full data to file ...', Sys.time(), sep=' '))
    write.table(top, file=paste0(dir.loc, "diffExpr_szVctl_full.txt"),
		sep='\t', row.names=FALSE, quote=FALSE)
                                        # Subset via FDR
    print(paste('Subsetting for FDR 0.05 ...', Sys.time(), sep=' '))
    top.fdr <- top[top$adj.P.Val<=0.05,]
    print(paste('There are:', dim(top.fdr)[1], 'DE features!', sep=' '))
    print(paste('Writing FDR DE features to file ...', Sys.time(), sep=' '))
    write.table(top.fdr, file=paste0(dir.loc, "diffExpr_szVctl_FDR05.txt"),
		sep='\t', row.names=FALSE, quote=FALSE)
    print(paste('Finished with differential expression ...', Sys.time(),
                sep=' '))
}

plot_graphs <- function(dir.loc)
{
                                        # Load Top Table Data
    load(paste0(dir.loc, 'topTable_allData.RData'))
                                        # Plot volcano plot
    pdf(file=paste0(dir.loc, "volcanoPlot.pdf"), 8, 6)
    with(top, plot(logFC, -log10(P.Value), pch=20, cex=0.6))
    with(subset(top, adj.P.Val<=0.05), points(logFC, -log10(P.Value),
                                              pch=20, col='red', cex=0.6))
    with(subset(top, abs(logFC)>0.50), points(logFC, -log10(P.Value),
                                              pch=20, col='orange', cex=0.6))
    with(subset(top, adj.P.Val<=0.05 & abs(logFC)>0.50),
         points(logFC, -log10(P.Value), pch=20, col='green', cex=0.6))
    dev.off()
                                        # Plot MA plot
    pdf(file=paste0(dir.loc, "MAplot.pdf"), 8, 6)
    with(top, plot(AveExpr, logFC, pch=20, cex=0.5))
    with(subset(top, adj.P.Val<0.05), 
         points(AveExpr, logFC, col="red", pch=20, cex=0.5))
    dev.off()
}

diffExpr_splice <- function(dir.loc)
{
    print(paste('Starting differential splicing analysis ...', Sys.time(), sep=' '))
    load(file=paste0(dir.loc, 'fitSVA.RData'))
    load(file=paste0(dir.loc, 'modelVariables.RData'))
    
    options(width=200)
    sp     <- diffSplice(fit, geneid="Geneid", exonid="Exonid")
    topSp <- topSplice(sp, coef=1, number=Inf)
                                        # Save differential splicing variables
    save(topSp, file=paste0(dir.loc, 'differentialSplice_ALL.RData'))
    print(paste('Writing full diffsplice data to file ...', Sys.time(), sep=' '))
    write.table(topSp, file=paste0(dir.loc, "diffSplice_szVctl_full.txt"),
		sep='\t', row.names=FALSE, quote=FALSE)
                                        # Subset via FDR
    print(paste('Subsetting for FDR 0.05 ...', Sys.time(), sep=' '))
    topSp.fdr <- topSp[topSp$FDR<=0.05,]
    print(paste('There are:', dim(topSp.fdr)[1], 'DS features!', sep=' '))
    print(paste('Writing FDR diffsplice to file ...', Sys.time(), sep=' '))
    write.table(topSp.fdr, file=paste0(dir.loc, "diffSplice_szVctl_FDR05.txt"),
                sep='\t', row.names=FALSE, quote=FALSE)
    print(paste('Finished with differential splicing analysis ...', Sys.time(),
                sep=' '))
}

opt = arguments$opt
counts_file = arguments$args[1]
dir_loc = arguments$args[2]

prepare_data(counts_file, dir_loc, opt$feature)
diff_expr(dir_loc, opt$feature)
plot_graphs(dir_loc)

if(opt$feature == "exons")
{
    diffExpr_splice(dir_loc)
}    

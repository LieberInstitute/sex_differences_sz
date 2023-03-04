## WGCNA analysis
filter_outliers <- function(expression, z_threshold = 2.5){
                                        # Input: an expression matrix
                                        # Output: an expression matrix with
                                        #         outliers removed
                                        # Remove samples with z normalized total
                                        # distance from other samples > z_threshold
    sample_distance = dist(expression)
    dist_z = scale(colSums(as.matrix(sample_distance)))
    stopifnot(all(rownames(dist_z) == rownames(expression)))
    keepSamples = dist_z < z_threshold
    new_expression = expression[keepSamples,]
    new_expression
}

prepare_data <- function(){
    suppressMessages(library(dplyr))
                                        # Load sample data
    load("../../../../differential_expression/hippocampus/_m/genes/voomSVA.RData")
    sample_table = v$design %>% as.data.frame %>% select(-Intercept) %>%
        rename("Sex"="Male") #%>% select(Ancestry, Sex, Age)
                                        # Load residualized expression
    fn = paste0("../../../../differential_expression/hippocampus/_m/",
                "genes/residualized_expression.tsv")
    vsd <- data.table::fread(fn) %>% replace(is.na(.), "") %>%
        tibble::column_to_rownames("V1")
    print(dim(vsd))
                                        # Keep only the columns and rows that
                                        # are present in both the sample table
                                        # and vsd file
    samples = intersect(colnames(vsd), rownames(sample_table))
    vsd = vsd[,samples]
    sample_table = sample_table[samples,]
                                        # WGCNA data import
    suppressMessages(library(WGCNA))
    options(stringsAsFactors = FALSE)
    datExpr0 = t(vsd)
                                        # Remove offending genes and samples
                                        # from the data
    gsg = goodSamplesGenes(datExpr0, verbose = 3);
    if (!gsg$allOK){
        datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
    }
    datExpr=datExpr0
                                        # Remove outliers
    datExpr = filter_outliers(datExpr0, z_threshold = 2.5)
    rm(datExpr0)
                                        # Clean data
    samples = intersect(rownames(datExpr), rownames(sample_table))
    sample_table = sample_table[samples,]
    datExpr = datExpr[samples,]
    print(dim(datExpr))
    save(datExpr, sample_table, file = '00.RData')
}

prepare_traits <- function(){
    lnames = load('00.RData')
                                        # Associate traits with samples
    traitRows = match(rownames(datExpr), rownames(sample_table))
    datTraits = sample_table[traitRows,]
                                        # Diagnostic plot: Sample dendrogram
                                        # and trait heatmap
    pdf(file='sample_dendrogram_and_trait_heatmap.pdf',height=22,width = 26)
    sampleTree2 = hclust(dist(datExpr), method = "average")
                                        # Convert traits to a color
                                        # representation: white means low, red
                                        # means high, grey means missing entry
    traitColors = numbers2colors(traitRows, signed=FALSE);
                                        # Plot the sample dendrogram and the
                                        # colors underneath.
    plotDendroAndColors(sampleTree2, traitColors, groupLabels="Avg. Counts",
                        main = "Sample dendrogram and trait heatmap",
                        cex.dendroLabels=0.7)
    dev.off()
    save(datExpr, sample_table, datTraits, file = "01.RData")
}

plot_power_parameter <- function(datExpr, plot_filename){
                                        # Choose a set of soft-thresholding powers
    powers = seq(from = 1, to=30, by=1)
                                        # Call the network topology analysis
                                        # function
    sft = pickSoftThreshold(datExpr, networkType = PARAM_NETWORK_TYPE,
                            powerVector = powers, verbose = 5)
                                        # Plot the results:
    pdf(file=plot_filename)
    par(mfcol = c(2,2)); par(mar = c(4.2, 4.5 , 2.2, 0.5),oma=c(0,0,2,0))
    cex1 = 0.7;
                                        # Scale-free topology fit index as a
                                        # function of the soft-thresholding power
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",
         ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="blue");
                                        # this line corresponds to using an R^2
                                        # cut-off of h
    abline(h=0.80,col="red")
                                        # Mean connectivity as a function of the
                                        # soft-thresholding power
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
         xlab="Soft Threshold (power)", ylab="Mean Connectivity",
         type="n", main = paste("Mean connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
         cex=cex1,col="blue")
    #####
    plot(sft$fitIndices[,1], sft$fitIndices[,6],
         xlab="Soft Threshold (power)", ylab="Median Connectivity",
         type="n", main = paste("Median connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,6], labels=powers, cex=cex1,
         col="blue")
    #####
    plot(sft$fitIndices[,1], sft$fitIndices[,7],
         xlab="Soft Threshold (power)", ylab="Max Connectivity",
         type="n", main = paste("Max connectivity"))
    text(sft$fitIndices[,1], sft$fitIndices[,7], labels=powers, cex=cex1,
         col="blue")
    dev.off()
}

figure_out_power_parameter <- function(){
    library(WGCNA)
    options(stringsAsFactors = FALSE);
    lnames = load(file = '01.RData')
    plot_power_parameter(datExpr, 'power_parameter_selection.pdf')
}

construct_network <- function(softpower){
  library(WGCNA)
  options(stringsAsFactors = FALSE);
  enableWGCNAThreads(nThreads=16)
  lnames = load(file = "01.RData")
                                        # softPower value from previous plot
                                        # power_parameter_selection.pdf
  softPower = softpower
                                        # ALWAYS choose a value equal or above
                                        # (better) 0.8
  cor <- WGCNA::cor
  net = blockwiseModules(datExpr, maxBlockSize=30000,
                         power = softPower, deepSplit = 2,
                         networkType = PARAM_NETWORK_TYPE,
                         TOMType = PARAM_NETWORK_TYPE,
                         numericLabels = TRUE, corType = "bicor",
                         saveTOMs = TRUE, mergeCutHeight = 0.15,
                         saveTOMFileBase = "TOM", minModuleSize = 50,
                         verbose = 3)
  moduleLabels = net$colors
  moduleColors = labels2colors(net$colors)
  MEs = net$MEs;
  geneTree = net$dendrograms[[1]];
  save(net, MEs, moduleLabels, moduleColors, geneTree, softPower,
       file="02.RData")
}

plot_cluster_dendrogram <- function(){
    library(WGCNA)
    options(stringsAsFactors = FALSE);
    enableWGCNAThreads(nThreads=16)
    load(file = "02.RData")
    pdf(file="cluster_dendrogram.pdf",height=16,width = 22)
    mergedColors = labels2colors(net$colors)
    plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                        "Module Colors", dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05, cex.dendroLabels=0.3)
    dev.off()
}

correlate_with_traits <- function(){
    library(WGCNA)
    options(stringsAsFactors = FALSE)
    enableWGCNAThreads(nThreads=16)
    lnames = load(file = "01.RData")
    lnames = load(file = "02.RData")
    # Define numbers of genes and samples
    nGenes = ncol(datExpr);
    nSamples = nrow(datExpr);
    # Recalculate MEs with color labels
    MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
    MEs = orderMEs(MEs0)
    moduleTraitCor = cor(MEs, datTraits, use = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
    # Plot
    pdf(file="module_trait_relationships.pdf", height=16,width = 22)
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "");
    dim(textMatrix) = dim(moduleTraitCor)
    par(mar = c(6, 8.5, 3, 3));
    # Display the correlation values within a heatmap plot
    labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(datTraits),
                   yLabels = names(MEs), ySymbols = names(MEs),
                   colorLabels = FALSE, naColor = "grey",
                   colors = blueWhiteRed(50), textMatrix = textMatrix,
                   setStdMargins = FALSE, cex.text = 0.9, zlim = c(-1,1),
                   main = paste("Module kME-Trait Correlation"))
    dev.off()
}

export_eigengene_tables <- function(){
    library(WGCNA)
    options(stringsAsFactors = FALSE)
    lnames = load(file = "01.RData")
    lnames = load(file = "02.RData")
                                        # Define numbers of genes and samples
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
                                        # Recalculate MEs with color labels
    MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
    rownames(MEs0) = rownames(datExpr)
    write.csv(MEs0, 'eigengenes.csv')
                                        # Write modules
    modules = data.frame(row.names=colnames(datExpr), module=moduleColors)
    write.csv(modules, 'modules.csv')
    save(datExpr,softPower,moduleColors, file = "cytoscapenetwork.Rdata")
}

## Run analysis
PARAM_NETWORK_TYPE = 'signed'
doParallel::registerDoParallel(cores=16)
prepare_data()
                                        # 1 - Sample dendrogram and trait heatmap
prepare_traits()
                                        # 2 - Scale Free Topology Model Fit
figure_out_power_parameter()
softPower = 8
construct_network(softPower)
                                        # 3 - TOM Dendogram
plot_cluster_dendrogram()
                                        # 4 - Module Eigenvalue Correlation with
                                        #     sample's traits
correlate_with_traits()
export_eigengene_tables()

## Reproducibility
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

## WGCNA analysis

construct_network <- function(softpower){
  library(WGCNA)
  options(stringsAsFactors = FALSE);
  enableWGCNAThreads(nThreads=4)
  lnames = load(file = "01.RData")
                                        # softPower value from previous plot
                                        # power_parameter_selection.pdf
  softPower = softpower; #Based on the hippocampus
                                        # ALWAYS choose a value equal or above
                                        # (better) 0.8
  cor <- WGCNA::cor
  net = blockwiseModules(datExpr, maxBlockSize=30000,
                         power = softPower, deepSplit = 3,
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
    enableWGCNAThreads(nThreads=4)
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
    enableWGCNAThreads(nThreads=4)
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
doParallel::registerDoParallel(cores=4)
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

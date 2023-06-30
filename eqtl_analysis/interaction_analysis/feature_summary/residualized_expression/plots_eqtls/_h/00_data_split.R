#### Script to subset expression
suppressMessages({
    library(magrittr)
})

#### MAIN analysis
                                        # Set variables for filtering data
feature <- "Gene"

                                        # Generate output directory
outdir  <- paste("out", sep="/")
dir.create(outdir, recursive=TRUE)

                                        # Extract interacting eQTL
infile <- "../../../_m/BrainSeq_sexGenotypes_4features_3regions.txt.gz"
eqtls0 <- read.table(gzfile(infile), comment.char="", header=TRUE, sep="\t")
cc <- dplyr::filter(eqtls0, feature_type == feature, region == "Caudate") %>%
    dplyr::arrange(gene_id, lfsr) %>% dplyr::group_by(gene_id) %>%
    dplyr::slice(1) %>% dplyr::arrange(lfsr) %>% as.data.frame()

dd <- dplyr::filter(eqtls0, feature_type == feature, region == "DLPFC") %>%
    dplyr::arrange(gene_id, lfsr) %>% dplyr::group_by(gene_id) %>%
    dplyr::slice(1) %>% dplyr::arrange(lfsr) %>% as.data.frame()

hh <- dplyr::filter(eqtls0, feature_type == feature, region == "Hippocampus") %>%
    dplyr::arrange(gene_id, lfsr) %>% dplyr::group_by(gene_id) %>%
    dplyr::slice(1) %>% dplyr::arrange(lfsr) %>% as.data.frame()

geneids <- setdiff(setdiff(cc$gene_id, dd$gene_id), hh$gene_id)
outfile <- paste0(outdir, "/eqtls_top25.tsv")
eqtls   <- dplyr::filter(eqtls0, feature_type == feature,
                         gene_id %in% geneids) %>%
    dplyr::arrange(gene_id, lfsr) %>% dplyr::group_by(gene_id) %>%
    dplyr::slice(1) %>% dplyr::arrange(lfsr) %>% as.data.frame() %>%
    head(25)

eqtls %>% data.table::fwrite(outfile, sep="\t")

                                        # Load residuals
infile  <- "../../_m/genes_residualized_expression.csv"
dat     <- read.table(gzfile(infile), comment.char="", header=TRUE, sep=",") %>%
    tidyr::separate(gene_id, c("gene_name", "gene_id"), "\\|") %>%
    dplyr::filter(gene_id %in% eqtls$gene_id)

                                        # Save data
write.table(dat, paste0(outdir,"/p"), col.names=TRUE,
            row.names=FALSE, quote=FALSE, sep="\t")

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

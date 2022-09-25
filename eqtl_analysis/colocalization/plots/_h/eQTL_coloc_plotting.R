suppressPackageStartupMessages({
    library(dplyr)
    library(eQTpLot)
    library(argparse)
})

## Create parser object
parser <- ArgumentParser()
parser$add_argument("-t", "--tissue", type="character",
                    help="tissue to run analysis on")
parser$add_argument("-f", "--feature", type="character",
                    help="feature name used in eQTL analysis")
parser$add_argument("-c", "--chrom", type="integer",
                    help="feature chromosome number")
parser$add_argument("-s", "--start", type="integer",
                    help="feature start site")
parser$add_argument("-e", "--end", type="integer",
                    help="feature end site")
parser$add_argument("-p", "--perm_pval", type="double", default=0.0001,
                    help="eQTL p-value threshold [default: %default]")
args <- parser$parse_args()

## Functions
get_gwas <- function(feature){
    fn <- paste0('gwas_pgc3_', feature, ".txt")
    return(data.table::fread(fn))
}

get_eqtl <- function(feature){
    fn <- paste0('eqtl_', feature, '.txt')
    return(data.table::fread(fn))
}

get_gene_df <- function(chrom, start, end, feature){
    return(data.frame('CHR'=chrom, 'Start'=start,
                      'Stop'=end, 'Gene'= feature,
                      'Build'= 'hg38'))
}

get_ld <- function(feature){
    fn <- paste0("shared_snps_", feature, ".ld")
    return(data.table::fread(fn, drop="V8"))
}

save_ggplots <- function(p, fn, w=6, h=6){
    for(ext in c('.svg', '.png', '.pdf')){
        ggsave(p, filename=paste0(fn, ext), width=w, height=h)
    }
}
## Main
                                        # Load data
gwas_df  <- get_gwas(args$feature)
eqtl_df  <- get_eqtl(args$feature)
genes_df <- get_gene_df(args$chrom, args$start,
                        args$end, args$feature)
ld_df    <- get_ld(args$feature)

                                        # eQTpLot
p = eQTpLot(GWAS.df=gwas_df, eQTL.df=eqtl_df, Genes.df=genes_df,
            getplot=FALSE, LD.df=ld_df, LDmin=10, R2min=0.25,
            LDcolor='black', gene=args$feature, trait='SCZD',
            gbuild='hg38', tissue=args$tissue,CollapseMethod="min",
            sigpvalue_eQTL=args$perm_pval,congruence=FALSE)

## Reproducibility information
Sys.time()
proc.time()
options(width=120)
sessioninfo::session_info()

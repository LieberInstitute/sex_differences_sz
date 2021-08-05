#!/usr/bin/R

library(dplyr)
library(argparse)

load_eqtl <- function(){
    ## Load eqtls, calculated FDR, convert to matrix eQTL format, and save
    eqtl_dt = data.table::fread("../../_m/Brainseq_LIBD.allpairs.txt.gz") %>%
        group_by(gene_id) %>% mutate(fdr=p.adjust(pval_nominal, method="fdr")) %>%
        mutate(statistic=slope / sqrt(slope_se/ma_samples + slope_se/ma_count)) %>%
        select(variant_id, gene_id, slope, statistic, pval_nominal, fdr) %>%
        rename("p-value"="pval_nominal")## %>%
        ## tidyr::separate(variant_id, c("chrom", "pos", "alleles"),
        ##                 sep=":", remove=FALSE, extra="merge")
    return(eqtl_dt)
}

memEQTL <- memoise::memoise(load_eqtl)

prep_genotype_matrix <- function(chrom){
    filename = paste0("/dcs04/lieber/ds2b/users/kynon/v4_phase3_paper/inputs/genotypes/",
                      "byChrom/a_transpose/_m/LIBD_Brain_TopMed.", chrom,
                      ".traw")
    geno_dt = data.table::fread(filename) %>% select(SNP, starts_with("Br"))
    colnames(geno_dt) <- gsub("_.*", "", colnames(geno_dt))
    ## Drop duplicated samples and save
    geno_dt %>% subset(., select = which(!duplicated(names(.)))) %>%
        data.table::fwrite(paste("genotypes", chrom, "txt.gz", sep='.'),
                           sep='\t')
}

extract_gen_positions <- function(chr_num){
    filename = paste0("/dcs04/lieber/ds2b/users/kynon/v4_phase3_paper/inputs/genotypes/",
                      "byChrom/a_transpose/_m/LIBD_Brain_TopMed.", chr_num,
                      ".traw")
    geno_dt = data.table::fread(filename) %>% select(SNP, CHR, POS)
    geno_dt %>%
        data.table::fwrite(paste("gen.position", chr_num, "txt", sep='.'),
                           sep='\t')
    return(geno_dt$SNP)
}

split_eqtl_by_chrom <- function(chr_num, snps){
    memEQTL() %>% filter(variant_id %in% snps) %>% as.data.frame %>%
        data.table::fwrite(paste("cis.eqtls", chr_num, "txt.gz", sep="."),
                           sep='\t')
}

## Create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option
parser$add_argument("-c", "--chrom", type="character",
                    help="Chromosome to run [required]")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

if(length(args$chrom) == 0){
    print("Missing chromosome!!")
} else {
    prep_genotype_matrix(args$chrom)
    snps = extract_gen_positions(args$chrom)
    split_eqtl_by_chrom(args$chrom, snps)
}

## ## Run file prep by chromosome
## for(chr_num in c(seq(1,22), "X")){
##     prep_genotype_matrix(chr_num)
##     snps = extract_gen_positions(chr_num)
##     split_eqtl_by_chrom(chr_num, snps)
## }

## Reproducibility Information
print(Sys.time())
print(proc.time())
options(width = 120)
print(sessioninfo::session_info())

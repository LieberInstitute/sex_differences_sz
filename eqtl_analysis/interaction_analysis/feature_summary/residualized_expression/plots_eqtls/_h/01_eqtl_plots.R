## Generating eQTL boxplots
suppressPackageStartupMessages({
    library(here)
    library(dplyr)
    library(ggpubr)
})

save_ggplots <- function(fn, p, w, h){
    for(ext in c('.pdf', '.png')){
        ggsave(paste0(fn, ext), plot=p, width=w, height=h)
    }
}

plot_eqtl <- function(variant_id, gene_id, gene_name, eqtl_annot){
    y0  <- min(big_dat[, gene_id]) - 0.2
    y1  <- max(big_dat[, gene_id]) + 0.2
    variant_name <- gsub(":", "\\.", variant_id)
    titlename <- paste(gene_name, gene_id, eqtl_annot, sep="\n")
    bxp <- big_dat %>% select(all_of(c(gene_id, variant_name, "Sex",
                                       "Region"))) %>%
        tidyr::drop_na() %>%
        ggboxplot(x=variant_name, y=gene_id, color="Sex", add="jitter",
                  ylab="Residualized Expression", xlab=variant_id,
                  add.params=list(alpha=0.5), alpha=0.4, legend="bottom",
                  palette="npg", ylim=c(y0,y1), facet.by="Region",
                  panel.labs.font=list(face="bold"), outlier.shape=NA,
                  ggtheme=theme_pubr(base_size=20, border=TRUE)) +
        font("xy.title", face="bold") + ggtitle(titlename) +
        theme(plot.title = element_text(hjust=0.5, face="bold"))
    return(bxp)
}

#### MAIN analysis
                                        # Run analysis
feature <- "genes"

                                        # Read residualized expression
outdir  <- paste("out", sep="/")
pfile   <- paste0(outdir, "/p")
p       <- read.table(pfile, comment.char="", header=TRUE)

                                        # gene meta
efile <- paste0(outdir, "/eqtls_top25.tsv")
e     <- read.table(efile, comment.char="", sep="\t", header=TRUE)

variants <- e[,"variant_id"]

vfile <- paste0(pfile, ".variants")
write.table(data.frame(variants), file=vfile, quote=FALSE, 
            col.names=FALSE, row.names=FALSE, sep=",")

                                        # get genotypes
plinkfile <- here("input/genotypes/subset_by_sex/shared_snps",
                  "_m/LIBD_Brain_TopMed")
outfile   <- paste0(pfile, "_temp")
command   <- paste("plink2 --bfile", plinkfile, "--extract", vfile,
                   "--export A", "--out", outfile)
system(command)

                                        # read sample phenotypes
sfile <- here('input/phenotypes/_m/phenotypes.csv')
s     <- read.table(sfile, comment.char="", header=TRUE, sep=",") %>%
    distinct(RNum, .keep_all=TRUE) %>%
    tibble::column_to_rownames("RNum")
snames <- select(s, BrNum)
                                        # select samples
gene_annot <- p %>% select(gene_id, gene_name) %>%
    tibble::column_to_rownames("gene_id")
p    <- p %>% tibble::column_to_rownames("gene_id") %>%
    select(-gene_name) %>% t()
idx  <- intersect(rownames(p), rownames(s))
p    <- p[is.element(rownames(p), idx),]
s    <- s[is.element(rownames(s), idx),]
                                        # keep samples with genotypes
pfam <- here("input/genotypes/subset_by_sex/shared_snps",
             "_m/LIBD_Brain_TopMed.fam")
ind  <- read.table(pfam, comment.char="", header=FALSE)
id   <- intersect(s$BrNum, ind[,1])
s    <- filter(s, BrNum %in% id)
idz  <- intersect(rownames(p), rownames(s))
p    <- p[is.element(rownames(p), idz),]

                                        # align samples phenotypes
idx  <- intersect(rownames(s), rownames(p))
p    <- p[match(idx, rownames(p)),]
s    <- s[match(idx, rownames(s)),]
geno <- read.table(paste0(outfile,".raw"), header=TRUE)
geno <- geno[match(id, geno$FID),] %>%
    select(-c("IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
names(geno) <- gsub("_[[:upper:]]+$", "", names(geno))

                                        # prep data
big_dat <- merge(s[,c("BrNum", "Age", "Sex", "Region", "Dx")], p,
                 by.x=0, by.y=0) %>%
    inner_join(geno, by=c("BrNum"="FID")) %>%
    distinct(Row.names, .keep_all=TRUE)

                                        # plot eQTL
plotdir <- paste("plots", sep="/")
dir.create(plotdir)

for(j in 1:dim(e)[1]){
    cat(j,"\n")        
    eqtl       <- e[j,]; gene_id <- eqtl$gene_id
    variant_id <- eqtl$variant_id; gene_name <- gene_annot[gene_id,]
    eqtl_annot <- paste("eQTL lfsr:", signif(eqtl$lfsr, 2))
                                        # Boxplot
    plotfile   <- paste0(plotdir, "/eqtl_", j)
    bxp        <- plot_eqtl(variant_id, gene_id, gene_name, eqtl_annot)
    save_ggplots(plotfile, bxp, 10, 6)
}

                                        # clean up
outfile <- paste0(pfile, "_temp.raw")
system(paste("rm", outfile, sep=" "))
outfile <- paste0(pfile, "_temp.log")
system(paste("rm", outfile, sep=" "))

#### Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

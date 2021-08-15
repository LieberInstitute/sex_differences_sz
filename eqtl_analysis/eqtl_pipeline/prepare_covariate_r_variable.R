## Gnerate covariates for eQTL and TWAS analysis

suppressMessages({library(SummarizedExperiment)
                  library(tidyverse)
                  library(jaffelab)
                  library(sva)})

## Load data
base_loc = "/ceph/projects/v4_phase3_paper/inputs/counts/_m"
                                        # Replace these files with correct tissue!
load(paste(base_loc,"caudate_brainseq_phase3_hg38_rseTx_merged_n464.rda"), sep='/')
load(paste(base_loc,"caudate_brainseq_phase3_hg38_rseJxn_merged_n464.rda"), sep='/')
load(paste(base_loc,"caudate_brainseq_phase3_hg38_rseExon_merged_n464.rda"), sep='/')
load(paste(base_loc,"caudate_brainseq_phase3_hg38_rseGene_merged_n464.rda"), sep='/')

## Filter based on gene expression
rse_gene = rse_gene[rowData(rse_gene)$meanExprs > 0.2,]
rse_exon = rse_exon[rowData(rse_exon)$meanExprs > 0.2,]
rowRanges(rse_jxn)$Length <- 100
jRp10m = recount::getRPKM(rse_jxn, 'Length')
rse_jxn = rse_jxn[rowMeans(jRp10m) > 0.4,]
rse_tx = rse_tx[rowMeans(assays(rse_tx)$tpm) > 0.4,]

## Sample selection
keepInd = which(colData(rse_gene)$Age > 13)
rse_gene = rse_gene[,keepInd]
rse_exon = rse_exon[,keepInd]
rse_jxn = rse_jxn[,keepInd]
rse_tx = rse_tx[,keepInd]

## Extract phenotypes and RPKMs
pd = colData(rse_gene)
geneRpkm = recount::getRPKM(rse_gene, "Length")
exonRpkm = recount::getRPKM(rse_exon, "Length")
jxnRp10m = recount::getRPKM(rse_jxn, 'Length')
txTpm = assays(rse_tx)$tpm

## Statistical model
### Add MDS to phenotyeps
mds = data.table::fread("../../_m/LIBD_Brain_TopMed.mds") %>%
    rename_at(.vars = vars(starts_with("C")),
              function(x){sub("C", "snpPC", x)}) %>%
    mutate_if(is.character, as.factor)

new_pd = pd %>% as.data.frame %>%
    inner_join(mds, by=c("BrNum"="FID")) %>%
    distinct(RNum, .keep_all = TRUE) %>% mutate(ids=RNum) %>%
    column_to_rownames("ids")

### Model
new_pd$Dx = factor(new_pd$Dx,levels = c("Control", "Schizo", "Bipolar"))
mod = model.matrix(~0+Dx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5,
                   data = new_pd)

## PCA
pcaGene = prcomp(t(log2(geneRpkm[, new_pd$RNum]+1)))
kGene = num.sv(log2(geneRpkm[, new_pd$RNum]+1), mod)
genePCs = pcaGene$x[,1:kGene]

pcaExon = prcomp(t(log2(exonRpkm[, new_pd$RNum]+1)))
kExon = num.sv(log2(exonRpkm[, new_pd$RNum]+1), mod, vfilter=50000)
exonPCs = pcaExon$x[,1:kExon]

pcaJxn = prcomp(t(log2(jxnRp10m[, new_pd$RNum]+1)))
kJxn = num.sv(log2(jxnRp10m[, new_pd$RNum]+1), mod, vfilter=50000)
jxnPCs = pcaJxn$x[,1:kJxn]

pcaTx = prcomp(t(log2(txTpm[, new_pd$RNum]+1)))
kTx = num.sv(log2(txTpm[, new_pd$RNum]+1), mod, vfilter=50000)
txPCs = pcaTx$x[,1:kTx]

save(genePCs, exonPCs, jxnPCs, txPCs,
     file="pcs_caudate_4features_filtered_over13.rda")

## Covariates
covsGene0 = t(cbind(mod[,-1],genePCs))
covsExon0 = t(cbind(mod[,-1],exonPCs))
covsJxn0 = t(cbind(mod[,-1],jxnPCs))
covsTx0 = t(cbind(mod[,-1],txPCs))

save(covsGene0, covsExon0, covsJxn0, covsTx0, file="covariates.rda")

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

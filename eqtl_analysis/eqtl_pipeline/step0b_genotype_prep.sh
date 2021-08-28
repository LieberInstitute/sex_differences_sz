#!/bin/sh
set -o errexit -o pipefail

BFILE="/PATH/TO/GENOTYPES"

Rscript ../_h/extract_shared_snps.R

plink2 --bfile $BFILE/LIBD_Brain_TopMed \
       --extract shared_snps.tsv \
       --make-bed --out LIBD_Brain_TopMed

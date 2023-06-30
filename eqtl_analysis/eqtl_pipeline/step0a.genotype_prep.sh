#!/bin/bash

## This script assumes PLINK2 is in your PATH

PHENO="/PATH/TO/PHENOTYPE/INFORMATION"
BFILE="/PATH/TO/GENOTYPE/BFILES"

echo "Starting ... $(date +'%Y/%m/%d %H:%M:%S')"

awk -F',' '$2=="Male"{print $18}' $PHENO | \
    grep -f - $BFILE/LIBD_Brain_TopMed.fam > keepMale.fam

plink2 --bfile $BFILE/LIBD_Brain_TopMed \
       --keep keepMale.fam --maf 0.05 \
       --mind 0.1 --geno 0.1 --hwe 1e-10 \
       --min-alleles 1 --make-pgen \
       --out LIBD_TOPMed_male

awk -F',' '$2=="Female"{print $18}' $PHENO | \
    grep -f - $BFILE/LIBD_Brain_TopMed.fam > keepFemale.fam

plink2 --bfile $BFILE/LIBD_Brain_TopMed \
       --keep keepFemale.fam --maf 0.05 \
       --mind 0.1 --geno 0.1 --hwe 1e-10 \
       --make-pgen --min-alleles 1 \
       --out LIBD_TOPMed_female

echo "Ending ... $(date +'%Y/%m/%d %H:%M:%S')"

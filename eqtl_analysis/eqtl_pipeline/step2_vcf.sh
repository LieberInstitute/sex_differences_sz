#!/bin/bash

## This assumes plink2 is in PATH.
## It is recommended to only run this once per tissue
## and think the VCF and indexed VCF file to other
## features (ln -s).

GENOTYPE_DIR='/PATH/TO/GENOTYPES'

plink2 --bfile ${GENOTYPE_DIR}/LIBD_Brain_TopMed \
       --keep ../../_m/keepFam.txt \
       --export vcf id-paste=fid vcf-dosage=DS \
       --out genotypes

perl -pe 's/^(##contig=<ID=)/$1chr/; s/^([^#])/chr$1/' genotypes.vcf \
     > genotypes_chr.vcf

bgzip genotypes_chr.vcf && tabix -p vcf genotypes_chr.vcf.gz

rm -f genotypes.vcf

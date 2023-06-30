#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=18.0G,h_vmem=20G,h_fsize=50G
#$ -N 'subset_genotypes'
#$ -o ./summary.out
#$ -e ./summary.out

GENOTYPE_DIR='../../../../../input/genotypes/subset_by_sex/shared_snps/_m'

module load plink/2.0
plink2 --bfile $GENOTYPE_DIR/LIBD_Brain_TopMed \
       --keep ../../_m/keepFam.txt \
       --make-bed --out genotypes

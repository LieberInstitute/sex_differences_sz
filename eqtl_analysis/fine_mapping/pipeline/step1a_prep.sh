#!/bin/bash
#$ -cwd -V
#$ -l h_fsize=50G
#$ -N cc_dapg_prep
#$ -e processing.log
#$ -o processing.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module list

echo "**** Run DAPG prep ****"
FEATURE="genes"; TISSUE="caudate"
VCF="../../../vcf/_m/genotypes_chr.vcf.gz"
EXPRESSION="../../../normalize_expression/_m/${FEATURE}.expression.bed.gz"
COVARIATES="../../../covariates/_m/${FEATURE}.combined_covariates.txt"

## Converting to sbams format
### Run processing script
perl ../_h/process.pl -e $EXPRESSION -g $VCF -c $COVARIATES -t $TISSUE

### Edit software location
sed -i 's|assemble.pl|../_h/assemble.pl|' ${TISSUE}.assemble.cmd

echo "**** Job ends ****"
date

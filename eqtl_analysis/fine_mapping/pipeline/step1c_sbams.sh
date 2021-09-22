#!/bin/bash

FEATURE="genes"
TISSUE="hippocampus"
VCF="../../../../vcf/_m/genotypes_chr.vcf.gz"
EXPRESSION="../../../_m/${FEATURE}.expression.bed.gz"
COVARIATES="../../../../covariates/_m/${FEATURE}.combined_covariates.txt"

## Converting to sbams format
### Run processing script
perl /ceph/opt/dap/gtex_v8_analysis/process.pl \
     -e $EXPRESSION -g $VCF -c $COVARIATES -t $TISSUE
### Edit software location
sed -i 's|assemble.pl|../_h/assemble.pl|' ${TISSUE}.assemble.cmd
### Execute batch assemble to get SBAMS
cat ${TISSUE}.assemble.cmd | parallel --jobs 32 {}

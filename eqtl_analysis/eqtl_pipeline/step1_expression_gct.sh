#!/bin/bash

FEATURE="gene" ## change based on feature to process
TISSUE="caudate" ## brain region to process
FAM="/ceph/projects/v4_phase3_paper/inputs/genotypes/_m/LIBD_Brain_TopMed.fam"
PHENO="/ceph/projects/v4_phase3_paper/inputs/phenotypes/_m/merged_phenotypes.csv"
TPM="/ceph/projects/v4_phase3_paper/inputs/counts/text_files_counts/tpm/_m/${TISSUE}/${FEATURE}/tpm.csv"

## Match counts with feature short hand
if [[ $FEATURE == 'transcript' ]]; then
    VAR="tx"
elif [[ $FEATURE == "junction" ]]; then
    VAR="jxn"
else
    VAR=$FEATURE
fi
COUNTS="/ceph/projects/v4_phase3_paper/inputs/counts/text_files_counts/_m/${TISSUE}/${VAR}_counts.txt"

python ../_h/prepare_gct.py --fam_file $FAM --pheno_file $PHENO \
       --tpm_file $TPM --counts_file $COUNTS --tissue $TISSUE

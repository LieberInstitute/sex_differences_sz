#!/bin/bash

FEATURE="gene" ## change based on feature to process
FAM="/location/of/FAM/genotypes.fam"
PHENO="/location/of/phenotype/file.csv"
TPM="/location/of/TPM/file.csv"

## Match counts with feature short hand
if [[ $FEATURE == 'transcript' ]]; then
    VAR="tx"
elif [[ $FEATURE == "junction" ]]; then
    VAR="jxn"
else
    VAR=$FEATURE
fi
COUNTS="/location/of/counts/file/${VAR}_counts.txt"

python3 ../_h/01_prepare_gct.py --fam_file $FAM --pheno_file $PHENO \
       --tpm_file $TPM --counts_file $COUNTS --tissue $TISSUE

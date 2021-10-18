#!/bin/bash

FEATURE="genes"
FASTQTL_DIR="/ceph/opt/enhanced_fastqtl"

python ../_h/run_FastQTL_threaded.py \
       --covariates ../../../covariates/_m/${FEATURE}.combined_covariates.txt \
       --permute 1000 10000 --ma_sample_threshold 10 --threads 32 --chunks 250 \
       --maf_threshold 0.01 --window 0.5e6 --fastqtl_dir $FASTQTL_DIR \
       ../../../vcf/_m/genotypes_chr.vcf.gz \
       ../../_m/${FEATURE}.expression.bed.gz \
       Brainseq_LIBD

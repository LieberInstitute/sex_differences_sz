#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N gene_covariates
#$ -o ./summary.out
#$ -e ./summary.out

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load conda_R/4.2.x
module load gcc/9.1.0
module load pandoc

module list

echo "**** Run covariates prep ****"
FEATURE="genes"; TISSUE="caudate"

Rscript ../_h/03_generate_covs.R --feature $FEATURE --region $TISSUE

echo "**** Job ends ****"
date

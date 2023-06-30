#!/bin/bash
#$ -cwd
#$ -l mem_free=2G,h_vmem=2G,h_fsize=10G
#$ -N eQTLplots_caudate
#$ -o ./logs/summary.log
#$ -e ./logs/summary.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load gcc/9.1.0
module load plink/2.0
module load pandoc
module load R

module list

## Edit with your job command
TISSUE="caudate"

echo "**** Run combine files ****"
Rscript ../_h/01_eqtl_plots.R --tissue $TISSUE

echo "**** Job ends ****"
date

#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N sz_coloc_comp
#$ -o ./summary.log
#$ -e ./summary.log

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

echo "**** Run comparison ****"
Rscript ../_h/coloc_sz.R

echo "**** Job ends ****"
date

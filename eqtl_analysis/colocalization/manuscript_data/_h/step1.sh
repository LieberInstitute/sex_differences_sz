#!/bin/bash
#$ -cwd
#$ -R y
#$ -pe local 5
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N coloc_supp
#$ -o ./summary.log
#$ -e ./summary.log
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load R/4.0.3
module list

echo "**** Run generate supplement ****"
Rscript ../_h/generate_supp.R

echo "**** Job ends ****"
date

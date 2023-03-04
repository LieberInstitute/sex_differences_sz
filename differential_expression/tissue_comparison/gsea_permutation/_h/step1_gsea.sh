#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=20.0G,h_vmem=20G,h_fsize=50G
#$ -N gsea_analysis
#$ -o ./summary.out
#$ -e ./summary.out
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load conda_R
module list

## Edit with your job command
echo "**** Run GOseq ****"
Rscript ../_h/perm_gsea.R

echo "**** Job ends ****"
date

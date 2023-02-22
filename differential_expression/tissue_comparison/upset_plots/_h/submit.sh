#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N upset_plot
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
module load gcc/9.1.0
module load R/4.0.3
module load pandoc

module list

## Edit with your job command
echo "**** Run upset plot ****"
Rscript ../_h/upsetplot.R
Rscript ../_h/upsetplot_cmc.R

echo "**** Job ends ****"
date

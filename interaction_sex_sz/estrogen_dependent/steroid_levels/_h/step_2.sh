#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N estrogen_dependent
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
module load conda_R/4.2.x
module load pandoc

module list

## Edit with your job command
echo "**** Explore sex hormone expression ****"

Rscript ../_h/plot_data.R

echo "**** Job ends ****"
date

#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N gsea_similarity
#$ -o ./similarity.log
#$ -e ./similarity.log

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

## Edit with your job command
echo "**** Run GO semantic similarity ****"
Rscript ../_h/go_similarity.R

echo "**** Job ends ****"
date

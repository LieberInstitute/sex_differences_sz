#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=25G,h_vmem=25G,h_fsize=100G
#$ -N extract_model
#$ -o ./mash.log
#$ -e ./mash.log

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

## Job command
echo "**** Run mashr prep ****"

Rscript ../_h/03_extract_mash.R

echo "**** Job ends ****"
date

#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=15.0G,h_vmem=15G,h_fsize=50G
#$ -N 'softpower_caudate'
#$ -o ./summary.out
#$ -e ./summary.out

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

#### Load modules
module load gcc/9.1.0
module load conda_R/4.2.x
module load pandoc

module list

echo "**** Run network analysis ****"
Rscript  ../_h/01_estimate_softpower.R

echo "**** Job ends ****"
date -Is

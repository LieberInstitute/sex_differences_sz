#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=2.0G,h_vmem=20G,h_fsize=50G
#$ -N 'resid_interaction_eqtl'
#$ -o ./summary.out
#$ -e ./summary.out

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
module load gcc/9.1.0
module load pandoc

module list

echo "**** Run residualization ****"
Rscript ../_h/generate_res.R

echo "**** Job ends ****"
date -Is

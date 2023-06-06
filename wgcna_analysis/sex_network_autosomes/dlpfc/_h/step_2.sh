#!/bin/bash
#$ -cwd
#$ -R y
#$ -pe local 15
#$ -l mem_free=8G,h_vmem=8G,h_fsize=100G
#$ -N 'network_dlpfc'
#$ -o ./network.out
#$ -e ./network.out

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

echo "**** Run network analysis ****"
Rscript ../_h/02_coexpression_network.R 

echo "**** Job ends ****"
date -Is


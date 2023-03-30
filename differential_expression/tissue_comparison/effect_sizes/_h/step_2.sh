#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=10.0G,h_vmem=10G,h_fsize=50G
#$ -N 'correlation'
#$ -o ./summary.out
#$ -e ./summary.out

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

#### List modules
module load conda_R/4.2.x
module load gcc/9.1.0
module load pandoc

module list

echo "**** Run summarize results ****"
mkdir sig/
python  ../_h/corr_N_sig.py

echo "**** Job ends ****"
date -Is

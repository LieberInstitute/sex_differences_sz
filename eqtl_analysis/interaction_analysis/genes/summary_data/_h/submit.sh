#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=100.0G,h_vmem=100G,h_fsize=50G
#$ -N 'summarize_eqtl'
#$ -o ./summary.out
#$ -e ./summary.out

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module list

echo "**** Run summarize results ****"
python3  ../_h/feature_summarize.py

echo "**** Job ends ****"
date -Is

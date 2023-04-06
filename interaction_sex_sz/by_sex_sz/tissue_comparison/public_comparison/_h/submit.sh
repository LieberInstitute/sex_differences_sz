#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=5G,h_vmem=5G,h_fsize=50G
#$ -N 'overlap_qin'
#$ -o ./summary.out
#$ -e ./summary.out

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load pandoc
module list

echo "**** Run summarize results ****"
python3  ../_h/generate_summary.py

echo "**** Job ends ****"
date -Is

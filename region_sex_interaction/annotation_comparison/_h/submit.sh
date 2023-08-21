#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=10.0G,h_vmem=10G,h_fsize=50G
#$ -N 'correlation'
#$ -o ./summary.log
#$ -e ./summary.log

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

#### List modules
module load pandoc

module list

echo "**** Run summarize results ****"
python3 ../_h/compare_versions.py

echo "**** Job ends ****"
date -Is

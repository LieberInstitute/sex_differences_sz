#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=20G,h_vmem=20G,h_fsize=50G
#$ -N 'enrichment_fishers'
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

echo "**** Run enrichment analysis ****"
python3  ../_h/de_enrichment_fishers.py

echo "**** Job ends ****"
date -Is

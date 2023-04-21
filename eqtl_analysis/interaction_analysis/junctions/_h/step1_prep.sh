#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=200G,h_vmem=200G,h_fsize=50G
#$ -N prepare_mashr
#$ -o ./summary.out
#$ -e ./summary.out

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module list

## Job command
FEATURE="junctions"

echo "**** Run mashr prep ****"
python3 ../_h/01_prepare_files.py --feature $FEATURE

echo "**** Job ends ****"
date

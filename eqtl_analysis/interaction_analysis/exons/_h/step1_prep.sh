#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=90G,h_vmem=90G,h_fsize=50G
#$ -N prepare_mashr
#$ -o ./summary.log
#$ -e ./summary.log

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
FEATURE="exons"

echo "**** Run mashr prep ****"
python3 ../_h/01_prepare_files.py --feature $FEATURE

echo "**** Job ends ****"
date

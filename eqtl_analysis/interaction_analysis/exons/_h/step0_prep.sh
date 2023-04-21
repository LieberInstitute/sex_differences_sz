#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N extract_eGenes
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
FEATURE="exons"

echo "**** Run mashr prep ****"
python3 ../_h/00_extract_top_eqtls.py --feature $FEATURE

echo "**** Job ends ****"
date

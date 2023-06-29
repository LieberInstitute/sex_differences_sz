#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N combine_data
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
module load pandoc

module list

## Edit with your job command
echo "**** Generate combined adata ****"
python3 ../_h/combine_data.py

echo "**** Job ends ****"
date

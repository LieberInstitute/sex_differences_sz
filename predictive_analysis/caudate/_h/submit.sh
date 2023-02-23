#!/bin/bash
#$ -cwd
#$ -l mem_free=10.0G,h_vmem=10G,h_fsize=50G
#$ -N dRFE_caudate
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
module load gcc/9.1.0
module load R/4.0.3
module load pandoc
module load python

module list

echo "**** Run dRFEtools for predictive analysis ****"
Rscript ../_h/differential_analysis.R --feature $FEATURE

echo "**** Job ends ****"
date

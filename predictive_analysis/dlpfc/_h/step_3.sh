#!/bin/bash
#$ -cwd
#$ -l mem_free=30.0G,h_vmem=30G,h_fsize=50G
#$ -N dRFE_exon_dlpfc
#$ -o ./summary_exon.log
#$ -e ./summary_exon.log

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
FEATURE="exons"; TISSUE="dlpfc"

python ../_h/dRFE_rf.py --feature $FEATURE --tissue $TISSUE

echo "**** Job ends ****"
date

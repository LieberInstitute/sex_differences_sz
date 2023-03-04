#!/bin/bash
#$ -cwd
#$ -l mem_free=20.0G,h_vmem=20G,h_fsize=50G
#$ -N dRFE_tx_hippocampus
#$ -o ./summary_tx.log
#$ -e ./summary_tx.log

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
FEATURE="transcripts"; TISSUE="hippocampus"

python ../_h/dRFE_rf.py --feature $FEATURE --tissue $TISSUE

echo "**** Job ends ****"
date

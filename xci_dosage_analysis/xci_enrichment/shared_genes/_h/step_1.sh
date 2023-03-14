#!/bin/bash
#$ -cwd
#$ -l mem_free=5.0G,h_vmem=5G,h_fsize=50G
#$ -N xci_enrichment
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
module load pandoc

module list

echo "**** Run XCI enrichment analysis ****"
python3 ../_h/01_examine_sharing.py

echo "**** Job ends ****"
date

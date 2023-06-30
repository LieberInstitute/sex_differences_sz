#!/bin/bash
#$ -cwd
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N combine_parquet
#$ -o ./combine.log
#$ -e ./combine.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load samtools
module load htslib
module list

## Edit with your job command
echo "**** Combine parquet files ****"
python3 ../_h/02_combine_parquet.py
bgzip -f BrainSEQ_TOPMed.interaction.txt
tabix -f BrainSEQ_TOPMed.interaction.txt.gz

echo "**** Job ends ****"
date

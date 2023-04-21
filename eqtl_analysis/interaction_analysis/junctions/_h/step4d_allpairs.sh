#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N mash_allpairs_clean
#$ -o ./clean_mash.log
#$ -e ./clean_mash.log

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
echo "**** Run combine mash results ****"

gzip posterior_mean.sex_interaction.txt
gzip lfsr.sex_interaction.txt

echo "**** Job ends ****"
date

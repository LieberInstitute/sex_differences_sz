#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N mash_allpairs_chunking
#$ -o ./chunking.log
#$ -e ./chunking.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load conda_R/4.2.x
module load gcc/9.1.0
module load pandoc

module list

## Job command
echo "**** Run mashr prep ****"
mkdir output
Rscript ../_h/04_generate_chunks.R --chunk_size 1250 --output output

echo "**** Job ends ****"
date

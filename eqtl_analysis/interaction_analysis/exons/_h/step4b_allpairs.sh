#!/bin/bash
#$ -cwd
#$ -R y
#$ -t 1-1250:1
#$ -tc 50
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N mash_allpairs_array
#$ -o ./output/allpairs_$TASK_ID.log
#$ -e ./output/allpairs_$TASK_ID.log

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

Rscript ../_h/05_all_association_mash.R \
	--chunk_num $SGE_TASK_ID --output output

echo "**** Job ends ****"
date

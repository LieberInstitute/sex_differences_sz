#!/bin/bash
#$ -cwd
#$ -R y
#$ -t 1-1000:1
#$ -tc 50
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N subsampling_sexBYsz
#$ -o ./logs/summary_$TASK_ID.log
#$ -e ./logs/summary_$TASK_ID.log

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
echo "**** Run DE subsampling permutation analysis ****"

Rscript ../_h/permutate_samples.R --perm_num $SGE_TASK_ID

echo "**** Job ends ****"
date

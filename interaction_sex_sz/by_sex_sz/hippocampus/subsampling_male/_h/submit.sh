#!/bin/bash
#$ -cwd
#$ -R y
#$ -t 1-1000:1
#$ -tc 50
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N subsampling_sexBYsz
#$ -o ./summary.log
#$ -e ./summary.log
#$ -m e -M jade.benjamin@libd.org

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load R
module list

## Job command
echo "**** Run DE subsampling permutation analysis ****"

Rscript ../_h/permutate_samples.R --perm_num $SGE_TASK_ID

echo "**** Job ends ****"
date

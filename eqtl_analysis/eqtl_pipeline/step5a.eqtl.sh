#!/bin/bash
#$ -cwd
#$ -l caracol,mem_free=50G,h_vmem=50G,h_fsize=100G
#$ -N gene_tensorqtl
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
module load python
module load R

module list

## Edit with your job command
export CUDA_VISIBLE_DEVICES=0

echo "**** Run tensorQTL ****"
FEATURE="genes"
python ../_h/01_eqtl_tensorqtl.py --feature $FEATURE

echo "**** Job ends ****"
date

#!/bin/bash
#$ -cwd
#$ -R y
#$ -t 1-2:1
#$ -l mem_free=5G,h_vmem=5G,h_fsize=100G
#$ -N mash_allpairs_combine
#$ -o ./allpairs_mash.log
#$ -e ./allpairs_mash.log

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
LABELS=("lfsr" "posterior_mean")
python3 ../_h/06_combine_results.py \
        --output output --label ${LABELS[$SGE_TASK_ID-1]}

echo "**** Job ends ****"
date

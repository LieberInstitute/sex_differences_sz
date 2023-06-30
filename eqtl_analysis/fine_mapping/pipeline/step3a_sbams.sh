#!/bin/bash
#$ -cwd -V
#$ -t 1-25979:1
#$ -N gen_sbams
#$ -e ./logs/sbams_$TASK_ID.out
#$ -o ./logs/sbams_$TASK_ID.out
#$ -tc 75

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load gcc/9.1.0
module load bcftools/1.9

module list

echo "**** Generate sbams ****"
TISSUE="caudate"
export OPENBLAS_NUM_THREADS=1

python3 ../_h/parallel_run_sbams.py \
	--task $SGE_TASK_ID --tissue $TISSUE

echo "**** Job ends ****"
date

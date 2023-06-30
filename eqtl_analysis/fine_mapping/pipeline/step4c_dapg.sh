#!/bin/bash
#$ -cwd -V
#$ -N hh_dapg
#$ -l mem_free=5G,h_vmem=5G,h_fsize=50G
#$ -e ./logs/dapg_$TASK_ID.out
#$ -o ./logs/dapg_$TASK_ID.out
#$ -t 1-22307:1
#$ -tc 50

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load gsl
module load gcc/9.1.0

module list

echo "**** Fine mapping with dap-g ****"
TISSUE="hippocampus"; PRIOR=(priors/*.prior)
export OPENBLAS_NUM_THREADS=2

## Running DAP-G
GENEID=`basename ${PRIOR[${SGE_TASK_ID}]} .prior`
# ID=`echo $GENEID | sed 's/\./_/'`
/dcs04/lieber/ds2b/opt/dap/dap_src/dap-g \
    -d ${TISSUE}/${GENEID}.sbams.dat \
    -p priors/${GENEID}.prior \
    -ld_control 0.5 --all -t 2 > dapg_out/${GENEID}.out

echo "**** Job ends ****"
date

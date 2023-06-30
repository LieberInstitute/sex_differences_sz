#!/bin/bash
#$ -cwd -V
#$ -l mem_free=30G,h_vmem=30G,h_fsize=50G
#$ -N torus_priors
#$ -e torus.out
#$ -o torus.out

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module list

echo "**** Run torus ****"
FEATURE="genes"; TISSUE="dlpfc"
EQTL="../../_m/BrainSEQ_TOPMed.interaction.txt.gz"

/dcs04/lieber/ds2b/opt/torus/src/torus.static \
    -d $EQTL --fastqtl -dump_prior priors

echo "**** Job ends ****"
date

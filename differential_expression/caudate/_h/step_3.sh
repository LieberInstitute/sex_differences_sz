#!/bin/bash
#$ -cwd
#$ -l mem_free=35.0G,h_vmem=35G,h_fsize=50G
#$ -N de_sex_exon_caudate
#$ -o ./summary_exon.log
#$ -e ./summary_exon.log

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
module load tex

module list

echo "**** Run exonral correlation ****"
FEATURE="exons"

Rscript ../_h/differential_analysis.R --feature $FEATURE

echo "**** Job ends ****"
date

#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N coloc_supp
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
module list

echo "**** Run generate supplement ****"
mkdir eQTpLot
for TISSUE in caudate dlpfc hippocampus; do
    mkdir eQTpLot/${TISSUE}
    for SEX in female male; do
	mkdir eQTpLot/${TISSUE}/${SEX}
	cp -v ../../plots/_m/${TISSUE}/${SEX}/*png \
	   eQTpLot/${TISSUE}/${SEX}/
    done
done

## Save file
tar -czvf BrainSeq_colocalization_eQTpLoTs.tar.gz eQTpLot

## Clean directory
rm -rf eQTpLot

echo "**** Job ends ****"
date

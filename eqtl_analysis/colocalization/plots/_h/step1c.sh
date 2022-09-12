#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=10G,h_vmem=10G,h_fsize=100G
#$ -N coloc_plots_hippocampus
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
module load plink/1.90b6.6
module load R/4.0.3
module list


echo "**** Run eQTpLot ****"
TISSUE="hippocampus"
BED="/dcs04/lieber/statsgen/jbenjami/projects/sex_differences_sz/input/counts/text_files_counts/_m/caudate"
PGC3_COLOC="../../../../prep_eqtl_analysis/${TISSUE}/genes/prepare_expression/fastqtl_nominal/dapg_finemap/fastENLOC/_m/pgc3.enloc.sig.out"

mkdir $TISSUE

for FEATURE in `awk '$6 > 0.5 {print $1}' $PGC3_COLOC | cut -f1 -d: | grep ENSG`; do
    grep $FEATURE $BED/gene.bed | while read -r PARAM; do
        GNAME=`echo $PARAM | awk '{print $1}'`
        CHR=`echo $PARAM | awk '{print $2}' | sed 's/chr//' -`
        START=`echo $PARAM | awk '{print $3}'`
        END=`echo $PARAM | awk '{print $4}'`
        echo $TISSUE $GNAME $CHR $START $END
        for SEX in female male; do
            mkdir $TISSUE/$SEX
            cd $TISSUE/$SEX
            python ../../../_h/prep_colocalization.py \
                   --chrom $CHR --start $START --end $END \
                   --feature $GNAME --tissue $TISSUE --sex $SEX \
                   --perm_pval 0.0001 ## This is too many eQTL
            cd ../../
        done
    done
done

echo "**** Job ends ****"
date

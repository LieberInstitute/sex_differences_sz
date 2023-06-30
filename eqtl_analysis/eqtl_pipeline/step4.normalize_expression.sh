#!/bin/bash
#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N norm_expression
#$ -o ./summary.out
#$ -e ./summary.out

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

## List current modules for reproducibility
module load htslib
module load samtools
module load python/3.9.10

module list

echo "**** Run normalized expression ****"
SOFTWARE="../_h"
BED="../../../../../input/counts/text_files_counts/_m/caudate/gene.bed"

ln -sfn ../_h/rnaseqnorm.py .

python3 $SOFTWARE/04_eqtl_prepare_expression.py \
        --feature gene --bed_file $BED -o . \
	../../_m/tpm.gct ../../_m/counts.gct \
        ../../_m/sample_id_to_brnum.tsv \
	../../_m/vcf_chr_list.txt genes

echo "**** Job ends ****"
date

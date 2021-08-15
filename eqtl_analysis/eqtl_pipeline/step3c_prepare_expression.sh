#!/bin/bash

SOFTWARE="../_h"
BED="/ceph/projects/v3_phase3_paper/inputs/phase3/text_files_counts/_m"
GTF="/ceph/genome/human/gencode25/gtf.CHR/_m/gencode.v25.annotation.gtf"

ln -sfn ../_h/rnaseqnorm.py .

python $SOFTWARE/eqtl_prepare_expression.py \
       --feature exon \
       --bed_file $BED/exon.bed \
       -o . ../../_m/tpm.gct \
       ../../_m/counts.gct \
       $GTF ../../_m/sample_id_to_brnum.tsv \
       ../../_m/vcf_chr_list.txt exons

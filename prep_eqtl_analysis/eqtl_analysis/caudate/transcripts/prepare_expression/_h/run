#!/bin/bash

SOFTWARE="../_h"
GTF="/ceph/genome/human/gencode25/gtf.CHR/_m/gencode.v25.annotation.gtf"

ln -sfn ../_h/rnaseqnorm.py .

python $SOFTWARE/eqtl_prepare_expression.py \
       --feature transcript -o . ../../_m/tpm.gct \
       ../../_m/counts.gct \
       $GTF ../../_m/sample_id_to_brnum.tsv \
       ../../_m/vcf_chr_list.txt transcripts

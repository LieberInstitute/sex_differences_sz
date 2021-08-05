#!/bin/bash

SOFTWARE="/dcl02/lieber/apaquola/opt/eigenMT"

qsub -V -N genepos -o log.out -e log.out -cwd \
     -b y "module load R; Rscript ../_h/prep_genepos.R"

for CHR in `seq 1 22` X; do
    qsub -V -l mem_free=18.0G,h_vmem=20G -N eigen_${CHR}_dt -o log_${CHR}.out -e log_${CHR}.out -cwd \
	 -b y "module load R; \
	       Rscript ../_h/prep_eigenMT.R --chrom ${CHR}; \
               python $SOFTWARE/eigenMT.py \
               	      --CHROM ${CHR} --cis_dist 0.5e6 \
                      --QTL cis.eqtls.${CHR}.txt.gz \
                      --GEN genotypes.${CHR}.txt.gz \
                      --GENPOS gen.position.${CHR}.txt \
                      --PHEPOS phe.position.txt \
                      --OUT Brainseq_LIBD.${CHR}.tsv"
done

# rm *gz *txt

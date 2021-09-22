#!/bin/bash

GWAS="/ceph/projects/v4_phase3_paper/inputs/sz_gwas/pgc2_clozuk/map_phase3/zscore/tabix/fine_mapped_loci/torus_pip"
SOFTWARE="/ceph/opt/fastenloc/src"
VCF="../../../../../vcf/_m"

## Prepare probabilistic eQTL annotation
$SOFTWARE/summarize_dap2enloc.pl \
    -dir ../../_m/dapg_out/ -vcf $VCF/genotypes_chr.vcf.gz \
    -tissue hippocampus | gzip > fastenloc.eqtl.annotation.vcf.gz

## Run fastENLOC
$SOFTWARE/fastenloc.static -eqtl fastenloc.eqtl.annotation.vcf.gz \
                           -gwas $GWAS/_m/snps_in_gwas_loci.pip.gz \
                           -total_variants 7341682 -t hippocampus -thread 32 \
                           -prefix fastenloc > fastenloc.log 2> fastenloc.stderr

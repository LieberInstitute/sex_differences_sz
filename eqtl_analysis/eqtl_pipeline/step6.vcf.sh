#!/bin/bash
#$ -cwd
#$ -R y
#$ -l mem_free=18.0G,h_vmem=20G,h_fsize=50G
#$ -N 'subset_genotypes'
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
module load plink/2.0

module list

echo "**** Generate VCF ****"
GENOTYPE_DIR='../../../../../input/genotypes/subset_by_sex/shared_snps/_m'

plink2 --bfile ${GENOTYPE_DIR}/LIBD_Brain_TopMed \
       --keep ../../_m/keepFam.txt \
       --export vcf id-paste=fid vcf-dosage=DS \
       --out genotypes

perl -pe 's/^(##contig=<ID=)/$1chr/; s/^([^#])/chr$1/' genotypes.vcf \
     > genotypes_chr.vcf

bgzip genotypes_chr.vcf && tabix -p vcf genotypes_chr.vcf.gz

rm -f genotypes.vcf

echo "**** Job ends ****"
date

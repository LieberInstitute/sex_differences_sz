#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N 'auto_go_enrichment'
#$ -o ./summary.out
#$ -e ./summary.out

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module list

echo "**** Run GO plotting ****"
python3  ../_h/01_gene_ontology.py

echo "**** Job ends ****"
date -Is

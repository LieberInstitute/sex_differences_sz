#$ -cwd
#$ -R y
#$ -l mem_free=50G,h_vmem=50G,h_fsize=50G
#$ -N 'pi1_analysis'
#$ -o ./summary.out
#$ -e ./summary.out

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load conda_R/4.2.x
module load gcc/9.1.0
module load pandoc

module list

echo "**** Run pi1 estimation ****"
Rscript  ../_h/pi1_estimate.R

echo "**** Job ends ****"
date -Is

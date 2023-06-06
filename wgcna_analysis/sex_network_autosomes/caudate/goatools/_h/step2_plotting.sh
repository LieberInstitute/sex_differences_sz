#$ -cwd
#$ -R y
#$ -l h_fsize=50G
#$ -N 'go_plotting'
#$ -o ./plotting.out
#$ -e ./plotting.out

echo "**** Job starts ****"
date -Is

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"

module load R/3.6.1
module list

echo "**** Run GO plotting ****"
python3  ../_h/02_plot_go.py

echo "**** Job ends ****"
date -Is

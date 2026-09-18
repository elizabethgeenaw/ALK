#!/bin/bash
#SBATCH --job-name=figs_ctspec_bar
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/figs_ctspec_bar_%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/figs_ctspec_bar_%j.err

cd /project/gilad/geena/publication_clean/code/treatment/differential_expression/dreamlet
mkdir -p /project/gilad/geena/publication_clean/code/treatment/differential_expression/logs

module load gsl/2.7
module load R/4.4.1
export LD_LIBRARY_PATH=/software/hdf5-1.12.0-el8-x86_64/lib:/software/gsl-2.7-el8-x86_64/lib:$LD_LIBRARY_PATH

echo "[$(date)] Starting figures_ctspecific_bar.R on $(hostname)"
Rscript figures_ctspecific_bar.R
EXIT_CODE=$?
echo "[$(date)] figures_ctspecific_bar.R finished with exit code ${EXIT_CODE}"
exit $EXIT_CODE

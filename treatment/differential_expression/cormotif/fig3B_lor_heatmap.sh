#!/bin/bash
#SBATCH --job-name=fig3B_lor_hm
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/fig3B_lor_hm_%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/fig3B_lor_hm_%j.err

cd /project/gilad/geena/publication_clean/code/treatment/differential_expression/cormotif
mkdir -p /project/gilad/geena/publication_clean/code/treatment/differential_expression/logs

module load gsl/2.7
module load R/4.4.1
export LD_LIBRARY_PATH=/software/hdf5-1.12.0-el8-x86_64/lib:/software/gsl-2.7-el8-x86_64/lib:$LD_LIBRARY_PATH

echo "[$(date)] Starting fig3B_lor_motif_gene_heatmap.R on $(hostname)"
Rscript fig3B_lor_motif_gene_heatmap.R
EXIT_CODE=$?
echo "[$(date)] fig3B_lor_motif_gene_heatmap.R finished with exit code ${EXIT_CODE}"
exit $EXIT_CODE

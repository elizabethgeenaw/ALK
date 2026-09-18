#!/bin/bash
#SBATCH --job-name=mat_int
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/matrixeqtl/slurm_output/interaction_chr_%A_%a.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/matrixeqtl/slurm_output/interaction_chr_%A_%a.err
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --array=1-22
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G
#SBATCH --time=36:00:00

module load R/4.4.1

# UMI_LABEL can be passed as an environment variable:
#   sbatch --export=UMI_LABEL=medianUMI10 4_matrixeqtl_interaction.sh
# Default is medianUMI1.
UMI_LABEL=${UMI_LABEL:-medianUMI1}

CHR_RANGE="${SLURM_ARRAY_TASK_ID}-${SLURM_ARRAY_TASK_ID}"
echo "Processing chromosome: $CHR_RANGE  UMI label: $UMI_LABEL"

mkdir -p /project/gilad/geena/publication_clean/code/treatment/matrixeqtl/slurm_output
mkdir -p /project/gilad/geena/publication_clean/results/treatment/matrixeqtl/interaction/tmp

Rscript /project/gilad/geena/ALK_paper_code/treatment/matrixeqtl/04_matrixeqtl_interaction.R \
  --chr_range "$CHR_RANGE" \
  --umi_label "$UMI_LABEL" \
  --nc 1

#! /bin/bash
#SBATCH --job-name=t_mat  # Job name
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/matrixeqtl/slurm_output/medianUMI1_eqtl_chr_%A_%a.out  # Standard output log
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/matrixeqtl/slurm_output/medianUMI1_eqtl_chr_%A_%a.err   # Error log
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --array=1-22                # Array job with 8 tasks
#SBATCH --ntasks=1                 # Run a single task
#SBATCH --ntasks-per-node=1
#SBATCH --mem=30G                  # Memory per job
#SBATCH --time=36:00:00            # Time limit hrs:min:sec

module load R/4.4.1


# Define chromosome based on the array task ID
CHR_RANGE="${SLURM_ARRAY_TASK_ID}-${SLURM_ARRAY_TASK_ID}"

# Print the chromosome being processed
echo "Processing chromosome: $CHR_RANGE"

# Run R script with the chromosome range argument
Rscript /project/gilad/geena/ALK_paper_code/treatment/matrixeqtl/02_matrixeqtl_treatment.R --chr_range "$CHR_RANGE" --nc 1
#!/bin/bash
#SBATCH --job-name=counts_power_clean
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=1:00:00
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/counts_power_clean-%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/counts_power_clean-%j.err

module load gsl/2.7
module load R/4.4.1
export LD_LIBRARY_PATH=/software/hdf5-1.12.0-el8-x86_64/lib:/software/gsl-2.7-el8-x86_64/lib:$LD_LIBRARY_PATH

Rscript /project/gilad/geena/ALK_paper_code/treatment/differential_expression/dreamlet/20_supp_S3_counts_power.R

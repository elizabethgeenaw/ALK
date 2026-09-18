#!/bin/bash
#SBATCH --job-name=supp_S14
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=2:00:00
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/matrixeqtl/logs/supp_S14-%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/matrixeqtl/logs/supp_S14-%j.err

module load R/4.4.1
Rscript /project/gilad/geena/ALK_paper_code/treatment/matrixeqtl/supp_S14_interaction_concordance.R

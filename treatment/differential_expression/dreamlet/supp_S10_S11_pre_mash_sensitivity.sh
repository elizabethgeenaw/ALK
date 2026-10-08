#!/bin/bash
#SBATCH --job-name=pre_mash_sensitivity
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/pre_mash_sensitivity-%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/pre_mash_sensitivity-%j.err

module load R/4.4.1

Rscript /project/gilad/geena/ALK_paper_code/treatment/differential_expression/dreamlet/supp_S10_S11_pre_mash_sensitivity.R

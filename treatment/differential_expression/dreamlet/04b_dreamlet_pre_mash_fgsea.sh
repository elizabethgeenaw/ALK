#!/bin/bash
#SBATCH --job-name=dreamlet_pre_mash_fgsea
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/dreamlet_pre_mash_fgsea-%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/dreamlet_pre_mash_fgsea-%j.err

module load R/4.4.1

Rscript /project/gilad/geena/ALK_paper_code/treatment/differential_expression/dreamlet/04b_dreamlet_pre_mash_fgsea.R

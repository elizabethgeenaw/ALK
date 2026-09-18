#!/bin/bash
#SBATCH --job-name=fig2C_hallmark
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=30:00
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/fig2C_hallmark-%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/fig2C_hallmark-%j.err

module load R/4.4.1
Rscript /project/gilad/geena/ALK_paper_code/treatment/differential_expression/dreamlet/fig2C_hallmark_heatmap.R

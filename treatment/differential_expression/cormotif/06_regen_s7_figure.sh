#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=regen_s7_figure
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/regen_s7_figure-%j.out

module load R/4.2.0

Rscript /project/gilad/geena/ALK_paper_code/treatment/differential_expression/cormotif/06_regen_s7_figure.R

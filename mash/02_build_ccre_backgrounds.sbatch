#!/bin/bash
#SBATCH --job-name=build_ccre_bg
#SBATCH --output=/project/gilad/geena/publication_clean/code/mash/logs/build_ccre_backgrounds_%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/mash/logs/build_ccre_backgrounds_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad-hm


module load R/4.4.1

Rscript /project/gilad/geena/ALK_paper_code/mash/02_build_ccre_backgrounds.R

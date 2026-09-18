#!/bin/bash
#SBATCH --job-name=supp_S5fig_S6
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=1:00:00
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/supp_S5fig_S6-%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/supp_S5fig_S6-%j.err

module load R/4.4.1
Rscript /project/gilad/geena/ALK_paper_code/treatment/differential_expression/dreamlet/supp_S5fig_S6_gsea_ct_patterns.R

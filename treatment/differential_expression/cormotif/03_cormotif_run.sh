#!/bin/bash
#SBATCH --time=80:00:00
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200G
#SBATCH --job-name=cormotif_new_filters
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/cormotif_new_filters-%j.out


module load R/4.2.0

Rscript /project/gilad/geena/ALK_paper_code/treatment/differential_expression/cormotif/03_cormotif_run_S7fig_S10_S14tables.R

#!/bin/bash
#SBATCH --job-name=drug_upset
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad-hm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/drug_specific_upset-%j.out

module load R/4.4.1

Rscript /project/gilad/geena/ALK_paper_code/treatment/differential_expression/cormotif/20_supp_S8_drug_specific_upset.R

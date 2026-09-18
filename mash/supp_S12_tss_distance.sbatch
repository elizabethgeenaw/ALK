#!/bin/bash
#SBATCH --job-name=tss_multipage
#SBATCH --output=/project/gilad/geena/publication_clean/code/mash/logs/tss_multipage_%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/mash/logs/tss_multipage_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad-hm

module load gsl/2.7
module load R/4.4.1
export R_LIBS_USER=/project/gilad/geena/Rlibs/R-4.4.1

# Usage: sbatch plot_tss_distance_multipage.sbatch [lfc]
# lfc defaults to 3 (paper's response-eQTL threshold).
LFC=${1:-3}

Rscript /project/gilad/geena/ALK_paper_code/mash/13_supp_S12_tss_distance.R "$LFC"

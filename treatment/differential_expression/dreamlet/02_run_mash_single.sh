#!/bin/bash
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --mem=32G
#SBATCH --time=04:00:00

# Runs one per-drug mash fit
#
# This script is meant to be called by 02_run_mash_launcher.sh, which handles
# submitting the 5 per-drug jobs in parallel. You can also sbatch it directly
# for a one-off run on a single drug.
#
# Most relevant code in 02_run_mash.R; this file just sets up the
# environment and passes the drug name straight through.
#
# Usage (direct):
#   sbatch 02_run_mash_single.sh --drug <Drug>

module load R/4.4.1
module load gsl/2.6

cd /project/gilad/geena/publication_clean/code/treatment/differential_expression/dreamlet
# OPENBLAS_MAIN_FREE=1 stops OpenBLAS from forking off worker threads that
# collide with mash's own multiprocessing on shared nodes.
export OPENBLAS_MAIN_FREE=1

if [ "$1" != "--drug" ] || [ -z "$2" ]; then
  echo "Usage: sbatch 02_run_mash_single.sh --drug <Drug>"
  exit 1
fi

Rscript 02_run_mash.R --drug "$2"

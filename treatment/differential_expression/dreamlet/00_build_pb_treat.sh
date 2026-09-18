#!/bin/bash
#SBATCH --job-name=build_pb_treat
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad-hm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200G
#SBATCH --time=12:00:00
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/build_pb_treat-%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/differential_expression/logs/build_pb_treat-%j.err

module load gsl/2.7
module load R/4.4.1
export LD_LIBRARY_PATH=/software/hdf5-1.12.0-el8-x86_64/lib:/software/gsl-2.7-el8-x86_64/lib:$LD_LIBRARY_PATH

Rscript /project/gilad/geena/ALK_paper_code/treatment/differential_expression/dreamlet/00_build_pb_treat.R

#!/bin/bash
#SBATCH --job-name=sketch
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/preprocessing/logs/sketch_umap_%j.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/preprocessing/logs/sketch_umap_%j.err
#SBATCH --time=80:00:00
#SBATCH --account=pi-gilad
#SBATCH --partition=gilad
#SBATCH --nodelist=midway-l16b-31
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=40G

module load R/4.2.0
module load geos/3.7.0
module load hdf5_hl
module load hdf5
module load cmake

Rscript /project/gilad/geena/ALK_paper_code/treatment/preprocessing/05_fig1B_sketch_umap.R

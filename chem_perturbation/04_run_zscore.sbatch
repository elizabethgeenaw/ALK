#!/bin/bash
#SBATCH --job-name=zscore
#SBATCH --partition=caslake
#SBATCH --account=pi-gilad
#SBATCH --time=00:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=/project/gilad/geena/chem_perturbation_clean/logs/zscore_%j.out
#SBATCH --error=/project/gilad/geena/chem_perturbation_clean/logs/zscore_%j.err

source /software/modules/init/bash
module load python/miniforge-25.3.0
source activate /project/gilad/geena/envs/crisp_env

python /project/gilad/geena/ALK_paper_code/chem_perturbation/04_zscore.py

#!/bin/bash
#SBATCH --job-name=crisp_rg
#SBATCH --output=/project/gilad/geena/chem_perturbation_clean/logs/step1_2_radial_glia.out
#SBATCH --error=/project/gilad/geena/chem_perturbation_clean/logs/step1_2_radial_glia.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --account=pi-gilad

source /software/modules/init/bash
module load python/miniforge-25.3.0
source activate /project/gilad/geena/envs/crisp_env

cd /project/gilad/geena/chem_perturbation_clean

export SINGLE_CT=Radial_glia
export CRISP_BATCH_SIZE=100000

echo "Step 1: Preprocess + scGPT embed"
python 1_prep_input.py

echo "Step 2: CRISP predictions"
python 2_run_crisp.py

echo "Done: Radial_glia"

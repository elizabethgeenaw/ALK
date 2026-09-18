#!/bin/bash
#SBATCH --job-name=crisp_highmem
#SBATCH --output=/project/gilad/geena/chem_perturbation_clean/logs/step1_2_highmem_%A_%a.out
#SBATCH --error=/project/gilad/geena/chem_perturbation_clean/logs/step1_2_highmem_%A_%a.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --account=pi-gilad
#SBATCH --array=0-3%2

source /software/modules/init/bash
module load python/miniforge-25.3.0
source activate /project/gilad/geena/envs/crisp_env

cd /project/gilad/geena/chem_perturbation_clean

# Cell types that OOM at 60G — run with 120G and smaller batch size
CTS=(
    Neuronal_IPC
    Radial_glia
    Stromal_cells
    iPSCs
)

CT=${CTS[$SLURM_ARRAY_TASK_ID]}
echo "Cell type: $CT (task $SLURM_ARRAY_TASK_ID)"

export SINGLE_CT=$CT
export CRISP_BATCH_SIZE=30000

echo "Step 1: Preprocess + scGPT embed"
python 1_prep_input.py

echo "Step 2: CRISP predictions"
python 2_run_crisp.py

echo "Done: $CT"

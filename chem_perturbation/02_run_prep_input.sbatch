#!/bin/bash
#SBATCH --job-name=crisp_full
#SBATCH --output=/project/gilad/geena/chem_perturbation_clean/logs/step1_2_%A_%a.out
#SBATCH --error=/project/gilad/geena/chem_perturbation_clean/logs/step1_2_%A_%a.err
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH --account=pi-gilad
#SBATCH --array=0-16%2

source /software/modules/init/bash
module load python/miniforge-25.3.0
source activate /project/gilad/geena/envs/crisp_env

cd /project/gilad/geena/chem_perturbation_clean

# High-memory cell types (Neuronal_IPC, Radial_glia, Stromal_cells, iPSCs) are in run_step1_2_highmem.sbatch
CTS=(
    Bronchiolar_and_alveolar_epithelial_cells
    Cardiac_mesoderm_cells
    Ciliated_epithelial_cells
    Early_mesoderm_cells
    ENS_glia
    Epicardial_fat_cells
    Glioblast
    Goblet_cells
    Hepatoblasts
    Metanephric_cells
    Neuroblast
    Neuron
    Retinal_neurons
    Retinal_pigment_cells
    Retinal_progenitors_and_Muller_glia
    Squamous_epithelial_cells
    Vascular_endothelial_cells
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

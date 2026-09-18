#!/bin/bash
#SBATCH --job-name=fix_normstats
#SBATCH --output=/project/gilad/geena/chem_perturbation_clean/logs/fix_normstats.out
#SBATCH --error=/project/gilad/geena/chem_perturbation_clean/logs/fix_normstats.err
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --partition=caslake
#SBATCH --account=pi-gilad

source /software/modules/init/bash
module load python/miniforge-25.3.0
source activate /project/gilad/geena/envs/crisp_env

cd /project/gilad/geena/chem_perturbation_clean

echo "=== Step 1: Regenerate normstats from raw descriptors ==="
python make_rdkit2d_normstats.py

echo "=== Step 2: Validate smile_df against reference parquet ==="
python validate_smile_df.py

#!/bin/bash
#SBATCH --job-name=crisp_val
#SBATCH --account=pi-gilad
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=/project/gilad/geena/chem_perturbation_clean/validation/logs/validation.out
#SBATCH --error=/project/gilad/geena/chem_perturbation_clean/validation/logs/validation.err

PYTHON=/project/gilad/geena/envs/crisp_env/bin/python
VALDIR=/project/gilad/geena/chem_perturbation_clean/validation

echo "=== 1. Positive control ==="
$PYTHON $VALDIR/1_positive_control.py

echo ""
echo "=== 2. ALK biology ==="
$PYTHON $VALDIR/2_alk_biology.py

echo ""
echo "=== 3. Drug structure-response ==="
$PYTHON $VALDIR/3_drug_structure_response.py

echo ""
echo "=== 4. Cross-cell-type consistency ==="
$PYTHON $VALDIR/4_crosscelltype_consistency.py

echo ""
echo "=== All validations complete ==="

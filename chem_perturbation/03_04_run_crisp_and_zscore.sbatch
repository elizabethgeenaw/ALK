#!/bin/bash
#SBATCH --job-name=zscore_fgsea
#SBATCH --output=/project/gilad/geena/chem_perturbation_clean/logs/step3_4_%j.out
#SBATCH --error=/project/gilad/geena/chem_perturbation_clean/logs/step3_4_%j.err
#SBATCH --partition=caslake
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --account=pi-gilad

source /software/modules/init/bash
module load python/miniforge-25.3.0
source activate /project/gilad/geena/envs/crisp_env
module load R/4.4.1

cd /project/gilad/geena/chem_perturbation_clean

echo "Z-scores + GSEA"
python 3_zscore_gsea.py

echo "Convert parquet to TSV for FGSEA"
ZSCORE_DIR="results/predictions/zscore_analysis"
python -c "import pandas as pd; pd.read_parquet('${ZSCORE_DIR}/gene_zscores_all.parquet').to_csv('${ZSCORE_DIR}/gene_zscores_all.tsv', sep='\t', index=False)"

echo "FGSEA"
Rscript 4_fgsea.R


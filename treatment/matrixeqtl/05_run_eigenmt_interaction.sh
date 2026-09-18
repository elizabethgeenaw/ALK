#!/bin/bash
#SBATCH --job-name=eigenmt_interact
#SBATCH --output=/project/gilad/geena/publication_clean/code/treatment/matrixeqtl/slurm_output/eigenmt_interact_%A_%a.out
#SBATCH --error=/project/gilad/geena/publication_clean/code/treatment/matrixeqtl/slurm_output/eigenmt_interact_%A_%a.err
#SBATCH --array=1-22
#SBATCH --time=30:00:00
#SBATCH --account=pi-gilad
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G

# eigenMT local multiple-testing correction for the G×T interaction model results.
# Run after 4_matrixeqtl_interaction.R finishes.
# UMI_LABEL controls which result set to correct (default: medianUMI1).
# To run on medianUMI10: UMI_LABEL=medianUMI10 sbatch 9b_run_eigenmt_interaction.sh

module load python

UMI_LABEL=${UMI_LABEL:-medianUMI1}

BASE_DIR="/project/gilad/geena/publication_clean/results/treatment/matrixeqtl/interaction/${UMI_LABEL}"
GENO_DIR="/project/gilad/geena/ALK_analysis/matrixeqtl_interaction_input/geno/pseudoreplicate"
PHEPOS="/project/gilad/geena/ALK_analysis/matrixeqtl_treatment_input/genopos.txt"
EIGENMT="/project/gilad/geena/ALK_analysis/eigenMT/eigenMT.py"

mkdir -p "${BASE_DIR}/eigenMT"

chr=chr${SLURM_ARRAY_TASK_ID}

declare -a arr=()
while IFS= read -r line; do
    arr+=("$line")
done < <(ls "${BASE_DIR}" | grep "_interact_cis\.txt$" | sed -E 's/_chr[0-9]+_interact_cis\.txt//' | sort | uniq)

for treatment_celltype in "${arr[@]}"; do
    python "${EIGENMT}" --CHROM ${chr} \
        --QTL "${BASE_DIR}/${treatment_celltype}_${chr}_interact_cis.txt" \
        --GEN "${GENO_DIR}/${chr}_snp.txt" \
        --GENPOS "${GENO_DIR}/${chr}_snpsloc.txt" \
        --PHEPOS "${PHEPOS}" \
        --OUT "${BASE_DIR}/eigenMT/${treatment_celltype}_${chr}_interact_cis.txt"
done

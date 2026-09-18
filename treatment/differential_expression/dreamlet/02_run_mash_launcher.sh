#!/bin/bash
# Launcher for the 5 per-drug dreamlet::run_mash fits on the treatment 1v4 data.
#
# Each per-drug mash fit shrinks that drug's 30-cell-type coefficient vector
# using a mixture-of-covariance prior estimated from the drug's own
# coefficient matrix. Each is the input to every downstream mash-based
# figure and supplementary table (Fig 2C, 2D, Supp Figs S5/S6, Supp Tables
# S6/S7/S8/S9). Submits all 5 as parallel sbatch jobs (~30 min each on
# a standard node).
#
# Both modes call 02_run_mash_single.sh, which loads the R/gsl modules
# and calls 02_run_mash.R --drug <Drug>.
#
# Usage:
#   bash 02_run_mash_launcher.sh          # submit all 5 per-drug jobs
#   bash 02_run_mash_launcher.sh Cer      # submit only the Cer job

SCRIPT_DIR=/project/gilad/geena/publication_clean/code/treatment/differential_expression/dreamlet
mkdir -p "$SCRIPT_DIR/logs"

TARGET=${1:-}

DRUGS=("Ale" "Bri" "Cer" "Cri" "Lor")

# If a specific drug name was passed, override the full list.
if [ -n "$TARGET" ]; then
  DRUGS=("$TARGET")
fi

for drug in "${DRUGS[@]}"; do
  echo "Submitting per-drug mash job for: $drug"
  sbatch \
    --job-name="dreamlet_mash_by_drug_${drug}" \
    --output="${SCRIPT_DIR}/logs/dreamlet_mash_by_drug_${drug}.out" \
    --error="${SCRIPT_DIR}/logs/dreamlet_mash_by_drug_${drug}.err" \
    "$SCRIPT_DIR/02_run_mash_single.sh" --drug "$drug"
done
echo "Submitted ${#DRUGS[@]} per-drug job(s)."

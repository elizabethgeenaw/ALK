#!/bin/bash

base_dir="/work"

list_file="$base_dir/medianUMI1.list"
workflow_path="$base_dir/fastqtl_to_mash.ipynb"
data_folder="/project/gilad/geena/publication_clean/mash_input/medianUMI1_4.30.26"
log_file="$base_dir/LOG_medianUMI1_fastqtl2mash_4.30.26.txt"

mkdir -p "$data_folder"

# Change space to underscore in .out.gz file names
for file in /project/gilad/geena/publication_clean/results/combined/medianUMI1/fastqtl_like/*; do
    mv "$file" "${file// /_}"
done

# Initialize log file if missing
if [ ! -f "$log_file" ]; then
  echo "Starting fastqtl2mash processing..." > "$log_file"
  echo "Log generated on $(date)" >> "$log_file"
  echo "----------------------------------" >> "$log_file"
fi


echo "Processing mash input for joint analysis for all treatments and cell types..." | tee -a "$log_file"

# Run the mash script for the entire combined list
sos run "$workflow_path" \
  --data-list "$list_file" \
  -j "${SLURM_CPUS_PER_TASK:-1}" \
  --cwd "$data_folder/" \
  --cols 3 4 5 >> "$log_file" 2>&1

if [ $? -eq 0 ]; then
  echo "Completed processing for joint analysis." | tee -a "$log_file"
else
  echo "Error encountered during joint analysis." | tee -a "$log_file"
fi

echo "fastqtl2mash joint processing completed. Check the log file at $log_file for details."

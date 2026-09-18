import anndata as ad
ad.__version__ # 0.8.0
import pandas as pd

# Load the combined pre-QC data
KR=ad.read_h5ad('/project2/gilad/katie/ebQTL/highpass_combinedFiles/102andPilot_MetaAdded_noNorm.h5ad')
print("Dataset dimensions:", KR.shape) # 1825288 rows (cells),36601 columns (features)
print(KR.obs.columns) # check colnames

# Check the number of lines 
KR_donor=KR.obs.donor_id.unique()
print(f"Number of unique donors: {len(KR_donor)}") 

# Filter valid donors
valid_donors = [d for d in KR.obs['donor_id'].unique() if d not in ['doublet', 'unassigned', 'SCM-9', 'SCM-10', 'H28834']]
print(f"Valid donors: {len(valid_donors)}")

# Number of splits
n_splits = 8  
donor_splits = [
    valid_donors[i::n_splits] for i in range(n_splits)
]

# Save each split
output_dir = "/project/gilad/geena/ALK_analysis/QC_preprocessing_control"
for i, donors in enumerate(donor_splits):
    subset = KR[KR.obs['donor_id'].isin(donors), :]
    output_path = f"{output_dir}/split_by_donor_{i+1}.h5ad"
    subset.write(output_path)
    print(f"Saved split {i+1} to {output_path}")

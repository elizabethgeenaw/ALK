# ---- V4 outputs ----
# CRISP pipeline shared config (drug names, paths, thresholds)
# --------------------

"""
Config for CRISP zero-shot perturbation prediction pipeline.

Predicts gene expression changes for 5 ALK inhibitors (1 µM) across 18 cell types using all control cells (no subsampling).  
Drug doses are normalized to the NeurIPS training scale.  
Drug embeddings are computed with descriptastorus matching CRISP's rdkit2D column schema. 
scGPT_human checkpoint is the correct one for the released CRISP weights.
"""

# Inputs 
CONTROL_RDS = "/project/gilad/geena/publication_clean/results/control/celltype_annotation/combined_control_filtnorm_cellid.rds"
CELLTYPE_COL = "cellid_type"
VALID_CELLTYPES_FILE = "/project/gilad/geena/publication_clean/matrixeqtl_input/control/pseudobulk_15cells_1indiv_summary/valid_celltypes_20indiv.tsv"

# Per-celltype MTX data exported by 0_export_control_cells.R
DATA_DIR = "/project/gilad/geena/chem_perturbation_clean/data/control_h5ad"

# CRISP repo (needed on sys.path so we can `import CRISP`)
CRISP_REPO_DIR = "/project/gilad/geena/chem_perturbation/code/chem_perturbation/CRISP"

# Downloaded checkpoints and reference data (copied to data/)
CRISP_MODEL_PATH = "/project/gilad/geena/chem_perturbation_clean/data/model_weights/nips_split_model.pt"
SCGPT_CHECKPOINT_DIR = "/project/gilad/geena/chem_perturbation_clean/data/scgpt_checkpoint/scGPT_human"
RDKIT2D_REF_PARQUET = "/project/gilad/geena/chem_perturbation_clean/data/drug_embeddings/rdkit2D_embedding_lincs_nips.parquet"
NIPS_REF_H5AD = "/project/gilad/geena/chem_perturbation_clean/data/nips_pp_scFM_resplit.h5ad"

# Output
OUTPUT_DIR = "/project/gilad/geena/chem_perturbation_clean/results/predictions"

# Drugs
TREATMENTS = ["Ale", "Bri", "Cer", "Cri", "Lor"]

DRUG_SMILES = {
    "Ale": "CCC1=CC2=C(C=C1N3CCC(CC3)N4CCOCC4)C(C5=C(C2=O)C6=C(N5)C=C(C=C6)C#N)(C)C",
    "Bri": "CN1CCN(CC1)C2CCN(CC2)C3=CC(=C(C=C3)NC4=NC=C(C(=N4)NC5=CC=CC=C5P(=O)(C)C)Cl)OC",
    "Cer": "CC1=CC(=C(C=C1C2CCNCC2)OC(C)C)NC3=NC=C(C(=N3)NC4=CC=CC=C4S(=O)(=O)C(C)C)Cl",
    "Cri": "C[C@H](C1=C(C=CC(=C1Cl)F)Cl)OC2=C(N=CC(=C2)C3=CN(N=C3)C4CCNCC4)N",
    "Lor": "C[C@@H]1C2=C(C=CC(=C2)F)C(=O)N(CC3=NN(C(=C3C4=CC(=C(N=C4)N)O1)C#N)C)C",
}

# Doses normalized to NeurIPS training scale (as in the zero-shot tutorial)
# NeurIPS DMSO reference dose was 14.1 µM; training data divides by this -> 0.0709
# All ALK inhibitors were dosed at 1 µM in the experiment.
_NIPS_DOSE_NORM = 14.1
DRUG_DOSES = {t: 1.0 / _NIPS_DOSE_NORM for t in TREATMENTS}

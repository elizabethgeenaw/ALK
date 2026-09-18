# ---- V4 outputs ----
# CRISP validation script (sanity checks, not directly in V4)
# --------------------

"""
Drug structure-response: identify genes that differ most between ALK inhibitors.
Crizotinib (multi-kinase: ALK/MET/ROS1) vs Lorlatinib (ALK-selective) should
differ on MET/ROS1-related genes rather than core ALK targets.
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg

OUT_DIR = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(OUT_DIR, exist_ok=True)

PRED_DIR = cfg.OUTPUT_DIR

# Get all completed cell types
cts = set()
for f in os.listdir(PRED_DIR):
    if f.endswith("_control_mean.parquet"):
        cts.add(f.replace("_control_mean.parquet", ""))

# Collect logFC per drug per cell type
all_lfc = {}
for ct in sorted(cts):
    ctrl = pd.read_parquet(os.path.join(PRED_DIR, f"{ct}_control_mean.parquet"))
    gene_names = ctrl.index.tolist()
    ctrl_vals = ctrl.values.flatten()

    for treatment in cfg.TREATMENTS:
        pred_path = os.path.join(PRED_DIR, f"{ct}_{treatment}_predicted_mean.parquet")
        if not os.path.exists(pred_path):
            continue
        pred = pd.read_parquet(pred_path).values.flatten()
        lfc = pred - ctrl_vals
        all_lfc[(ct, treatment)] = pd.Series(lfc, index=gene_names)

if not all_lfc:
    print("No predictions found.")
    sys.exit(0)

# Average logFC across cell types for each drug
drug_mean_lfc = {}
for treatment in cfg.TREATMENTS:
    drug_lfcs = [all_lfc[(ct, treatment)] for ct in sorted(cts) if (ct, treatment) in all_lfc]
    if drug_lfcs:
        drug_mean_lfc[treatment] = pd.concat(drug_lfcs, axis=1).mean(axis=1)

lfc_df = pd.DataFrame(drug_mean_lfc)
lfc_df.to_csv(os.path.join(OUT_DIR, "drug_mean_lfc.csv"))

# Pairwise drug-drug correlation of logFC profiles
print("=== Pairwise Pearson r of mean logFC profiles ===")
corr = lfc_df.corr()
print(corr.round(4).to_string())
corr.to_csv(os.path.join(OUT_DIR, "drug_lfc_correlation.csv"))

# Genes that differ most between Crizotinib (multi-kinase) and Lorlatinib (ALK-selective)
if "Cri" in lfc_df.columns and "Lor" in lfc_df.columns:
    diff = lfc_df["Cri"] - lfc_df["Lor"]
    diff_df = pd.DataFrame({
        "gene": diff.index,
        "Cri_lfc": lfc_df["Cri"].values,
        "Lor_lfc": lfc_df["Lor"].values,
        "diff_Cri_minus_Lor": diff.values,
        "abs_diff": diff.abs().values,
    }).sort_values("abs_diff", ascending=False)
    diff_df.to_csv(os.path.join(OUT_DIR, "cri_vs_lor_diff.csv"), index=False)

    print("\n=== Top 20 genes differing between Crizotinib and Lorlatinib ===")
    print(diff_df.head(20)[["gene", "Cri_lfc", "Lor_lfc", "diff_Cri_minus_Lor"]].to_string(index=False))

    # Check if MET/ROS1 are among top differentiators
    met_ros = diff_df[diff_df["gene"].isin(["MET", "ROS1", "HGF"])]
    if len(met_ros) > 0:
        print("\n=== MET/ROS1 pathway genes ===")
        print(met_ros[["gene", "Cri_lfc", "Lor_lfc", "diff_Cri_minus_Lor"]].to_string(index=False))

# Variance across drugs per gene (genes most drug-specific)
lfc_df["drug_variance"] = lfc_df[cfg.TREATMENTS].var(axis=1)
top_var = lfc_df.nlargest(30, "drug_variance")
print("\n=== Top 30 genes with highest variance across drugs ===")
print(top_var[cfg.TREATMENTS + ["drug_variance"]].round(4).to_string())

top_var[cfg.TREATMENTS + ["drug_variance"]].to_csv(
    os.path.join(OUT_DIR, "top_drug_variable_genes.csv")
)

print(f"\nDone. Results saved to {OUT_DIR}")

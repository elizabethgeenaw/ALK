# ---- V4 outputs ----
# CRISP validation script (sanity checks, not directly in V4)
# --------------------

"""
Cross-cell-type consistency: check whether core drug response genes are
consistently perturbed across cell types, while allowing for cell-type-specific effects.

Reports:
- Per-drug rank correlation of gene logFC across all cell type pairs
- Genes with consistently high |logFC| across cell types (core response)
- Genes with high variance across cell types (cell-type-specific)
"""
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

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
cts = sorted(cts)

# Collect logFC: {drug: DataFrame with genes as rows, cell types as columns}
drug_lfc = {t: {} for t in cfg.TREATMENTS}

for ct in cts:
    ctrl = pd.read_parquet(os.path.join(PRED_DIR, f"{ct}_control_mean.parquet"))
    ctrl_vals = ctrl.values.flatten()
    gene_names = ctrl.index.tolist()

    for treatment in cfg.TREATMENTS:
        pred_path = os.path.join(PRED_DIR, f"{ct}_{treatment}_predicted_mean.parquet")
        if not os.path.exists(pred_path):
            continue
        pred = pd.read_parquet(pred_path).values.flatten()
        drug_lfc[treatment][ct] = pd.Series(pred - ctrl_vals, index=gene_names)

for treatment in cfg.TREATMENTS:
    df = pd.DataFrame(drug_lfc[treatment])
    if df.shape[1] < 2:
        print(f"{treatment}: only {df.shape[1]} cell types, skipping")
        continue

    # Pairwise Spearman rank correlation across cell types
    ct_list = df.columns.tolist()
    corrs = []
    for i in range(len(ct_list)):
        for j in range(i + 1, len(ct_list)):
            r, p = spearmanr(df[ct_list[i]], df[ct_list[j]])
            corrs.append({"ct1": ct_list[i], "ct2": ct_list[j], "spearman_r": r, "p": p})
    corr_df = pd.DataFrame(corrs)
    mean_r = corr_df["spearman_r"].mean()
    min_r = corr_df["spearman_r"].min()
    max_r = corr_df["spearman_r"].max()
    print(f"\n{treatment}: {len(ct_list)} cell types, {len(corrs)} pairs")
    print(f"  Spearman r: mean={mean_r:.4f}, min={min_r:.4f}, max={max_r:.4f}")

    corr_df.to_csv(os.path.join(OUT_DIR, f"crossct_{treatment}_pairwise_spearman.csv"), index=False)

    # Genes with consistent effects: high mean |logFC| and low CV across cell types
    gene_stats = pd.DataFrame({
        "mean_lfc": df.mean(axis=1),
        "median_lfc": df.median(axis=1),
        "std_lfc": df.std(axis=1),
        "mean_abs_lfc": df.abs().mean(axis=1),
        "n_positive": (df > 0).sum(axis=1),
        "n_negative": (df < 0).sum(axis=1),
        "n_celltypes": df.shape[1],
    })
    # Fraction of cell types with same sign
    gene_stats["sign_consistency"] = gene_stats[["n_positive", "n_negative"]].max(axis=1) / gene_stats["n_celltypes"]

    gene_stats.to_csv(os.path.join(OUT_DIR, f"crossct_{treatment}_gene_stats.csv"))

    # Core response: high mean |logFC| AND high sign consistency
    core = gene_stats[(gene_stats["mean_abs_lfc"] > 0.1) & (gene_stats["sign_consistency"] >= 0.8)]
    core_sorted = core.sort_values("mean_abs_lfc", ascending=False)
    print(f"  Core response genes (|logFC|>0.1, sign consistent >=80%): {len(core_sorted)}")
    print(f"  Top 15:")
    print(core_sorted.head(15)[["mean_lfc", "mean_abs_lfc", "sign_consistency"]].round(4).to_string())

    core_sorted.to_csv(os.path.join(OUT_DIR, f"crossct_{treatment}_core_genes.csv"))

    # Cell-type-specific: high variance relative to mean
    variable = gene_stats[gene_stats["mean_abs_lfc"] > 0.05].copy()
    variable["cv"] = variable["std_lfc"] / (variable["mean_abs_lfc"] + 1e-6)
    variable_sorted = variable.sort_values("cv", ascending=False)
    print(f"\n  Most cell-type-variable genes (top 15):")
    print(variable_sorted.head(15)[["mean_lfc", "std_lfc", "cv", "sign_consistency"]].round(4).to_string())

    variable_sorted.head(100).to_csv(os.path.join(OUT_DIR, f"crossct_{treatment}_variable_genes.csv"))

print(f"\nDone. Results saved to {OUT_DIR}")

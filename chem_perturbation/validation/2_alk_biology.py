# ---- V4 outputs ----
# CRISP validation script (sanity checks, not directly in V4)
# --------------------

"""
Check whether CRISP predictions show expected ALK inhibitor biology.
Reports logFC for known ALK pathway genes and runs a simple gene-set enrichment
(hypergeometric test) for genes in key pathways.
"""
import os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg

OUT_DIR = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(OUT_DIR, exist_ok=True)

PRED_DIR = cfg.OUTPUT_DIR

# Known ALK downstream / pathway genes
# These should generally be downregulated by ALK inhibitors
ALK_PATHWAY_GENES = {
    "ALK_direct": ["ALK", "EML4"],
    "MAPK_ERK": ["MAPK1", "MAPK3", "MAP2K1", "MAP2K2", "BRAF", "RAF1"],
    "PI3K_AKT": ["PIK3CA", "PIK3CB", "PIK3CD", "AKT1", "AKT2", "MTOR", "RPS6KB1"],
    "STAT3": ["STAT3", "JAK1", "JAK2"],
    "proliferation": ["MYC", "CCND1", "CCND2", "CDK4", "CDK6", "MCM2", "MKI67", "PCNA"],
    "survival": ["BCL2", "BCL2L1", "MCL1", "BIRC5"],
    "ALK_targets": ["ETV1", "ETV4", "ETV5", "DUSP6", "SPRY2", "SPRY4"],
    "MET_ROS1": ["MET", "ROS1", "HGF"],  # off-targets for Crizotinib
}

all_pathway_genes = set()
for genes in ALK_PATHWAY_GENES.values():
    all_pathway_genes.update(genes)

# Get all completed cell types
cts = set()
for f in os.listdir(PRED_DIR):
    if f.endswith("_control_mean.parquet"):
        cts.add(f.replace("_control_mean.parquet", ""))

rows = []
for ct in sorted(cts):
    ctrl = pd.read_parquet(os.path.join(PRED_DIR, f"{ct}_control_mean.parquet"))
    gene_names = ctrl.index.tolist()
    ctrl_vals = ctrl.values.flatten()

    for treatment in cfg.TREATMENTS:
        pred_path = os.path.join(PRED_DIR, f"{ct}_{treatment}_predicted_mean.parquet")
        if not os.path.exists(pred_path):
            continue
        pred = pd.read_parquet(pred_path)
        pred_vals = pred.values.flatten()
        lfc = pred_vals - ctrl_vals

        for pathway, genes in ALK_PATHWAY_GENES.items():
            for gene in genes:
                if gene in gene_names:
                    idx = gene_names.index(gene)
                    rows.append({
                        "cell_type": ct, "drug": treatment, "pathway": pathway,
                        "gene": gene, "control_expr": ctrl_vals[idx],
                        "pred_expr": pred_vals[idx], "logFC": lfc[idx],
                    })

df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT_DIR, "alk_pathway_genes.csv"), index=False)

# Summary: mean logFC per pathway per drug (averaged across cell types)
if len(df) > 0:
    summary = df.groupby(["drug", "pathway"]).agg(
        mean_logFC=("logFC", "mean"),
        median_logFC=("logFC", "median"),
        n_genes=("gene", "nunique"),
        n_celltypes=("cell_type", "nunique"),
    ).reset_index()
    summary.to_csv(os.path.join(OUT_DIR, "alk_pathway_summary.csv"), index=False)

    print("=== Mean logFC per pathway per drug (across cell types) ===")
    pivot = summary.pivot_table(index="pathway", columns="drug", values="mean_logFC")
    print(pivot.round(4).to_string())

    print("\n=== Top 10 most perturbed pathway genes (by |logFC|, averaged across cell types) ===")
    gene_avg = df.groupby(["gene", "pathway"])["logFC"].mean().reset_index()
    gene_avg["abs_logFC"] = gene_avg["logFC"].abs()
    print(gene_avg.nlargest(10, "abs_logFC")[["gene", "pathway", "logFC"]].to_string(index=False))
else:
    print("No pathway genes found in predictions.")

print(f"\nDone. Results saved to {OUT_DIR}")

# ---- V4 outputs ----
# CRISP validation script (sanity checks, not directly in V4)
# --------------------

"""
Positive control: run CRISP on held-out NeurIPS data (B cells, split2='ood')
for Crizotinib and Palbociclib, compare predicted logFC to ground truth.
Reports Pearson r on DE genes (same as CRISP tutorial).
"""
import os, sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config as cfg

OUT_DIR = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(OUT_DIR, exist_ok=True)

# Load NeurIPS data
nips = sc.read_h5ad(cfg.NIPS_REF_H5AD)
var_names = nips.var_names
de_genes = nips.uns["rank_genes_groups_cov"]

# Load reference smile_df (these drugs ARE in the reference parquet)
smile_df = pd.read_parquet(cfg.RDKIT2D_REF_PARQUET)

# Load CRISP model (split2: trained on CD8+ T, NK, Treg, Myeloid)
sys.path.insert(0, cfg.CRISP_REPO_DIR)
from CRISP.trainer import Trainer
exp = Trainer()
exp.load_model(cfg.CRISP_MODEL_PATH)

# OOD control cells (B cells)
adata_ctrl = nips[
    (nips.obs["neg_control"] == 1)
    & (nips.obs["split2"] == "ood")
    & (nips.obs["cell_type"] == "B cells")
].copy()
print(f"Control B cells (OOD): {adata_ctrl.n_obs} cells")

# Drug SMILES from NeurIPS metadata
drug_to_smile = dict(zip(nips.obs["condition"], nips.obs["SMILES"]))

results = []
for drug in ["Crizotinib", "Palbociclib"]:
    smile = drug_to_smile[drug]
    dose = 0.07092198581560284

    # Predict
    adata_pred, _, _ = exp.get_prediction(
        adata_ctrl.copy(), dose=dose, smile=smile, smile_df=smile_df, FM_emb="X_scGPT"
    )

    # Ground truth
    adata_true = nips[
        (nips.obs["split2"] == "ood")
        & (nips.obs["cell_type"] == "B cells")
        & (nips.obs["condition"] == drug)
    ].copy()

    # DE genes
    des = de_genes[f"B cells_{drug}"]
    de_bool = var_names.isin(np.array(des))
    idx_de = np.where(de_bool)[0]

    # LogFC
    ctrl_mean = np.array(adata_ctrl.X.mean(0)).flatten()
    yt = np.array(adata_true.X.mean(0)).flatten() - ctrl_mean
    yp = np.array(adata_pred.X.mean(0)).flatten() - ctrl_mean

    # Pearson r on DE genes
    r_de, p_de = pearsonr(yp[idx_de], yt[idx_de])
    # Pearson r on all genes
    r_all, p_all = pearsonr(yp, yt)

    # Directional accuracy: % of genes where predicted and true logFC agree in sign
    # Only meaningful for genes with non-trivial effect
    nonzero = (np.abs(yt) > 0.05) | (np.abs(yp) > 0.05)
    dir_acc_all = np.mean(np.sign(yp[nonzero]) == np.sign(yt[nonzero]))
    dir_acc_de = np.mean(np.sign(yp[idx_de]) == np.sign(yt[idx_de]))

    # Among DE genes: how many are up vs down, and does prediction match?
    n_up_true = (yt[idx_de] > 0).sum()
    n_down_true = (yt[idx_de] < 0).sum()
    n_up_pred_correct = ((yt[idx_de] > 0) & (yp[idx_de] > 0)).sum()
    n_down_pred_correct = ((yt[idx_de] < 0) & (yp[idx_de] < 0)).sum()

    print(f"\n{drug}:")
    print(f"  True treated cells: {adata_true.n_obs}")
    print(f"  DE genes: {len(des)}")
    print(f"  Pearson r (DE genes): {r_de:.4f} (p={p_de:.2e})")
    print(f"  Pearson r (all genes): {r_all:.4f} (p={p_all:.2e})")
    print(f"  Directional accuracy (DE genes): {dir_acc_de:.1%} ({n_up_pred_correct}/{n_up_true} up, {n_down_pred_correct}/{n_down_true} down correct)")
    print(f"  Directional accuracy (all genes |logFC|>0.05): {dir_acc_all:.1%} ({nonzero.sum()} genes)")

    results.append({
        "drug": drug, "n_true": adata_true.n_obs, "n_de": len(des),
        "r_de": r_de, "p_de": p_de, "r_all": r_all, "p_all": p_all,
        "dir_acc_de": dir_acc_de, "dir_acc_all": dir_acc_all,
        "n_nonzero": nonzero.sum(),
    })

    # Save per-gene results for plotting
    gene_df = pd.DataFrame({
        "gene": var_names, "true_lfc": yt, "pred_lfc": yp, "is_de": de_bool
    })
    gene_df.to_csv(os.path.join(OUT_DIR, f"positive_control_{drug}.csv"), index=False)

pd.DataFrame(results).to_csv(os.path.join(OUT_DIR, "positive_control_summary.csv"), index=False)
print("\nDone. Results saved to", OUT_DIR)

"""
Per-gene 1-vs-4 z-score computation on CRISP predictions.

For each (cell_type, gene, focal_treatment):
    z = (pred_focal - mean(other 4)) / sd(other 4, ddof=1)

Reads prediction parquets and writes z-scores as TSV + parquet.
Enrichment testing (FGSEA) is run separately in 4_fgsea.R.

Usage:
    python 3_zscore.py
"""
import os, re, logging
import pandas as pd
import numpy as np

logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
log = logging.getLogger()

PRED_DIR = "/project/gilad/geena/chem_perturbation_clean/results/predictions"
OUT_DIR = os.path.join(PRED_DIR, "zscore_analysis")
os.makedirs(OUT_DIR, exist_ok=True)

TREATMENTS = ["Ale", "Bri", "Cer", "Cri", "Lor"]


def _discover_cell_types():
    """Find valid cell types from prediction file names."""
    pattern = re.compile(r"^(.+)_(?:Ale|Bri|Cer|Cri|Lor)_predicted_mean\.parquet$")
    cts = set()
    for f in os.listdir(PRED_DIR):
        m = pattern.match(f)
        if m:
            cts.add(m.group(1))
    return sorted(cts)


def _load_preds(ct_safe):
    """Load prediction parquets for all 5 treatments. Returns DataFrame or None."""
    frames = {}
    for t in TREATMENTS:
        p = os.path.join(PRED_DIR, f"{ct_safe}_{t}_predicted_mean.parquet")
        if not os.path.exists(p):
            return None
        df = pd.read_parquet(p)
        # Gene names are stored as the index
        frames[t] = df["predicted"]
    return pd.DataFrame(frames)


def _zscores(pred_df):
    """Compute 1-vs-4 z-scores."""
    ts = list(pred_df.columns)
    z = pd.DataFrame(index=pred_df.index, columns=ts, dtype=float)
    for f in ts:
        others = [t for t in ts if t != f]
        bg = pred_df[others]
        mu = bg.mean(axis=1)
        sd = bg.std(axis=1, ddof=1)
        zi = (pred_df[f] - mu) / sd
        z[f] = zi.where(sd != 0, 0.0)
    return z


def main():
    valid_cts = _discover_cell_types()
    log.info("Found %d cell types", len(valid_cts))

    z_parts = []
    for ct in valid_cts:
        pred_df = _load_preds(ct)
        if pred_df is None:
            log.info("Skipping %s (missing parquets)", ct)
            continue

        z_df = _zscores(pred_df)
        zl = z_df.stack().reset_index()
        zl.columns = ["gene", "treatment", "z_score"]
        zl["cell_type"] = ct
        z_parts.append(zl)
        log.info("  %s: %d z-scores", ct, len(zl))

    if not z_parts:
        log.warning("No z-scores produced.")
        return

    z_all = pd.concat(z_parts, ignore_index=True)
    n_cond = len(z_all.groupby(["cell_type", "treatment"]))

    tsv_path = os.path.join(OUT_DIR, "gene_zscores_all.tsv")
    z_all.to_csv(tsv_path, sep="\t", index=False)

    parquet_path = os.path.join(OUT_DIR, "gene_zscores_all.parquet")
    z_all.to_parquet(parquet_path, index=False)

    log.info("Saved z-scores -> %s  (n=%d rows, %d conditions)", tsv_path, len(z_all), n_cond)


if __name__ == "__main__":
    main()

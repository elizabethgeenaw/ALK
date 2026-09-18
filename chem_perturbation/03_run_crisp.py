"""
Run CRISP zero-shot predictions on preprocessed cells.  For each cell type and treatment, save the mean predicted expression (parquet).

To use: python 2_run_crisp.py
If cell type is huge (ex: radial glia ~289k cells) and GPU out of memory and need to break up into chunks: 
CRISP_BATCH_SIZE=100000 python 2_run_crisp.py
where if >0 and n_cells > CRISP_BATCH_SIZE -> process in chunks
"""
import os
import sys
import logging
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, os.path.dirname(__file__))
import config as cfg

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
log = logging.getLogger(__name__)

VALID_CTS = [
    "Bronchiolar_and_alveolar_epithelial_cells",
    "Cardiac_mesoderm_cells",
    "Ciliated_epithelial_cells",
    "Early_mesoderm_cells",
    "ENS_glia",
    "Epicardial_fat_cells",
    "Glioblast",
    "Goblet_cells",
    "Hepatoblasts",
    "Metanephric_cells",
    "Neuroblast",
    "Neuron",
    "Neuronal_IPC",
    "Radial_glia",
    "Retinal_neurons",
    "Retinal_pigment_cells",
    "Retinal_progenitors_and_Muller_glia",
    "Squamous_epithelial_cells",
    "Stromal_cells",
    "Vascular_endothelial_cells",
    "iPSCs",
]


def _build_smile_df(smiles_dict):
    """Convert drug SMILES strings into z-scored RDKit2D descriptor vectors (CRISP expects as input).
    Same pipeline as CRISP's sciplex_dataprocess.ipynb: generate 200 RDKit2D descriptors with descriptastorus, drop constant columns, z-score normalize based on mean and SD of drugs used for training.
    Use training drug mean/SD from paper so our embeddings are on the same scale
    """
    from descriptastorus.descriptors.DescriptorGenerator import MakeGenerator

    # Get normalization stats precomputed from CRISP
    stats_path = os.path.join(os.path.dirname(__file__), "data", "rdkit2d_normstats.npz")
    stats = np.load(stats_path, allow_pickle=True)
    mu = pd.Series(stats["mu"], index=list(stats["kept_cols"]))
    sd = pd.Series(stats["sd"], index=list(stats["kept_cols"]))
    drop_cols = list(stats["drop_cols"]) # constant/uninformative descriptors removed during training
    ref_cols = list(stats["ref_cols"]) # final column order expected by CRISP

    # Generate 200 RDKit2D descriptors for each drug
    gen = MakeGenerator(("RDKit2D",))
    rows = {}
    for drug, smi in smiles_dict.items():
        result = gen.process(smi)
        if result is None or not result[0]:  # result[0] is a success flag
            raise ValueError(f"descriptastorus failed for {drug}")
        rows[smi] = np.asarray(result[1:], dtype=np.float64)  # result[1:] are the 200 descriptor values

    full = pd.DataFrame(
        np.vstack(list(rows.values())),
        index=list(rows.keys()),
        columns=[f"latent_{i+1}" for i in range(200)],
    )
    # Drop constant columns, z-score, reorder to match CRISP
    kept = full.drop(columns=drop_cols)
    zscored = (kept - mu) / sd
    zscored = zscored.fillna(0.0)  # NaN from 0/0 division -> 0
    smile_df = zscored[ref_cols]
    smile_df.index.name = "smiles"
    return smile_df


def _load_trainer():
    """Load the CRISP model from the repo and pre-trained weights."""
    sys.path.insert(0, cfg.CRISP_REPO_DIR)  # add CRISP repo to path so we can import it
    from CRISP.trainer import Trainer
    trainer = Trainer()
    trainer.load_model(cfg.CRISP_MODEL_PATH)
    log.info("Loaded CRISP model: %s", cfg.CRISP_MODEL_PATH)
    return trainer


def main():
    pred_dir = cfg.OUTPUT_DIR
    os.makedirs(pred_dir, exist_ok=True)
    prepped_dir = os.path.join(pred_dir, "_prepped")  # h5ad files from 1_prep_input

    trainer = _load_trainer()
    smile_df = _build_smile_df(cfg.DRUG_SMILES) # get drugs from config file

    cts_to_run = [os.environ["SINGLE_CT"]] if os.environ.get("SINGLE_CT") else VALID_CTS

    for ct_safe in cts_to_run:
        in_path = os.path.join(prepped_dir, f"{ct_safe}.h5ad")
        if not os.path.exists(in_path):
            log.warning("Missing prepped h5ad for %s, skipping.", ct_safe)
            continue

        log.info("--- %s ---", ct_safe)
        adata_ctrl = sc.read_h5ad(in_path)
        log.info("Loaded %d cells x %d genes", adata_ctrl.n_obs, adata_ctrl.n_vars)
        if "X_scGPT" not in adata_ctrl.obsm:
            log.error("X_scGPT missing — skip.")
            continue

        # Save control mean expression (average of input cells, NOT a CRISP prediction)
        # This is the raw observed expression, saved for reference 
        ctrl_mean = pd.Series(
            np.array(adata_ctrl.X.mean(axis=0)).flatten(),
            index=adata_ctrl.var_names,
            name="control_mean",
        )
        ctrl_mean.to_frame().to_parquet(
            os.path.join(pred_dir, f"{ct_safe}_control_mean.parquet")
        )

        for treatment in cfg.TREATMENTS:
            mean_out = os.path.join(pred_dir, f"{ct_safe}_{treatment}_predicted_mean.parquet")
            # Skip if already predicted (allows restarting without redoing work)
            if os.path.exists(mean_out):
                log.info("  %s x %s exists, skipping", ct_safe, treatment)
                continue

            smi = cfg.DRUG_SMILES[treatment]
            dose = cfg.DRUG_DOSES[treatment]
            log.info("  Predict %s x %s (%d cells, dose=%.5f) ...",
                     ct_safe, treatment, adata_ctrl.n_obs, dose)

            batch_size = int(os.environ.get("CRISP_BATCH_SIZE", 0))
            try:
                if batch_size > 0 and adata_ctrl.n_obs > batch_size:
                    # Batching: process cells in chunks to avoid GPU out-of-memory errors.
                    log.info("Batching %d cells into chunks of %d",
                             adata_ctrl.n_obs, batch_size)
                    pred_sum = np.zeros(adata_ctrl.n_vars, dtype=np.float64)
                    for start in range(0, adata_ctrl.n_obs, batch_size):
                        end = min(start + batch_size, adata_ctrl.n_obs)
                        chunk = adata_ctrl[start:end].copy()
                        adata_chunk, _, _ = trainer.get_prediction(
                            chunk,
                            dose=dose,
                            smile=smi,
                            smile_df=smile_df,
                            FM_emb="X_scGPT",  # tells CRISP to use scGPT embeddings from obsm
                        )
                        pred_sum += np.array(adata_chunk.X.sum(axis=0)).flatten()
                        log.info("    Batch %d-%d done", start, end)
                    mean_vals = pred_sum / adata_ctrl.n_obs
                else:
                    # Predict all cells at once, take mean
                    adata_pred, _, _ = trainer.get_prediction(
                        adata_ctrl.copy(),
                        dose=dose,
                        smile=smi,
                        smile_df=smile_df,
                        FM_emb="X_scGPT",
                    )
                    mean_vals = np.array(adata_pred.X.mean(axis=0)).flatten()
            except Exception as e:
                log.error("FAILED %s x %s: %s", ct_safe, treatment, e)
                continue

            # Save only the mean prediction per gene (not per-cell) to save disk space
            mean_pred = pd.Series(
                mean_vals,
                index=adata_ctrl.var_names,
                name="predicted",
            )
            mean_pred.to_frame().to_parquet(mean_out)
            log.info("Saved mean prediction -> %s", mean_out)

    log.info("Done.")


if __name__ == "__main__":
    main()

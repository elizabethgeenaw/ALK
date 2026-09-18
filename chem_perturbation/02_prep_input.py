"""
Preprocess control cells for CRISP inference: normalize, align genes to the
CRISP training set, and compute scGPT cell embeddings.
To run: python 1_prep_input.py
To run for single cell type: SINGLE_CT=Glioblast python 1_prep_input.py
"""
import os
import sys
import logging
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.io import mmread
from scipy.sparse import csr_matrix

sys.path.insert(0, os.path.dirname(__file__)) # add script's directory to Python's import search path to find config file
import config as cfg

logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s") # get python log to print messages (ex: INFO Loading CRISP reference....)
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


def _load_celltype_dir(ct_dir: str) -> ad.AnnData:
    """Load a cell type matrix info (from 0_export_control_cells.R) into AnnData.
    The R export writes counts as genes x cells (writeMM default), 
    so need to transpose to cells x genes for scanpy/AnnData convention.
    """
    counts = mmread(os.path.join(ct_dir, "counts.mtx")).tocsr() # convert to compressed sparse row (csr) for efficient row-wise operations
    barcodes = [l.strip() for l in open(os.path.join(ct_dir, "barcodes.tsv"))]
    genes = [l.strip() for l in open(os.path.join(ct_dir, "genes.tsv"))]
    obs = pd.read_csv(os.path.join(ct_dir, "obs.tsv"), sep="\t", index_col=0)
    X = counts.T.tocsr()  # transpose: genes x cells -> cells x genes
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
    adata.obs_names = barcodes
    return adata


def _align_genes_to_crisp(adata: ad.AnnData, crisp_genes: list) -> ad.AnnData:
    """Reorder genes to match CRISP's training gene vocabulary (18,211 genes).
    Our data has ~30k genes/transcripts but CRISP only has 18,211. 
    Create a new matrix with exactly those genes in the same order. 
    Genes present in our data are copied over and genes missing from our data are filled with 0.
    """
    n_cells = adata.n_obs
    n_genes = len(crisp_genes)
    in_idx = {g: i for i, g in enumerate(adata.var_names)}
    overlap = [g for g in crisp_genes if g in in_idx]
    log.info("Gene overlap with CRISP training set: %d / %d", len(overlap), n_genes)
    X_in = adata.X.tocsr() if hasattr(adata.X, "tocsr") else csr_matrix(adata.X)
    out = np.zeros((n_cells, n_genes), dtype=np.float32)
    for j, g in enumerate(crisp_genes):
        if g in in_idx:
            out[:, j] = X_in[:, in_idx[g]].toarray().flatten()
    new = ad.AnnData(X=out, obs=adata.obs.copy(), var=pd.DataFrame(index=crisp_genes))
    new.obs_names = adata.obs_names
    return new


def _embed_with_scgpt(adata: ad.AnnData, scgpt_dir: str):
    """Compute scGPT cell embeddings and store in adata.obsm['X_scGPT'].
    scGPT is a foundation model for single-cell data. It reads each cell's
    expression profile and produces a 512-dim embedding capturing the cell's
    transcriptomic state. CRISP uses these embeddings as its representation of cell identity.
    """
    from scgpt.tasks import embed_data
    try:
        embedded = embed_data(
            adata,
            model_dir=scgpt_dir,
            gene_col="index", # gene names are in AnnData index (var_names), not column in var
            batch_size=64, # scGPT default internal batch size for GPU inference
            return_new_adata=True, # return new AnnData where .X contains embeddings
        )
        adata.obsm["X_scGPT"] = np.array(embedded.X, dtype=np.float32)
        log.info("scGPT embedding shape: %s", adata.obsm["X_scGPT"].shape)
    except Exception as e:
        log.error("scGPT embed_data failed: %s", e)
        import traceback
        traceback.print_exc()


def main():
    prepped_dir = os.path.join(cfg.OUTPUT_DIR, "_prepped")
    os.makedirs(prepped_dir, exist_ok=True)

    # Load the CRISP training reference to get the exact gene list and order.
    # This is a large file (~3.6 GB) but we only need var_names, so we delete
    # the object immediately after extracting the gene list.
    log.info("Loading CRISP reference h5ad to get training gene list ...")
    nips_ref = sc.read_h5ad(cfg.NIPS_REF_H5AD)
    crisp_genes = nips_ref.var_names.tolist()
    log.info("CRISP training gene list: %d genes", len(crisp_genes))
    del nips_ref

    # SINGLE_CT env var allows running one cell type at a time (for SLURM array jobs). 
    # If not set, process all 18 cell types sequentially.
    cts_to_run = [os.environ["SINGLE_CT"]] if os.environ.get("SINGLE_CT") else VALID_CTS

    for ct_safe in cts_to_run:
        ct_dir = os.path.join(cfg.DATA_DIR, ct_safe)
        out_path = os.path.join(prepped_dir, f"{ct_safe}.h5ad")

        if not os.path.isdir(ct_dir):
            log.warning("Skipping %s: no exported dir", ct_safe)
            continue
        # Skip if already prepped (allows restarting without redoing work)
        if os.path.exists(out_path):
            log.info("%s already prepped, skipping.", ct_safe)
            continue

        log.info("--- %s ---", ct_safe)
        adata = _load_celltype_dir(ct_dir)
        log.info("Loaded: %d cells x %d genes", adata.n_obs, adata.n_vars)

        # Standard scRNA-seq normalization to match CRISP's expected input (counts per 10k, then log1p)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        # Reorder genes to match CRISP's vocabulary (missing genes -> 0)
        adata = _align_genes_to_crisp(adata, crisp_genes)

        # Compute scGPT embeddings (this is the slowest step and requires GPU)
        log.info("  Computing scGPT embeddings for %d cells ...", adata.n_obs)
        _embed_with_scgpt(adata, cfg.SCGPT_CHECKPOINT_DIR)

        if "X_scGPT" not in adata.obsm:
            log.error("  X_scGPT missing after embedding — skipping %s", ct_safe)
            continue

        # Set columns to str to avoid h5py errors when writing h5ad (h5py can't handle mixed-type object arrays)
        for col in adata.obs.columns:
            if adata.obs[col].dtype == object:
                adata.obs[col] = adata.obs[col].astype(str)

        # Save preprocessed AnnData with scGPT embeddings for step 2
        adata.write_h5ad(out_path)
        log.info("  Saved -> %s", out_path)

    log.info("Done.")


if __name__ == "__main__":
    main()

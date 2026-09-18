# ---- V4 outputs ----
# CRISP pipeline utility: validates SMILES input
# --------------------

"""
Validate that _build_smile_df produces embeddings consistent with
CRISP's reference rdkit2D parquet.

Picks a few SMILES from the reference parquet, computes embeddings
via _build_smile_df, and compares against the reference values.

Usage: python validate_smile_df.py
"""
import sys
import os
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))
import config as cfg

# reuse the same function from the main pipeline
from importlib.util import spec_from_file_location, module_from_spec
_spec = spec_from_file_location("run_crisp", os.path.join(os.path.dirname(__file__), "2_run_crisp.py"))
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)
_build_smile_df = _mod._build_smile_df

N_TEST = 5  # number of reference drugs to test


def main():
    ref_path = cfg.RDKIT2D_REF_PARQUET
    if not os.path.exists(ref_path):
        print(f"Reference parquet not found: {ref_path}")
        sys.exit(1)

    ref_df = pd.read_parquet(ref_path)
    print(f"Reference parquet: {ref_df.shape[0]} drugs, {ref_df.shape[1]} columns")

    # Pick a few SMILES from the reference parquet as test cases
    test_smiles = {f"ref_drug_{i}": smi for i, smi in enumerate(ref_df.index[:N_TEST])}
    print(f"Testing {len(test_smiles)} reference SMILES through _build_smile_df...")

    smile_df = _build_smile_df(test_smiles)

    shared_cols = smile_df.columns.intersection(ref_df.columns)
    print(f"Shared columns: {len(shared_cols)} / {len(smile_df.columns)} (ours) vs {ref_df.shape[1]} (ref)")

    if len(shared_cols) == 0:
        print("ERROR: No shared columns — column naming mismatch.")
        sys.exit(1)

    all_ok = True
    for smi in smile_df.index:
        ours = smile_df.loc[smi, shared_cols].values.astype(np.float64)
        ref = ref_df.loc[smi, shared_cols].values.astype(np.float64)
        max_diff = np.nanmax(np.abs(ours - ref))
        mean_diff = np.nanmean(np.abs(ours - ref))
        if max_diff > 0.01:
            print(f"  MISMATCH  {smi[:60]}  max={max_diff:.4f}  mean={mean_diff:.6f}")
            all_ok = False
        else:
            print(f"  OK        {smi[:60]}  max={max_diff:.6f}  mean={mean_diff:.8f}")

    if all_ok:
        print("\nAll test SMILES validated successfully.")
    else:
        print("\nSome SMILES had mismatches — check normalization stats.")
        sys.exit(1)


if __name__ == "__main__":
    main()

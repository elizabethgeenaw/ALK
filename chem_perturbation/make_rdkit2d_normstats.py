# ---- V4 outputs ----
# CRISP pipeline utility: computes RDKit 2D normalization stats
# --------------------

"""
make_rdkit2d_normstats.py
Compute normalization stats from RAW RDKit2D descriptors of CRISP's training drugs.

The reference parquet (rdkit2D_embedding_lincs_nips.parquet) is already z-scored,
so we cannot derive raw mean/SD from it. Instead, this script:
1. Reads the reference parquet to get the list of training drug SMILES
2. Regenerates raw RDKit2D descriptors for all training drugs using descriptastorus
3. Drops latent_0 (the RDKit2D_calculated success flag, not a real descriptor)
4. Drops columns with std <= 0.01 (matching CRISP's sciplex_dataprocess.ipynb)
5. Computes mean and SD of the remaining raw columns
6. Saves to data/rdkit2d_normstats.npz

The resulting .npz is used by 2_run_crisp.py to z-score our ALK inhibitors'
descriptors on the same scale the CRISP model was trained on.

To run: python make_rdkit2d_normstats.py
(Takes ~10 min for ~17k drugs; can run on a compute node.)
"""
import os
import numpy as np
import pandas as pd
from descriptastorus.descriptors.DescriptorGenerator import MakeGenerator

RDKIT2D_PARQUET = "/project/gilad/geena/chem_perturbation_clean/data/drug_embeddings/rdkit2D_embedding_lincs_nips.parquet"
OUT_PATH = os.path.join(os.path.dirname(__file__), "data", "rdkit2d_normstats.npz")

STD_THRESHOLD = 0.01  # CRISP drops columns with std <= this


def main():
    # Get training drug SMILES from the reference parquet
    ref_df = pd.read_parquet(RDKIT2D_PARQUET)
    smiles_list = list(ref_df.index)
    print(f"Reference parquet: {len(smiles_list)} drugs, {ref_df.shape[1]} columns")

    # Regenerate raw RDKit2D descriptors for all training drugs
    gen = MakeGenerator(("RDKit2D",))
    print(f"Generating raw RDKit2D descriptors for {len(smiles_list)} drugs...")
    rows = []
    failed = 0
    for i, smi in enumerate(smiles_list):
        result = gen.process(smi)
        if result is None or not result[0]:
            # Match CRISP's handling: fill with zeros for failed drugs
            rows.append(np.zeros(200, dtype=np.float64))
            failed += 1
        else:
            rows.append(np.asarray(result[1:], dtype=np.float64))  # skip result[0] (success flag)
        if (i + 1) % 2000 == 0:
            print(f"  {i+1}/{len(smiles_list)} done")

    print(f"Done. {failed} drugs failed descriptor generation.")

    # Build dataframe: 200 descriptors named latent_1 .. latent_200
    raw_df = pd.DataFrame(
        np.vstack(rows),
        index=smiles_list,
        columns=[f"latent_{i+1}" for i in range(200)],
    )

    # Replace NaN/Inf with 0 (matching CRISP's preprocessing)
    raw_df.replace([np.inf, -np.inf], 0, inplace=True)
    raw_df.fillna(0, inplace=True)

    # Drop columns with std <= threshold (matching CRISP's sciplex_dataprocess.ipynb)
    sds = raw_df.std(axis=0)
    drop_cols = sds[sds <= STD_THRESHOLD].index.tolist()
    print(f"Columns with std <= {STD_THRESHOLD}: {len(drop_cols)} -> {drop_cols}")

    kept = raw_df.drop(columns=drop_cols)
    kept_cols = kept.columns.tolist()
    print(f"Kept columns: {len(kept_cols)}")

    # Verify column count matches reference parquet
    if len(kept_cols) != ref_df.shape[1]:
        print(f"WARNING: kept {len(kept_cols)} columns but reference has {ref_df.shape[1]}.")
        print(f"  Reference columns: {list(ref_df.columns)}")
        print(f"  Our columns: {kept_cols}")

    # Compute raw mean and SD (these are the correct z-score stats)
    mu = kept.mean(axis=0).values
    sd = kept.std(axis=0).values

    # Save
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    np.savez(
        OUT_PATH,
        mu=mu,
        sd=sd,
        kept_cols=np.array(kept_cols, dtype=object),
        drop_cols=np.array(drop_cols, dtype=object),
        ref_cols=np.array(kept_cols, dtype=object),
    )
    print(f"Saved -> {OUT_PATH}")


if __name__ == "__main__":
    main()

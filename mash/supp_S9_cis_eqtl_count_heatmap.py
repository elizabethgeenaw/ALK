# ---- V4 outputs ----
# Supp Fig S9: /project/gilad/geena/ALK_draft/V4/Supplementary Figures/S9_fig_cis_eqtl_count_heatmap.pdf
# --------------------

"""
Cis-eQTL count heatmap: unique gene_snp pairs (lfsr < 0.05) per
cell type x condition (5 drugs + Control).

Reads:
  /project/gilad/geena/publication_clean/results/mash/medianUMI1/supp_table_all_significant_eqtls.tsv

Writes:
  /project/gilad/geena/publication_clean/results/mash/medianUMI1/plots/
    fig_cis_eqtl_count_heatmap.png
    fig_cis_eqtl_count_heatmap.pdf
"""

import os
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['font.family'] = 'DejaVu Sans'

SIG_TSV = "/project/gilad/geena/publication_clean/results/mash/medianUMI1/supp_table_all_significant_eqtls.tsv"
OUT_DIR = "/project/gilad/geena/publication_clean/results/mash/medianUMI1/plots"

CONDITION_ORDER = ['Ale', 'Bri', 'Cer', 'Cri', 'Lor', 'Control']


def make_pivot(df, condition_order):
    pivot = (
        df.groupby(['celltype', 'treatment'])['gene_snp']
        .nunique()
        .unstack('treatment')
        .reindex(columns=condition_order)
    )
    pivot = pivot.fillna(0).astype(int)
    pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).index]
    return pivot


def make_heatmap(pivot, out_path_png, out_path_pdf, title, cbar_label):
    n_rows, n_cols = pivot.shape

    fig_h = max(5.0, min(10.0, 1.5 + n_rows * 0.38))
    fig_w = 6 + 0.4 * max(0, n_cols - 5)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    values = pivot.values.astype(float)
    im = ax.imshow(values, cmap='YlOrRd', aspect='auto',
                   vmin=0, vmax=values.max())

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(list(pivot.columns), fontsize=10)
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(pivot.index, fontsize=8)

    vmax = values.max()
    for i in range(n_rows):
        for j in range(n_cols):
            v = int(values[i, j])
            ax.text(j, i, v, ha='center', va='center', fontsize=7,
                    color='white' if vmax > 0 and v > vmax * 0.6 else 'black')

    ax.set_title(title, fontsize=11, pad=8)
    plt.colorbar(im, ax=ax, shrink=0.6, label=cbar_label)
    plt.tight_layout()
    plt.savefig(out_path_png, dpi=200, bbox_inches='tight')
    plt.savefig(out_path_pdf, bbox_inches='tight')
    plt.close(fig)
    print(f"  wrote {out_path_png}")
    print(f"  wrote {out_path_pdf}")


if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)

    print(f"Reading {SIG_TSV} ...")
    sig = pd.read_csv(SIG_TSV, sep='\t')

    pivot = make_pivot(sig, CONDITION_ORDER)
    print(pivot)

    make_heatmap(
        pivot,
        f"{OUT_DIR}/fig_cis_eqtl_count_heatmap.png",
        f"{OUT_DIR}/fig_cis_eqtl_count_heatmap.pdf",
        "Significant cis-eQTLs (lfsr < 0.05) per cell type × condition",
        "# cis-eQTLs (gene–SNP pairs)",
    )

    print("\nDone.")

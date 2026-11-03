# ---- V4 outputs ----
# Supp Fig S17: /project/gilad/geena/ALK_draft/V4/Supplementary Figures/S17_crisp_methods_schematic.pdf
# --------------------

"""
Methods schematic for the z-score + fGSEA CRISP analysis.

Six-panel summary figure (A-F) describing the canonical pipeline:
  A. CRISP-predicted expression per gene (3 real example genes, focal = Bri)
  B. 1-vs-4 z-score formula + worked example (using panel-A high-z gene)
  C. Z-score distribution histogram (cardiac mesoderm x Bri, with response-eGene rug)
  D. Ranked gene list + fGSEA running ES curve
  E. fGSEA enrichment test description
  F. -log10(nominal p) heatmap across 21 cell types x 5 treatments,
     with BH-corrected significant cells marked.

Reads (does not modify):
  results/predictions/{ct}_{tr}_predicted_mean.parquet
  results/predictions/zscore_analysis/gene_zscores_all.tsv
  results/predictions/zscore_analysis/fgsea_results.tsv
  publication_clean/results/mash/medianUMI1/lfc3/treatment_responsive_eqtls_signfixed.tsv
  publication_clean/results/mash/medianUMI1/combined_mash_results.tsv

Writes:
  results/predictions/zscore_analysis/slides/fig6_methods_schematic.png
  results/predictions/zscore_analysis/slides/fig6_methods_schematic.pdf

Usage:
  source activate /project/gilad/geena/envs/crisp_env
  python 06_supp_S17_methods_schematic.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec

# ── Paths ─────────────────────────────────────────────────────────────────────
PRED_DIR   = "/project/gilad/geena/chem_perturbation_clean/results/predictions"
ZSCORE_DIR = os.path.join(PRED_DIR, "zscore_analysis")
OUT_DIR    = os.path.join(ZSCORE_DIR, "slides")
os.makedirs(OUT_DIR, exist_ok=True)

RESP_PATH  = "/project/gilad/geena/publication_clean/results/mash/medianUMI1/lfc3/treatment_responsive_eqtls_signfixed.tsv"
MASH_PATH  = "/project/gilad/geena/publication_clean/results/mash/medianUMI1/combined_mash_results.tsv"
FGSEA_PATH = os.path.join(ZSCORE_DIR, "fgsea_results.tsv")
ZSCORE_PATH = os.path.join(ZSCORE_DIR, "gene_zscores_all.tsv")

# ── Constants ─────────────────────────────────────────────────────────────────
TREATMENTS = ["Ale", "Bri", "Cer", "Cri", "Lor"]
DRUG_NAMES = {
    "Ale": "Alectinib", "Bri": "Brigatinib", "Cer": "Ceritinib",
    "Cri": "Crizotinib", "Lor": "Lorlatinib",
}
DRUG_COLORS = {
    "Ale": "#4C72B0", "Bri": "#DD8452", "Cer": "#55A868",
    "Cri": "#C44E52", "Lor": "#8172B3",
}
CT_LABELS = {
    "Bronchiolar_and_alveolar_epithelial_cells": "Bronch./alv. epithelial",
    "Cardiac_mesoderm_cells":                    "Cardiac mesoderm",
    "Ciliated_epithelial_cells":                 "Ciliated epithelial",
    "ENS_glia":                                  "ENS glia",
    "Early_mesoderm_cells":                      "Early mesoderm",
    "Epicardial_fat_cells":                      "Epicardial fat",
    "Glioblast":                                 "Glioblast",
    "Goblet_cells":                              "Goblet",
    "Hepatoblasts":                              "Hepatoblasts",
    "Metanephric_cells":                         "Metanephric",
    "Neuroblast":                                "Neuroblast",
    "Neuron":                                    "Neuron",
    "Neuronal_IPC":                              "Neuronal IPC",
    "Radial_glia":                               "Radial glia",
    "Retinal_neurons":                           "Retinal neurons",
    "Retinal_pigment_cells":                     "Retinal pigment",
    "Retinal_progenitors_and_Muller_glia":       "Retinal prog./Muller",
    "Squamous_epithelial_cells":                 "Squamous epithelial",
    "Stromal_cells":                             "Stromal",
    "Vascular_endothelial_cells":                "Vascular endothelial",
    "iPSCs":                                     "iPSCs",
}

# Example condition for panels A-D
EX_CT   = "Cardiac_mesoderm_cells"
EX_TR   = "Bri"

# Panel A genes (precomputed from the data; see implementation notes)
PANEL_A_GENES = ["TGFBR3", "RPL13A", "ANK3"]

OTHER_COLOR  = "#A8C5DA"
FOCAL_COLOR  = "#D62728"
RESP_COLOR   = "#E63946"
ES_COLOR     = "#2166ac"

# Self-contained figure caption (no supp methods section in the manuscript).
CAPTION = (
    "Figure SX. The z-score + fGSEA pipeline for testing whether CRISP-predicted drug "
    "responses are enriched for response eQTLs. Panels A–C use cardiac mesoderm × brigatinib "
    "as the worked example. "
    "(A) CRISP-predicted mean expression for three example genes spanning the z-score range "
    "(TGFBR3, z = +7.34; RPL13A, z = −0.07; ANK3, z = −13.39). Red = brigatinib (focal "
    "drug); blue = the other four ALK inhibitors; dashed line = mean of the other four. "
    "(B) Per-gene 1-vs-4 z-score: for each gene and focal drug, z = (focal drug's predicted "
    "expression − mean across the other four drugs) / sample SD across the other four drugs "
    "(Bessel n−1 correction). No threshold; signed z is the fGSEA ranking statistic. "
    "(C) fGSEA running enrichment score across genes ranked by z-score (high → low). "
    "Negative NES indicates response eGenes cluster at the brigatinib-suppressed tail of the "
    "ranking. Red ticks below = response-eGene positions in the ranked list; vertical dashed "
    "line = running-ES peak. NES = −1.48, nominal p = 0.003. "
    "(D) −log₁₀(nominal p) for all 21 cell types × 5 ALK inhibitors. Per condition, genes "
    "are ranked by z-score and tested for enrichment of that condition's response eGenes "
    "via fgseaMultilevel (minSize = 5, eps = 0). Response eGenes are defined by lfsr < 0.05 "
    "in the mash joint analysis AND ≥ 3-fold change in mash effect size vs. the matched "
    "control cell type. White outlines mark conditions with nominal p < 0.05 (7 of 105). "
    "No condition survives Benjamini–Hochberg correction across the 105-condition family "
    "(Storey π̂₀ ≈ 0.65)."
)


# ── Shared loaders ────────────────────────────────────────────────────────────
def load_predictions(ct):
    """Load all 5 treatments' predictions as one DataFrame indexed by gene."""
    frames = {}
    for t in TREATMENTS:
        p = os.path.join(PRED_DIR, f"{ct}_{t}_predicted_mean.parquet")
        frames[t] = pd.read_parquet(p)["predicted"]
    return pd.DataFrame(frames)


def load_zscores_for(ct, tr, mash_universe):
    """Per-gene z-scores for one condition, restricted to mash universe."""
    z = pd.read_csv(ZSCORE_PATH, sep="\t")
    sub = z[(z["cell_type"] == ct) & (z["treatment"] == tr)].copy()
    sub = sub[sub["gene"].isin(mash_universe)]
    sub = sub.drop_duplicates("gene")
    return sub


def load_mash_universe():
    mash = pd.read_csv(MASH_PATH, sep="\t", usecols=["gene_snp"])
    return set(mash["gene_snp"].str.split("_").str[0])


# ── Panel A ───────────────────────────────────────────────────────────────────
def draw_panel_A(outer_ax, pred_df, z_for_panelB):
    """Three real example genes, side-by-side mini bar charts."""
    outer_ax.set_axis_off()
    # Three rows: title, legend, bar charts — stacked to avoid overlap.
    inner_gs = gridspec.GridSpecFromSubplotSpec(
        3, 3, subplot_spec=outer_ax.get_subplotspec(),
        height_ratios=[0.14, 0.14, 1.0], wspace=0.70, hspace=0.70,
    )

    # Row 0: title (full width)
    ax_title = outer_ax.figure.add_subplot(inner_gs[0, :])
    ax_title.set_axis_off()
    ax_title.text(0.0, 0.5,
                  "A.  Predicted expression per gene",
                  transform=ax_title.transAxes, ha="left", va="center",
                  fontsize=10, fontweight="bold")

    # Row 1: legend (full width, centered, horizontal)
    ax_leg = outer_ax.figure.add_subplot(inner_gs[1, :])
    ax_leg.set_axis_off()
    legend_patches = [
        mpatches.Patch(color=FOCAL_COLOR, label="Brigatinib (focal)"),
        mpatches.Patch(color=OTHER_COLOR, label="Other 4 drugs"),
    ]
    ax_leg.legend(handles=legend_patches, fontsize=7.5, ncol=2,
                  loc="center", frameon=False)

    for i, gene in enumerate(PANEL_A_GENES):
        ax = outer_ax.figure.add_subplot(inner_gs[2, i])
        vals = pred_df.loc[gene]
        others = [t for t in TREATMENTS if t != EX_TR]
        mu = vals[others].mean()
        sd = vals[others].std(ddof=1)
        z  = (vals[EX_TR] - mu) / sd

        colors = [FOCAL_COLOR if t == EX_TR else OTHER_COLOR for t in TREATMENTS]
        ax.bar(TREATMENTS, [vals[t] for t in TREATMENTS],
               color=colors, edgecolor="white", linewidth=0.6, width=0.7)
        ax.axhline(mu, color="#444", linestyle="--", linewidth=0.9, zorder=4)
        ymax = max(vals) * 1.30
        ax.set_ylim(0, ymax)
        ax.set_title(f"{gene}\n$z$ = {z:+.2f}", fontsize=8, fontweight="bold")
        ax.tick_params(axis="x", labelsize=6.5, rotation=45)
        ax.tick_params(axis="y", labelsize=6.5)
        for tick in ax.get_xticklabels():
            tick.set_horizontalalignment("right")
        if i == 0:
            ax.set_ylabel("Predicted expression\n(logCPM)", fontsize=8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Store the example gene's numbers for panel B (must match panel B + fig 7)
    gene_hi = "TGFBR3"
    vals = pred_df.loc[gene_hi]
    others = [t for t in TREATMENTS if t != EX_TR]
    mu = vals[others].mean()
    sd = vals[others].std(ddof=1)
    z  = (vals[EX_TR] - mu) / sd
    z_for_panelB.update(dict(
        gene=gene_hi, focal=vals[EX_TR], others=vals[others].tolist(),
        mu=mu, sd=sd, z=z,
    ))


# ── Panel B ───────────────────────────────────────────────────────────────────
def draw_panel_B(ax, ex):
    ax.set_axis_off()
    ax.set_title("B.  1-vs-4 z-score",
                 fontsize=10, fontweight="bold", loc="left")

    formula = (r"For each gene $g$, focal drug $f$:" "\n\n"
               r"$z_{g,f}=\dfrac{x_{g,f}\;-\;\overline{x}_{g,\,\mathrm{other\,4}}}"
               r"{\mathrm{SD}\!\left(x_{g,\,\mathrm{other\,4}}\right)}$")
    ax.text(0.02, 0.95, formula, transform=ax.transAxes,
            ha="left", va="top", fontsize=9)

    worked = (rf"Example: {ex['gene']}, focal = Bri" "\n\n"
              rf"$x_{{\mathrm{{Bri}}}} = {ex['focal']:.3f}$" "\n"
              rf"$\overline{{x}}_{{\mathrm{{others}}}} = {ex['mu']:.3f}$" "\n"
              rf"$\mathrm{{SD}}_{{\mathrm{{others}}}} = {ex['sd']:.4f}$" "\n\n"
              rf"$z = \dfrac{{{ex['focal']:.3f} - {ex['mu']:.3f}}}{{{ex['sd']:.4f}}}"
              rf" = {ex['z']:+.2f}$")
    ax.text(0.02, 0.50, worked, transform=ax.transAxes,
            ha="left", va="top", fontsize=8.5, linespacing=1.5,
            bbox=dict(boxstyle="round,pad=0.8", facecolor="#fff8f0",
                      edgecolor="#cccccc", alpha=0.95))


# ── Panel C ───────────────────────────────────────────────────────────────────
def draw_panel_C(outer_ax, zsub, resp_genes, fgsea_row):
    """Running ES curve + response-eGene tick rug below."""
    outer_ax.set_axis_off()
    inner_gs = gridspec.GridSpecFromSubplotSpec(
        3, 1, subplot_spec=outer_ax.get_subplotspec(),
        height_ratios=[0.10, 4, 1], hspace=0.05,
    )
    ax_head = outer_ax.figure.add_subplot(inner_gs[0, 0])
    ax_head.set_axis_off()
    ax_head.text(0.0, 0.5,
                 "C.  fGSEA running ES",
                 transform=ax_head.transAxes, ha="left", va="center",
                 fontsize=10, fontweight="bold")
    ax_es   = outer_ax.figure.add_subplot(inner_gs[1, 0])
    ax_rug  = outer_ax.figure.add_subplot(inner_gs[2, 0], sharex=ax_es)

    sub_sorted = zsub.sort_values("z_score", ascending=False).reset_index(drop=True)
    is_resp = sub_sorted["gene"].isin(resp_genes).values
    N  = len(sub_sorted)
    nh = is_resp.sum()
    p_hit, p_miss = 1.0 / nh, 1.0 / (N - nh)
    incr = np.where(is_resp, p_hit, -p_miss)
    running = np.cumsum(incr)
    peak_idx = int(np.argmin(running))

    ax_es.plot(np.arange(N), running, color=ES_COLOR, linewidth=1.4, zorder=3)
    ax_es.fill_between(np.arange(N), running, 0,
                       where=running < 0, alpha=0.18, color=ES_COLOR)
    ax_es.axhline(0, color="gray", linewidth=0.6, zorder=2)
    ax_es.axvline(peak_idx, color=RESP_COLOR, linestyle="--",
                  linewidth=0.9, alpha=0.7)
    ax_es.set_ylabel("Enrichment score", fontsize=8.5)
    ax_es.tick_params(axis="x", labelbottom=False)
    ax_es.tick_params(axis="y", labelsize=7)
    ax_es.spines["top"].set_visible(False)
    ax_es.spines["right"].set_visible(False)

    hit_positions = np.where(is_resp)[0]
    ax_rug.vlines(hit_positions, 0, 1, color=RESP_COLOR,
                  linewidth=0.5, alpha=0.7)
    ax_rug.set_xlim(0, N)
    ax_rug.set_ylim(0, 1)
    ax_rug.set_yticks([])
    ax_rug.set_xlabel(
        "Genes ranked by z-score (high → low)",
        fontsize=8.5,
    )
    ax_rug.tick_params(axis="x", labelsize=7)
    ax_rug.text(-0.06, 0.5, "Response\neGenes",
                ha="right", va="center", fontsize=7,
                transform=ax_rug.transAxes, color=RESP_COLOR)
    ax_rug.spines["top"].set_visible(False)
    ax_rug.spines["right"].set_visible(False)


# ── Panel E ───────────────────────────────────────────────────────────────────
def draw_panel_E(ax):
    ax.set_axis_off()
    ax.set_title("E.  fGSEA enrichment test",
                 fontsize=10, fontweight="bold", loc="left")

    body = (
        "Per (cell type, treatment):\n"
        "  • Rank genes by z-score\n"
        "  • Gene set = response\n"
        "      eGenes for the condition\n"
        "  • $\\mathrm{fgseaMultilevel}$ →\n"
        "      NES, nominal $p$")
    ax.text(0.02, 0.92, body, transform=ax.transAxes,
            ha="left", va="top", fontsize=9,
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f7f4",
                      edgecolor="#cccccc", alpha=0.95))


# ── Panel D ───────────────────────────────────────────────────────────────────
def draw_panel_D(ax, fgsea):
    """21 x 5 -log10(nominal p) heatmap, nominal-p<0.05 cells outlined."""
    ct_order = sorted(fgsea["cell_type"].unique(), key=lambda x: CT_LABELS.get(x, x))
    mat   = np.zeros((len(ct_order), len(TREATMENTS)))
    nomsig = np.zeros_like(mat, dtype=bool)
    bhsig  = np.zeros_like(mat, dtype=bool)
    for i, ct in enumerate(ct_order):
        for j, tr in enumerate(TREATMENTS):
            row = fgsea[(fgsea["cell_type"] == ct) & (fgsea["treatment"] == tr)]
            if len(row):
                p   = float(row["p_value"].iloc[0])
                adj = float(row["adj_p_value"].iloc[0])
                mat[i, j]    = -np.log10(max(p, 1e-10))
                nomsig[i, j] = p   < 0.05
                bhsig[i, j]  = adj < 0.05

    im = ax.imshow(mat, aspect="auto", cmap="viridis",
                   vmin=0, vmax=max(3.0, mat.max()))
    ax.set_xticks(range(len(TREATMENTS)))
    ax.set_xticklabels([DRUG_NAMES[t] for t in TREATMENTS],
                       rotation=0, fontsize=7)
    ax.set_yticks(range(len(ct_order)))
    ax.set_yticklabels([CT_LABELS.get(c, c) for c in ct_order], fontsize=6.5)
    ax.tick_params(axis="both", length=0)

    # Outline cells: white outline scales with significance level.
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if bhsig[i, j]:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           fill=False, edgecolor="white", linewidth=2.5))
            elif nomsig[i, j]:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1,
                                           fill=False, edgecolor="white", linewidth=1.4))

    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label(r"$-\log_{10}(p_{\mathrm{nominal}})$", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    ax.set_title("D.  fGSEA across 21 cell types × 5 drugs",
                 fontsize=10, fontweight="bold", loc="left", pad=6)

    # Outline legend + test-setup blurb, packed below the heatmap.
    legend_handles = [
        mpatches.Patch(facecolor="#3b3b58", edgecolor="white", linewidth=2.5,
                       label="BH-adjusted $p < 0.05$"),
        mpatches.Patch(facecolor="#3b3b58", edgecolor="white", linewidth=1.4,
                       label="Nominal $p < 0.05$"),
    ]
    ax.legend(handles=legend_handles, fontsize=7, loc="upper left",
              bbox_to_anchor=(0.0, -0.07), ncol=2, framealpha=0.95)

    test_text = (
        "Per condition, genes ranked by z-score; gene set = condition's response eGenes; "
        "enrichment tested with fgseaMultilevel (minSize = 5, eps = 0). "
        "BH applied across all 105 conditions."
    )
    ax.text(0.0, -0.20, test_text, transform=ax.transAxes,
            ha="left", va="top", fontsize=7, wrap=True)


# ── Rationale sketch (separate output file) ───────────────────────────────────
def make_rationale_sketch(pred_df, gene=None):
    """Standalone single-panel sketch: focal vs. Gaussian background of other 4.

    Mirrors the worked example in panel B as a graphic. Saves
    fig7_zscore_rationale.png / .pdf.
    """
    if gene is None:
        gene = "TGFBR3"  # same example as panel B
    vals = pred_df.loc[gene]
    others = [t for t in TREATMENTS if t != EX_TR]
    mu = vals[others].mean()
    sd = vals[others].std(ddof=1)
    z  = (vals[EX_TR] - mu) / sd

    # x-range: extend out to include focal point and ±3 SD around mu
    x_min = min(mu - 3 * sd, vals[EX_TR] - 0.5 * abs(vals[EX_TR] - mu))
    x_max = max(mu + 3 * sd, vals[EX_TR] + 0.5 * abs(vals[EX_TR] - mu))
    pad = 0.05 * (x_max - x_min)
    x_min -= pad
    x_max += pad
    xs = np.linspace(x_min, x_max, 600)
    density = (1.0 / (sd * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((xs - mu) / sd) ** 2)

    fig, ax = plt.subplots(figsize=(11, 4.2))
    ax.fill_between(xs, density, color="#4C72B0", alpha=0.18, zorder=2)
    ax.plot(xs, density, color="#4C72B0", linewidth=1.8, zorder=3,
            label="Background (other 4 treatments)")
    ax.axvline(mu, color="#4C72B0", linestyle=":", linewidth=1.0, alpha=0.7, zorder=2)

    # |z| = 2 cutoff (in original x-axis units)
    for cut in (mu - 2 * sd, mu + 2 * sd):
        ax.axvline(cut, color="#444444", linestyle="--", linewidth=1.0, zorder=2)
    ax.plot([], [], color="#444444", linestyle="--", linewidth=1.0, label="|z| = 2 cutoff")

    # Other 4 treatments as small blue dots just under the x-axis
    y_dot = -0.06 * density.max()
    for t in others:
        ax.plot(vals[t], y_dot, "o", color="#4C72B0", markersize=7,
                markeredgecolor="white", markeredgewidth=0.8, zorder=5)
        ax.annotate(DRUG_NAMES[t],
                    xy=(vals[t], y_dot),
                    xytext=(vals[t], y_dot - 0.18 * density.max()),
                    ha="center", va="top", fontsize=9, color="#4C72B0",
                    rotation=35, rotation_mode="anchor",
                    arrowprops=dict(arrowstyle="-", color="#4C72B0",
                                    lw=0.6, alpha=0.5))

    # Focal treatment as a big red dot, label below
    ax.plot(vals[EX_TR], y_dot, "o", color=FOCAL_COLOR, markersize=14,
            markeredgecolor="white", markeredgewidth=1.2, zorder=6,
            label=f"{DRUG_NAMES[EX_TR]} (target)")
    ax.annotate(DRUG_NAMES[EX_TR],
                xy=(vals[EX_TR], y_dot),
                xytext=(vals[EX_TR], y_dot - 0.20 * density.max()),
                ha="center", va="top", fontsize=10.5, fontweight="bold",
                color=FOCAL_COLOR, rotation=35, rotation_mode="anchor",
                arrowprops=dict(arrowstyle="-", color=FOCAL_COLOR, lw=0.8))

    ax.set_ylabel("Density", fontsize=11)
    ax.set_xlim(x_min, x_max)
    ymin_padded = -0.32 * density.max()
    ax.set_ylim(ymin_padded, density.max() * 1.18)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(fontsize=9.5, loc="upper right", framealpha=0.9)

    title = (f"{DRUG_NAMES[EX_TR]} as target  |  $z$ = {z:+.2f}  |  "
             f"({vals[EX_TR]:.3f} − {mu:.3f}) / {sd:.3f}")
    ax.set_title(title, fontsize=12, pad=10)
    fig.suptitle(
        f"z-score rationale  ·  {gene}  ·  Cardiac mesoderm  ·  "
        f"focal = {DRUG_NAMES[EX_TR]}",
        fontsize=13, fontweight="bold", y=1.00,
    )

    plt.tight_layout()
    png = os.path.join(OUT_DIR, "fig7_zscore_rationale.png")
    pdf = os.path.join(OUT_DIR, "fig7_zscore_rationale.pdf")
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {png}")
    print(f"Saved {pdf}")


# ── Assemble ──────────────────────────────────────────────────────────────────
def main():
    print("Loading data…")
    mash_universe = load_mash_universe()
    pred_df = load_predictions(EX_CT)
    zsub    = load_zscores_for(EX_CT, EX_TR, mash_universe)

    resp = pd.read_csv(RESP_PATH, sep="\t")
    resp_genes = set(resp[(resp["celltype"] == EX_CT) &
                          (resp["treatment"] == EX_TR)]["gene"])
    resp_genes &= mash_universe
    print(f"  zsub: {len(zsub)} genes; response eGenes in cardiac×Bri: {len(resp_genes)}")

    fgsea = pd.read_csv(FGSEA_PATH, sep="\t")
    fgsea_row = fgsea[(fgsea["cell_type"] == EX_CT) &
                      (fgsea["treatment"] == EX_TR)].iloc[0]
    print(f"  fGSEA cardiac×Bri: NES={fgsea_row['NES']:.3f}, "
          f"p={fgsea_row['p_value']:.4f}, BH={fgsea_row['adj_p_value']:.3f}")

    fig = plt.figure(figsize=(11, 8.5))
    fig.suptitle(
        "z-score + fGSEA: are CRISP-predicted drug responses enriched for response eQTLs?",
        fontsize=11, fontweight="bold", y=0.975,
    )
    gs = gridspec.GridSpec(
        2, 3, figure=fig,
        width_ratios=[1.45, 0.85, 1.45], height_ratios=[1.0, 1.10],
        left=0.07, right=0.94, top=0.93, bottom=0.12,
        hspace=0.45, wspace=0.65,
    )

    ax_A = fig.add_subplot(gs[0, 0:2])
    ax_B = fig.add_subplot(gs[0, 2])
    ax_C = fig.add_subplot(gs[1, 0])
    ax_D = fig.add_subplot(gs[1, 1:])

    panelB_ex = {}
    draw_panel_A(ax_A, pred_df, panelB_ex)
    draw_panel_B(ax_B, panelB_ex)
    draw_panel_C(ax_C, zsub, resp_genes, fgsea_row)
    draw_panel_D(ax_D, fgsea)

    # Caption is supplied separately in the manuscript text — not drawn here.
    # The CAPTION string above is kept as a single source of truth for that text.

    png_path = os.path.join(OUT_DIR, "fig6_methods_schematic.png")
    pdf_path = os.path.join(OUT_DIR, "fig6_methods_schematic.pdf")
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Saved {png_path}")
    print(f"Saved {pdf_path}")

    # Also produce the standalone z-score rationale sketch
    make_rationale_sketch(pred_df)


if __name__ == "__main__":
    main()

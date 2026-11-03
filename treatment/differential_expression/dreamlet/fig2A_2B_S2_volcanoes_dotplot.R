# Main Fig 2A: 2A_fig_deg_dotplot_fraction_preLOO_horizontal.pdf
# Main Fig 2B: 2B_fig_volcano_featured_preLOO.pdf (Cardiac mesoderm + Hepatoblasts)
# Supp Fig S2 (5 panels): S2_fig_volcano_supp_allct_preLOO_p{1..5}.png (all 30 cell types)

# Combined figure script for the two main-text DE panels (Fig 2A and 2B) and
# the supplementary volcano gallery (Fig S2). All three read the same dreamlet
# 1v4 DE output, so they share data loading and germ-layer color setup.

# Fig 2A: horizontal dotplot with drugs on the y-axis, 30 cell types on the
# x-axis, dot size = % of tested genes that are DE (FDR < 0.05) in that
# (drug, cell type), and dot fill = signed median log2FC among the DE genes.
# Top marginal bar shows pooled DEG % per cell type; right marginal bar shows
# pooled DEG % per drug.

# Fig 2B: featured volcanoes for two cell types (Cardiac mesoderm cells,
# Hepatoblasts) x 5 drugs, laid out as one row per cell type.

# Fig S2: same volcano panel layout as Fig 2B but for ALL 30 cell types,
# broken across 5 pages so each page fits on letter paper.

# Inputs:
#   results/treatment/dreamlet_1v4_DE/top_all_1v4.csv
#   results/treatment/dreamlet_1v4_DE_summaries/deg_counts_by_assay_drug.csv
#   results/treatment/dreamlet_1v4_DE_summaries/deg_fraction_by_assay_drug.csv

# Outputs (all to results/treatment/dreamlet_1v4_DE_summaries/plots/):
#   deg_dotplot_fraction.{pdf,png}   -- Fig 2A
#   volcano_featured.{pdf,png}       -- Fig 2B
#   volcano_supp_allct.pdf + volcano_supp_allct_p{1..5}.png -- Supp Fig S2

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(scales)

source("/project/gilad/geena/ALK_paper_code/shared/figure_defs.R")

# Parameters
fdr_threshold <- 0.05
drug_order <- c("Ale", "Bri", "Cer", "Cri", "Lor")
drug_colors <- c(Ale = "#0072B2", Bri = "#D55E00", Cer = "#009E73",
                 Cri = "#CC79A7", Lor = "#E69F00")
nonsig_color <- "grey75"

# Volcano layout parameters (used by Fig 2B and Supp S2)
featured_cts <- c("Cardiac mesoderm cells", "Hepatoblasts")
n_label <- 6      # top-N genes to text-label per panel
cts_per_page <- 6 # cell types per supp figure page
strip_h <- 2.5    # inches per cell-type row
strip_w <- 12     # total figure width in inches

# Germ-layer colors: used for the top marginal in Fig 2A and for the row
# label backgrounds in the volcano strips.
germ_colors <- c(
  Epithelium = "#4292C6", Ectoderm = "#F16913",
  Mesoderm = "#44AA99", Endoderm = "#9970AB", iPSC = "#CC6677"
)

# Cell type -> germ layer mapping. Order here defines the cell-type ordering
# used across all three figures (germ layer, then alphabetical within layer).
cell_type_germ_layer <- c(
  "Goblet cells"                              = "Epithelium",
  "Ciliated epithelial cells"                 = "Epithelium",
  "Intestinal epithelial cells"               = "Epithelium",
  "Squamous epithelial cells"                 = "Epithelium",
  "Bronchiolar and alveolar epithelial cells" = "Epithelium",
  "ENS glia"                                  = "Ectoderm",
  "Glioblast"                                 = "Ectoderm",
  "Neuroblast"                                = "Ectoderm",
  "Neuron"                                    = "Ectoderm",
  "Neuronal IPC"                              = "Ectoderm",
  "Radial glia"                               = "Ectoderm",
  "Retinal neurons"                           = "Ectoderm",
  "Retinal pigment cells"                     = "Ectoderm",
  "Retinal progenitors and Muller glia"       = "Ectoderm",
  "Early mesoderm cells"                      = "Mesoderm",
  "Cardiac mesoderm cells"                    = "Mesoderm",
  "Cardiomyocytes"                            = "Mesoderm",
  "Epicardial fat cells"                      = "Mesoderm",
  "Smooth muscle cells"                       = "Mesoderm",
  "Vascular endothelial cells"                = "Mesoderm",
  "Stromal cells"                             = "Mesoderm",
  "Metanephric cells"                         = "Mesoderm",
  "Skeletal muscle cells"                     = "Mesoderm",
  "Lymphoid cells"                            = "Mesoderm",
  "Erythroblasts"                             = "Mesoderm",
  "Megakaryocytes"                            = "Mesoderm",
  "Myeloid cells"                             = "Mesoderm",
  "Hepatoblasts"                              = "Endoderm",
  "Stellate cells"                            = "Endoderm",
  "iPSCs"                                     = "iPSC"
)

# Paths
base_dir <- "/project/gilad/geena/publication_clean"
in_dir <- file.path(base_dir, "results/treatment/dreamlet_1v4_DE")
results_dir <- file.path(base_dir, "results/treatment/dreamlet_1v4_DE_summaries")
plots_dir <- file.path(results_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

ts <- function() format(Sys.time(), "[%H:%M:%S]")

# Helper: save a ggplot as both PDF (vector, cairo for good fonts) and PNG
save_gg <- function(p, name, w, h) {
  ggsave(file.path(plots_dir, paste0(name, ".pdf")), p,
         width = w, height = h, device = cairo_pdf)
  ggsave(file.path(plots_dir, paste0(name, ".png")), p,
         width = w, height = h, dpi = 200)
}

# Load the dreamlet 1v4 DE output. Every row is one gene tested in one cell
# type x drug contrast; columns include logFC, adj.P.Val, and P.Value.
message(ts(), " Loading top_all_1v4.csv...")
top_all <- read.csv(file.path(in_dir, "top_all_1v4.csv"),
                    stringsAsFactors = FALSE)
message(ts(), " ", nrow(top_all), " rows, ",
        n_distinct(top_all$assay), " cell types, ",
        n_distinct(top_all$drug), " drugs")


# Fig 2A: DEG fraction dotplot

# Read the pre-computed per-(cell type, drug) DEG counts and fractions
# (produced by supp_S3_counts_power.R).
message(ts(), " Loading pre-computed count/fraction tables (Fig 2A)...")
counts_by_assay <- read.csv(file.path(results_dir, "deg_counts_by_assay_drug.csv"),
                            stringsAsFactors = FALSE)
fraction_by_assay <- read.csv(file.path(results_dir, "deg_fraction_by_assay_drug.csv"),
                              stringsAsFactors = FALSE)

# Cell type ordering: germ layer, then alphabetical within layer. This is the
# x-axis order for Fig 2A and the row order for Fig 2B / Supp S2.
ct_layer_idx <- match(cell_type_germ_layer, names(germ_colors))
ordered_ct <- names(cell_type_germ_layer)[order(ct_layer_idx,
                                                names(cell_type_germ_layer))]
col_order <- intersect(ordered_ct, unique(counts_by_assay$assay))

# Signed median log2FC per (cell type, drug), computed over the DE subset.
# Drives the dot fill color in the Fig 2A dotplot.
sig <- filter(top_all, adj.P.Val < fdr_threshold)
effsz_signed_sig <- sig %>%
  group_by(assay, drug) %>%
  summarise(median_logFC = median(logFC), .groups = "drop")
color_lim <- max(abs(effsz_signed_sig$median_logFC), na.rm = TRUE)

# Marginal helpers: pool DEG counts across cell types or drugs to produce a
# single "% of tested genes DE" number per group. Used for the top and right
# marginal bars in Fig 2A.
marginal_pooled_pct_by_drug <- function(counts_tbl, ddg_col) {
  counts_tbl %>%
    group_by(drug) %>%
    summarise(total = 100 * sum(.data[[ddg_col]]) / sum(n_tested),
              .groups = "drop")
}
marginal_pooled_pct_by_ct <- function(counts_tbl, ddg_col) {
  counts_tbl %>%
    group_by(assay) %>%
    summarise(total = 100 * sum(.data[[ddg_col]]) / sum(n_tested),
              .groups = "drop")
}
pct_label_fn <- function(x) paste0(round(x, 1), "%")

# Horizontal dotplot builder: takes the metric table (with size_val and
# median_logFC columns), the two marginal tables, and returns a patchwork
# assembly with the main dotplot, top marginal (per-cell-type pooled %), and
# right marginal (per-drug pooled %).
make_dotplot_horizontal <- function(metric_tbl, top_marginal, right_marginal,
                                    size_name, size_breaks, size_trans,
                                    size_label_fn, marginal_label_fn,
                                    marginal_axis_label) {
  # Restrict to the ordered cell types and drugs, drop zero-size dots
  plot_dat <- metric_tbl %>%
    filter(assay %in% col_order, drug %in% drug_order, size_val > 0) %>%
    mutate(assay = factor(assay, levels = col_order),
           drug = factor(drug, levels = rev(drug_order)))

  # Main dot plot: cell type on x, drug on y (Ale on top since we reversed).
  # Dots are filled by signed median log2FC (symmetric divergent scale).
  p_main <- ggplot(plot_dat, aes(x = assay, y = drug)) +
    geom_point(aes(size = size_val, fill = median_logFC),
               shape = 21, color = "grey30", stroke = 0.2) +
    scale_size_continuous(trans = size_trans, range = c(0.6, 7),
                          breaks = size_breaks,
                          labels = size_label_fn(size_breaks),
                          name = size_name) +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         midpoint = 0, limits = c(-color_lim, color_lim),
                         na.value = "grey85",
                         name = expression(median~log[2]~FC)) +
    scale_x_discrete(labels = relabel_celltypes) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 9) +
    theme(
      # Color x-axis tick labels by germ layer so the ordering is visible
      axis.text.x = element_text(
        angle = 90, hjust = 1, vjust = 0.5,
        color = germ_colors[cell_type_germ_layer[levels(plot_dat$assay)]]
      ),
      # Color y-axis tick labels by drug (matches the dot colors readers see)
      axis.text.y = element_text(
        face = "bold",
        color = drug_colors[levels(plot_dat$drug)]
      ),
      panel.grid.major = element_line(color = "grey93", linewidth = 0.3),
      legend.position = "right",
      legend.box = "vertical"
    )

  # Top marginal: per-cell-type pooled DEG %, one bar per cell type, filled by
  # germ layer.
  top_tot <- top_marginal %>%
    filter(assay %in% col_order) %>%
    mutate(assay = factor(assay, levels = col_order),
           germ = cell_type_germ_layer[as.character(assay)])
  p_top <- ggplot(top_tot, aes(x = assay, y = total, fill = germ)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = germ_colors, guide = "none") +
    scale_y_continuous(labels = marginal_label_fn,
                       expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = marginal_axis_label) +
    theme_classic(base_size = 8) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(2, 2, 0, 2))

  # Right marginal: per-drug pooled DEG %, one bar per drug (same order as
  # main y-axis). Small enough to text-label each bar.
  right_tot <- right_marginal %>%
    mutate(drug = factor(drug, levels = rev(drug_order)))
  p_right <- ggplot(right_tot, aes(x = total, y = drug, fill = drug)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = marginal_label_fn(total)), hjust = -0.1, size = 2.6) +
    scale_fill_manual(values = drug_colors, guide = "none") +
    scale_x_continuous(labels = marginal_label_fn,
                       expand = expansion(mult = c(0, 0.22)),
                       position = "top") +
    labs(x = marginal_axis_label, y = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.margin = margin(2, 2, 2, 0))

  # Assemble: top marginal spans across the top of the main dotplot; right
  # marginal to the right of the main dotplot; legends collected on the right.
  (p_top | patchwork::plot_spacer()) /
    (p_main | p_right) +
    plot_layout(heights = c(2, 3), widths = c(6, 1), guides = "collect") &
    theme(legend.position = "right", legend.box = "vertical")
}

# Build and save Fig 2A
p_dot <- make_dotplot_horizontal(
  metric_tbl = fraction_by_assay %>%
    left_join(effsz_signed_sig, by = c("assay", "drug")) %>%
    rename(size_val = pct_sig),
  top_marginal = marginal_pooled_pct_by_ct(fraction_by_assay, "n_sig"),
  right_marginal = marginal_pooled_pct_by_drug(fraction_by_assay, "n_sig"),
  size_name = "DEGs (FDR < 0.05)\n(% of tested)",
  size_breaks = c(0.1, 1, 5, 15, 40),
  size_trans = "sqrt",
  size_label_fn = function(x) paste0(x, "%"),
  marginal_label_fn = pct_label_fn,
  marginal_axis_label = "Pooled %"
)

ggsave(file.path(plots_dir, "deg_dotplot_fraction.pdf"), p_dot,
       width = 14, height = 7, device = cairo_pdf)
ggsave(file.path(plots_dir, "deg_dotplot_fraction.png"), p_dot,
       width = 14, height = 7, dpi = 200)
message(ts(), " Written: deg_dotplot_fraction.{pdf,png}")


# Fig 2B and Supp S2: volcano panels

# Split top_all by (cell type, drug) for O(1) panel lookup. The two-byte
# separator is arbitrary; it just needs to not appear in cell-type or drug
# names.
top_split <- split(top_all, paste(top_all$assay, top_all$drug, sep = "\x01"))

# Cell-type row order for the volcano strips (same as Fig 2A column order).
all_cts <- unique(top_all$assay)
germ_idx <- match(cell_type_germ_layer[all_cts], names(germ_colors))
ct_order <- all_cts[order(germ_idx, all_cts, na.last = TRUE)]

# Per-drug x-axis limit: 99.9th percentile of |logFC| across all cell types.
# Ensures panels of the same drug share an axis so they're visually comparable.
message(ts(), " Computing per-drug x-axis limits (Fig 2B / S2)...")
xlim_per_drug <- setNames(
  lapply(drug_order, function(d) {
    fc <- top_all$logFC[top_all$drug == d]
    ceiling(quantile(abs(fc), 0.999, na.rm = TRUE) * 10) / 10
  }),
  drug_order
)

# Build a single volcano panel: one cell type x one drug. Colors sig genes
# (FDR < 0.05) with the drug color, non-sig grey. Labels the top-N genes by
# raw p-value with gene symbols.
make_volcano <- function(assay_name, drug) {
  key <- paste(assay_name, drug, sep = "\x01")
  df <- top_split[[key]]
  if (is.null(df) || nrow(df) == 0) return(ggplot() + theme_void())

  xlim_val <- xlim_per_drug[[drug]]

  df <- df %>% mutate(sig = adj.P.Val < fdr_threshold)

  labeled_ids <- df %>%
    filter(sig) %>%
    arrange(P.Value) %>%
    slice_head(n = n_label) %>%
    pull(ID)
  df <- df %>% mutate(pt_label = ifelse(ID %in% labeled_ids, ID, NA_character_))

  # Dashed line at the FDR boundary (the p-value of the least-significant DE
  # gene). Skipped when no genes are significant.
  fdr_line <- if (any(df$sig, na.rm = TRUE)) {
    max(df$P.Value[df$sig], na.rm = TRUE)
  } else NA_real_

  ggplot(df, aes(x = logFC, y = -log10(P.Value))) +
    geom_point(data = filter(df, !sig),
               color = nonsig_color, size = 0.4, alpha = 0.5) +
    geom_point(data = filter(df, sig),
               color = drug_colors[[drug]], size = 0.7, alpha = 0.75) +
    { if (!is.na(fdr_line))
        geom_hline(yintercept = -log10(fdr_line),
                   linetype = "dashed", color = "grey40", linewidth = 0.3) } +
    geom_text_repel(
      aes(label = pt_label), size = 2.0, fontface = "italic",
      max.overlaps = Inf, min.segment.length = 0,
      segment.size = 0.25, box.padding = 0.6, point.padding = 0.2,
      force = 6, force_pull = 0.5, max.iter = 20000, seed = 42,
      color = "black", na.rm = TRUE, show.legend = FALSE
    ) +
    scale_x_continuous(limits = c(-xlim_val, xlim_val), oob = scales::squish) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(title = paste0(drug, " vs. others"),
         x = expression(log[2]~"FC"),
         y = expression(-log[10]~italic(P))) +
    theme_classic(base_size = 7) +
    theme(plot.title = element_text(size = 7, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 6),
          axis.text = element_text(size = 5),
          plot.margin = margin(3, 3, 3, 3))
}

# Germ-layer row label used at the left of each volcano strip. The rectangle
# background is colored by germ layer so cell types cluster visually.
make_row_label <- function(label, fill_color = "grey80") {
  wrapped <- paste(strwrap(label, width = 22), collapse = "\n")
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = wrapped,
             angle = 90, size = 2.5, fontface = "bold", hjust = 0.5,
             color = "white", lineheight = 0.9) +
    theme_void() +
    theme(plot.background = element_rect(fill = fill_color, color = NA),
          plot.margin = margin(0, 2, 0, 2))
}

# One volcano-strip row = one cell type across all 5 drugs, with a colored
# germ-layer label at the left.
make_strip <- function(ct) {
  panels <- lapply(drug_order, function(d) make_volcano(ct, d))
  germ <- cell_type_germ_layer[ct]
  color <- if (!is.na(germ)) germ_colors[[germ]] else "grey80"
  make_row_label(relabel_celltypes(ct), fill_color = color) +
    wrap_plots(panels, nrow = 1) +
    plot_layout(widths = c(0.25, 5))
}

# Fig 2B: featured cell types (2 rows).
message(ts(), " Generating featured volcanoes (Fig 2B)...")
strips_featured <- lapply(featured_cts, make_strip)
p_feat <- wrap_plots(strips_featured, ncol = 1)
save_gg(p_feat, "volcano_featured",
        w = strip_w, h = length(featured_cts) * strip_h)
message(ts(), " Written: volcano_featured.{pdf,png}")

# Supp Fig S2: all 30 cell types, paginated so each page holds cts_per_page
# rows. One multi-page PDF plus one PNG per page for easy in-line inclusion.
message(ts(), " Generating supplement volcanoes (Fig S2), ",
        length(ct_order), " cell types...")
n_pages <- ceiling(length(ct_order) / cts_per_page)
pdf_path <- file.path(plots_dir, "volcano_supp_allct.pdf")
cairo_pdf(pdf_path, width = strip_w, height = cts_per_page * strip_h,
          onefile = TRUE)
for (pg in seq_len(n_pages)) {
  idx <- seq((pg - 1) * cts_per_page + 1,
             min(pg * cts_per_page, length(ct_order)))
  ct_pg <- ct_order[idx]
  strips <- lapply(ct_pg, make_strip)
  page_p <- wrap_plots(strips, ncol = 1)
  print(page_p)
  message(ts(), " PDF page ", pg, "/", n_pages, " done")

  # Also save the same page as a standalone PNG (one per page).
  png_path <- file.path(plots_dir, paste0("volcano_supp_allct_p", pg, ".png"))
  ggsave(png_path, page_p,
         width = strip_w, height = length(ct_pg) * strip_h, dpi = 200)
}
dev.off()
message(ts(), " Written: volcano_supp_allct.pdf + per-page PNGs")

message(ts(), " fig2A_2B_S2_volcanoes_dotplot.R complete.")

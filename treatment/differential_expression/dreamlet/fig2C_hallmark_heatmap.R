# Main Fig 2C: 2C_DE1v4_mash_fig_hallmark_heatmap_UNBIASED.pdf

# Build a Hallmark-pathway NES heatmap summarizing the mash-shrunken 1v4 DE
# results across drugs. For each drug we take the top-N Hallmark pathways
# ranked by mean absolute NES across cell types (an objective, data-driven
# selection with no manual curation), take the union across drugs as the
# heatmap rows, and plot the drug-by-pathway mean-NES matrix with the
# rows clustered. Stars mark pathways significant in a large majority of the
# 30 tested cell types (single star: 15/30, double star: 25/30).

# Input:
#   results/treatment/dreamlet_mash/gsea/dreamlet_mash_gsea.tsv

# Output:
#   results/treatment/dreamlet_mash/gsea_plots/hallmark_heatmap.{pdf,png}

library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(tibble)
library(data.table)

# Parameters
fdr_threshold <- 0.05
drug_order <- c("Ale", "Bri", "Cer", "Cri", "Lor")

# Number of Hallmark pathways to take from each drug (per-drug top-N by
# |mean NES|). Their union becomes the heatmap rows.
top_n_hallmark_per_drug <- 8

# Minimum number of cell types where a pathway must be significant before it
# is even considered for selection. Filters out pathways that are noise-only.
min_n_sig_for_selection <- 3

# Cell-type breadth thresholds used to draw stars on top of heatmap cells.
# Interpretation: "sig in majority" (>=15/30 CTs) gets one star; "sig in
# almost all" (>=25/30 CTs) gets two stars.
min_n_sig_for_star <- 15
min_n_sig_for_star_strong <- 25

# Paths
base_dir <- "/project/gilad/geena/publication_clean"
gsea_dir <- file.path(base_dir, "results/treatment/dreamlet_mash/gsea")
plots_dir <- file.path(base_dir, "results/treatment/dreamlet_mash/gsea_plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

ts <- function() format(Sys.time(), "[%H:%M:%S]")

# Load fgsea output (one row per drug x cell type x pathway) and restrict
# to the Hallmark collection.
message(ts(), " Loading dreamlet_mash_gsea.tsv...")
gsea_all <- fread(file.path(gsea_dir, "dreamlet_mash_gsea.tsv"),
                  sep = "\t", data.table = FALSE)

hallmark_res <- filter(gsea_all, collection == "Hallmark")

# Per (drug, pathway) summary: mean and median NES across the 30 cell types,
# and how many cell types are significant. Drives selection + star matrix.
hallmark_summary <- hallmark_res %>%
  group_by(drug, pathway) %>%
  summarise(mean_nes = mean(NES, na.rm = TRUE),
            median_nes = median(NES, na.rm = TRUE),
            n_tested = n(),
            n_sig = sum(padj < fdr_threshold, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(drug = factor(drug, levels = drug_order))

# Objective top-N pathway selection per drug. Pathways must be significant in
# at least min_n_sig_for_selection cell types to qualify; among qualifying
# pathways, keep the top-N by |mean NES|. Union across drugs becomes rows.
message(ts(), " Computing top-N Hallmark selection...")
objective_selection <- hallmark_summary %>%
  filter(n_sig >= min_n_sig_for_selection) %>%
  group_by(drug) %>%
  slice_max(order_by = abs(mean_nes), n = top_n_hallmark_per_drug,
            with_ties = FALSE) %>%
  ungroup()
selected_pathways <- unique(objective_selection$pathway)

# Build the drug-by-pathway mean-NES matrix that goes into the heatmap.
build_hallmark_matrix <- function(pathway_set) {
  mat <- hallmark_summary %>%
    filter(pathway %in% pathway_set) %>%
    select(drug, pathway, mean_nes) %>%
    pivot_wider(names_from = drug, values_from = mean_nes) %>%
    column_to_rownames("pathway") %>%
    as.matrix()
  mat[, drug_order, drop = FALSE]
}
mat_h <- build_hallmark_matrix(selected_pathways)

# Turn HALLMARK_TERM_LIKE_THIS into a readable "Term Like This" row label.
row_labels <- rownames(mat_h) %>%
  sub("^HALLMARK_", "", .) %>%
  gsub("_", " ", .) %>%
  tolower() %>%
  tools::toTitleCase()

# Star matrix: for every displayed pathway in every drug, count cell types
# passing the significance threshold and map to "", "*", or "**". This is
# drawn on top of each heatmap cell via cell_fun below.
star_wide <- hallmark_res %>%
  filter(pathway %in% rownames(mat_h)) %>%
  group_by(drug, pathway) %>%
  summarise(n_sig = sum(padj < fdr_threshold, na.rm = TRUE), .groups = "drop") %>%
  mutate(star = case_when(
    n_sig >= min_n_sig_for_star_strong ~ "**",
    n_sig >= min_n_sig_for_star ~ "*",
    TRUE ~ ""
  )) %>%
  pivot_wider(id_cols = pathway, names_from = drug, values_from = star,
              values_fill = "") %>%
  column_to_rownames("pathway")
star_mat <- as.matrix(star_wide[rownames(mat_h), drug_order])

# Symmetric diverging color scale around 0, so up- and down-regulated pathways
# are equally saturated.
max_val <- max(abs(mat_h), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_val, 0, max_val),
                      c("#D55E00", "white", "#0072B2"))

ht <- Heatmap(
  mat_h, col = col_fun, name = "Mean NES", row_labels = row_labels,
  cluster_rows = TRUE, cluster_columns = FALSE,
  column_names_gp = gpar(fontsize = 9), row_names_gp = gpar(fontsize = 7.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (nchar(star_mat[i, j]) > 0)
      grid.text(star_mat[i, j], x, y,
                gp = gpar(fontsize = 11, fontface = "bold"))
  },
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                              labels_gp = gpar(fontsize = 7)),
  column_title = "Hallmark",
  column_title_gp = gpar(fontsize = 10, fontface = "italic")
)

# Write PDF + PNG side by side (same heatmap, two formats for the manuscript).
for (ext in c("pdf", "png")) {
  f <- file.path(plots_dir, paste0("hallmark_heatmap.", ext))
  if (ext == "pdf") pdf(f, width = 6, height = 9)
  else png(f, width = 6, height = 9, units = "in", res = 200)
  draw(ht, heatmap_legend_side = "right")
  dev.off()
}
message(ts(), " Written: hallmark_heatmap.{pdf,png}")

message(ts(), " fig2C_hallmark_heatmap.R complete.")

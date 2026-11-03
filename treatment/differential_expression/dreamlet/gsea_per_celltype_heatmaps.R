# Supp Fig S5: S5_fig_supp_full_hallmark_ct_heatmap_all_drugs_stacked.pdf
# Supp Fig S7 (5 drugs): S7_fig_supp_top_gobp_ct_heatmap_{Ale,Bri,Cer,Cri,Lor}.pdf
# (The pre-mash GO:BP companion is Supp Fig S6, built by supp_S5fig_S6_premash.R.)

# Two supplementary heatmaps derived from the mash-shrunken 1v4 fgsea results.
# Both use the same 30-cell-type row order (germ layer, then alphabetical
# within layer), the same significance rule (padj < 0.05, marked with "*"),
# and the same signed NES coloring (blue = down, white = 0, red = up).

# Supp Fig S5: all-drug stacked Hallmark heatmap. 5 drug panels x 30 cell
# types = 150 rows; 50 Hallmark pathways as columns. Pathway order is fixed
# across all 5 panels (clustered once on the drug-averaged NES), and the
# color scale is shared across drugs so a reader can scan any pathway column
# top-to-bottom and directly compare drugs.

# Supp Fig S7: one heatmap per drug of top-N GO:BP pathways, split into
# up-regulated (top-N by max positive NES) and down-regulated (top-N by
# min negative NES) blocks. Gives the strongest positive and negative
# GO:BP signals per drug without cherry-picking.

# Inputs:
#   results/treatment/dreamlet_mash/gsea/dreamlet_mash_gsea.tsv
#   reference_sets/celltype_germlayer.csv

# Outputs (all to results/treatment/dreamlet_mash/ct_patterns_plots/):
#   fig_supp_full_hallmark_ct_heatmap_all_drugs_stacked.{pdf,png}
#   fig_supp_top_gobp_ct_heatmap_{Ale,Bri,Cer,Cri,Lor}.{pdf,png}

library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

source("/project/gilad/geena/ALK_paper_code/shared/figure_defs.R")

# Parameters
fdr_threshold <- 0.05
drug_order <- c("Ale", "Bri", "Cer", "Cri", "Lor")
drug_full_name <- c(Ale = "Alectinib", Bri = "Brigatinib",
                    Cer = "Ceritinib", Cri = "Crizotinib",
                    Lor = "Lorlatinib")

# Number of pathways per direction (up + down) per drug for Supp Fig S6.
top_n_gobp_per_direction <- 15L

germ_layer_order <- c("Epithelium", "Ectoderm", "Mesoderm",
                      "Endoderm", "iPSCs")
germ_layer_colors <- c(Epithelium = "#4292C6",
                       Ectoderm = "#F16913",
                       Mesoderm = "#44AA99",
                       Endoderm = "#9970AB",
                       iPSCs = "#CC6677")

# Paths
base_dir <- "/project/gilad/geena/publication_clean"
gsea_dir <- file.path(base_dir, "results/treatment/dreamlet_mash/gsea")
out_dir <- file.path(base_dir, "results/treatment/dreamlet_mash/ct_patterns_plots")
ref_dir <- file.path(base_dir, "reference_sets")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ts <- function() format(Sys.time(), "[%H:%M:%S]")

# Load the mash fgsea output (one row per drug x cell type x pathway) and
# attach germ-layer labels so we can split cell types into germ-layer blocks
# on the heatmaps.
message(ts(), " Loading GSEA + germ-layer references ...")
gsea_all <- fread(file.path(gsea_dir, "dreamlet_mash_gsea.tsv"),
                  sep = "\t", data.table = FALSE)
gl_map <- read_csv(file.path(ref_dir, "celltype_germlayer.csv"),
                   show_col_types = FALSE) %>%
  select(assay = cellid_type, germlayer)

gsea_all <- gsea_all %>%
  inner_join(gl_map, by = "assay") %>%
  mutate(germlayer = factor(germlayer, levels = germ_layer_order),
         drug = factor(drug, levels = drug_order))

hallmark_res <- filter(gsea_all, collection == "Hallmark")
gobp_res <- filter(gsea_all, collection == "GOBP")

# Fixed cell-type row order used across all heatmaps: germ layer block, then
# alphabetical within layer.
celltype_order_df <- gl_map %>%
  filter(assay %in% unique(gsea_all$assay)) %>%
  mutate(germlayer = factor(germlayer, levels = germ_layer_order)) %>%
  arrange(germlayer, assay)
celltype_order <- celltype_order_df$assay
celltype_layers <- setNames(as.character(celltype_order_df$germlayer),
                            celltype_order_df$assay)

# Turn HALLMARK_TERM / GOBP_TERM_LIKE_THIS into a readable "Term Like This"
# column label.
readable_pathway_name <- function(x) {
  sub("^HALLMARK_", "", x) %>% sub("^GOBP_", "", .) %>%
    gsub("_", " ", .) %>% tolower() %>% tools::toTitleCase()
}

# Build a cell-type x pathway NES matrix, ordering rows by celltype_order
# and columns by the pathway_set argument.
build_nes_matrix <- function(df, pathway_set, celltype_order_use = celltype_order) {
  mat <- df %>%
    filter(pathway %in% pathway_set) %>%
    select(assay, pathway, NES) %>%
    pivot_wider(names_from = pathway, values_from = NES) %>%
    column_to_rownames("assay") %>%
    as.matrix()
  mat <- mat[intersect(celltype_order_use, rownames(mat)), , drop = FALSE]
  mat <- mat[, intersect(pathway_set, colnames(mat)), drop = FALSE]
  mat
}

# Same shape as build_nes_matrix but the values are "*" (padj < fdr_threshold)
# or "" (not sig). Drawn as text on top of each heatmap cell.
build_sig_matrix <- function(df, pathway_set, celltype_order_use = celltype_order) {
  sig <- df %>%
    filter(pathway %in% pathway_set) %>%
    mutate(sig = ifelse(padj < fdr_threshold, "*", "")) %>%
    select(assay, pathway, sig) %>%
    pivot_wider(names_from = pathway, values_from = sig, values_fill = "") %>%
    column_to_rownames("assay") %>%
    as.matrix()
  sig <- sig[intersect(celltype_order_use, rownames(sig)), , drop = FALSE]
  sig <- sig[, intersect(pathway_set, colnames(sig)), drop = FALSE]
  sig
}

# Draw a single cell-type x pathway heatmap with a germ-layer row annotation,
# significance stars in each cell, and PDF + PNG output. Used for Supp Fig S6.
draw_ct_heatmap <- function(mat, sig_mat, layers_vec, title, out_stem,
                            width, height, cluster_cols = TRUE,
                            fixed_color_scale = NULL,
                            column_split = NULL) {
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    message(ts(), " skipping ", basename(out_stem), " (empty matrix)")
    return(invisible(NULL))
  }
  # Symmetric diverging color scale around 0; use fixed_color_scale to keep
  # multiple heatmaps directly comparable.
  max_val <- if (is.null(fixed_color_scale)) max(abs(mat), na.rm = TRUE)
             else fixed_color_scale
  if (!is.finite(max_val) || max_val == 0) max_val <- 1
  col_fun <- colorRamp2(c(-max_val, 0, max_val),
                        c("#D55E00", "white", "#0072B2"))
  row_layers <- factor(layers_vec[rownames(mat)], levels = germ_layer_order)

  row_ann <- rowAnnotation(
    `Germ layer` = row_layers,
    col = list(`Germ layer` = germ_layer_colors),
    annotation_name_gp = gpar(fontsize = 8),
    annotation_name_side = "top",
    annotation_legend_param = list(title_gp = gpar(fontsize = 8),
                                   labels_gp = gpar(fontsize = 7))
  )

  col_labels <- readable_pathway_name(colnames(mat))

  ht_args <- list(
    matrix = mat, col = col_fun, name = "NES", na_col = "grey90",
    cluster_rows = FALSE, cluster_columns = cluster_cols,
    row_split = row_layers, row_gap = unit(1.2, "mm"),
    row_title_gp = gpar(fontsize = 8, fontface = "bold"),
    column_labels = col_labels,
    row_labels = relabel_celltypes(rownames(mat)),
    row_names_gp = gpar(fontsize = 6.5),
    column_names_gp = gpar(fontsize = 7),
    column_names_rot = 45,
    left_annotation = row_ann,
    cell_fun = function(j, i, x, y, width, height, fill) {
      s <- sig_mat[i, j]
      if (!is.na(s) && nchar(s) > 0)
        grid.text(s, x, y, gp = gpar(fontsize = 8, fontface = "bold"))
    },
    heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                                labels_gp = gpar(fontsize = 7)),
    column_title = title,
    column_title_gp = gpar(fontsize = 10, fontface = "bold")
  )
  if (!is.null(column_split)) {
    ht_args$column_split <- column_split
    ht_args$cluster_column_slices <- FALSE
    ht_args$column_gap <- unit(3, "mm")
  }
  ht <- do.call(Heatmap, ht_args)

  # Generous bottom/left padding so long rotated column labels don't clip.
  pad <- unit(c(4, 4, 0.5, 0.5), "cm")
  for (ext in c("pdf", "png")) {
    f <- paste0(out_stem, ".", ext)
    if (ext == "pdf") pdf(f, width = width, height = height)
    else png(f, width = width, height = height, units = "in", res = 200)
    draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
         merge_legend = TRUE, padding = pad)
    dev.off()
  }
  message(ts(), " Written: ", basename(out_stem), ".{pdf,png}")
}


# Supp Fig S5: all-drug stacked Hallmark heatmap

# Cluster Hallmark pathways once on the drug-averaged NES matrix so the same
# pathway order is used across all 5 drug panels. Shared color scale so
# blue-vs-red intensity is directly comparable across drugs.
message(ts(), " Building all-drug stacked Hallmark heatmap (Supp Fig S5) ...")

all_hallmark_pathways <- unique(hallmark_res$pathway)
hallmark_scale <- max(abs(hallmark_res$NES), na.rm = TRUE)

avg_nes_mat <- hallmark_res %>%
  group_by(assay, pathway) %>%
  summarise(NES = mean(NES, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = pathway, values_from = NES) %>%
  column_to_rownames("assay") %>%
  as.matrix()
pw_order <- if (ncol(avg_nes_mat) > 1) {
  colnames(avg_nes_mat)[hclust(dist(t(avg_nes_mat)))$order]
} else colnames(avg_nes_mat)

# Stack the 5 per-drug NES matrices vertically, using row names "drug|CT"
# to keep track of which drug each row belongs to. Same for sig stars.
combined_mat_v <- do.call(rbind, lapply(drug_order, function(d) {
  m <- build_nes_matrix(hallmark_res %>% filter(drug == d), pw_order)[, pw_order,
                                                                     drop = FALSE]
  rownames(m) <- paste(d, rownames(m), sep = "|")
  m
}))
combined_sig_v <- do.call(rbind, lapply(drug_order, function(d) {
  s <- build_sig_matrix(hallmark_res %>% filter(drug == d), pw_order)[, pw_order,
                                                                     drop = FALSE]
  rownames(s) <- paste(d, rownames(s), sep = "|")
  s
}))

# Recover drug + cell type from the composite row names for annotation.
row_ct <- sub("^[^|]+\\|", "", rownames(combined_mat_v))
row_drug <- factor(sub("\\|.*$", "", rownames(combined_mat_v)),
                   levels = drug_order)
row_layers_v <- factor(celltype_layers[row_ct], levels = germ_layer_order)

col_fun_v <- colorRamp2(c(-hallmark_scale, 0, hallmark_scale),
                        c("#D55E00", "white", "#0072B2"))
row_ann_v <- rowAnnotation(
  `Germ layer` = row_layers_v,
  col = list(`Germ layer` = germ_layer_colors),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_name_side = "top",
  annotation_legend_param = list(title_gp = gpar(fontsize = 8),
                                 labels_gp = gpar(fontsize = 7))
)

ht_v <- Heatmap(
  matrix = combined_mat_v, col = col_fun_v,
  name = "NES", na_col = "grey90",
  cluster_rows = FALSE, cluster_columns = FALSE,
  cluster_row_slices = FALSE,
  row_split = row_drug,
  row_gap = unit(3, "mm"),
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  column_labels = readable_pathway_name(colnames(combined_mat_v)),
  row_labels = relabel_celltypes(row_ct),
  row_names_gp = gpar(fontsize = 6.5),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  left_annotation = row_ann_v,
  cell_fun = function(j, i, x, y, width, height, fill) {
    s <- combined_sig_v[i, j]
    if (!is.na(s) && nchar(s) > 0)
      grid.text(s, x, y, gp = gpar(fontsize = 8, fontface = "bold"))
  },
  heatmap_legend_param = list(title_gp = gpar(fontsize = 8),
                              labels_gp = gpar(fontsize = 7)),
  column_title = "All drugs: Hallmark NES by cell type (5 drug panels stacked, shared pathway order + color scale)",
  column_title_gp = gpar(fontsize = 10, fontface = "bold")
)

pad_v <- unit(c(4, 4, 0.5, 0.5), "cm")
for (ext in c("pdf", "png")) {
  f <- file.path(out_dir,
                 paste0("fig_supp_full_hallmark_ct_heatmap_all_drugs_stacked.", ext))
  if (ext == "pdf") pdf(f, width = 14, height = 26)
  else png(f, width = 14, height = 26, units = "in", res = 200)
  draw(ht_v, heatmap_legend_side = "right",
       annotation_legend_side = "right",
       merge_legend = TRUE, padding = pad_v)
  dev.off()
}
message(ts(), " Written: fig_supp_full_hallmark_ct_heatmap_all_drugs_stacked.{pdf,png}")


# Supp Fig S6: per-drug top-N GO:BP heatmaps

# For each drug, pick the top-N GO:BP pathways going UP (by max positive
# NES among cell types) and the top-N going DOWN (by min negative NES).
# Require at least one significant cell type in the same direction so we
# aren't picking noise. Split columns into Up / Down blocks in the heatmap.
message(ts(), " Building per-drug top-", top_n_gobp_per_direction,
        " GO:BP heatmaps (Supp Fig S6) ...")

gobp_dir_summary <- gobp_res %>%
  group_by(drug, pathway) %>%
  summarise(max_NES = max(NES, na.rm = TRUE),
            min_NES = min(NES, na.rm = TRUE),
            n_sig_up = sum(padj < fdr_threshold & NES > 0, na.rm = TRUE),
            n_sig_dn = sum(padj < fdr_threshold & NES < 0, na.rm = TRUE),
            .groups = "drop")

for (d in drug_order) {
  # Top-N pathways going up: sort by max positive NES; require sig in >=1 CT.
  top_up <- gobp_dir_summary %>%
    filter(drug == d, n_sig_up >= 1L, max_NES > 0) %>%
    slice_max(order_by = max_NES,
              n = top_n_gobp_per_direction, with_ties = FALSE) %>%
    pull(pathway)
  # Top-N pathways going down: sort by min negative NES.
  top_dn <- gobp_dir_summary %>%
    filter(drug == d, n_sig_dn >= 1L, min_NES < 0) %>%
    slice_min(order_by = min_NES,
              n = top_n_gobp_per_direction, with_ties = FALSE) %>%
    pull(pathway)

  if (length(top_up) == 0 && length(top_dn) == 0) {
    message(ts(), " no directional top GO:BP for ", d, "; skipping")
    next
  }

  top_paths <- c(top_up, top_dn)
  col_split <- factor(c(rep("Up (top by max NES)", length(top_up)),
                        rep("Down (top by min NES)", length(top_dn))),
                      levels = c("Up (top by max NES)",
                                 "Down (top by min NES)"))

  # Build the matrix by hand rather than via build_nes_matrix so the top_up
  # and top_dn column blocks stay in order (build_nes_matrix reorders by set
  # intersection which would lose the up/down grouping).
  sub_g <- gobp_res %>% filter(drug == d)
  mat <- sub_g %>% filter(pathway %in% top_paths) %>%
    select(assay, pathway, NES) %>%
    pivot_wider(names_from = pathway, values_from = NES) %>%
    column_to_rownames("assay") %>% as.matrix()
  mat <- mat[intersect(celltype_order, rownames(mat)), top_paths, drop = FALSE]
  sig <- sub_g %>% filter(pathway %in% top_paths) %>%
    mutate(sig = ifelse(padj < fdr_threshold, "*", "")) %>%
    select(assay, pathway, sig) %>%
    pivot_wider(names_from = pathway, values_from = sig, values_fill = "") %>%
    column_to_rownames("assay") %>% as.matrix()
  sig <- sig[intersect(celltype_order, rownames(sig)), top_paths, drop = FALSE]

  draw_ct_heatmap(
    mat, sig, celltype_layers,
    title = paste0(drug_full_name[d], " (", d,
                   "): top-", top_n_gobp_per_direction,
                   " up + top-", top_n_gobp_per_direction,
                   " down GO:BP by cell type"),
    out_stem = file.path(out_dir, paste0("fig_supp_top_gobp_ct_heatmap_", d)),
    width = 16, height = 14,
    # Keep up/down blocks intact and sorted within block; don't cluster cols.
    cluster_cols = FALSE,
    column_split = col_split
  )
}

message(ts(), " supp_S5fig_S6_gsea_ct_patterns.R complete.")

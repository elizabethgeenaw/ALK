#!/usr/bin/env Rscript
# Companion to Fig 3B: supplementary z-scored logFC heatmap (row z-scored
# relative to the 5-drug row) for lorlatinib-specific motif genes,
# partitioned into two blocks per cell type:
#   1. Shared: genes in BOTH the hepatoblast Lor motif (n=671) and the
#      cardiac-mesoderm Lor motif (n=2225).
#   2. Unique to that cell type.
# Rows are hierarchically clustered *within* each block so the block boundary
# is preserved.
#
# This panel was previously main Fig 3B but was demoted because the z-score
# is true by construction: genes were selected for the heatmap because they
# differ in Lor, so the row ordering carries no independent evidence for the
# biological claim. Main Fig 3B is now the drug x cell-type occurrence matrix
# (`fig3B_drug_celltype_matrix.pdf`, from cormotif_figures.R).
#
# Inputs:
#   Lor-specific motif gene lists:
#     Hepatoblasts_motif5_genes_K2plus.txt              (n=671)
#     Cardiac mesoderm cells_motif3_genes_K2plus.txt    (n=2225)
#   DE table with per-(drug, cell type) logFC:
#     top_all_1v4.csv
#
# Outputs (all under COR_DIR):
#   figS_lor_motif_heatmap_partitioned.{pdf,png}
#   figS_lor_motif_heatmap_partitioned.rds
#   figS_lor_motif_heatmap_partitioned_data.csv
#
# Run on gilad-hm with R/4.4.1.

library(dplyr)
library(readr)
library(tidyr)
library(pheatmap)
library(grid)

# Paths.
base_dir <- "/project/gilad/geena/publication_clean"
cor_dir  <- file.path(base_dir,
  "results/treatment/cormotif/proteincoding_newfilter_15cells_4arm")
hep_list <- file.path(cor_dir,
  "Hepatoblasts/Hepatoblasts_motif5_genes_K2plus.txt")
cm_list  <- file.path(cor_dir,
  "Cardiac mesoderm cells/Cardiac mesoderm cells_motif3_genes_K2plus.txt")
de_tab   <- file.path(base_dir,
  "results/treatment/dreamlet_clean_proteincoding/top_all_1v4.csv")
s2_xlsx  <- file.path(cor_dir, "S2_gene_lists.xlsx")

out_dir <- cor_dir

# Manuscript-reported values (used as consistency checks below).
drug_order      <- c("Ale", "Bri", "Cer", "Cri", "Lor")
expected_shared <- 173

# Load motif gene lists and cross-check against manuscript counts.
hep_genes <- trimws(readLines(hep_list))
hep_genes <- hep_genes[nchar(hep_genes) > 0]
cm_genes  <- trimws(readLines(cm_list))
cm_genes  <- cm_genes[nchar(cm_genes) > 0]

cat("Hepatoblasts Lor-motif genes: ", length(hep_genes),
    "  (expected 671)\n", sep = "")
cat("Cardiac mesoderm Lor-motif genes: ", length(cm_genes),
    " (expected 2225)\n", sep = "")

shared_genes <- intersect(hep_genes, cm_genes)
cat("Shared (intersection): ", length(shared_genes),
    "  (manuscript claim: ", expected_shared, ")\n", sep = "")

# Hard-fail if the shared count disagrees with the manuscript. Silently
# fabricating a heatmap that contradicts a printed number would be worse
# than breaking the pipeline.
if (length(shared_genes) != expected_shared) {
  cat("MOTIF-OVERLAP-MISMATCH: got ", length(shared_genes),
      ", expected ", expected_shared, "\n", sep = "")
  cat("Stopping - do not fabricate the figure with a mismatched count.\n")
  quit(status = 1)
}

# Optional cross-check against S2 workbook (sheet names only).
if (requireNamespace("readxl", quietly = TRUE) && file.exists(s2_xlsx)) {
  sheets <- readxl::excel_sheets(s2_xlsx)
  cat("S2_gene_lists.xlsx sheets: ", paste(sheets, collapse = ", "),
      "\n", sep = "")
} else {
  cat("(readxl not available or S2 workbook missing - skipping cross-check.)\n")
}

# Load DE table and build per-cell-type gene x drug logFC matrices. Rows with
# any NA drug are dropped because row-wise z-scoring requires a defined SD.
de <- read_csv(de_tab, show_col_types = FALSE) %>%
  filter(drug %in% drug_order) %>%
  select(assay, ID, drug, logFC)

build_mat <- function(ct, gene_set) {
  m <- de %>%
    filter(assay == ct, ID %in% gene_set) %>%
    pivot_wider(names_from = drug, values_from = logFC,
                values_fill = NA_real_) %>%
    as.data.frame()
  rownames(m) <- m$ID; m$ID <- NULL; m$assay <- NULL
  m <- as.matrix(m[, intersect(drug_order, colnames(m)), drop = FALSE])
  m <- m[complete.cases(m), , drop = FALSE]
  m
}

mat_hep <- build_mat("Hepatoblasts", hep_genes)
mat_cm  <- build_mat("Cardiac mesoderm cells", cm_genes)

cat("\nHepatoblasts matrix (after dropping NA rows): ",
    nrow(mat_hep), " x ", ncol(mat_hep), "\n", sep = "")
cat("Cardiac mesoderm matrix (after dropping NA rows): ",
    nrow(mat_cm),  " x ", ncol(mat_cm),  "\n", sep = "")

# Row z-score across the 5 drugs. If a row has zero SD (perfectly flat)
# we set the divisor to 1 rather than propagating NaN.
zscore_rows <- function(m) {
  m_c <- m - rowMeans(m)
  sds <- apply(m, 1, sd)
  m_c / ifelse(sds == 0, 1, sds)
}

z_hep <- zscore_rows(mat_hep)
z_cm  <- zscore_rows(mat_cm)

# Partition rows into "shared" vs "unique" blocks and hierarchically cluster
# within each block. Row order is shared-cluster then unique-cluster, which
# preserves the block boundary as a visible transition in the heatmap.
cluster_within <- function(z_mat, gene_set, shared_set) {
  present   <- rownames(z_mat)
  in_shared <- present %in% shared_set
  in_unique <- !in_shared

  z_shared <- z_mat[in_shared, , drop = FALSE]
  z_unique <- z_mat[in_unique, , drop = FALSE]

  order_within <- function(m) {
    if (nrow(m) <= 1) return(rownames(m))
    d <- dist(m, method = "euclidean")
    h <- hclust(d, method = "average")
    rownames(m)[h$order]
  }
  ord_shared  <- order_within(z_shared)
  ord_unique  <- order_within(z_unique)
  final_order <- c(ord_shared, ord_unique)

  list(z_ord   = z_mat[final_order, , drop = FALSE],
       shared  = ord_shared,
       unique_ = ord_unique)
}

hep_p <- cluster_within(z_hep, hep_genes, shared_genes)
cm_p  <- cluster_within(z_cm,  cm_genes,  shared_genes)

cat("\nHepatoblasts blocks: shared=", length(hep_p$shared),
    "  unique=", length(hep_p$unique_), "\n", sep = "")
cat("Cardiac mesoderm blocks: shared=", length(cm_p$shared),
    "  unique=", length(cm_p$unique_), "\n", sep = "")

# Row annotation column identifying shared vs unique block.
make_ann <- function(shared_ids, unique_ids) {
  ann <- data.frame(
    block = c(rep("shared", length(shared_ids)),
              rep("unique", length(unique_ids))),
    row.names = c(shared_ids, unique_ids),
    stringsAsFactors = FALSE)
  ann$block <- factor(ann$block, levels = c("shared", "unique"))
  ann
}
ann_hep <- make_ann(hep_p$shared, hep_p$unique_)
ann_cm  <- make_ann(cm_p$shared,  cm_p$unique_)

# Shared block emphasised so the overlap reads at a glance.
ann_colors <- list(block = c(shared = "#F1B430", unique = "#BEBEBE"))

# Shared color scale across the two heatmaps for direct visual comparison.
# Cap at 3 to prevent extreme outliers from dominating the palette.
abs_max <- max(abs(c(hep_p$z_ord, cm_p$z_ord)), na.rm = TRUE)
abs_max <- min(abs_max, 3)
breaks  <- seq(-abs_max, abs_max, length.out = 101)
colors  <- colorRampPalette(c("#0072B2", "white", "#D55E00"))(100)

# Save one heatmap per cell type. Height scales with row count so the tall
# CM heatmap doesn't crush its labels; capped at 20 in to keep the file usable.
save_hm <- function(z_mat, ann, name) {
  n_rows <- nrow(z_mat)
  h <- max(4, min(20, 0.006 * n_rows + 4))  # 671 -> ~8, 2225 -> ~17

  ph <- pheatmap(
    z_mat,
    cluster_rows       = FALSE,   # order already fixed within blocks above
    cluster_cols       = FALSE,
    annotation_row     = ann,
    annotation_colors  = ann_colors,
    annotation_names_row = FALSE,
    show_rownames      = FALSE,
    show_colnames      = TRUE,
    color              = colors,
    breaks             = breaks,
    border_color       = NA,
    fontsize           = 8,
    fontsize_col       = 10,
    main               = paste0(name, " - Lor-specific motif genes (z-scored logFC)"),
    silent             = TRUE)

  for (ext in c("pdf", "png")) {
    f <- file.path(out_dir,
                   paste0("figS_lor_motif_heatmap_",
                          gsub("[ ]+", "_", name), ".", ext))
    if (ext == "pdf") pdf(f, width = 5, height = h)
    else              png(f, width = 5, height = h, units = "in",
                          res = 300, type = "cairo")
    grid.newpage(); grid.draw(ph$gtable); dev.off()
    cat("  wrote: ", f, "\n", sep = "")
  }
  ph
}

cat("\nRendering heatmaps...\n")
ph_hep <- save_hm(hep_p$z_ord, ann_hep, "Hepatoblasts")
ph_cm  <- save_hm(cm_p$z_ord,  ann_cm,  "Cardiac mesoderm cells")

# Bundled RDS + long-format row-order CSV for downstream inspection.
saveRDS(list(hep = hep_p, cm = cm_p, shared = shared_genes,
             ann_hep = ann_hep, ann_cm = ann_cm),
        file.path(out_dir, "figS_lor_motif_heatmap_partitioned.rds"))

row_meta <- bind_rows(
  tibble::tibble(cell_type = "Hepatoblasts",
                 gene      = c(hep_p$shared, hep_p$unique_),
                 block     = c(rep("shared", length(hep_p$shared)),
                               rep("unique", length(hep_p$unique_))),
                 row_order = seq_len(nrow(hep_p$z_ord))),
  tibble::tibble(cell_type = "Cardiac mesoderm cells",
                 gene      = c(cm_p$shared, cm_p$unique_),
                 block     = c(rep("shared", length(cm_p$shared)),
                               rep("unique", length(cm_p$unique_))),
                 row_order = seq_len(nrow(cm_p$z_ord))))
write_csv(row_meta,
          file.path(out_dir, "figS_lor_motif_heatmap_partitioned_data.csv"))

cat("\nDone.\n")

# Export all control cells from Seurat RDS to MatrixMarket format 
# https://jr-leary7.github.io/quarto-site/tutorials/Seurat_AnnData_Conversion.html
# Running CRISP in python but had some issues with reticulate package so doing this instead
# One directory per cell type
# Need to run with R/4.2.0 (old node) bc R/4.4.1 doesn't have BPCells

library(Seurat)
library(Matrix)
library(BPCells)

control_rds <- "/project/gilad/geena/publication_clean/results/control/celltype_annotation/combined_control_filtnorm_cellid.rds"
out_dir <- "/project/gilad/geena/chem_perturbation_clean/data/control_h5ad"
celltype_col <- "cellid_type"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

valid_cts <- c(
  "Bronchiolar and alveolar epithelial cells",
  "Cardiac mesoderm cells",
  "Ciliated epithelial cells",
  "Early mesoderm cells",
  "ENS glia",
  "Epicardial fat cells",
  "Glioblast",
  "Goblet cells",
  "Hepatoblasts",
  "Metanephric cells",
  "Neuroblast",
  "Neuron",
  "Neuronal IPC",
  "Radial glia",
  "Retinal neurons",
  "Retinal pigment cells",
  "Retinal progenitors and Muller glia",
  "Squamous epithelial cells",
  "Stromal cells",
  "Vascular endothelial cells",
  "iPSCs"
)
cat("Valid cell types:", length(valid_cts), "\n")


cat("Loading", control_rds, "...\n")
obj <- readRDS(control_rds)
cat("Loaded:", ncol(obj), "cells x", nrow(obj), "genes\n")

counts <- GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
meta <- obj@meta.data
ct_vec <- as.character(meta[[celltype_col]])

for (ct in valid_cts) {
  print(ct)
  # Make directory
  ct_safe <- gsub("[^A-Za-z0-9_-]", "_", ct)
  ct_dir <- file.path(out_dir, ct_safe)
  dir.create(ct_dir, showWarnings = FALSE)

  # Skip if already exported
  if (file.exists(file.path(ct_dir, "barcodes.tsv"))) {
    n_existing <- length(readLines(file.path(ct_dir, "barcodes.tsv")))
    cat(sprintf("  %s -> already exported (%d cells), skipping\n", ct, n_existing))
    next
  }

  # Find cells belonging to this cell type
  idx <- which(ct_vec == ct)
  if (length(idx) == 0) {
    cat(sprintf("  %s -> 0 cells, skipping\n", ct))
    next
  }

  # Subset counts and metadata; convert to sparse matrix for writeMM
  sub_counts <- counts[, idx, drop = FALSE]
  sub_counts <- as(sub_counts, "dgCMatrix")
  sub_meta <- meta[idx, , drop = FALSE]
  sub_meta$original_celltype <- ct
  sub_meta$cellid_type <- ct

  # Write MatrixMarket counts, cell barcodes, gene names, and metadata
  writeMM(sub_counts, file.path(ct_dir, "counts.mtx"))
  writeLines(colnames(sub_counts), file.path(ct_dir, "barcodes.tsv"))
  writeLines(rownames(sub_counts), file.path(ct_dir, "genes.tsv"))
  write.table(sub_meta, file.path(ct_dir, "obs.tsv"),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  cat(sprintf("  %s -> %d cells written to %s\n",
              ct, ncol(sub_counts), ct_dir))
}

cat("Done.\n")

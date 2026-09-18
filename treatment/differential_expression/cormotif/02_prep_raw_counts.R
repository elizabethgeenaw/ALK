# Get per-drug raw pseudobulk count matrices from the dreamlet pseudobulk SCE
# Writes them as per-drug RDS files for load_and_merge_pseudobulk() helper

library(dreamlet)
library(SingleCellExperiment)

pb_treat_path <- "/project/gilad/geena/publication_clean/results/treatment/dreamlet_clean/pb_treat"
arm_csv <- "/project/gilad/geena/publication_clean/results/treatment/cormotif/arm_counts_cormotif.csv"
out_dir <- "/project/gilad/geena/publication_clean/results/treatment/cormotif/pseudobulk_15cells_4indiv"
drugs <- c("Ale", "Bri", "Cer", "Cri", "Lor")
min_cells <- 15
min_indiv <- 4

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load pseudobulk SCE
# pb_treat stores each cell type as a separate assay (assayNames() returns cell type names)
# assay(x, ct) is the genes x samples count matrix for that cell type
cat("Loading pb_treat ...\n")
pb_treat <- readRDS(pb_treat_path)
cat("Class:", class(pb_treat), "\n")
cat("N assays (cell types):", length(assayNames(pb_treat)), "\n")
cat("N samples:", ncol(pb_treat), "\n\n")

cc  <- cellCounts(pb_treat)
arm <- substr(rownames(cc), 1, 3)
ind <- sub(paste0("^(", paste(drugs, collapse = "|"), ")"), "", rownames(cc))

# Exclude "unassigned" because motifs derived from unlabeled cells are not
# biologically interpretable.
arm_counts <- read.csv(arm_csv, stringsAsFactors = FALSE)
exclude_cts  <- "unassigned"

keep_cts <- arm_counts$cell_type[
  apply(arm_counts[, drugs, drop = FALSE], 1, min) >= min_indiv &
  !arm_counts$cell_type %in% exclude_cts
]
cat("Cell types to extract:", length(keep_cts), "\n")

cat("Assay names (first 5):", paste(head(assayNames(pb_treat), 5), collapse = ", "), "\n\n")

# Sanity check: any eligible cell type that's not actually an assay in the
# SCE gets dropped with a warning rather than silently failing.
missing_cts <- setdiff(keep_cts, assayNames(pb_treat))
if (length(missing_cts) > 0) {
  cat("WARNING: these cell types are in arm_counts but missing from pb_treat assays:\n")
  cat(paste0("  ", missing_cts, "\n"))
  keep_cts <- intersect(keep_cts, assayNames(pb_treat))
  cat("Proceeding with", length(keep_cts), "cell types.\n")
}

# Extract per-drug count matrices. Structure: drug_data[[drug]][[cell_type]]
# is a genes x individuals matrix containing only the samples that pass the
# >= min_cells threshold for that (drug, cell type) arm.
drug_data <- setNames(vector("list", length(drugs)), drugs)
for (d in drugs) drug_data[[d]] <- list()

for (ct in keep_cts) {
  cat("Extracting:", ct, "\n")

  cnt_mat   <- assay(pb_treat, ct)
  smp_names <- colnames(cnt_mat)

  for (d in drugs) {
    # Pseudobulk samples for this drug arm that meet the cell count threshold.
    # intersect() with actual column names just to double check
    pass_smp <- rownames(cc)[arm == d & cc[, ct] >= min_cells]
    pass_smp <- intersect(pass_smp, smp_names)

    if (length(pass_smp) == 0) {
      cat("  ", d, ": 0 samples with >=", min_cells, "cells - skipping arm\n")
      next
    }

    mat <- cnt_mat[, pass_smp, drop = FALSE]
    colnames(mat) <- sub(paste0("^", d), "", colnames(mat))

    drug_data[[d]][[ct]] <- mat
    cat("  ", d, ":", ncol(mat), "individuals,", nrow(mat), "genes\n")
  }
}

# Write one RDS per drug
for (d in drugs) {
  out_path <- file.path(out_dir, paste0("pseudo_bulk_", d, ".rds"))
  saveRDS(drug_data[[d]], out_path)
  cat("Saved", d, "->", out_path,
      "(", length(drug_data[[d]]), "cell types)\n")
}

# Verification summary: cell-type count and per-cell-type sample counts.
cat("\n--- Verification ---\n")
for (d in drugs) {
  cts_d <- names(drug_data[[d]])
  cat(d, ": ", length(cts_d), " cell types | sample counts: ",
      paste(sapply(cts_d, function(ct) ncol(drug_data[[d]][[ct]])), collapse = ", "),
      "\n", sep = "")
}

cat("\nDone. Outputs in:", out_dir, "\n")
cat("Next: sbatch cormotif_new_filters.sh\n")

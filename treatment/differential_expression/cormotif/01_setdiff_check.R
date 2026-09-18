# Upstream of Fig 3, Supp Fig S7/S8, Supp Tables S10-S14: decides which cell
# types are eligible for the Cormotif refit. Cormotif fits its own per-arm
# regressions independently of the eQTL pipeline, so the >=20-individuals gate
# used for eQTL is too strict here. This script relaxes the gate to >=4
# individuals per drug arm (matching dreamlet's own min.samples default), then
# reports which cell types are gained relative to the current cormotif output.
#
# Output: arm_counts_cormotif.csv, used by 02_prep_raw_counts.R to decide
# which cell types to extract.

library(dreamlet)
library(dplyr)

pb_treat_path <- "/project/gilad/geena/publication_clean/results/treatment/dreamlet_clean/pb_treat"
cormotif_dir  <- "/project/gilad/geena/publication_clean/results/treatment/cormotif/protein_coding_only_15cells"
out_dir       <- "/project/gilad/geena/publication_clean/results/treatment/cormotif"
drugs         <- c("Ale", "Bri", "Cer", "Cri", "Lor")
min_cells     <- 15L
min_indiv     <- 4L

cat("Loading pb_treat ...\n")
pb_treat <- readRDS(pb_treat_path)
cat("Class:", class(pb_treat), "\n")

# cellCounts returns samples x cell_types; each entry is the number of cells
# from that pseudobulk sample assigned to that cell type. Sample names are
# formatted "<drug><individual_id>", so we split off the first 3 chars for arm.
cc  <- cellCounts(pb_treat)
arm <- substr(rownames(cc), 1, 3)
ind <- sub(paste0("^(", paste(drugs, collapse = "|"), ")"), "", rownames(cc))

cat("Samples in pb_treat:", nrow(cc), "\n")
cat("Cell types in pb_treat:", ncol(cc), "\n")

# Per-arm counts: for each cell type, how many individuals in each drug arm
# have >= min_cells cells of that type. This is the actual filter used to
# decide whether the arm is fittable.
arm_counts <- do.call(rbind, lapply(colnames(cc), function(ct) {
  keep <- cc[, ct] >= min_cells
  tab  <- table(factor(arm[keep], levels = drugs))
  data.frame(cell_type = ct, t(as.integer(tab)), check.names = FALSE,
             stringsAsFactors = FALSE)
}))
colnames(arm_counts) <- c("cell_type", drugs)

write.csv(arm_counts, file.path(out_dir, "arm_counts_cormotif.csv"), row.names = FALSE)
cat("Saved arm_counts_cormotif.csv\n")

# Pairwise intersection counts (diagnostic): how many individuals contribute
# to both arms of each pairwise comparison with >= min_cells cells.
comp_pairs <- combn(drugs, 2, simplify = FALSE)

pair_counts <- do.call(rbind, lapply(colnames(cc), function(ct) {
  row <- list(cell_type = ct)
  for (p in comp_pairs) {
    d1 <- p[[1]]; d2 <- p[[2]]
    pass1 <- ind[arm == d1 & cc[, ct] >= min_cells]
    pass2 <- ind[arm == d2 & cc[, ct] >= min_cells]
    row[[paste0(d1, "_", d2)]] <- length(intersect(pass1, pass2))
  }
  as.data.frame(row, stringsAsFactors = FALSE)
}))

write.csv(pair_counts, file.path(out_dir, "pairwise_counts_cormotif.csv"), row.names = FALSE)
cat("Saved pairwise_counts_cormotif.csv\n")

# Apply the per-arm gate: keep cell types where ALL 5 arms have >= min_indiv
# individuals with >= min_cells cells. Then compare to the existing cormotif
# output directory to see which cell types are added / dropped.
keep_cts <- arm_counts$cell_type[
  apply(arm_counts[, drugs, drop = FALSE], 1, min) >= min_indiv
]

current_cts <- list.dirs(cormotif_dir, recursive = FALSE, full.names = FALSE)

added   <- setdiff(keep_cts, current_cts)
dropped <- setdiff(current_cts, keep_cts)

# Marginal cell types: pass the gate but with any arm at only 4-6 individuals.
# Candidates for a seed-stability check because low-N fits are more sensitive
# to random initialisation.
keep_arm <- arm_counts[arm_counts$cell_type %in% keep_cts, drugs, drop = FALSE]
arm_min  <- apply(keep_arm, 1, min)
marginal <- arm_counts$cell_type[arm_counts$cell_type %in% keep_cts][arm_min <= 6]

# Assemble a human-readable markdown report so the decision can be reviewed
# before launching the long-running cormotif refit.
fmt_list <- function(x) if (length(x) == 0) "  (none)" else paste0("  - ", x, collapse = "\n")

report <- c(
  "# Cormotif Cell-Type Set Difference",
  paste0("Date: ", Sys.Date()),
  "",
  paste0("Filter: all 5 per-arm counts >= ", min_indiv,
         " individuals with >= ", min_cells, " cells"),
  paste0("Cormotif-eligible set size:  ", length(keep_cts)),
  paste0("Current Cormotif set size:   ", length(current_cts),
         "  (eQTL-gated >= 20/arm)"),
  "",
  "ADDED vs current (switching to >= 4/arm threshold):",
  fmt_list(added),
  "",
  "DROPPED vs current (should be none; >= 20 subset of >= 4):",
  fmt_list(dropped),
  "",
  "Marginal cell types (min per-arm 4-6):",
  fmt_list(marginal),
  "",
  "Pairwise intersection counts saved to: pairwise_counts_cormotif.csv",
  "",
  "--- Arm counts for added cell types ---"
)

if (length(added) > 0) {
  added_rows <- arm_counts[arm_counts$cell_type %in% added, , drop = FALSE]
  for (i in seq_len(nrow(added_rows))) {
    report <- c(report, paste0(
      "  ", added_rows$cell_type[i], ": ",
      paste(drugs, as.integer(added_rows[i, drugs]), sep = "=", collapse = ", ")
    ))
  }
}

writeLines(report, file.path(out_dir, "celltype_set_diff.md"))
cat("Saved celltype_set_diff.md\n")

cat("\n--- SUMMARY ---\n")
cat("Cell types added:  ", length(added),   "\n")
cat("Cell types dropped:", length(dropped), "\n")
cat("Marginal:          ", length(marginal), "\n")
if (length(added)   > 0) { cat("Added:\n");   cat(paste0("  ", added,   "\n")) }
if (length(dropped) > 0) { cat("Dropped:\n"); cat(paste0("  ", dropped, "\n")) }

cat("\n=== Review celltype_set_diff.md before proceeding to 02_prep_raw_counts.R ===\n")

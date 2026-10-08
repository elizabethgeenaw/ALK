#!/usr/bin/env Rscript

# 1v4 differential expression run for the ALK inhibitor paper.
# For each drug in turn, we contrast it against the mean of the other four
# ALK inhibitors ("one-vs-others"), fit a dreamlet model with individual
# as random effect, and write out topTable results per drug plus a combined
# top_all_1v4.csv used by nearly every downstream DE analysis. Every gene
# retained by processAssays() is tested (no biotype / MT- / RP[SL] filter).
# Cell types with any drug arm below the min-individual threshold are dropped.
# Also writes a per-arm QC table (DE_QC_powertier.csv) and the subset of
# rows that fall below the exclusion threshold (DE_excluded_conditions.csv).
# Writes per-cell-type gene universes (universes/<cell type>.txt) for Cormotif.

library(BPCells)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(dreamlet)
library(zenith)
library(variancePartition)
options(future.globals.maxSize = 50e+09)

packageVersion("dreamlet")

# Parameters
fdr_threshold <- 0.05
min_cells <- 15          # pseudobulk cell-count floor per (donor, cell type)
min_arm_exclude <- 4     # drop cell types with any drug arm below this many donors
drug_order <- c("Ale", "Bri", "Cer", "Cri", "Lor")

# The five 1v4 contrasts: each drug minus the mean of the other four. Written
# as strings so dreamlet::dreamlet parses them into a proper contrast matrix.
contrasts <- c(
  "Ale_vs_others = StimStatusAle - (StimStatusBri + StimStatusCer + StimStatusCri + StimStatusLor)/4",
  "Bri_vs_others = StimStatusBri - (StimStatusAle + StimStatusCer + StimStatusCri + StimStatusLor)/4",
  "Cer_vs_others = StimStatusCer - (StimStatusBri + StimStatusAle + StimStatusCri + StimStatusLor)/4",
  "Cri_vs_others = StimStatusCri - (StimStatusBri + StimStatusCer + StimStatusAle + StimStatusLor)/4",
  "Lor_vs_others = StimStatusLor - (StimStatusBri + StimStatusCer + StimStatusCri + StimStatusAle)/4"
)

# Paths
base_dir <- "/project/gilad/geena/publication_clean"
results_dir <- file.path(base_dir, "results/treatment/dreamlet_1v4_DE")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Load the pseudobulk object. No pre-filter on biotype or gene name is
# applied here; processAssays applies its per-cell-type expression filter
# and the surviving genes are what dreamlet tests.
message("Loading pb_treat...")
pb_treat <- readRDS(file.path(base_dir, "results/treatment/dreamlet_clean/pb_treat"))
message("Loaded pb_treat: ", nrow(pb_treat), " genes, ", ncol(pb_treat), " samples")

# Get per-(donor, cell type) cell counts. `arm` is the drug label extracted
# from the pseudobulk sample name (first 3 characters). Both cc and arm are
# used both here (arm-based cell-type filter) and later (per-arm QC table).
cc <- cellCounts(pb_treat)
arm <- substr(rownames(cc), 1, 3)

# Cell-type arm filter: for each cell type, count the donors in each drug arm
# that pass the min_cells threshold. Drop cell types where any arm has fewer
# than min_arm_exclude donors -- those contrasts can't be fit reliably.
arm_counts_prefilt <- do.call(rbind, lapply(colnames(cc), function(ct) {
  keep <- cc[, ct] >= min_cells
  tab <- table(factor(arm[keep], levels = drug_order))
  data.frame(cell_type = ct, t(as.integer(tab)))
}))
colnames(arm_counts_prefilt) <- c("cell_type", drug_order)

excluded_cts <- arm_counts_prefilt$cell_type[
  apply(arm_counts_prefilt[, drug_order], 1, min) < min_arm_exclude
]
cat("Excluding ", length(excluded_cts),
    " cell types with any arm < ", min_arm_exclude, " individuals:\n", sep = "")
cat(paste(" -", excluded_cts, collapse = "\n"), "\n")

keep_cts <- setdiff(assayNames(pb_treat), excluded_cts)

# Normalize + compute voom precision weights.
# processAssays() applies dreamlet's per-cell-type gene filter and produces
# the input to dreamlet(). Saved to disk because it's the expensive step.
message("Running processAssays...")
res.proc.indivRE_treat <- processAssays(
  pb_treat,
  ~ StimStatus + (1|individual),
  min.cells = min_cells,
  min.count = 5
)
saveRDS(res.proc.indivRE_treat, file.path(results_dir, "resprocindivRE_treat"))
message("Saved: resprocindivRE_treat")

# Write per-cell-type gene universes (plain text, one gene per line) for
# Cormotif. Cormotif reads these to restrict its input to the exact gene
# set dreamlet tested in that cell type, keeping the two analyses on the
# same universe. File names use the raw cell type with "/" replaced to "_"
# so there are no path-separator surprises.
universe_dir <- file.path(results_dir, "universes")
dir.create(universe_dir, recursive = TRUE, showWarnings = FALSE)
for (ct in assayNames(res.proc.indivRE_treat)) {
  genes <- rownames(assay(res.proc.indivRE_treat, ct))
  safe_ct <- gsub("/", "_", ct, fixed = TRUE)
  writeLines(genes, file.path(universe_dir, paste0(safe_ct, ".txt")))
}
message("Wrote per-cell-type gene universes to ", universe_dir,
        " (", length(assayNames(res.proc.indivRE_treat)), " files)")

# Cell types that survived both the arm filter above AND processAssays'
# gene filter.
keep_assays <- setdiff(intersect(keep_cts, assayNames(res.proc.indivRE_treat)),
                       "unassigned")

# Fit dreamlet with individual as random effect and the five 1v4 contrasts.
message("Running dreamlet (", length(keep_assays), " assays)...")
res.dl2_treat <- dreamlet(res.proc.indivRE_treat,
                          ~ 0 + StimStatus + (1|individual),
                          contrasts = contrasts,
                          assays = keep_assays)
saveRDS(res.dl2_treat, file.path(results_dir, "res.dl2_treat"))
message("Saved: res.dl2_treat")

# Extract per-drug top tables (one row per gene x cell type). flatten_topTable
# drops any list columns that dreamlet may attach so the result is writable
# as CSV.
flatten_topTable <- function(df, contrast_name) {
  df <- as.data.frame(df)
  df <- df[, sapply(df, function(x) !is.list(x))]
  df$contrast <- contrast_name
  df
}

message("Extracting topTable results...")
top_Ale <- topTable(res.dl2_treat, coef = contrasts[1], number = Inf)
top_Bri <- topTable(res.dl2_treat, coef = contrasts[2], number = Inf)
top_Cer <- topTable(res.dl2_treat, coef = contrasts[3], number = Inf)
top_Cri <- topTable(res.dl2_treat, coef = contrasts[4], number = Inf)
top_Lor <- topTable(res.dl2_treat, coef = contrasts[5], number = Inf)

write.csv(top_Ale, file.path(results_dir, "top_Ale_1v4.csv"))
write.csv(top_Bri, file.path(results_dir, "top_Bri_1v4.csv"))
write.csv(top_Cer, file.path(results_dir, "top_Cer_1v4.csv"))
write.csv(top_Cri, file.path(results_dir, "top_Cri_1v4.csv"))
write.csv(top_Lor, file.path(results_dir, "top_Lor_1v4.csv"))

# Bind the five per-drug tables into one long table (top_all_1v4.csv) that
# downstream scripts read.
top_all <- bind_rows(
  flatten_topTable(top_Ale, "Ale_vs_others"),
  flatten_topTable(top_Bri, "Bri_vs_others"),
  flatten_topTable(top_Cer, "Cer_vs_others"),
  flatten_topTable(top_Cri, "Cri_vs_others"),
  flatten_topTable(top_Lor, "Lor_vs_others")
) %>% mutate(drug = sub("_vs_others", "", contrast))

write.csv(top_all, file.path(results_dir, "top_all_1v4.csv"), row.names = FALSE)
message("Written: top_all_1v4.csv (", nrow(top_all), " rows)")

# Per-arm QC table. For every (cell type, drug) we record:
#   n_focal_arm: donors in the focal drug arm that pass min_cells
#   n_other_min: minimum donors across the other four drug arms
#   n_DEG, median_absLFC: DE burden from top_all above
# Rows are flagged "exclude" if any relevant arm is below min_arm_exclude.
# Wrapped in tryCatch so a tidyeval glitch here never kills the whole DE
# pipeline; top_all_1v4.csv + resprocindivRE_treat + res.dl2_treat are the
# real outputs and already written above.
message("Computing per-arm QC table...")
tryCatch({
  arm_counts <- do.call(rbind, lapply(colnames(cc), function(ct) {
    keep <- cc[, ct] >= min_cells
    tab <- table(factor(arm[keep], levels = drug_order))
    data.frame(cell_type = ct, t(as.integer(tab)))
  }))
  colnames(arm_counts) <- c("cell_type", drug_order)

  # Long-form arm count table built with base R (avoids tidyr edge cases).
  arm_long <- do.call(rbind, lapply(drug_order, function(d) {
    data.frame(cell_type = arm_counts$cell_type,
               drug = d,
               n_focal_arm = arm_counts[[d]],
               stringsAsFactors = FALSE)
  }))

  other_min_mat <- do.call(rbind, lapply(seq_len(nrow(arm_counts)), function(i) {
    sapply(drug_order, function(d) min(arm_counts[i, setdiff(drug_order, d)]))
  }))
  other_min_long <- do.call(rbind, lapply(drug_order, function(d) {
    data.frame(cell_type = arm_counts$cell_type,
               drug = d,
               n_other_min = other_min_mat[, d],
               stringsAsFactors = FALSE)
  }))
  arm_long <- merge(arm_long, other_min_long,
                    by = c("cell_type", "drug"), all.x = TRUE)

  deg_summary <- top_all %>%
    filter(adj.P.Val < fdr_threshold) %>%
    group_by(assay, drug) %>%
    summarise(n_DEG = n(), median_absLFC = median(abs(logFC)), .groups = "drop") %>%
    as.data.frame()

  qc <- merge(arm_long, deg_summary,
              by.x = c("cell_type", "drug"),
              by.y = c("assay", "drug"),
              all.x = TRUE)
  qc$n_DEG[is.na(qc$n_DEG)] <- 0L
  qc$min_arm <- pmin(qc$n_focal_arm, qc$n_other_min)
  qc$power_tier <- ifelse(
    qc$n_focal_arm < min_arm_exclude | qc$n_other_min < min_arm_exclude,
    "exclude", "pass"
  )
  names(qc)[names(qc) == "cell_type"] <- "assay"
  qc <- qc[order(qc$min_arm, -qc$median_absLFC, na.last = TRUE), ]

  write.csv(qc, file.path(results_dir, "DE_QC_powertier.csv"), row.names = FALSE)
  message("Written: DE_QC_powertier.csv")

  cat("\nPower tier distribution:\n")
  print(table(qc$power_tier))

  excluded <- qc[qc$power_tier == "exclude", ]
  excluded <- excluded[order(excluded$assay, excluded$drug), ]
  cat("\nExcluded conditions (n = ", nrow(excluded), "):\n", sep = "")
  print(excluded[, c("assay", "drug", "n_focal_arm", "n_other_min",
                     "n_DEG", "median_absLFC")])

  write.csv(excluded, file.path(results_dir, "DE_excluded_conditions.csv"),
            row.names = FALSE)
  message("Written: DE_excluded_conditions.csv")
}, error = function(e) {
  message("QC table generation failed (non-fatal): ", conditionMessage(e))
  message("Main DE outputs (top_all_1v4.csv, res.dl2_treat, universes/) are already written.")
})

message("01_dreamlet_1v4_DE.R complete.")

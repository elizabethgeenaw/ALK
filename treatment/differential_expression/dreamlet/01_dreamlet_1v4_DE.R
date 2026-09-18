#!/usr/bin/env Rscript

# 1v4 differential expression run for the ALK inhibitor paper.
# For each drug in turn, we contrast it against the mean of the other four
# ALK inhibitors ("one-vs-others"), fit a dreamlet model with individual
# as random effect, and write out topTable results per drug plus a combined
# top_all_1v4.csv used by nearly every downstream DE analysis. Genes are
# restricted to protein-coding; cell types with any drug arm below the
# min-individual threshold are dropped.
# Also writes a per-arm QC table (DE_QC_powertier.csv) and the subset of
# rows that fall below the exclusion threshold (DE_excluded_conditions.csv).

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
results_dir <- file.path(base_dir, "results/treatment/dreamlet_clean_proteincoding")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Load protein-coding gene list and restrict the pseudobulk object to those
# genes only. This filter happens once here; processAssays and dreamlet
# operate on the filtered object.
pc_genes <- readLines(file.path(base_dir, "reference_sets/HgProteinCodingGenes.tsv"))
message("Loaded ", length(pc_genes), " protein-coding genes")

message("Loading pb_treat...")
pb_treat_all <- readRDS(file.path(base_dir, "results/treatment/dreamlet_clean/pb_treat"))

keep_genes <- rownames(pb_treat_all) %in% pc_genes
pb_treat <- pb_treat_all[keep_genes, ]
message("Genes after protein-coding filter: ", nrow(pb_treat),
        " (from ", nrow(pb_treat_all), ")")
rm(pb_treat_all)

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
message("Computing per-arm QC table...")

arm_counts <- do.call(rbind, lapply(colnames(cc), function(ct) {
  keep <- cc[, ct] >= min_cells
  tab <- table(factor(arm[keep], levels = drug_order))
  data.frame(cell_type = ct, t(as.integer(tab)))
}))
colnames(arm_counts) <- c("cell_type", drug_order)

arm_long <- arm_counts %>%
  pivot_longer(all_of(drug_order), names_to = "drug", values_to = "n_focal_arm")

other_min_mat <- do.call(rbind, lapply(seq_len(nrow(arm_counts)), function(i) {
  sapply(drug_order, function(d) min(arm_counts[i, setdiff(drug_order, d)]))
}))
other_min_df <- data.frame(cell_type = arm_counts$cell_type, other_min_mat) %>%
  pivot_longer(cols = all_of(drug_order),
               names_to = "drug", values_to = "n_other_min")

arm_long <- left_join(arm_long, other_min_df, by = c("cell_type", "drug"))

deg_summary <- top_all %>%
  filter(adj.P.Val < fdr_threshold) %>%
  group_by(assay, drug) %>%
  summarise(n_DEG = n(), median_absLFC = median(abs(logFC)), .groups = "drop")

qc <- arm_long %>%
  left_join(deg_summary, by = c("cell_type" = "assay", "drug")) %>%
  replace_na(list(n_DEG = 0L, median_absLFC = NA_real_)) %>%
  mutate(
    min_arm = pmin(n_focal_arm, n_other_min),
    power_tier = ifelse(
      n_focal_arm < min_arm_exclude | n_other_min < min_arm_exclude,
      "exclude", "pass"
    )
  ) %>%
  rename(assay = cell_type) %>%
  arrange(min_arm, desc(median_absLFC))

write.csv(qc, file.path(results_dir, "DE_QC_powertier.csv"), row.names = FALSE)
message("Written: DE_QC_powertier.csv")

cat("\nPower tier distribution:\n")
print(table(qc$power_tier))

excluded <- qc %>% filter(power_tier == "exclude") %>% arrange(assay, drug)
cat("\nExcluded conditions (n = ", nrow(excluded), "):\n", sep = "")
print(excluded %>%
        select(assay, drug, n_focal_arm, n_other_min, n_DEG, median_absLFC) %>%
        as.data.frame())

write.csv(excluded, file.path(results_dir, "DE_excluded_conditions.csv"),
          row.names = FALSE)
message("Written: DE_excluded_conditions.csv")

message("01_dreamlet_1v4_DE.R complete.")

# Supp Table S6: S6_DE1v4_supp_table_de_gene_breadth_premash.csv
# Supp Table S7: S7_DE1v4_supp_table_de_gene_breadth_postmash_with_ntested.csv
# Supp Table S8: S8_DE1v4_mash_supp_table_ct_specific_with_ntested.csv

# Three supplementary tables summarizing per-(drug, gene) DE breadth from the
# 1v4 dreamlet run. All three read the same upstream data (dreamlet output +
# mash-shrunken PM/lfsr per drug), so they share loading and are built together.

# S6 (premash breadth): per (drug, gene) counts of cell types where the gene
# is DE in the raw dreamlet output (adj.P.Val < 0.05). One row per
# (drug, gene) that is DE in at least one cell type; columns include the
# breadth count, up/down splits, effect-size summaries, and a semicolon-
# separated "details" string listing each DE cell type with its adj.P and
# logFC.

# S7 (postmash breadth): same shape as S6 but computed from the mash-shrunken
# posterior (lfsr < 0.05, effect size = PM). Reflects mash's cross-cell-type
# borrowing of strength.

# S8 (cell-type-specific): every (drug, cell type, gene) where the gene is
# significant in that cell type AND in exactly one cell type total for the
# drug (breadth_K == 1). Includes n_ct_tested so a "K=1 out of 4 tested"
# call (narrow expression) can be distinguished from a "K=1 out of 30 tested"
# call (real cell-type-specific effect).

# All three tables are written as both .csv (machine-readable) and .xlsx
# (with the gene column forced to text so Excel doesn't auto-convert
# MARCH1/SEPT1/DEC1/OCT1 to dates).

# Inputs:
#   results/treatment/dreamlet_clean_proteincoding/top_all_1v4.csv
#   results/treatment/dreamlet_mash/by_drug_{Ale,Bri,Cer,Cri,Lor}/{mash_pm.csv, mash_lfsr.csv, logFC_original.csv}

# Outputs (written to results/treatment/dreamlet_mash/regen_v4/):
#   supp_table_de_gene_breadth_premash.{csv,xlsx}
#   supp_table_de_gene_breadth_postmash_with_ntested.{csv,xlsx}
#   supp_table_ct_specific_with_ntested.{csv,xlsx}

.libPaths(c("/project/gilad/geena/Rlibs/R-4.4.1", .libPaths()))
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(openxlsx)

# Paths
results_dir <- "/project/gilad/geena/publication_clean/results/treatment"
dreamlet_dir <- file.path(results_dir, "dreamlet_clean_proteincoding")
mash_dir <- file.path(results_dir, "dreamlet_mash")

# Write into a regen_v4 subdir so re-runs don't clobber the V4-referenced
# tables in the mash_dir top level.
out_dir <- file.path(mash_dir, "regen_v4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Constants
drugs <- c("Ale", "Bri", "Cer", "Cri", "Lor")

# adj.P.Val threshold for the premash (raw dreamlet) DE call
fdr_threshold <- 0.05
# lfsr threshold for the postmash (mash posterior) DE call
lfsr_threshold <- 0.05

ts <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")

# Small formatters for the de_ct_details string column.
# fmt_p: p-values in mixed sci/fixed with 3 sig figs (so 1e-40 stays sci and
# 0.03 stays fixed).
# fmt_fc: log2FC / posterior means at 4 decimal places (aesthetic; matches V4).
fmt_p <- function(x) format(x, digits = 3)
fmt_fc <- function(x) sprintf("%.4f", x)

# Helper: write a data frame as both a .csv (for programmatic consumers) and
# an .xlsx with the gene column force-formatted as text. The text format
# stops Excel from turning MARCH1 -> "March-1" etc. when the file is opened
# by double-click.
write_both <- function(df, csv_path, gene_col = "gene") {
  write_csv(df, csv_path)
  xlsx_path <- sub("\\.csv$", ".xlsx", csv_path)
  wb <- createWorkbook()
  sheet <- "table"
  addWorksheet(wb, sheet)
  writeData(wb, sheet, df, startRow = 1, startCol = 1)
  gene_col_idx <- which(names(df) == gene_col)
  if (length(gene_col_idx) == 1) {
    addStyle(wb, sheet,
             style = createStyle(numFmt = "TEXT"),
             rows = seq_len(nrow(df) + 1),
             cols = gene_col_idx,
             gridExpand = TRUE)
  }
  freezePane(wb, sheet, firstRow = TRUE)
  saveWorkbook(wb, xlsx_path, overwrite = TRUE)
  cat("  also wrote:", basename(xlsx_path), "\n")
}


# S6: premash DE breadth (from raw dreamlet 1v4 output)

# Load the dreamlet 1v4 DE table (one row per gene x cell type x drug).
# Rename to the S6/S7/S8 schema: dreamlet's "assay" column is cell type;
# dreamlet's "ID" column is the gene symbol.
cat(ts(), " Reading top_all_1v4.csv ...\n")
top <- read.csv(file.path(dreamlet_dir, "top_all_1v4.csv"),
                stringsAsFactors = FALSE)
cat(ts(), " top rows:", nrow(top),
    " drugs:", paste(unique(top$drug), collapse = ","), "\n")
top <- top %>% rename(gene = ID, cell_type = assay)

# For each (drug, gene), compute breadth = number of cell types where the
# gene is DE (adj.P.Val < 0.05), plus up/down splits, effect-size summaries
# across the DE subset, and a formatted "details" string listing each DE
# cell type ordered by adj.P.Val.
build_s6 <- function(df) {
  df %>%
    group_by(drug, gene) %>%
    summarize(
      n_ct_tested = n(),
      n_ct_de_premash = sum(adj.P.Val < fdr_threshold, na.rm = TRUE),
      n_ct_up = sum(adj.P.Val < fdr_threshold & logFC > 0, na.rm = TRUE),
      n_ct_down = sum(adj.P.Val < fdr_threshold & logFC < 0, na.rm = TRUE),
      # Effect-size summaries computed only over the DE cell types (not
      # all tested cell types), so mean_abs_logFC reflects "how strong the
      # DE calls are" rather than being diluted by non-significant cells.
      mean_abs_logFC = if (any(adj.P.Val < fdr_threshold))
                         mean(abs(logFC[adj.P.Val < fdr_threshold])) else NA_real_,
      max_abs_logFC = if (any(adj.P.Val < fdr_threshold))
                        max(abs(logFC[adj.P.Val < fdr_threshold])) else NA_real_,
      # top_ct = the cell type in which the gene has the largest |logFC| among
      # the DE cell types. Not always the most-significant cell type.
      top_ct = if (any(adj.P.Val < fdr_threshold)) {
                 sig <- adj.P.Val < fdr_threshold
                 cell_type[sig][which.max(abs(logFC[sig]))]
               } else NA_character_,
      top_logFC = if (any(adj.P.Val < fdr_threshold)) {
                    sig <- adj.P.Val < fdr_threshold
                    logFC[sig][which.max(abs(logFC[sig]))]
                  } else NA_real_,
      # Details string: for each DE cell type (ordered by adj.P.Val ascending),
      # print "cell type (adj.P=..., logFC=...)".
      de_ct_details = {
        sig <- adj.P.Val < fdr_threshold
        if (!any(sig)) "" else {
          ord <- order(adj.P.Val[sig])
          paste(sprintf("%s (adj.P=%s, logFC=%s)",
                        cell_type[sig][ord],
                        fmt_p(adj.P.Val[sig][ord]),
                        fmt_fc(logFC[sig][ord])),
                collapse = "; ")
        }
      },
      .groups = "drop"
    ) %>%
    # Drop (drug, gene) rows with no DE call anywhere; sort broadest first,
    # then largest effect.
    filter(n_ct_de_premash > 0) %>%
    arrange(drug, desc(n_ct_de_premash), desc(max_abs_logFC))
}

s6 <- build_s6(top)
s6_out <- file.path(out_dir, "supp_table_de_gene_breadth_premash.csv")
write_both(s6, s6_out)
cat(ts(), " Wrote S6:", s6_out, " (", nrow(s6), "rows)\n")


# S7 and S8: mash posterior tables

# Load the per-drug mash output. Each drug has three matrices in its
# by_drug_<Drug>/ folder: mash_pm (posterior means), mash_lfsr (local false
# sign rates), logFC_original (raw dreamlet log2FC pre-shrinkage). Rows are
# genes, columns are cell types; NA means the gene was not tested in that
# cell type.
cat(ts(), " Loading mash per-drug outputs ...\n")

load_mash_drug <- function(drug) {
  d <- file.path(mash_dir, paste0("by_drug_", drug))
  pm <- read.csv(file.path(d, "mash_pm.csv"), row.names = 1, check.names = FALSE)
  lfsr <- read.csv(file.path(d, "mash_lfsr.csv"), row.names = 1, check.names = FALSE)
  lfc0 <- read.csv(file.path(d, "logFC_original.csv"), row.names = 1, check.names = FALSE)
  # All three matrices should share the same gene rows and cell-type columns.
  stopifnot(identical(rownames(pm), rownames(lfsr)))
  stopifnot(identical(rownames(pm), rownames(lfc0)))
  stopifnot(identical(colnames(pm), colnames(lfsr)))
  stopifnot(identical(colnames(pm), colnames(lfc0)))

  # Pivot each matrix to long form (one row per gene x cell type) and inner-join
  # them on (gene, cell type). Drop rows with NA lfsr (untested).
  pm_l <- as.data.frame(as.table(as.matrix(pm)), stringsAsFactors = FALSE)
  lfsr_l <- as.data.frame(as.table(as.matrix(lfsr)), stringsAsFactors = FALSE)
  lfc_l <- as.data.frame(as.table(as.matrix(lfc0)), stringsAsFactors = FALSE)
  names(pm_l) <- c("gene", "cell_type", "PM")
  names(lfsr_l) <- c("gene", "cell_type", "lfsr")
  names(lfc_l) <- c("gene", "cell_type", "logFC_original")

  pm_l %>%
    inner_join(lfsr_l, by = c("gene", "cell_type")) %>%
    inner_join(lfc_l, by = c("gene", "cell_type")) %>%
    filter(!is.na(lfsr)) %>%
    mutate(drug = drug,
           gene = as.character(gene),
           cell_type = as.character(cell_type))
}

de_long <- map_dfr(drugs, load_mash_drug)
cat(ts(), " de_long rows (non-NA lfsr):", nrow(de_long), "\n")


# S7: postmash DE breadth per (drug, gene). Same shape as S6 but sig calls
# come from lfsr < 0.05 on the mash posterior, effect sizes are PM. Sorted
# alphabetically (drug, gene) to match V4.
build_s7 <- function(df) {
  df %>%
    group_by(drug, gene) %>%
    summarize(
      n_ct_tested = n(),
      n_ct_de_postmash = sum(lfsr < lfsr_threshold, na.rm = TRUE),
      n_ct_up = sum(lfsr < lfsr_threshold & PM > 0, na.rm = TRUE),
      n_ct_down = sum(lfsr < lfsr_threshold & PM < 0, na.rm = TRUE),
      mean_abs_PM = if (any(lfsr < lfsr_threshold))
                      mean(abs(PM[lfsr < lfsr_threshold])) else NA_real_,
      max_abs_PM = if (any(lfsr < lfsr_threshold))
                     max(abs(PM[lfsr < lfsr_threshold])) else NA_real_,
      top_ct = if (any(lfsr < lfsr_threshold)) {
                 sig <- lfsr < lfsr_threshold
                 cell_type[sig][which.max(abs(PM[sig]))]
               } else NA_character_,
      top_PM = if (any(lfsr < lfsr_threshold)) {
                 sig <- lfsr < lfsr_threshold
                 PM[sig][which.max(abs(PM[sig]))]
               } else NA_real_,
      de_ct_details = {
        sig <- lfsr < lfsr_threshold
        if (!any(sig)) "" else {
          ord <- order(lfsr[sig])
          paste(sprintf("%s (lfsr=%s, PM=%s)",
                        cell_type[sig][ord],
                        fmt_p(lfsr[sig][ord]),
                        fmt_fc(PM[sig][ord])),
                collapse = "; ")
        }
      },
      .groups = "drop"
    ) %>%
    filter(n_ct_de_postmash > 0) %>%
    arrange(drug, gene)
}

s7 <- build_s7(de_long)
s7_out <- file.path(out_dir, "supp_table_de_gene_breadth_postmash_with_ntested.csv")
write_both(s7, s7_out)
cat(ts(), " Wrote S7:", s7_out, " (", nrow(s7), "rows)\n")


# S8: cell-type-specific mash DE calls (breadth_K == 1). For each drug, keep
# rows where the gene is significant in exactly one cell type. Include
# n_ct_tested (from de_long) so readers can distinguish narrow expression
# (K=1 out of few tested) from real cell-type-specificity (K=1 out of many).
s8 <- de_long %>%
  filter(lfsr < lfsr_threshold) %>%
  group_by(drug, gene) %>%
  mutate(breadth_K = n()) %>%
  ungroup() %>%
  filter(breadth_K == 1) %>%
  select(drug, cell_type, gene, PM, logFC_original, lfsr, breadth_K)

ntested <- de_long %>%
  group_by(drug, gene) %>%
  summarize(n_ct_tested = n(), .groups = "drop")

s8 <- s8 %>%
  left_join(ntested, by = c("drug", "gene")) %>%
  arrange(drug, cell_type, desc(abs(PM)))

s8_out <- file.path(out_dir, "supp_table_ct_specific_with_ntested.csv")
write_both(s8, s8_out)
cat(ts(), " Wrote S8:", s8_out, " (", nrow(s8), "rows)\n")

cat(ts(), " 05_supp_S6_S7_S8_breadth_tables.R complete.\n")

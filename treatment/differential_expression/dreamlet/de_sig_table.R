# ---- V4 outputs ----
# Supp Table S5: /project/gilad/geena/ALK_draft/V4/Supplementary Tables/S5_DE1v4_top_all_1v4_sig.xlsx
# --------------------

#!/usr/bin/env Rscript
# 22_supp_S5_de_sig_table.R
# Pre-LOO significant 1v4 DE supplementary Excel.
# Sheet 1: flat significant rows (adj.P.Val < 0.05).
# Sheet 2: drug x cell-type count matrix.

library(dplyr)
library(tidyr)
library(openxlsx)

results_dir <- "/project/gilad/geena/publication_clean/results/treatment/dreamlet_1v4_DE"
out_path    <- file.path(results_dir, "supp_table_DE_sig.xlsx")

message("Reading top_all_1v4.csv...")
top <- read.csv(file.path(results_dir, "top_all_1v4.csv"))
message("  Rows: ", nrow(top))

message("Filtering to adj.P.Val < 0.05...")
sig <- top %>%
  filter(adj.P.Val < 0.05) %>%
  select(gene = ID, cell_type = assay, drug, contrast,
         logFC, AveExpr, t, P.Value, adj.P.Val, B) %>%
  arrange(cell_type, drug, adj.P.Val)
message("  Significant rows: ", nrow(sig))

message("Building drug x cell-type count matrix...")
counts <- sig %>%
  count(cell_type, drug) %>%
  pivot_wider(names_from = drug, values_from = n, values_fill = 0) %>%
  arrange(cell_type)

drug_cols <- setdiff(names(counts), "cell_type")
counts$Total <- rowSums(counts[, drug_cols, drop = FALSE])

total_row <- as.data.frame(
  c(list(cell_type = "Total"),
    as.list(colSums(counts[, c(drug_cols, "Total"), drop = FALSE]))),
  stringsAsFactors = FALSE
)
counts <- bind_rows(counts, total_row)

stopifnot(
  sum(counts[counts$cell_type != "Total", drug_cols]) == nrow(sig)
)

message("Cell types: ", sum(counts$cell_type != "Total"),
        "; drugs: ", length(drug_cols))

message("Writing workbook...")
wb <- createWorkbook()

addWorksheet(wb, "DE_sig_1v4")
writeData(wb, "DE_sig_1v4", sig)
addStyle(wb, "DE_sig_1v4",
         style      = createStyle(numFmt = "TEXT"),
         rows       = seq_len(nrow(sig) + 1),
         cols       = 1,
         gridExpand = TRUE)
freezePane(wb, "DE_sig_1v4", firstRow = TRUE)
setColWidths(wb, "DE_sig_1v4", cols = seq_along(sig), widths = "auto")

addWorksheet(wb, "counts_drug_by_celltype")
writeData(wb, "counts_drug_by_celltype", counts)
freezePane(wb, "counts_drug_by_celltype",
           firstActiveRow = 2, firstActiveCol = 2)
setColWidths(wb, "counts_drug_by_celltype",
             cols = seq_along(counts), widths = "auto")

saveWorkbook(wb, out_path, overwrite = TRUE)
message("Saved: ", out_path)

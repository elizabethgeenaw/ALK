# Utility script: re-derives the motif pattern table from the RDS files
# already written by cormotif, without refitting cormotif itself. Useful
# after a bug fix in the classification logic (e.g. drug_pattern) or when
# regenerating the xlsx with additional columns.
#
# For each cell type directory, reads the fitted motif RDS and any GO/gene
# outputs, then classifies each motif by which pairwise comparisons it is
# active in and whether that pattern matches "all comparisons involving drug X
# and no others" (drug-specific), all 10 comparisons (pan), none (null), or
# something else (mixed).

library(dplyr)
library(openxlsx)

# The 10 pairwise comparisons between the 5 drugs, in a stable order that
# the cormotif fit uses. comparison_drug1/2 give the two drugs on each side
# of each comparison so we can look up "which comparisons involve drug X".
comparison_labels <- c(
  "Ale vs Bri", "Ale vs Cer", "Ale vs Cri", "Ale vs Lor",
  "Bri vs Cer", "Bri vs Cri", "Bri vs Lor",
  "Cer vs Cri", "Cer vs Lor", "Cri vs Lor"
)

drug_colors <- c(Ale = "#0072B2", Bri = "#D55E00", Cer = "#009E73",
                 Cri = "#CC79A7", Lor = "#E69F00")
comparison_drug1 <- c("Ale","Ale","Ale","Ale","Bri","Bri","Bri","Cer","Cer","Cri")
comparison_drug2 <- c("Bri","Cer","Cri","Lor","Cer","Cri","Lor","Cri","Lor","Lor")

# Walk every cell-type subdirectory under output_base_dir, read the fitted
# motif RDS, and emit one row per (cell_type, motif). motif_q_threshold is
# the posterior probability cutoff above which a comparison is considered
# "active" for a given motif.
generate_pattern_table <- function(output_base_dir, motif_q_threshold = 0.5) {
  drugs <- names(drug_colors)
  # For each drug, which comparison indices involve that drug on either side.
  drug_comp_idx <- setNames(lapply(drugs, function(d) {
    which(comparison_drug1 == d | comparison_drug2 == d)
  }), drugs)

  cell_type_dirs <- list.dirs(output_base_dir, recursive = FALSE)
  rows <- list()

  for (ct_dir in cell_type_dirs) {
    cell_type <- basename(ct_dir)
    rds_file  <- file.path(ct_dir, paste0(cell_type, "_motif_fittedK2plus.rds"))
    if (!file.exists(rds_file)) next

    motif.fitted    <- readRDS(rds_file)
    motif_q         <- motif.fitted$bestmotif$motif.q       # motifs x comparisons
    posterior_probs <- motif.fitted$bestmotif$clustlike     # genes x motifs
    K           <- nrow(motif_q)
    total_genes <- nrow(posterior_probs)

    for (m in seq_len(K)) {
      q_row   <- motif_q[m, ]
      active  <- q_row > motif_q_threshold
      # Genes assigned to this motif under a soft threshold on cluster posterior.
      n_genes <- sum(posterior_probs[, m] > 0.5)
      is_null <- !any(active)

      active_comparisons <- paste(comparison_labels[active], collapse = "; ")

      # Strict drug-specific classification: active in ALL comparisons
      # involving drug X AND inactive in ALL comparisons that exclude drug X.
      # If more than one drug passes the test we call it a combined pattern.
      drug_flags <- sapply(drugs, function(d) {
        involved     <- drug_comp_idx[[d]]
        not_involved <- setdiff(seq_along(comparison_labels), involved)
        all(active[involved]) && !any(active[not_involved])
      })
      active_drugs <- drugs[drug_flags]
      drug_label <- if (length(active_drugs) == 0) {
        if (is_null) "null" else "mixed"
      } else if (length(active_drugs) == 1) {
        paste0(active_drugs, "-specific")
      } else {
        paste(active_drugs, collapse = "+")
      }

      # If a GO:BP enrichment CSV exists for this motif, take the top 3 terms
      # by elimFDR (topGO's elim-corrected Fisher p) and paste them into a
      # single cell for readability in the xlsx.
      go_file   <- file.path(ct_dir, paste0(cell_type, "_motif", m, "_GO_BP_K2plus.csv"))
      top_terms <- ""
      top_fdrs  <- ""
      if (file.exists(go_file)) {
        go_df  <- read.csv(go_file)
        go_sig <- go_df[!is.na(go_df$elimFDR) & go_df$elimFDR < 0.05, ]
        if (nrow(go_sig) > 0) {
          go_top    <- head(go_sig[order(go_sig$elimFDR), ], 3)
          top_terms <- paste(go_top$Term, collapse = " | ")
          top_fdrs  <- paste(signif(go_top$elimFDR, 2), collapse = " | ")
        }
      }

      # Genes list is stored as one gene per line; concat for the xlsx column.
      gene_file <- file.path(ct_dir, paste0(cell_type, "_motif", m, "_genes_K2plus.txt"))
      gene_list_str <- ""
      if (file.exists(gene_file)) {
        gene_list_str <- paste(readLines(gene_file), collapse = "; ")
      }

      rows[[paste(cell_type, m)]] <- data.frame(
        cell_type          = cell_type,
        motif              = m,
        n_genes            = n_genes,
        pct_genes          = round(100 * n_genes / total_genes, 1),
        is_null            = is_null,
        drug_pattern       = drug_label,
        active_comparisons = active_comparisons,
        top_BP_terms       = top_terms,
        top_BP_FDR         = top_fdrs,
        gene_list          = gene_list_str,
        stringsAsFactors   = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    cat("No results found in:", output_base_dir, "\n")
    return(invisible(NULL))
  }

  out      <- bind_rows(rows)
  out_file <- file.path(output_base_dir, "motif_pattern_table.xlsx")
  write.xlsx(out, out_file, rowNames = FALSE)
  cat("Saved motif_pattern_table.xlsx:", nrow(out), "rows\n")
  invisible(out)
}

# Rerun for the 15-cell / protein-coding-only cormotif output that feeds
# the manuscript figures.
output_base_dir_pc <- "/project/gilad/geena/publication_clean/results/treatment/cormotif/protein_coding_only_15cells/"
generate_pattern_table(output_base_dir_pc)

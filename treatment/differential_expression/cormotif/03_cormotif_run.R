# Cormotif refit + downstream outputs.
# Produces Supp Fig S7 (GO_BP_elim multi-page combined plots) and
# Supp Tables S10-S14 (motif pattern table + gene lists + significant GO
# terms for BP/MF/CC).
#
# Pipeline:
#   1. filter_genes_cormotif(): protein-coding restriction, drop MT-/RPS-/RPL-
#      genes, keep only genes with avg log2CPM > 0. Returns RAW counts.
#      voom + TMM happen downstream inside limmafit.
#   2. run_pairwise_limma(): fits limma for all 10 pairwise drug comparisons
#      with global filter only (no per-comparison filterByExpr, which would
#      drop drug-specific genes from the comparisons that would detect them).
#   3. Cormotif runtype = "limmafits" (bypasses a bug in the installed
#      CormotifCounts::limmafit.counts).
#   4. extract_genes_from_motif() uses rownames(p.post) as the gene universe
#      so it stays consistent with what Cormotif actually saw.
#   5. Post-hoc: topGO enrichment, motif pattern table, mixed-motif
#      asymmetry, per-cell-type combined GO plots, and multi-cell-type Supp
#      Fig S7 assembly.
#
# Cell-type inclusion: all 5 per-arm counts >= 4 individuals with >= 15 cells

library(dplyr)
library(ggplot2)
library(Cormotif)
library(edgeR)
library(CormotifCounts)
library(topGO)
library(org.Hs.eg.db)
library(openxlsx)

# Constants

comparison_labels <- c(
  "Ale vs Bri", "Ale vs Cer", "Ale vs Cri", "Ale vs Lor",
  "Bri vs Cer", "Bri vs Cri", "Bri vs Lor",
  "Cer vs Cri", "Cer vs Lor", "Cri vs Lor"
)

drug_colors <- c(Ale = "#0072B2",
                 Bri = "#D55E00",
                 Cer = "#009E73",
                 Cri = "#CC79A7",
                 Lor = "#E69F00")
comparison_drug1 <- c("Ale","Ale","Ale","Ale","Bri","Bri","Bri","Cer","Cer","Cri")
comparison_drug2 <- c("Bri","Cer","Cri","Lor","Cer","Cri","Lor","Cri","Lor","Lor")
drugs            <- names(drug_colors)

# Germ layer colors — exact values from 5_sketch_umap_visualization.R (germ_colors)
germ_layer_colors <- c(
  Epithelium = "#4292C6",
  Ectoderm   = "#F16913",
  Mesoderm   = "#44AA99",
  Endoderm   = "#9970AB",
  iPSC       = "#CC6677"
)

# Germ layer assignment for the 30 DE cell types.
# Derived from the ectoderm/mesoderm/endoderm/epithelium vectors in
# 5_sketch_umap_visualization.R (lines 18–21).
# Order matches the UMAP legend: Epithelium → Ectoderm → Mesoderm → Endoderm → iPSCs.
cell_type_germ_layer <- c(
  # Epithelium
  "Goblet cells"                            = "Epithelium",
  "Ciliated epithelial cells"               = "Epithelium",
  "Intestinal epithelial cells"             = "Epithelium",
  "Squamous epithelial cells"               = "Epithelium",
  "Bronchiolar and alveolar epithelial cells" = "Epithelium",
  # Ectoderm
  "ENS glia"                                = "Ectoderm",
  "Glioblast"                               = "Ectoderm",
  "Neuroblast"                              = "Ectoderm",
  "Neuron"                                  = "Ectoderm",
  "Neuronal IPC"                            = "Ectoderm",
  "Radial glia"                             = "Ectoderm",
  "Retinal neurons"                         = "Ectoderm",
  "Retinal pigment cells"                   = "Ectoderm",
  "Retinal progenitors and Muller glia"     = "Ectoderm",
  # Mesoderm
  "Early mesoderm cells"                    = "Mesoderm",
  "Cardiac mesoderm cells"                  = "Mesoderm",
  "Cardiomyocytes"                          = "Mesoderm",
  "Epicardial fat cells"                    = "Mesoderm",
  "Smooth muscle cells"                     = "Mesoderm",
  "Vascular endothelial cells"              = "Mesoderm",
  "Stromal cells"                           = "Mesoderm",
  "Metanephric cells"                       = "Mesoderm",
  "Skeletal muscle cells"                   = "Mesoderm",
  "Lymphoid cells"                          = "Mesoderm",
  "Erythroblasts"                           = "Mesoderm",
  "Megakaryocytes"                          = "Mesoderm",
  "Myeloid cells"                           = "Mesoderm",
  # Endoderm
  "Hepatoblasts"                            = "Endoderm",
  "Stellate cells"                          = "Endoderm",
  # iPSCs last
  "iPSCs"                                   = "iPSC"
)

# Motif specificity helpers

# Returns the drug name if the motif is exactly "active in all comparisons
# involving that drug, inactive in all others"; otherwise NA.
classify_motif <- function(motif_row_q, threshold = 0.5) {
  active <- motif_row_q >= threshold
  for (drug in drugs) {
    x_idx     <- which(comparison_drug1 == drug | comparison_drug2 == drug)
    non_x_idx <- setdiff(seq_along(motif_row_q), x_idx)
    if (all(active[x_idx]) && !any(active[non_x_idx])) return(drug)
  }
  NA_character_
}

motif_specificity <- function(motif.fitted, threshold = 0.5) {
  apply(motif.fitted$bestmotif$motif.q, 1, function(q) classify_motif(q, threshold = threshold))
}


# Computes motif classification metadata for a single motif's q-vector.
# Returns a list with four fields used by generate_pattern_table(),
# generate_supplementary_tables(), and compute_mixed_motif_asymmetry().
#
# Uses >= threshold (matching classify_motif()) — this fixes the > vs >=
# inconsistency that existed in generate_pattern_table().
#
# motif_type distinguishes pan-comparison (all 10 active) from generic "mixed":
#   "null"          — no comparisons active
#   "pan-comparison"— all 10 comparisons active
#   "{drug}-specific"— exactly the 4 comparisons involving that drug are active
#   "mixed"         — active, but not all-10 and not drug-specific
#
# drug_pattern is the raw label (same as drug_label in old generate_pattern_table):
#   "null" / "mixed" / "{drug}-specific" / "{drug1}+{drug2}" (unreachable; see below)

compute_motif_meta <- function(q_row, threshold = 0.5) {
  active   <- q_row >= threshold
  n_active <- sum(active)
  is_null  <- n_active == 0

  drug_comp_idx <- setNames(lapply(drugs, function(d)
    which(comparison_drug1 == d | comparison_drug2 == d)), drugs)

  drug_flags <- sapply(drugs, function(d) {
    inv <- drug_comp_idx[[d]]
    non <- setdiff(seq_along(q_row), inv)
    all(active[inv]) && !any(active[non])
  })
  active_drugs <- drugs[drug_flags]

  # drug_pattern: raw pattern label
  drug_pattern <- if (length(active_drugs) == 0) {
    if (is_null) "null" else "mixed"
  } else if (length(active_drugs) == 1) {
    paste0(active_drugs, "-specific")
  } else {
    # Two drugs simultaneously satisfying drug_flags is mathematically impossible:
    # each drug's "not_involved" set overlaps the other's "involved" set, so
    # the required active/inactive constraints contradict each other. This branch
    # is unreachable in practice.
    paste(active_drugs, collapse = "+")
  }

  # motif_type: reader-facing category
  motif_type <- if (is_null) {
    "null"
  } else if (n_active == length(q_row)) {
    "pan-comparison"
  } else if (length(active_drugs) == 1) {
    paste0(active_drugs, "-specific")
  } else {
    "mixed"
  }

  list(
    is_null              = is_null,
    n_active_comparisons = n_active,
    drug_pattern         = drug_pattern,
    motif_type           = motif_type,
    active_comparisons   = paste(comparison_labels[active], collapse = "; ")
  )
}


# Excel-safe writer: forces character columns to text format so Excel cannot
# auto-convert gene names (SEPT2, MARCH1, etc.) to dates.
write_xlsx_safe <- function(df, path) {
  # Excel hard limit: 32,767 characters per cell. Long gene-list strings exceed
  # this; truncate before writeData so openxlsx does not warn or corrupt the
  # sharedStrings table. Full strings are preserved in the companion CSV.
  excel_cell_max <- 32767L
  df <- as.data.frame(lapply(df, function(col) {
    if (!is.character(col)) return(col)
    too_long <- !is.na(col) & nchar(col) > excel_cell_max
    col[too_long] <- paste0(substr(col[too_long], 1L, excel_cell_max - 14L),
                            "...[truncated]")
    col
  }), stringsAsFactors = FALSE)

  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "Sheet1")
  openxlsx::writeData(wb, "Sheet1", df)
  text_style <- openxlsx::createStyle(numFmt = "TEXT")
  char_cols  <- which(sapply(df, is.character))
  for (col in char_cols) {
    openxlsx::addStyle(wb, "Sheet1", style = text_style,
                       rows = 2:(nrow(df) + 1), cols = col,
                       gridExpand = TRUE)
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  cat("Saved:", path, "\n")
}


# Gene filter (replaces normalize_matrix).
# Returns filtered RAW COUNTS; voom + TMM happen inside limmafit.counts.
# Side effect: appends a row to gene_filter_log in the calling environment
# via <<- when cell_type is non-NULL.
filter_genes_cormotif <- function(celltype_counts, protein_coding_genes = NULL,
                                   cell_type = NULL) {
  m <- celltype_counts
  n_input <- nrow(m)

  if (!is.null(protein_coding_genes))
    m <- m[rownames(m) %in% protein_coding_genes, , drop = FALSE]
  n_pc <- nrow(m)

  mt_rp <- grepl("^MT-", rownames(m)) | grepl("^RP[SL]", rownames(m))
  m <- m[!mt_rp, , drop = FALSE]
  n_mtrp <- nrow(m)

  dge <- edgeR::DGEList(counts = m)
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  avglogcpm <- rowMeans(edgeR::cpm(dge, log = TRUE, prior.count = 1))
  m <- m[avglogcpm > 0, , drop = FALSE]
  n_cpm <- nrow(m)

  ct_label <- if (!is.null(cell_type)) cell_type else "?"
  cat("  ", ct_label, ": ", n_input,
      " -> ", n_pc,    " (protein-coding)",
      " -> ", n_mtrp,  " (MT/RP)",
      " -> ", n_cpm,   " (avgCPM>0)\n",
      sep = "")

  if (!is.null(cell_type) && exists("gene_filter_log", envir = .GlobalEnv)) {
    gene_filter_log[[cell_type]] <<- data.frame(
      cell_type              = cell_type,
      n_input                = n_input,
      n_after_protein_coding = n_pc,
      n_after_mt_rp          = n_mtrp,
      n_after_avglogcpm      = n_cpm,
      stringsAsFactors       = FALSE
    )
  }

  m  # raw counts
}


# Plotting helpers

plot_motif_ggplot <- function(motif.fitted, celltype_name, output_dir) {
  u  <- 1:dim(motif.fitted$bestmotif$motif.q)[2]
  v  <- 1:dim(motif.fitted$bestmotif$motif.q)[1]
  Ng <- nrow(motif.fitted$bestmotif$p.post)

  png(file.path(output_dir, paste0(celltype_name, "_motif_plot_K2plus.png")),
      width = 10, height = 9, units = "in", res = 300)
  old_par <- par(mar = c(9, 4, 4, 2))
  on.exit({ par(old_par); dev.off() })

  layout(matrix(1:2, ncol = 2), widths = c(2, 1))

  image(u, v, t(motif.fitted$bestmotif$motif.q),
        col  = gray(seq(from = 1, to = 0, by = -0.1)),
        xaxt = "n", yaxt = "n",
        xlab = "", ylab = "Corr. Motifs",
        main = paste(celltype_name, "pattern"))
  axis(2, at = seq_along(v))
  axis(1, at = u, labels = FALSE)

  usr    <- par("usr")
  yrange <- diff(usr[3:4])
  sh     <- yrange * 0.07

  text(x = u[1] - 0.6, y = usr[3] - sh / 2,
       labels = "Drug 1", adj = c(1, 0.5), xpd = TRUE, cex = 0.7)
  text(x = u[1] - 0.6, y = usr[3] - 1.5 * sh,
       labels = "Drug 2", adj = c(1, 0.5), xpd = TRUE, cex = 0.7)

  rect(u - 0.5, usr[3] - sh,      u + 0.5, usr[3],
       col = drug_colors[comparison_drug1], border = NA, xpd = TRUE)
  text(x = u, y = usr[3] - sh / 2,
       labels = comparison_drug1, col = "white",
       srt = 0, adj = c(0.5, 0.5), xpd = TRUE, cex = 0.65, font = 2)

  rect(u - 0.5, usr[3] - 2 * sh,  u + 0.5, usr[3] - sh,
       col = drug_colors[comparison_drug2], border = NA, xpd = TRUE)
  text(x = u, y = usr[3] - 1.5 * sh,
       labels = comparison_drug2, col = "white",
       srt = 0, adj = c(0.5, 0.5), xpd = TRUE, cex = 0.65, font = 2)

  text(x = u, y = usr[3] - 2 * sh - yrange * 0.02,
       labels = comparison_labels,
       srt = 45, adj = 1, xpd = TRUE, cex = 0.9)

  for (i in 1:(length(u) + 1)) abline(v = (i - 0.5))
  for (i in 1:(length(v) + 1)) abline(h = (i - 0.5))

  par(mar = c(9, 1, 4, 2))
  NK          <- nrow(motif.fitted$bestmotif$motif.q)
  genecount   <- colSums(motif.fitted$bestmotif$clustlike > 0.5)
  motif_p_post <- genecount / Ng
  plot(0, 0.7, pch = ".", xlim = c(0, 1.2), ylim = c(0.75, NK + 0.25),
       frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "",
       main = "Number of genes in motif")
  segments(0, 0.7, motif_p_post[1], 0.7)
  rect(0, 1:NK - 0.3, motif_p_post, 1:NK + 0.3, col = "dark grey")
  mtext(1:NK, at = 1:NK, side = 2, cex = 0.8)
  text(motif_p_post + 0.15, 1:NK, labels = genecount)
}


# Internal: heatmap panel with drug-specific color overlay and Roman numeral labels.
# Ported from replot_motifs_combined_specificity.R.
plot_motif_heatmap_panel <- function(motif.fitted, celltype_name, specificity) {
  u       <- 1:dim(motif.fitted$bestmotif$motif.q)[2]
  v       <- 1:dim(motif.fitted$bestmotif$motif.q)[1]
  motif_q <- motif.fitted$bestmotif$motif.q

  image(u, v, t((motif_q >= 0.5) * 1L),
        col    = c("white", "black"),
        breaks = c(-0.5, 0.5, 1.5),
        xaxt = "n", yaxt = "n",
        xlab = "", ylab = "Corr. Motifs",
        main = celltype_name,
        cex.lab = 1.8, cex.main = 2.0)

  # Recolor active cells for drug-specific motifs
  for (i in seq_along(specificity)) {
    if (!is.na(specificity[i])) {
      active_cols <- which(motif_q[i, ] >= 0.5)
      for (j in active_cols)
        rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5,
             col = drug_colors[specificity[i]], border = NA)
    }
  }

  axis(2, at = seq_along(v), labels = as.roman(seq_along(v)),
       las = 1, cex.axis = 1.8)

  usr    <- par("usr")
  yrange <- diff(usr[3:4])
  sh     <- yrange * 0.07

  text(x = u[1] - 0.6, y = usr[3] - sh / 2,
       labels = "Drug 1", adj = c(1, 0.5), xpd = TRUE, cex = 1.4)
  text(x = u[1] - 0.6, y = usr[3] - 1.5 * sh,
       labels = "Drug 2", adj = c(1, 0.5), xpd = TRUE, cex = 1.4)

  rect(u - 0.5, usr[3] - sh, u + 0.5, usr[3],
       col = drug_colors[comparison_drug1], border = NA, xpd = TRUE)
  text(x = u, y = usr[3] - sh / 2, labels = comparison_drug1, col = "white",
       adj = c(0.5, 0.5), xpd = TRUE, cex = 1.3, font = 2)

  rect(u - 0.5, usr[3] - 2 * sh, u + 0.5, usr[3] - sh,
       col = drug_colors[comparison_drug2], border = NA, xpd = TRUE)
  text(x = u, y = usr[3] - 1.5 * sh, labels = comparison_drug2, col = "white",
       adj = c(0.5, 0.5), xpd = TRUE, cex = 1.3, font = 2)

  for (i in 1:(length(u) + 1)) abline(v = (i - 0.5))
  for (i in 1:(length(v) + 1)) abline(h = (i - 0.5))
}

# Internal: gene-count bar panel with drug-specific annotation.
# Ported from replot_motifs_combined_specificity.R.
plot_barcount_panel <- function(motif.fitted, specificity) {
  NK        <- nrow(motif.fitted$bestmotif$motif.q)
  genecount <- colSums(motif.fitted$bestmotif$clustlike > 0.5)

  # Motif I is the non-DE background pool (active in no pairwise comparison).
  # On a linear scale it compresses the drug-specific motifs to invisibility
  # (e.g. 16,048 vs 342 vs 671 vs 2,225 in Hepatoblasts/Cardiac mesoderm).
  # Render it as an inline text annotation and rescale the DE motif bars so
  # they are comparable to each other.
  bar_max    <- 0.9
  bar_widths <- rep(0, NK)
  if (NK >= 2) {
    de_counts <- genecount[2:NK]
    denom     <- max(de_counts, na.rm = TRUE)
    if (denom > 0) bar_widths[2:NK] <- de_counts / denom * bar_max
  }

  plot(0, 0.7, pch = ".", xlim = c(0, 2.3), ylim = c(0.75, NK + 0.25),
       frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "",
       main = "Number of genes in motif", cex.main = 1.8)

  if (NK >= 2)
    rect(0, (2:NK) - 0.3, bar_widths[2:NK], (2:NK) + 0.3, col = "dark grey")

  mtext(as.character(as.roman(1:NK)), at = 1:NK, side = 2, las = 1, cex = 1.7)

  text(0.02, 1,
       labels = paste0(format(genecount[1], big.mark = ","), " (non-DE)"),
       cex = 1.3, adj = c(0, 0.5), font = 3, col = "grey30")

  if (NK >= 2)
    text(bar_widths[2:NK] + 0.03, 2:NK,
         labels = format(genecount[2:NK], big.mark = ","),
         cex = 1.5, adj = c(0, 0.5))

  for (i in seq_len(NK)) {
    if (!is.na(specificity[i]))
      text(1.2, i, labels = paste0(specificity[i], "-specific"),
           col = drug_colors[specificity[i]], cex = 1.6, font = 2,
           adj = c(0, 0.5), xpd = TRUE)
  }
}

plot_motif_publication <- function(motif.fitted, celltype_name, output_dir) {
  spec    <- motif_specificity(motif.fitted)
  out_pdf <- file.path(output_dir, paste0(celltype_name, "_motif_plot_publication.pdf"))

  cairo_pdf(out_pdf, width = 11, height = 7, family = "Helvetica")
  old_par <- par(mar = c(7, 4, 4, 1))
  on.exit({ par(old_par); dev.off() })

  layout(matrix(1:2, ncol = 2), widths = c(2, 1.5))
  plot_motif_heatmap_panel(motif.fitted, celltype_name, spec)
  par(mar = c(7, 1, 4, 3))
  plot_barcount_panel(motif.fitted, spec)

  cat("Saved publication motif plot:", out_pdf, "\n")
}


plot_go_panel <- function(cell_type, motif_number, ontology = "BP",
                          algorithm = "elim",
                          top_n = 8, output_dir = NULL,
                          go_df = NULL, output_base_dir = NULL,
                          save_plot = FALSE,
                          base_size  = 12,   # ggplot2 base font size in pt
                          wrap_width = 60) { # strwrap width for GO term labels
  # Filter and rank on raw p-value (elim or classic Fisher).
  # topGO authors advise against BH correction on GO p-values because tests are
  # positively correlated via the GO hierarchy; elim additionally conditions on
  # children, making BH over all terms extremely conservative.
  p_col <- if (algorithm == "elim") "elimFisher" else "classicFisher"

  if (is.null(go_df)) {
    go_file <- file.path(output_base_dir, cell_type,
                         paste0(cell_type, "_motif", motif_number,
                                "_GO_", ontology, "_K2plus.csv"))
    if (!file.exists(go_file)) {
      cat("GO file not found:", go_file, "\n")
      return(invisible(NULL))
    }
    go_df <- read.csv(go_file)
  }

  go_sig <- go_df[!is.na(go_df[[p_col]]) & go_df[[p_col]] < 0.05, ]
  if (nrow(go_sig) == 0) {
    cat("No significant GO terms for", cell_type, "motif", motif_number, ontology,
        "(", algorithm, ")\n")
    return(invisible(NULL))
  }

  go_top <- head(go_sig[order(go_sig[[p_col]]), ], top_n)
  go_top$neg_log10_p <- -log10(go_top[[p_col]])
  go_top$Term_wrapped <- sapply(go_top$Term, function(x)
    paste(strwrap(x, width = wrap_width), collapse = "\n"))
  go_top$Term_wrapped <- factor(go_top$Term_wrapped,
                                levels = rev(go_top$Term_wrapped))

  p <- ggplot(go_top, aes(x = neg_log10_p, y = Term_wrapped)) +
    geom_col(fill = "#0072B2", width = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "gray50", linewidth = 0.5) +
    labs(x    = expression(-log[10](italic(p))),
         y    = NULL,
         title = paste0(cell_type, "  |  Motif ", motif_number,
                        " | ", ontology, " (", algorithm, ")")) +
    theme_classic(base_size = base_size) +
    theme(plot.title  = element_text(hjust = 0.5, face = "bold", size = base_size),
          axis.text.y = element_text(size = base_size * 0.83),
          axis.text.x = element_text(size = base_size * 0.83))

  if (save_plot) {
    if (is.null(output_dir)) output_dir <- file.path(output_base_dir, cell_type)
    out_file <- file.path(output_dir,
                          paste0(cell_type, "_motif", motif_number,
                                 "_GO_", ontology, "_", algorithm, "_panel.png"))
    ggsave(out_file, p,
           width = 8, height = 2 + top_n * 0.4, dpi = 300, limitsize = FALSE)
    cat("Saved GO panel:", out_file, "\n")
  }
  invisible(p)
}


# Pairwise limma fits (TMM + voom, no per-comparison filterByExpr).
# Per-comparison filterByExpr would exclude genes expressed specifically in one
# drug arm; we use the global avglogcpm > 0 filter in filter_genes_cormotif
# instead. voom's mean-variance weighting handles low-expression genes within
# specific comparisons by assigning them near-zero t-statistics (effectively
# neutral in the EM). This matches the lab's CormotifCounts_v2 approach.
#
# The installed CormotifCounts::limmafit.counts has a separate coding bug
# (not.nas %in% keep with a logical vector collapses genes to index 1). We
# bypass it by running the 10 comparisons here and using runtype = "limmafits".

run_pairwise_limma <- function(count_mat, group_labels_numeric, compid,
                               norm.method = "TMM") {
  compnum <- nrow(compid)
  comp_names <- c("Ale_Bri","Ale_Cer","Ale_Cri","Ale_Lor",
                  "Bri_Cer","Bri_Cri","Bri_Lor","Cer_Cri","Cer_Lor","Cri_Lor")
  fitlist <- setNames(vector("list", compnum), comp_names)

  for (i in seq_len(compnum)) {
    message("Running limma for comparison ", i, " / ", compnum)
    selid1 <- which(group_labels_numeric == compid[i, 1])
    selid2 <- which(group_labels_numeric == compid[i, 2])
    counts_sub <- cbind(count_mat[, selid1, drop = FALSE],
                        count_mat[, selid2, drop = FALSE])
    g1num <- length(selid1)
    g2num <- length(selid2)

    designmat <- cbind(base  = rep(1L, g1num + g2num),
                       delta = c(rep(0L, g1num), rep(1L, g2num)))

    # No per-comparison filterByExpr: genes expressed specifically in one drug arm
    # (near-zero counts in the comparison arm) would fail filterByExpr for the
    # comparisons against that arm, making them invisible to the EM. The global
    # avglogcpm > 0 filter in filter_genes_cormotif is the correct inclusion
    # boundary; voom assigns near-zero t-statistics (effectively neutral) to
    # lowly expressed genes within specific comparisons without discarding them.
    selected <- which(!apply(counts_sub, 1, anyNA))

    if (length(selected) == 0) {
      cat("  WARNING: 0 genes pass filters for comparison", i, "— skipping\n")
      next
    }

    d   <- edgeR::DGEList(counts_sub[selected, , drop = FALSE])
    d   <- edgeR::calcNormFactors(d, method = norm.method)
    y   <- limma::voom(d, designmat)
    fit <- limma::lmFit(y, designmat)
    fit <- limma::eBayes(fit)

    fitlist[[i]] <- fit
  }
  fitlist
}


# Cormotif fit

run_cormotif <- function(count_data, celltype_name, output_dir) {
  print("Running cormotif...")
  sample_names <- colnames(count_data)

  # Column names are {drug}NA{individual} (e.g. AleNA18486)
  group_labels <- sapply(sample_names, function(x) strsplit(x, "NA")[[1]][1])
  group_labels_numeric <- as.numeric(
    factor(group_labels, levels = c("Ale", "Bri", "Cer", "Cri", "Lor"))
  )

  compid <- matrix(c(
    1, 2,  # Ale vs Bri
    1, 3,  # Ale vs Cer
    1, 4,  # Ale vs Cri
    1, 5,  # Ale vs Lor
    2, 3,  # Bri vs Cer
    2, 4,  # Bri vs Cri
    2, 5,  # Bri vs Lor
    3, 4,  # Cer vs Cri
    3, 5,  # Cer vs Lor
    4, 5   # Cri vs Lor
  ), ncol = 2, byrow = TRUE)

  fitlist <- run_pairwise_limma(as.matrix(count_data), group_labels_numeric, compid)
  saveRDS(fitlist, file.path(output_dir, paste0(celltype_name, "_fitlist.rds")))

  set.seed(1)
  motif.counts <- CormotifCounts::cormotiffit(
    fitlist,
    runtype  = "limmafits",   # bypasses buggy limmafit.counts; uses limmafit.list
    K        = 2:8,
    max.iter = 1000,
    BIC      = TRUE
  )

  # CormotifCounts with runtype="limmafits" does not propagate gene names to
  # p.post/clustlike. Recover them from fitlist: limmafit.list builds the
  # t-stat matrix by appending genes in first-appearance order across comparisons.
  if (is.null(rownames(motif.counts$bestmotif$p.post))) {
    gene_order <- character(0)
    for (fit in fitlist) {
      if (!is.null(fit)) {
        new_genes  <- rownames(fit$coefficients)
        gene_order <- c(gene_order, setdiff(new_genes, gene_order))
      }
    }
    n <- nrow(motif.counts$bestmotif$p.post)
    if (length(gene_order) == n) {
      rownames(motif.counts$bestmotif$p.post)    <- gene_order
      rownames(motif.counts$bestmotif$clustlike) <- gene_order
    } else {
      warning("Gene name recovery failed: p.post has ", n,
              " rows but fitlist has ", length(gene_order), " unique genes.")
    }
  }

  saveRDS(motif.counts,
          file.path(output_dir, paste0(celltype_name, "_motif_fittedK2plus.rds")))


  plot_motif_ggplot(motif.counts, celltype_name, output_dir)
  plot_motif_publication(motif.counts, celltype_name, output_dir)
  return(motif.counts)
}


# Gene extraction (uses p.post rownames as gene universe).
# filterByExpr may drop genes in specific pairwise comparisons (they get NA in
# the t-matrix and are excluded from p.post). Use p.post rownames as the
# authoritative gene universe rather than the input count matrix rownames.
extract_genes_from_motif <- function(motif.fitted, motif_number, threshold = 0.5) {
  posterior_probs <- motif.fitted$bestmotif$clustlike
  gene_indices    <- which(posterior_probs[, motif_number] > threshold)
  rownames(motif.fitted$bestmotif$p.post)[gene_indices]
}


# topGO enrichment

run_topGO <- function(gene_symbols, all_genes_list, ontology = "BP") {
  print(paste("Running GO:", ontology))

  # allGenes is a named factor(0/1): 1 = gene is in the motif gene set.
  # geneSel = function(p) p == 1 selects factor entries equal to 1 (the motif genes).
  # This works by duck-typed R factor comparison; would need updating if topGO
  # changes its allGenes contract in a future version.
  GOdata <- new("topGOdata",
                ontology = ontology,
                allGenes = all_genes_list,
                geneSel  = function(p) p == 1,
                annot    = annFUN.org,
                mapping  = "org.Hs.eg.db",
                ID       = "symbol")

  result_classic <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  result_elim    <- runTest(GOdata, algorithm = "elim",    statistic = "fisher")

  # Retrieve ALL tested GO terms (not just top 50) so BH FDR is applied globally.
  # Previously topNodes = 50 meant FDR was corrected only within 50 terms, making
  # it invalid as a global correction.
  n_tested <- length(usedGO(GOdata))

  all_res <- GenTable(GOdata,
                      classicFisher = result_classic,
                      elimFisher    = result_elim,
                      orderBy       = "elimFisher",
                      topNodes      = n_tested,
                      numChar       = 200)

  parse_topgo_pval <- function(x) as.numeric(sub("^< *", "", as.character(x)))
  all_res$elimFisher    <- parse_topgo_pval(all_res$elimFisher)
  all_res$classicFisher <- parse_topgo_pval(all_res$classicFisher)

  return(all_res)
}


# Per-cell-type analysis

analyze_single_cell_type <- function(merged_matrix, cell_type, output_dir,
                                     protein_coding_genes) {
  cat("Processing:", cell_type, "\n")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Filter genes; returns raw counts. Gene filter log updated as side effect.
  raw_counts <- filter_genes_cormotif(merged_matrix, protein_coding_genes, cell_type)

  motif.fitted <- run_cormotif(raw_counts, cell_type, output_dir)

  # Gene universe from the motif object (accounts for filterByExpr on counts path)
  all_genes  <- rownames(motif.fitted$bestmotif$p.post)
  num_motifs <- ncol(motif.fitted$bestmotif$clustlike)
  go_results <- list()

  for (motif_num in 1:num_motifs) {
    cat("Analyzing motif:", motif_num, "\n")

    motif_genes <- extract_genes_from_motif(motif.fitted, motif_num)
    writeLines(motif_genes,
               file.path(output_dir,
                         paste0(cell_type, "_motif", motif_num, "_genes_K2plus.txt")))

    all_genes_list <- factor(as.integer(all_genes %in% motif_genes))
    names(all_genes_list) <- all_genes

    motif_key <- paste0("motif", motif_num)
    go_results[[motif_key]] <- list()

    for (ont in c("BP", "MF", "CC")) {
      go_df <- tryCatch(
        run_topGO(motif_genes, all_genes_list, ontology = ont),
        error = function(e) {
          cat("  topGO error for", ont, ":", conditionMessage(e), "\n")
          NULL
        }
      )
      go_results[[motif_key]][[ont]] <- go_df

      if (!is.null(go_df)) {
        write.csv(go_df,
                  file = file.path(output_dir,
                                   paste0(cell_type, "_motif", motif_num,
                                          "_GO_", ont, "_K2plus.csv")),
                  row.names = FALSE)
        # Individual per-motif panels are not saved here; replot_all() produces
        # combined per-cell-type PDFs (one per algorithm) from the CSVs above.
      }
    }
  }

  return(go_results)
}


# Load and merge pseudobulk

load_and_merge_pseudobulk <- function(pb_dir,
                                       drugs = c("Ale", "Bri", "Cer", "Cri", "Lor")) {
  cat("Loading per-drug pseudobulk matrices...\n")
  pb_list <- lapply(drugs, function(drug) {
    f <- file.path(pb_dir, paste0("pseudo_bulk_", drug, ".rds"))
    cat("  Reading:", f, "\n")
    readRDS(f)
  })
  names(pb_list) <- drugs

  cell_types <- Reduce(intersect, lapply(pb_list, names))
  cat("Cell types in all drugs:", length(cell_types), "\n")

  merged <- lapply(cell_types, function(ct) {
    mats <- lapply(drugs, function(drug) {
      m <- pb_list[[drug]][[ct]]
      colnames(m) <- paste0(drug, colnames(m))
      m
    })
    all_genes    <- Reduce(intersect, lapply(mats, rownames))
    mats_aligned <- lapply(mats, function(m) m[all_genes, , drop = FALSE])
    do.call(cbind, mats_aligned)
  })
  names(merged) <- cell_types
  return(merged)
}


# Run analysis for all cell types

analyze_all_cell_types <- function(merged_pb, output_base_dir, protein_coding_genes) {
  all_results <- list()

  for (cell_type in names(merged_pb)) {
    cat("Processing:", cell_type, "\n")

    tryCatch({
      output_dir <- file.path(output_base_dir, cell_type)
      go_results <- analyze_single_cell_type(merged_pb[[cell_type]], cell_type,
                                             output_dir, protein_coding_genes)
      all_results[[cell_type]] <- go_results
    }, error = function(e) {
      cat("Error in processing:", cell_type, "\n")
      print(e)
      all_results[[cell_type]] <<- list(error = e)
    })
  }

  return(all_results)
}


# Post-hoc helpers: re-run topGO from saved RDS, replot all figures,
# generate supplementary tables, generate motif pattern table.

rerun_topgo_from_saved <- function(output_base_dir) {
  cat("Re-running topGO from saved RDS + gene files in:", output_base_dir, "\n")
  cell_type_dirs <- list.dirs(output_base_dir, recursive = FALSE)

  for (ct_dir in cell_type_dirs) {
    cell_type <- basename(ct_dir)
    rds_file  <- file.path(ct_dir, paste0(cell_type, "_motif_fittedK2plus.rds"))
    if (!file.exists(rds_file)) {
      cat("  No RDS found, skipping:", cell_type, "\n")
      next
    }

    cat("Processing:", cell_type, "\n")
    motif.fitted <- readRDS(rds_file)
    all_genes    <- rownames(motif.fitted$bestmotif$p.post)
    K            <- ncol(motif.fitted$bestmotif$clustlike)

    for (m in seq_len(K)) {
      gene_file <- file.path(ct_dir, paste0(cell_type, "_motif", m, "_genes_K2plus.txt"))
      if (!file.exists(gene_file)) next
      motif_genes    <- readLines(gene_file)
      all_genes_list <- factor(as.integer(all_genes %in% motif_genes))
      names(all_genes_list) <- all_genes

      for (ont in c("BP", "MF", "CC")) {
        go_df <- tryCatch(
          run_topGO(motif_genes, all_genes_list, ontology = ont),
          error = function(e) {
            cat("  topGO error motif", m, ont, ":", conditionMessage(e), "\n")
            NULL
          }
        )
        if (!is.null(go_df)) {
          write.csv(go_df,
                    file      = file.path(ct_dir,
                                          paste0(cell_type, "_motif", m,
                                                 "_GO_", ont, "_K2plus.csv")),
                    row.names = FALSE)
          cat("  Saved:", cell_type, "motif", m, ont, "\n")
        }
      }
    }
  }
  cat("Done re-running topGO.\n")
}


replot_all <- function(output_base_dir, ontology = "BP",
                       algorithms = "elim", top_n = 8) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("patchwork is required for combined GO figures: install.packages('patchwork')")
  library(patchwork)

  cat("Regenerating plots from existing results in:", output_base_dir, "\n")
  cell_type_dirs <- list.dirs(output_base_dir, recursive = FALSE)

  for (ct_dir in cell_type_dirs) {
    cell_type <- basename(ct_dir)
    cat("Replotting:", cell_type, "\n")

    # Motif pattern heatmaps: once per cell type, not per algorithm.
    rds_file <- file.path(ct_dir, paste0(cell_type, "_motif_fittedK2plus.rds"))
    motif.fitted <- NULL
    if (file.exists(rds_file)) {
      tryCatch({
        motif.fitted <- readRDS(rds_file)
        plot_motif_ggplot(motif.fitted, cell_type, ct_dir)
        plot_motif_publication(motif.fitted, cell_type, ct_dir)
      }, error = function(e) cat("  Pattern plot error:", conditionMessage(e), "\n"))
    }

    # Discover GO CSVs once; reuse across both algorithms.
    go_csvs <- list.files(ct_dir,
                          pattern = paste0("_GO_", ontology, "_K2plus\\.csv$"),
                          full.names = TRUE)
    if (length(go_csvs) == 0) {
      cat("  No", ontology, "GO CSVs found, skipping GO figures.\n")
      next
    }

    motif_nums <- as.integer(sub(paste0(".*_motif(\\d+)_GO_", ontology, ".*"), "\\1",
                                 basename(go_csvs)))
    go_csvs    <- go_csvs[order(motif_nums)]
    motif_nums <- sort(motif_nums)

    # Pre-read GO dataframes and motif titles once, share across algorithms
    go_dfs <- lapply(go_csvs, function(f)
      tryCatch(read.csv(f), error = function(e) NULL))

    motif_metas <- lapply(motif_nums, function(m) {
      if (!is.null(motif.fitted))
        compute_motif_meta(motif.fitted$bestmotif$motif.q[m, ])
      else
        NULL
    })

    panel_titles <- lapply(seq_along(motif_nums), function(i) {
      m    <- motif_nums[i]
      meta <- motif_metas[[i]]
      if (!is.null(meta)) {
        n_genes <- sum(motif.fitted$bestmotif$clustlike[, m] > 0.5)
        paste0("Motif ", m, " | ", meta$motif_type, " | ", n_genes, " genes")
      } else {
        paste0("Motif ", m)
      }
    })

    # One combined PDF per algorithm.
    for (algorithm in algorithms) {
      panels <- list()
      for (i in seq_along(motif_nums)) {
        if (is.null(go_dfs[[i]])) next
        m <- motif_nums[i]

        # Skip null motifs — no active comparisons, GO enrichment reflects
        # non-responsive background rather than drug biology
        if (!is.null(motif_metas[[i]]) && motif_metas[[i]]$motif_type == "null") next

        p_base <- tryCatch(
          plot_go_panel(cell_type    = cell_type,
                        motif_number = m,
                        ontology     = ontology,
                        algorithm    = algorithm,
                        top_n        = top_n,
                        go_df        = go_dfs[[i]],
                        save_plot    = FALSE),
          error = function(e) {
            cat("  GO panel error motif", m, "(", algorithm, ") -",
                conditionMessage(e), "\n")
            NULL
          }
        )
        if (!is.null(p_base))
          panels[[length(panels) + 1]] <- p_base +
            ggtitle(panel_titles[[i]]) +
            theme(plot.title = element_text(size = 9, face = "plain"))
      }

      if (length(panels) == 0) {
        cat("  No significant GO terms (", algorithm, ") — skipping.\n")
        next
      }

      ncols    <- min(length(panels), 3L)
      nrows    <- ceiling(length(panels) / ncols)
      combined <- wrap_plots(panels, ncol = ncols) +
        plot_annotation(title    = cell_type,
                        subtitle = paste0(ontology, " enrichment  |  ", algorithm,
                                          " algorithm  |  p < 0.05"),
                        theme    = theme(plot.title    = element_text(face = "bold",
                                                                      size = 12),
                                         plot.subtitle = element_text(size = 9)))

      out_pdf <- file.path(ct_dir,
                           paste0(cell_type, "_GO_", ontology, "_",
                                  algorithm, "_all_motifs.pdf"))
      pdf_w <- ncols * 6
      pdf_h <- nrows * (2.5 + top_n * 0.35)
      tryCatch({
        cairo_pdf(out_pdf, width = pdf_w, height = pdf_h, family = "Helvetica")
        print(combined)
        dev.off()
        cat("  Saved:", basename(out_pdf), "\n")
      }, error = function(e) cat("  Combined figure error:", conditionMessage(e), "\n"))
    }
  }
  cat("Done replotting.\n")
}


generate_supplementary_tables <- function(output_base_dir) {
  cat("Generating supplementary tables for:", output_base_dir, "\n")

  cell_type_dirs <- list.dirs(output_base_dir, recursive = FALSE)
  s1_rows <- list()
  s2_rows <- list()
  go_rows  <- list(BP = list(), MF = list(), CC = list())

  for (ct_dir in cell_type_dirs) {
    cell_type <- basename(ct_dir)
    rds_file  <- file.path(ct_dir, paste0(cell_type, "_motif_fittedK2plus.rds"))
    if (!file.exists(rds_file)) next

    motif.fitted    <- readRDS(rds_file)
    posterior_probs <- motif.fitted$bestmotif$clustlike
    motif_q         <- motif.fitted$bestmotif$motif.q
    K           <- ncol(posterior_probs)
    total_genes <- nrow(posterior_probs)

    s1_row <- data.frame(cell_type = cell_type, K = K, total_genes = total_genes)
    for (m in seq_len(K)) {
      n <- sum(posterior_probs[, m] > 0.5)
      s1_row[[paste0("motif", m, "_n_genes")]]  <- n
      s1_row[[paste0("motif", m, "_pct_genes")]] <- round(100 * n / total_genes, 1)
    }
    s1_rows[[cell_type]] <- s1_row

    for (m in seq_len(K)) {
      gene_file <- file.path(ct_dir, paste0(cell_type, "_motif", m, "_genes_K2plus.txt"))
      if (file.exists(gene_file)) {
        genes <- readLines(gene_file)
        if (length(genes) > 0) {
          s2_rows[[paste0(cell_type, "_motif", m)]] <- data.frame(
            cell_type = cell_type, motif = m, gene = genes
          )
        }
      }

      meta <- compute_motif_meta(motif_q[m, ])

      for (ont in c("BP", "MF", "CC")) {
        go_file <- file.path(ct_dir,
                             paste0(cell_type, "_motif", m, "_GO_", ont, "_K2plus.csv"))
        if (!file.exists(go_file)) next
        go_df  <- read.csv(go_file)
        go_sig <- go_df[!is.na(go_df$elimFisher) & go_df$elimFisher < 0.05, ]
        if (nrow(go_sig) > 0) {
          go_sig$cell_type            <- cell_type
          go_sig$motif                <- m
          go_sig$motif_type           <- meta$motif_type
          go_sig$n_active_comparisons <- meta$n_active_comparisons
          go_sig$drug_pattern         <- meta$drug_pattern
          go_sig$active_comparisons   <- meta$active_comparisons
          leading <- c("cell_type", "motif", "motif_type", "n_active_comparisons",
                       "drug_pattern", "active_comparisons")
          drop   <- intersect(c("classicFisher", "elimFDR", "classicFDR"), colnames(go_sig))
          go_sig <- go_sig[, c(leading, setdiff(colnames(go_sig), c(leading, drop)))]
          go_rows[[ont]][[paste0(cell_type, "_motif", m)]] <- go_sig
        }
      }
    }
  }

  save_table <- function(df, stem) {
    write.csv(df, file.path(output_base_dir, paste0(stem, ".csv")), row.names = FALSE)
    cat("  Saved", paste0(stem, ".csv"), "\n")
    write_xlsx_safe(df, file.path(output_base_dir, paste0(stem, ".xlsx")))
  }

  if (length(s1_rows) > 0) save_table(bind_rows(s1_rows), "S1_motif_summary")
  if (length(s2_rows) > 0) save_table(bind_rows(s2_rows), "S2_gene_lists")

  ont_stems <- c(BP = "S3_GO_BP_significant",
                 MF = "S4_GO_MF_significant",
                 CC = "S5_GO_CC_significant")
  for (ont in c("BP", "MF", "CC")) {
    if (length(go_rows[[ont]]) > 0)
      save_table(bind_rows(go_rows[[ont]]), ont_stems[[ont]])
  }

  if (length(go_rows[["BP"]]) > 0) {
    bp_all <- bind_rows(go_rows[["BP"]])
    bp_compact <- do.call(rbind, lapply(
      split(bp_all, paste(bp_all$cell_type, bp_all$motif)), function(x)
        head(x[order(x$elimFisher), ], 3)
    ))
    rownames(bp_compact) <- NULL
    save_table(bp_compact, "S6_GO_BP_top3_per_motif")
  }
}


generate_pattern_table <- function(output_base_dir, motif_q_threshold = 0.5) {
  cell_type_dirs <- list.dirs(output_base_dir, recursive = FALSE)
  rows <- list()

  for (ct_dir in cell_type_dirs) {
    cell_type <- basename(ct_dir)
    rds_file  <- file.path(ct_dir, paste0(cell_type, "_motif_fittedK2plus.rds"))
    if (!file.exists(rds_file)) next

    motif.fitted    <- readRDS(rds_file)
    motif_q         <- motif.fitted$bestmotif$motif.q
    posterior_probs <- motif.fitted$bestmotif$clustlike
    K           <- nrow(motif_q)
    total_genes <- nrow(posterior_probs)

    for (m in seq_len(K)) {
      meta    <- compute_motif_meta(motif_q[m, ], threshold = motif_q_threshold)
      n_genes <- sum(posterior_probs[, m] > 0.5)

      read_top_go <- function(ont) {
        f <- file.path(ct_dir, paste0(cell_type, "_motif", m, "_GO_", ont, "_K2plus.csv"))
        if (!file.exists(f)) return(list(terms = "", elim_p = ""))
        go_df  <- read.csv(f)
        go_sig <- go_df[!is.na(go_df$elimFisher) & go_df$elimFisher < 0.05, ]
        if (nrow(go_sig) == 0) return(list(terms = "", elim_p = ""))
        go_top <- head(go_sig[order(go_sig$elimFisher), ], 3)
        list(terms  = paste(go_top$Term, collapse = " | "),
             elim_p = paste(signif(go_top$elimFisher, 2), collapse = " | "))
      }

      bp <- read_top_go("BP")
      mf <- read_top_go("MF")
      cc <- read_top_go("CC")

      gene_file <- file.path(ct_dir,
                             paste0(cell_type, "_motif", m, "_genes_K2plus.txt"))
      gene_list_str <- ""
      if (!meta$is_null && file.exists(gene_file))
        gene_list_str <- paste(readLines(gene_file), collapse = "; ")

      rows[[paste(cell_type, m)]] <- data.frame(
        cell_type            = cell_type,
        motif                = m,
        n_genes              = n_genes,
        pct_genes            = round(100 * n_genes / total_genes, 1),
        is_null              = meta$is_null,
        motif_type           = meta$motif_type,
        n_active_comparisons = meta$n_active_comparisons,
        drug_pattern         = meta$drug_pattern,
        active_comparisons   = meta$active_comparisons,
        top_BP_terms         = bp$terms,
        top_BP_elim_p        = bp$elim_p,
        top_MF_terms         = mf$terms,
        top_MF_elim_p        = mf$elim_p,
        top_CC_terms         = cc$terms,
        top_CC_elim_p        = cc$elim_p,
        gene_list            = gene_list_str,
        stringsAsFactors     = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    cat("No results found in:", output_base_dir, "\n")
    return(invisible(NULL))
  }

  out <- bind_rows(rows)
  write.csv(out, file.path(output_base_dir, "motif_pattern_table.csv"), row.names = FALSE)
  cat("Saved motif_pattern_table.csv:", nrow(out), "rows\n")
  write_xlsx_safe(out, file.path(output_base_dir, "motif_pattern_table.xlsx"))
  invisible(out)
}


# Mixed-motif investigation tooling.
# Characterizes non-drug-specific motifs (mixed AND pan-comparison) by asking:
# for each drug, do the 4 comparisons involving that drug tend to have higher
# motif.q values than the 6 comparisons that don't?
#
# Included motifs: null and drug-specific are skipped. Mixed AND pan-comparison
# are retained. Pan-comparison (all 10 active) is intentionally included because
# even an all-active motif can be asymmetric - one drug's comparisons may be
# consistently more confidently active (higher raw q) than others.
#
# Outputs:
#   mixed_motif_raw_q.{csv,xlsx}  - one row per motif: raw q for all 10
#                                   comparisons, for interactive inspection

compute_mixed_motif_asymmetry <- function(output_base_dir, motif_q_threshold = 0.5) {
  drug_comp_idx  <- setNames(lapply(drugs, function(d)
    which(comparison_drug1 == d | comparison_drug2 == d)), drugs)
  comp_col_names <- paste0("q_", gsub(" vs ", "_", comparison_labels))

  asym_rows <- list()
  raw_rows  <- list()

  for (ct_dir in list.dirs(output_base_dir, recursive = FALSE)) {
    cell_type <- basename(ct_dir)
    rds_file  <- file.path(ct_dir, paste0(cell_type, "_motif_fittedK2plus.rds"))
    if (!file.exists(rds_file)) next

    motif.fitted <- readRDS(rds_file)
    motif_q      <- motif.fitted$bestmotif$motif.q
    posterior    <- motif.fitted$bestmotif$clustlike

    for (m in seq_len(nrow(motif_q))) {
      q_row <- motif_q[m, ]
      meta  <- compute_motif_meta(q_row, threshold = motif_q_threshold)

      # Skip null motifs (nothing to analyze) and drug-specific motifs (their
      # asymmetry is fully captured by which drug the motif belongs to).
      # Keep mixed AND pan-comparison (drug_pattern == "mixed" for both).
      # Note: a "multi-drug-combo" drug_pattern (e.g. "Ale+Bri") is mathematically
      # impossible — see compute_motif_meta() — so this filter is exhaustive.
      if (meta$drug_pattern %in% c("null") || endsWith(meta$drug_pattern, "-specific")) next

      n_genes <- sum(posterior[, m] > 0.5)

      # Raw q table: one row per motif with all 10 raw q-values.
      # Useful for visualizing the full q-vector alongside the asymmetry summary.
      raw_row <- setNames(as.data.frame(t(q_row)), comp_col_names)
      raw_rows[[paste(cell_type, m)]] <- cbind(
        data.frame(cell_type = cell_type, motif = m, n_genes = n_genes,
                   motif_type = meta$motif_type,
                   stringsAsFactors = FALSE),
        raw_row)

      # Asymmetry table: for each drug D, compare mean q across:
      #   - the 4 comparisons involving D  (e.g. Ale vs Bri, Ale vs Cer, Ale vs Cri,
      #     Ale vs Lor for D = Ale)
      #   - the 6 comparisons NOT involving D
      # asymmetry_score = mean_q_involving - mean_q_not_involving.
      # A positive score for drug D means D's comparisons are more confidently active
      # in this motif. This does NOT imply drug-specificity — other comparisons are
      # also active (otherwise the motif would be drug-specific, not mixed).
      # q-values are EM probabilities (0–1), not p-values; two motifs with the same
      # binary active/inactive pattern can differ in asymmetry.
      for (d in drugs) {
        inv  <- drug_comp_idx[[d]]
        non  <- setdiff(seq_along(q_row), inv)
        asym_rows[[paste(cell_type, m, d)]] <- data.frame(
          cell_type            = cell_type,
          motif                = m,
          motif_type           = meta$motif_type,
          drug                 = d,
          n_genes              = n_genes,
          mean_q_involving     = round(mean(q_row[inv]), 4),
          mean_q_not_involving = round(mean(q_row[non]), 4),
          asymmetry_score      = round(mean(q_row[inv]) - mean(q_row[non]), 4),
          stringsAsFactors     = FALSE
        )
      }
    }
  }

  n_mixed <- length(raw_rows)
  if (n_mixed == 0) {
    cat("No mixed motifs found in:", output_base_dir, "\n")
    return(invisible(NULL))
  }

  asym_df <- bind_rows(asym_rows)
  raw_df  <- bind_rows(raw_rows)

  cat("Mixed motifs found:", n_mixed,
      "across cell types:", paste(unique(raw_df$cell_type), collapse = ", "),
      "\n")

  save_mixed <- function(df, stem) {
    write.csv(df, file.path(output_base_dir, paste0(stem, ".csv")), row.names = FALSE)
    write_xlsx_safe(df, file.path(output_base_dir, paste0(stem, ".xlsx")))
  }
  save_mixed(raw_df,  "mixed_motif_raw_q")

  invisible(list(asymmetry = asym_df, raw_q = raw_df))
}

# Supplementary figure assembly (Supp Fig S7).
# Builds multi-cell-type PNG pages from per-cell-type GO CSVs using ggplot2 +
# patchwork (no image-processing dependencies).
#
# Layout: colored vertical band (cell type label, rotated 90 deg) on the left
# of each cell type's panel block; cell types stacked vertically; n per page.
# Large cell types (> solo_threshold motifs) get their own page; smaller ones
# are paired up.
#
# Outputs: <output_base_dir>/GO_BP_elim_combined/supp_GO_BP_elim_page01.png, ...

assemble_go_supplementary_pages <- function(
    output_base_dir,
    output_subdir   = "GO_BP_elim_combined",
    n_per_page      = 2,       # max cell types per page (for small cell types)
    solo_threshold  = 4,       # cell types with > this many motifs get their own page
    dpi             = 200,
    ontology        = "BP",
    algorithm       = "elim",
    top_n           = 8,
    page_width        = 10,    # inches — each panel gets ~4.8" at 2 cols
    page_height       = 9.5,   # inches minimum
    panel_base_size   = 9,     # pt — smaller than standalone (12) to fit narrow panels
    panel_wrap_width  = 50,    # strwrap chars — at 9pt in ~3.4" panels, 50 chars fits one line
    panel_row_height  = 2.8,   # inches per row of 2 panels; page height scales with this
    cell_type_order   = NULL,  # optional character vector to reorder cell types
    group_colors      = NULL   # optional named vector: cell_type -> hex color
                               # e.g. setNames(germ_layer_colors[cell_type_germ_layer],
                               #               names(cell_type_germ_layer))
) {
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("patchwork is required: install.packages('patchwork')")
  library(patchwork)

  output_dir <- file.path(output_base_dir, output_subdir)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Discover cell types that have GO CSV results
  ct_dirs    <- list.dirs(output_base_dir, recursive = FALSE)
  has_go     <- vapply(ct_dirs, function(d)
    length(list.files(d, pattern = paste0("_GO_", ontology, "_K2plus\\.csv$"))) > 0,
    logical(1))
  ct_dirs    <- ct_dirs[has_go]
  cell_types <- basename(ct_dirs)

  if (length(cell_types) == 0)
    stop("No GO CSVs found in: ", output_base_dir,
         "\nRun analyze_single_cell_type() or rerun_topgo_from_saved() first.")

  # Apply ordering if provided
  if (!is.null(cell_type_order)) {
    ordered    <- cell_type_order[cell_type_order %in% cell_types]
    cell_types <- c(ordered, setdiff(cell_types, ordered))
  }
  ct_dirs <- file.path(output_base_dir, cell_types)

  # Band color per cell type
  if (!is.null(group_colors)) {
    missing_ct <- setdiff(cell_types, names(group_colors))
    if (length(missing_ct) > 0)
      warning("group_colors missing for: ", paste(missing_ct, collapse = ", "),
              " — using grey50.")
    band_colors <- setNames(
      ifelse(cell_types %in% names(group_colors), group_colors[cell_types], "grey50"),
      cell_types)
  } else {
    band_colors <- setNames(
      rep_len(c("#4E79A7","#F28E2B","#E15759","#76B7B2","#59A14F",
                "#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC"),
              length(cell_types)),
      cell_types)
  }

  # Pre-compute panel rows per cell type for dynamic page height.
  # Count GO CSVs (proxy for non-null motifs); ceiling(n / 2) rows at 2 cols.
  n_panel_rows <- vapply(ct_dirs, function(d) {
    n <- length(list.files(d, pattern = paste0("_GO_", ontology, "_K2plus\\.csv$")))
    ceiling(max(n, 1L) / 2L)
  }, numeric(1))

  # Colored label strip: annotate("rect") fills as a data layer so patchwork
  # cannot clear it; panel.background + plot.background both set for safety.
  make_strip <- function(label, color) {
    ggplot() +
      annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
               fill = color, color = NA) +
      annotate("text", x = 0.5, y = 0.5, label = label,
               angle = 90, color = "white", fontface = "bold",
               size = 5, hjust = 0.5, vjust = 0.5) +
      scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      theme_void() +
      theme(panel.background = element_rect(fill = color, color = NA),
            plot.background  = element_rect(fill = color, color = "white", linewidth = 1),
            plot.margin      = margin(0, 0, 0, 0))
  }

  # Build one cell type block:
  #   left: narrow colored strip (full height, rotated label)
  #   right: motif panels in a 3-column-wrapped grid (same layout as replot_all)
  build_ct_row <- function(ct_dir, cell_type, color) {
    rds_file     <- file.path(ct_dir, paste0(cell_type, "_motif_fittedK2plus.rds"))
    motif.fitted <- if (file.exists(rds_file))
      tryCatch(readRDS(rds_file), error = function(e) NULL) else NULL

    go_csvs <- list.files(ct_dir,
                          pattern = paste0("_GO_", ontology, "_K2plus\\.csv$"),
                          full.names = TRUE)
    if (length(go_csvs) == 0) return(NULL)

    motif_nums <- as.integer(sub(paste0(".*_motif(\\d+)_GO_", ontology, ".*"),
                                 "\\1", basename(go_csvs)))
    go_csvs    <- go_csvs[order(motif_nums)]
    motif_nums <- sort(motif_nums)

    panels <- list()
    for (i in seq_along(motif_nums)) {
      m     <- motif_nums[i]
      go_df <- tryCatch(read.csv(go_csvs[i]), error = function(e) NULL)
      if (is.null(go_df)) next

      if (!is.null(motif.fitted)) {
        meta <- compute_motif_meta(motif.fitted$bestmotif$motif.q[m, ])
        if (meta$motif_type == "null") next
        n_genes <- sum(motif.fitted$bestmotif$clustlike[, m] > 0.5)
        title <- paste0("Motif ", m, " | ", meta$motif_type, "\n", n_genes, " genes")
      } else {
        title <- paste0("Motif ", m)
      }

      p <- tryCatch(
        plot_go_panel(cell_type = cell_type, motif_number = m,
                      ontology = ontology, algorithm = algorithm,
                      top_n = top_n, go_df = go_df, save_plot = FALSE,
                      base_size  = panel_base_size,
                      wrap_width = panel_wrap_width),
        error = function(e) NULL
      )
      if (!is.null(p))
        panels[[length(panels) + 1]] <- p +
          ggtitle(title) +
          theme(plot.title = element_text(size = panel_base_size, face = "plain"))
    }
    if (length(panels) == 0) return(NULL)

    # Max 2 panels per row so each panel gets ~3.9 inches at page_width = 8.5
    ncols      <- min(length(panels), 2L)
    panel_grid <- wrap_plots(panels, ncol = ncols)
    strip <- make_strip(cell_type, color)
    wrap_plots(list(strip, panel_grid), nrow = 1, widths = c(3, 30))
  }

  # Dynamic page grouping:
  # - Large cell types (> solo_threshold motifs) → solo page
  # - Small cell types (≤ solo_threshold motifs) → paired across the full list,
  #   so a small cell type stranded next to a large one still gets a partner.
  # Pages are then sorted by their smallest cell-type index to preserve order.
  is_solo   <- (n_panel_rows * 2L) > solo_threshold
  solo_pages  <- as.list(which(is_solo))
  small_idx   <- which(!is_solo)
  small_pages <- if (length(small_idx) > 0)
    split(small_idx, ceiling(seq_along(small_idx) / n_per_page))
  else list()
  page_groups <- c(solo_pages, small_pages)
  page_groups <- page_groups[order(vapply(page_groups, min, numeric(1)))]

  cat("Assembling ", length(cell_types), " cell types into ",
      length(page_groups), " pages (",
      sum(lengths(page_groups) == 1L), " solo, ",
      sum(lengths(page_groups) > 1L), " paired)...\n", sep = "")

  for (pg in seq_along(page_groups)) {
    idx <- page_groups[[pg]]
    blocks <- list()

    block_rows <- integer(0)  # panel row count per built block (parallel to blocks)
    for (i in idx) {
      ct <- cell_types[i]
      b  <- tryCatch(
        build_ct_row(ct_dirs[i], ct, band_colors[[ct]]),
        error = function(e) { cat("  Error:", ct, "-", conditionMessage(e), "\n"); NULL }
      )
      if (!is.null(b)) {
        blocks[[length(blocks) + 1]] <- b
        block_rows <- c(block_rows, n_panel_rows[i])
      }
    }
    if (length(blocks) == 0) next

    # heights = block_rows ensures each cell type gets vertical space proportional
    # to its number of panel rows (not equal halves, which squishes dense cell types).
    page_plot <- wrap_plots(blocks, ncol = 1, heights = block_rows)

    pg_w <- page_width
    # Dynamic height: each panel row needs panel_row_height inches.
    # page_height is the minimum (caps short pages to a sensible floor).
    pg_h <- max(page_height, sum(block_rows) * panel_row_height)

    out_png <- file.path(output_dir,
                         sprintf("supp_GO_%s_%s_page%02d.png",
                                 ontology, algorithm, pg))
    tryCatch({
      png(out_png, width = pg_w, height = pg_h, units = "in", res = dpi)
      print(page_plot)
      dev.off()
      cat("  Saved:", basename(out_png), "\n")
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      cat("  Save error page", pg, ":", conditionMessage(e), "\n")
    })
  }
  cat("Done. PNGs written to:", output_dir, "\n")
  invisible(NULL)
}


# Main

pb_dir <- "/project/gilad/geena/publication_clean/results/treatment/cormotif/pseudobulk_15cells_4indiv"

hg_protein_coding_genes <- readLines(
  "/project/gilad/geena/publication_clean/reference_sets/hg_protein_coding_genes.tsv"
)

merged_pb <- load_and_merge_pseudobulk(pb_dir)

output_base_dir <- "/project/gilad/geena/publication_clean/results/treatment/cormotif/proteincoding_newfilter_15cells_4arm"
dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)

# gene_filter_log is populated as a side effect of filter_genes_cormotif() via <<-
# It must be initialized before any analyze_single_cell_type() call.
#gene_filter_log_single <- list()

# Test a single cell type
# analyze_single_cell_type(merged_pb[["Cardiac mesoderm cells"]], "Cardiac mesoderm cells",
#                             file.path(output_base_dir, "Cardiac mesoderm cells"), hg_protein_coding_genes)
#write.csv(dplyr::bind_rows(gene_filter_log_single),
#               file.path(output_base_dir, "gene_filter_log_single.csv"),
#               row.names = FALSE)
#     cat("Saved gene_filter_log_single.csv\n")

 
   
# Run for all cell types
#gene_filter_log <- list()
#all_results <- analyze_all_cell_types(merged_pb, output_base_dir, protein_coding_genes = hg_protein_coding_genes)
#write.csv(dplyr::bind_rows(gene_filter_log),
#            file.path(output_base_dir, "gene_filter_log.csv"),
#            row.names = FALSE)
#cat("Saved gene_filter_log.csv\n")

rerun_topgo_from_saved(output_base_dir)

generate_supplementary_tables(output_base_dir)
generate_pattern_table(output_base_dir)
compute_mixed_motif_asymmetry(output_base_dir)

replot_all(output_base_dir)

assemble_go_supplementary_pages(
  output_base_dir,
  cell_type_order = names(cell_type_germ_layer),
  group_colors    = setNames(germ_layer_colors[cell_type_germ_layer],
                             names(cell_type_germ_layer))
)

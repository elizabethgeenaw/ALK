# Main Fig 3A: side-by-side "motif pattern" and "genes per motif" panels for
# two exemplar cell types (Hepatoblasts and Cardiac mesoderm cells). The
# motif panel shows a K-motifs x 10-comparisons heatmap where each cell is 1
# if that comparison is active in that motif (q >= 0.5), 0 otherwise, and any
# motif classified as X-specific has its active cells recolored in drug X's
# color. The barcount panel next to it shows how many genes are assigned to
# each motif (posterior probability > 0.5) with a matching drug-specific
# annotation.

library(ggplot2)

output_base_dir <- "/project/gilad/geena/publication_clean/results/treatment/cormotif/cormotif_15cells_4arm"

# The 10 pairwise comparisons in the fit's canonical order, and the drug on
# each side of each comparison. drug_colors is the per-drug palette shared
# across all figures.
comparison_labels <- c(
  "Ale vs Bri", "Ale vs Cer", "Ale vs Cri", "Ale vs Lor",
  "Bri vs Cer", "Bri vs Cri", "Bri vs Lor",
  "Cer vs Cri", "Cer vs Lor", "Cri vs Lor"
)
comparison_drug1 <- c("Ale","Ale","Ale","Ale","Bri","Bri","Bri","Cer","Cer","Cri")
comparison_drug2 <- c("Bri","Cer","Cri","Lor","Cer","Cri","Lor","Cri","Lor","Lor")
drug_colors <- c(Ale = "#0072B2", Bri = "#D55E00", Cer = "#009E73",
                 Cri = "#CC79A7", Lor = "#E69F00")
drugs <- names(drug_colors)

# Classify a motif as drug X-specific iff its activity pattern across the 10
# comparisons is exactly "active in all X-containing comparisons, inactive
# elsewhere". Any deviation (e.g. active in one X and one non-X comparison)
# returns NA and is treated as mixed by downstream code.
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
  apply(motif.fitted$bestmotif$motif.q, 1, classify_motif, threshold = threshold)
}

# Render the K x 10 motif pattern panel using base graphics image(). When
# binary = TRUE we threshold q at 0.5 and use a 2-color palette; otherwise
# we render the continuous q values in grayscale. Drug-specific rows get
# their active cells recolored in the drug's color as an overlay.
plot_motif_panel <- function(motif.fitted, celltype_name, specificity, binary = FALSE) {
  u  <- 1:dim(motif.fitted$bestmotif$motif.q)[2]
  v  <- 1:dim(motif.fitted$bestmotif$motif.q)[1]

  motif_q <- motif.fitted$bestmotif$motif.q
  if (binary) {
    plot_q  <- (motif_q >= 0.5) * 1L
    col_vec <- c("white", "black")
  } else {
    plot_q  <- motif_q
    col_vec <- gray(seq(from = 1, to = 0, by = -0.1))
  }

  image(u, v, t(plot_q),
        col  = col_vec,
        xaxt = "n", yaxt = "n",
        xlab = "", ylab = "Corr. Motifs",
        main = paste(celltype_name, "pattern"),
        cex.lab = 1.8, cex.main = 2.0)

  # Recolor active cells of drug-specific motifs in the drug's color.
  for (i in seq_along(specificity)) {
    if (!is.na(specificity[i])) {
      drug_col    <- drug_colors[specificity[i]]
      active_cols <- which(motif_q[i, ] >= 0.5)
      for (j in active_cols) {
        rect(j - 0.5, i - 0.5, j + 0.5, i + 0.5,
             col = drug_col, border = NA)
      }
    }
  }

  # Roman numerals label the motif rows on the y-axis (matches the
  # cormotif publication convention).
  axis(2, at = seq_along(v), labels = as.roman(seq_along(v)),
       las = 1, cex.axis = 1.8)
  axis(1, at = u, labels = FALSE)

  # Draw the two-row Drug1/Drug2 color chip legend below the heatmap. usr
  # gives us the plot data-coordinate boundaries; sh (~7% of the y range)
  # is one chip height.
  usr    <- par("usr")
  yrange <- diff(usr[3:4])
  sh     <- yrange * 0.07

  text(x = u[1] - 0.6, y = usr[3] - sh / 2,
       labels = "Drug 1", adj = c(1, 0.5), xpd = TRUE, cex = 1.4)
  text(x = u[1] - 0.6, y = usr[3] - 1.5 * sh,
       labels = "Drug 2", adj = c(1, 0.5), xpd = TRUE, cex = 1.4)

  rect(u - 0.5, usr[3] - sh, u + 0.5, usr[3],
       col = drug_colors[comparison_drug1], border = NA, xpd = TRUE)
  text(x = u, y = usr[3] - sh / 2,
       labels = comparison_drug1, col = "white",
       srt = 0, adj = c(0.5, 0.5), xpd = TRUE, cex = 1.3, font = 2)

  rect(u - 0.5, usr[3] - 2 * sh, u + 0.5, usr[3] - sh,
       col = drug_colors[comparison_drug2], border = NA, xpd = TRUE)
  text(x = u, y = usr[3] - 1.5 * sh,
       labels = comparison_drug2, col = "white",
       srt = 0, adj = c(0.5, 0.5), xpd = TRUE, cex = 1.3, font = 2)

  # Grid lines separating cells.
  for (i in 1:(length(u) + 1)) abline(v = (i - 0.5))
  for (i in 1:(length(v) + 1)) abline(h = (i - 0.5))
}

# Horizontal bar chart panel: number of genes assigned to each motif
# (posterior probability > 0.5). Also draws the "X-specific" annotation
# at a constant x so labels line up cleanly across rows.
plot_barcount_panel <- function(motif.fitted, specificity) {
  Ng           <- nrow(motif.fitted$bestmotif$p.post)
  NK           <- nrow(motif.fitted$bestmotif$motif.q)
  genecount    <- colSums(motif.fitted$bestmotif$clustlike > 0.5)
  motif_p_post <- genecount / Ng

  # Wider xlim leaves room for the drug-specific label at spec_x = 1.15.
  plot(0, 0.7, pch = ".", xlim = c(0, 2.3), ylim = c(0.75, NK + 0.25),
       frame.plot = FALSE, axes = FALSE, xlab = "", ylab = "",
       main = "Number of genes in motif", cex.main = 1.8)
  segments(0, 0.7, motif_p_post[1], 0.7)
  rect(0, 1:NK - 0.3, motif_p_post, 1:NK + 0.3, col = "dark grey")
  mtext(as.character(as.roman(1:NK)), at = 1:NK, side = 2, las = 1, cex = 1.7)
  text(motif_p_post + 0.05, 1:NK, labels = format(genecount, big.mark = ","),
       cex = 1.7, adj = c(0, 0.5))

  spec_x <- 1.15
  for (i in seq_len(NK)) {
    if (!is.na(specificity[i])) {
      text(spec_x, i,
           labels = paste0(specificity[i], "-specific"),
           col = drug_colors[specificity[i]],
           cex = 2.0, font = 2, adj = c(0, 0.5))
    }
  }
}

# Load fitted motifs for the two exemplar cell types.
ct1_name <- "Hepatoblasts"
ct2_name <- "Cardiac mesoderm cells"
rds1 <- file.path(output_base_dir, ct1_name, paste0(ct1_name, "_motif_fittedK2plus.rds"))
rds2 <- file.path(output_base_dir, ct2_name, paste0(ct2_name, "_motif_fittedK2plus.rds"))

motif1 <- readRDS(rds1)
motif2 <- readRDS(rds2)

spec1 <- motif_specificity(motif1)
spec2 <- motif_specificity(motif2)

cat("Hepatoblasts motif specificity (I -> K):  ",
    paste(ifelse(is.na(spec1), "-", spec1), collapse = ", "), "\n")
cat("Cardiac mesoderm motif specificity (I -> K): ",
    paste(ifelse(is.na(spec2), "-", spec2), collapse = ", "), "\n")

# Compose the 4-panel figure: motif | barcount | motif | barcount, with
# custom margins so the drug1/drug2 chip legend below each motif panel
# has room.
draw_combined <- function(binary, device) {
  layout(matrix(1:4, nrow = 1), widths = c(2, 1.5, 2, 1.5))
  par(mar = c(7, 4, 4, 1)); plot_motif_panel(motif1, ct1_name, spec1, binary = binary)
  par(mar = c(7, 1, 4, 2)); plot_barcount_panel(motif1, spec1)
  par(mar = c(7, 4, 4, 1)); plot_motif_panel(motif2, ct2_name, spec2, binary = binary)
  par(mar = c(7, 1, 4, 2)); plot_barcount_panel(motif2, spec2)
}

# Emit binary-fill PNG + PDF versions (main figure uses the binary variant).
for (binary in c(TRUE)) {
  suffix <- if (binary) "_binary" else ""

  out_png <- file.path(output_base_dir,
                       paste0("cardiac_mesoderm_hepatoblasts_combined_spec", suffix, ".png"))
  out_pdf <- sub("\\.png$", ".pdf", out_png)

  png(out_png, width = 18, height = 7, units = "in", res = 300)
  draw_combined(binary)
  dev.off()

  pdf(out_pdf, width = 18, height = 7)
  draw_combined(binary)
  dev.off()

  cat("Wrote:", out_png, "\n       ", out_pdf, "\n")
}

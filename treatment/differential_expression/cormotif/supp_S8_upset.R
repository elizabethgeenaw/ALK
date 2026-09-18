# Supp Fig S8: per-drug UpSet plots of drug-specific motif genes, showing
# how much of the drug-specific gene set is shared across the cell types
# where that drug has a specific motif. One panel per drug that has at least
# 2 drug-specific motifs (Ale, Cer, Cri, Lor in the current fit; Bri has none).
#
# For each drug we collect the gene list of every cell type in which the drug
# is classified as drug-specific in motif_pattern_table.csv, build a binary
# membership matrix, and render the UpSet plot with the drug's color.

library(UpSetR)

# Paths and per-drug colors matching the rest of the manuscript.
cormotif_dir <- "/project/gilad/geena/publication_clean/results/treatment/cormotif/proteincoding_newfilter_15cells_4arm"
fig_dir      <- file.path(cormotif_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

drug_colors <- c(Ale = "#0072B2", Bri = "#D55E00", Cer = "#009E73",
                 Cri = "#CC79A7", Lor = "#E69F00")

# Pattern table is the source of truth for which (cell_type, motif) pairs
# are drug-specific. We derive the drug list from motif_type rather than
# hardcoding it, so an update to the fit that adds/removes drug-specific
# motifs propagates automatically.
pat <- read.csv(file.path(cormotif_dir, "motif_pattern_table.csv"),
                stringsAsFactors = FALSE)

drug_specific <- pat[grepl("-specific$", pat$motif_type), ]
drugs         <- sort(unique(sub("-specific$", "", drug_specific$motif_type)))

# For each drug: gather gene lists for every drug-specific motif, then draw
# the UpSet plot. Drugs with <2 cell types would produce a degenerate plot
# and are skipped.
for (drug in drugs) {
  rows <- drug_specific[drug_specific$motif_type == paste0(drug, "-specific"), ]

  if (nrow(rows) < 2) {
    message("Skipping ", drug, ": fewer than 2 cell types"); next
  }

  # Each row of `rows` points to a per-motif gene list on disk. We read them
  # all and drop any that failed to read.
  gene_lists <- lapply(seq_len(nrow(rows)), function(i) {
    ct  <- rows$cell_type[i]
    mot <- rows$motif[i]
    f   <- file.path(cormotif_dir, ct,
                     paste0(ct, "_motif", mot, "_genes_K2plus.txt"))
    if (!file.exists(f)) { warning("Not found: ", f); return(character(0)) }
    readLines(f)
  })
  names(gene_lists) <- rows$cell_type
  gene_lists <- gene_lists[lengths(gene_lists) > 0]

  if (length(gene_lists) < 2) { message("Skipping ", drug, ": < 2 files"); next }

  # Shorten " cells"/" cell" suffixes on cell-type names so the UpSet set
  # labels stay readable at plot size.
  names(gene_lists) <- gsub(" cells$| cell$", "", names(gene_lists))

  n_sets    <- length(gene_lists)
  n_genes   <- length(unique(unlist(gene_lists)))
  all_genes <- unique(unlist(gene_lists))

  # Binary presence matrix (rows = genes, columns = cell types).
  # check.names = FALSE keeps spaces in column names so the `sets =`
  # argument to upset() can match exactly.
  bin_mat  <- as.data.frame(
    lapply(gene_lists, function(g) as.integer(all_genes %in% g)),
    check.names = FALSE
  )
  ct_names <- colnames(bin_mat)

  # Sizing: cap the number of intersections shown, and scale width/height
  # with the number of sets so labels don't collide.
  nintersects <- min(40, 2^n_sets - 1)
  fig_w <- max(8, n_sets * 0.9)
  fig_h <- max(5, n_sets * 0.45 + 3)

  out_png <- file.path(fig_dir, paste0("drug_specific_upset_", drug, ".png"))

  # UpSetR uses base graphics. Wrapping in print() forces render in batch mode,
  # and PNG output avoids the blank first page that pdf() produces with UpSetR.
  png(out_png, width = fig_w, height = fig_h, units = "in", res = 200)
  print(upset(
    bin_mat,
    sets            = ct_names,
    nsets           = n_sets,
    nintersects     = nintersects,
    order.by        = "freq",
    mb.ratio        = c(0.6, 0.4),
    sets.bar.color  = drug_colors[drug],
    point.size      = 2.5,
    line.size       = 0.8,
    text.scale      = c(1.3, 1.2, 1.1, 1.0, 1.2, 1.0),
    mainbar.y.label = paste0(drug, "-specific genes in intersection"),
    sets.x.label    = paste0("Genes in motif  (", n_genes, " in union)")
  ))
  dev.off()

  cat("Saved:", out_png, "\n")

  # Per-drug summary to the log: how much of the union is unique to one cell type
  # (i.e. singleton bars in the UpSet plot).
  n_unique <- sum(sapply(ct_names, function(nm) {
    others <- unique(unlist(gene_lists[names(gene_lists) != nm]))
    sum(!gene_lists[[nm]] %in% others)
  }))
  cat(drug, "-specific: ", n_sets, " cell types, ", n_genes, " genes in union, ",
      n_unique, " (", round(100 * n_unique / n_genes, 1), "%) unique to one cell type\n",
      sep = "")
  for (nm in ct_names)
    cat("  ", nm, ":", length(gene_lists[[nm]]), "genes\n")
  cat("\n")
}

cat("Done. Figures in:", fig_dir, "\n")

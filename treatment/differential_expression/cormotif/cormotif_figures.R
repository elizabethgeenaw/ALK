# Sourced by fig3B_3C_refresh.R to render main Fig 3B, 3C, and 3D. The file
# is split into two sections: helper + panel-function definitions above the
# trailing `cat("===...` banner, and a runner block below it that calls the
# three fig3 panel functions. fig3B_3C_refresh.R uses that banner as a
# sentinel to eval only the definitions above (so it can skip running heavy
# panels like the rrvgo enrichment).

library(ggplot2)
library(dplyr)
library(tidyr)

# The 10 pairwise comparisons (Ale vs Bri, Ale vs Cer, ...) in the fixed
# order used by the cormotif fit. comparison_drug1/2 give the two drugs on
# each side; used to look up "which comparisons involve drug X".
comparison_labels <- c(
  "Ale vs Bri", "Ale vs Cer", "Ale vs Cri", "Ale vs Lor",
  "Bri vs Cer", "Bri vs Cri", "Bri vs Lor",
  "Cer vs Cri", "Cer vs Lor", "Cri vs Lor"
)
comparison_drug1 <- c("Ale","Ale","Ale","Ale","Bri","Bri","Bri","Cer","Cer","Cri")
comparison_drug2 <- c("Bri","Cer","Cri","Lor","Cer","Cri","Lor","Cri","Lor","Lor")
drugs <- c("Ale", "Bri", "Cer", "Cri", "Lor")

drug_colors <- c(Ale = "#0072B2", Bri = "#D55E00", Cer = "#009E73",
                 Cri = "#CC79A7", Lor = "#E69F00")

# Non-drug motif-type palette (used across panels for consistency).
# Drug-specific motifs get their own drug's color, added dynamically by
# motif_type_color() below. "group-specific" motifs (active in all 6
# between-group pairwise comparisons for a 2-drug grouping and in none of
# the 4 within-group comparisons) use a single purple so both the Cer+Lor
# and Bri+Cri variants render consistently.
motif_type_colors <- c(
  "pan-comparison" = "#2166AC",
  "group-specific" = "#5F3C8F",
  "mixed"          = "#999999",
  "null"           = "#E8E8E8"
)

output_base_dir <- "/project/gilad/geena/publication_clean/results/treatment/cormotif/cormotif_15cells_4arm"

# Germ-layer palette and cell-type -> germ-layer map. Order (Epithelium ->
# Ectoderm -> Mesoderm -> Endoderm -> iPSC) matches the UMAP legend, and
# this order is used to sort cell types along the row axis of the Fig 3B
# matrix.
germ_layer_colors <- c(
  Epithelium = "#4292C6",
  Ectoderm   = "#F16913",
  Mesoderm   = "#44AA99",
  Endoderm   = "#9970AB",
  iPSC       = "#CC6677"
)
cell_type_germ_layer <- c(
  "Goblet cells"                              = "Epithelium",
  "Ciliated epithelial cells"                 = "Epithelium",
  "Intestinal epithelial cells"               = "Epithelium",
  "Squamous epithelial cells"                 = "Epithelium",
  "Bronchiolar and alveolar epithelial cells" = "Epithelium",
  "ENS glia"                                  = "Ectoderm",
  "Glioblast"                                 = "Ectoderm",
  "Neuroblast"                                = "Ectoderm",
  "Neuron"                                    = "Ectoderm",
  "Neuronal IPC"                              = "Ectoderm",
  "Radial glia"                               = "Ectoderm",
  "Retinal neurons"                           = "Ectoderm",
  "Retinal pigment cells"                     = "Ectoderm",
  "Retinal progenitors and Muller glia"       = "Ectoderm",
  "Early mesoderm cells"                      = "Mesoderm",
  "Cardiac mesoderm cells"                    = "Mesoderm",
  "Cardiomyocytes"                            = "Mesoderm",
  "Epicardial fat cells"                      = "Mesoderm",
  "Smooth muscle cells"                       = "Mesoderm",
  "Vascular endothelial cells"                = "Mesoderm",
  "Stromal cells"                             = "Mesoderm",
  "Metanephric cells"                         = "Mesoderm",
  "Skeletal muscle cells"                     = "Mesoderm",
  "Lymphoid cells"                            = "Mesoderm",
  "Erythroblasts"                             = "Mesoderm",
  "Megakaryocytes"                            = "Mesoderm",
  "Myeloid cells"                             = "Mesoderm",
  "Hepatoblasts"                              = "Endoderm",
  "Stellate cells"                            = "Endoderm",
  "iPSCs"                                     = "iPSC"
)

# Helpers.

# figures_dir(): create-and-return the figures/ subdirectory inside the
# cormotif output tree. Every panel below writes here.
figures_dir <- function() {
  d <- file.path(output_base_dir, "figures")
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

# save_figure(): write a ggplot to both PDF and PNG under figures_dir().
# Standardises the DPI and file-naming across panels.
save_figure <- function(p, stem, w = 8, h = 6) {
  figs <- figures_dir()
  ggsave(file.path(figs, paste0(stem, ".pdf")), p, width = w, height = h)
  ggsave(file.path(figs, paste0(stem, ".png")), p, width = w, height = h, dpi = 300)
  cat("Saved:", stem, "\n")
  invisible(p)
}

# load_pattern_table(): read motif_pattern_table.csv and normalise the schema.
# Older versions of the table lacked a motif_type column; if it's missing we
# derive it from drug_pattern + active_comparisons so downstream code can
# always rely on motif_type being present.
load_pattern_table <- function() {
  f <- file.path(output_base_dir, "motif_pattern_table.csv")
  pt <- read.csv(f, stringsAsFactors = FALSE)

  if (!"motif_type" %in% colnames(pt)) {
    pt$n_active_comparisons <- ifelse(
      pt$is_null | pt$active_comparisons == "", 0L,
      lengths(regmatches(pt$active_comparisons,
                         gregexpr(";", pt$active_comparisons))) + 1L
    )
    pt$motif_type <- dplyr::case_when(
      pt$is_null                               ~ "null",
      pt$n_active_comparisons == 10L           ~ "pan-comparison",
      grepl("-specific$", pt$drug_pattern)     ~ pt$drug_pattern,
      TRUE                                     ~ "mixed"
    )
  }

  # Treat any NA motif_type as null so downstream filters behave.
  pt$motif_type[is.na(pt$motif_type)] <- "null"
  pt
}

# motif_type_color(): given a motif_type string, return the color to render
# it with. Drug-specific motifs get their drug's color; the remaining types
# come from motif_type_colors.
motif_type_color <- function(type) {
  drug_specific <- paste0(drugs, "-specific")
  drug_map      <- setNames(drug_colors, drug_specific)
  unname(c(motif_type_colors, drug_map)[type])
}


# Fig 3B: drug x cell-type occurrence matrix (main figure).
# Binary cell = 1 iff the drug has at least one drug-specific motif in that
# cell type, colored in the drug's palette color. Row order follows the
# germ-layer sequence so the neural clustering of Lor-specific motifs stays
# visually together, but germ layer is not encoded on the row labels
# themselves (the legend on Fig 3A already establishes the palette).
#
# include_cerlor_column (default FALSE) appends an extra "Cer+Lor vs others"
# column, filled with the group-specific palette color (#5F3C8F) for cell
# types that carry a group-specific motif whose active comparisons are the
# six between-group pairs for {Cer, Lor} vs {Ale, Bri, Cri}. Saves to a
# "_alt" filename so the base figure is unaffected.

fig3B_drug_celltype_matrix <- function(include_cerlor_column = FALSE) {
  pt <- load_pattern_table()

  ds      <- pt[grepl("-specific$", pt$motif_type), ]
  ds$drug <- sub("-specific$", "", ds$drug_pattern)

  # One row per (cell_type, drug); fill in missing combinations with 0 so
  # every cell in the matrix is defined.
  mat <- ds %>%
    dplyr::select(cell_type, drug) %>%
    dplyr::distinct() %>%
    dplyr::mutate(present = 1L) %>%
    tidyr::complete(cell_type, drug = drugs, fill = list(present = 0L))

  # Column-total sanity print (Bri should be empty in this fit).
  col_totals <- tapply(mat$present, mat$drug, sum)
  cat("=== Fig 3B drug-specific column totals ===\n"); print(col_totals)
  cat("  total drug-specific motifs =", sum(col_totals), "\n")
  if (col_totals["Bri"] != 0)
    warning("Expected Bri column empty; got: ", col_totals["Bri"])
  cat("\n")

  # Optionally add the "Cer+Lor vs others" group-specific overlay as an extra
  # column. One row per cell type: present = 1 if a group-specific motif in
  # that cell type matches the Cer+Lor-vs-Ale+Bri+Cri partition.
  cl_label <- "Cer+Lor vs others"
  col_order <- drugs
  if (include_cerlor_column) {
    cl_cts <- unique(pt$cell_type[pt$motif_type == "group-specific" &
                                  pt$drug_pattern == "Cer+Lor vs Ale+Bri+Cri"])
    cl_rows <- data.frame(
      cell_type = sort(unique(mat$cell_type)),
      drug      = cl_label,
      stringsAsFactors = FALSE
    )
    cl_rows$present <- as.integer(cl_rows$cell_type %in% cl_cts)
    mat <- dplyr::bind_rows(mat, cl_rows)
    col_order <- c(drugs, cl_label)
    cat("=== Fig 3B alt: Cer+Lor-vs-others column ===\n")
    cat("  cells filled: ", sum(cl_rows$present), " / ", nrow(cl_rows), "\n\n",
        sep = "")
  }

  # Row order: germ-layer sequence top to bottom. ggplot's discrete y-axis
  # renders factor levels bottom to top, so we rev() here.
  ct_top_to_bottom <- names(cell_type_germ_layer)[
    names(cell_type_germ_layer) %in% unique(mat$cell_type)]
  mat$cell_type <- factor(mat$cell_type, levels = rev(ct_top_to_bottom))
  mat$drug      <- factor(mat$drug, levels = col_order)

  # fill_drug is NA for absent cells so scale_fill_manual can render them
  # with na.value = "white".
  mat$fill_drug <- factor(ifelse(mat$present == 1L, as.character(mat$drug),
                                  NA_character_),
                          levels = col_order)

  fill_palette <- c(drug_colors, setNames("#5F3C8F", cl_label))
  xaxis_colors <- unname(fill_palette[col_order])

  p <- ggplot(mat, aes(x = drug, y = cell_type)) +
    geom_tile(aes(fill = fill_drug), color = "grey80", linewidth = 0.3) +
    scale_fill_manual(values = fill_palette, na.value = "white",
                      breaks = col_order, name = "Drug",
                      drop = FALSE, guide = "none") +
    labs(x = "Drug", y = "Cell type") +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_text(hjust = 0.5,
                                 color = xaxis_colors, face = "bold"),
      axis.text.y = element_text(size = 8, color = "black")
    )

  stem <- if (include_cerlor_column) "fig3B_drug_celltype_matrix_alt"
          else "fig3B_drug_celltype_matrix"
  w <- if (include_cerlor_column) 6.3 else 5.5
  save_figure(p, stem, w = w, h = 8)
}


# Fig 3C: GO enrichment for the three highlighted motif sets. Three
# sub-panels align with three biological claims in the manuscript text:
#   3C-i   Alectinib   -> splicing / spliceosome (12 Ale-specific motifs),
#                        plus an SR-factor membership strip
#   3C-ii  Ceritinib   -> ribosome biogenesis (3 Cer-specific motifs)
#   3C-iii Lorlatinib  -> glycolysis, concentrated in neural (ectoderm) types
#                        (all 14 Lor-specific motifs shown, incl. absences)
# Also emits fig3C_values.csv with the exact rows plotted, and prints
# manuscript-vs-data cross-checks without silently reconciling them.

# Load the pre-filtered GO:BP significant table.
.fig3C_load_go <- function() {
  read.csv(file.path(output_base_dir, "S3_GO_BP_significant.csv"),
           stringsAsFactors = FALSE)
}

# For a given drug and pathway keyword regex, take the top (lowest
# elimFisher) GO term per (cell_type, motif).
.fig3C_top_term <- function(go_df, drug, keyword_regex) {
  go_df %>%
    dplyr::filter(motif_type == paste0(drug, "-specific"),
                  grepl(keyword_regex, Term, ignore.case = TRUE),
                  !is.na(elimFisher)) %>%
    dplyr::group_by(cell_type, motif) %>%
    dplyr::slice_min(order_by = elimFisher, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
}

# All (cell_type, motif) pairs for a drug's drug-specific motifs. Used to
# render bars at 0 for cell types that have a motif but no matching-keyword
# GO hit (so absences show up as visible zero bars, not omitted rows).
.fig3C_all_motifs <- function(pt, drug) {
  pt[pt$motif_type == paste0(drug, "-specific"),
     c("cell_type", "motif")]
}

# Turn a vector of cell-type names into a factor ordered by the germ-layer
# sequence (reversed for ggplot's bottom-to-top axis).
.fig3C_ordered_ct <- function(cts) {
  ct_all <- names(cell_type_germ_layer)
  factor(cts, levels = rev(ct_all[ct_all %in% unique(as.character(cts))]))
}

# Build a single Fig 3C sub-panel. highlight_layer optionally recolors rows
# in a specific germ layer (used for Lor -> Ectoderm) while other rows go
# grey; without it, all bars use the drug's own color.
.fig3C_subpanel <- function(drug, keyword_regex, subpanel_title,
                            highlight_layer = NULL, pt, go_df) {
  all_pairs <- .fig3C_all_motifs(pt, drug)
  top_terms <- .fig3C_top_term(go_df, drug, keyword_regex)

  plot_df <- all_pairs %>%
    dplyr::left_join(top_terms, by = c("cell_type", "motif")) %>%
    dplyr::mutate(neg_log10_p = ifelse(is.na(elimFisher), 0, -log10(elimFisher)),
                  germ_layer  = unname(cell_type_germ_layer[cell_type]))

  if (!is.null(highlight_layer)) {
    plot_df$bar_group <- ifelse(plot_df$germ_layer == highlight_layer,
                                "highlight", "other")
    bar_col_map <- c("highlight" = unname(germ_layer_colors[highlight_layer]),
                     "other"     = "grey70")
  } else {
    plot_df$bar_group <- "drug"
    bar_col_map <- c("drug" = unname(drug_colors[drug]))
  }

  plot_df$cell_type <- .fig3C_ordered_ct(plot_df$cell_type)
  plot_df <- plot_df %>% dplyr::arrange(cell_type)

  p <- ggplot(plot_df, aes(x = neg_log10_p, y = cell_type, fill = bar_group)) +
    geom_col(width = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    scale_fill_manual(values = bar_col_map, guide = "none") +
    labs(x = expression(-log[10](italic(p))~"(elim Fisher)"),
         y = NULL, title = subpanel_title) +
    theme_classic(base_size = 10) +
    theme(plot.title  = element_text(hjust = 0, face = "bold", size = 11),
          axis.text.y = element_text(size = 8))

  list(plot = p,
       data = plot_df %>%
         dplyr::transmute(
           drug        = drug,
           cell_type   = as.character(cell_type),
           germ_layer  = germ_layer,
           motif_id    = motif,
           GO.ID       = GO.ID,
           Term        = Term,
           elimFisher  = elimFisher,
           Significant = Significant,
           Annotated   = Annotated,
           neg_log10_p = neg_log10_p))
}

# SR-factor membership strip for Fig 3C-i: for each Ale-specific motif's
# gene list, check which of a curated set of SR splicing factors is a
# member and render as a presence/absence tile matrix.
.fig3C_sr_strip <- function(pt) {
  sr_factors <- c("SRSF11", "SRSF2", "U2AF2", "SRRM2", "RBM25", "SRSF4", "DDX46")
  ale <- pt[pt$motif_type == "Ale-specific", c("cell_type", "motif")]

  rows <- do.call(rbind, lapply(seq_len(nrow(ale)), function(i) {
    ct <- ale$cell_type[i]; mot <- ale$motif[i]
    f  <- file.path(output_base_dir, ct,
                    paste0(ct, "_motif", mot, "_genes_K2plus.txt"))
    if (!file.exists(f)) { warning("Ale gene list missing: ", f); return(NULL) }
    genes <- readLines(f)
    data.frame(cell_type = ct,
               sr_factor = sr_factors,
               present   = sr_factors %in% genes,
               stringsAsFactors = FALSE)
  }))

  rows$sr_factor <- factor(rows$sr_factor, levels = sr_factors)
  rows$cell_type <- .fig3C_ordered_ct(rows$cell_type)

  totals <- tapply(rows$present, rows$sr_factor, sum)
  cat("=== Fig 3C-i SR-factor recurrence across ", length(unique(rows$cell_type)),
      " Ale-specific motifs ===\n", sep = "")
  print(totals)
  cat("(manuscript claim: SRSF11=11, SRSF2=6, U2AF2=6, SRRM2=6, RBM25=5, SRSF4=4)\n\n")

  p <- ggplot(rows, aes(x = sr_factor, y = cell_type, fill = present)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_manual(values = c(`TRUE`  = unname(drug_colors["Ale"]),
                                 `FALSE` = "grey92"),
                      breaks = c(TRUE, FALSE),
                      labels = c("Present", "Absent"),
                      name   = "In Ale motif") +
    labs(x = NULL, y = NULL, title = "SR-factor membership") +
    theme_classic(base_size = 10) +
    theme(plot.title  = element_text(hjust = 0, face = "bold", size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1,
                                     face = "italic", size = 8),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  list(plot = p, data = rows)
}

# Print manuscript-quoted enrichment p-values next to the values actually
# in the data so any drift becomes visible in the log rather than being
# silently reconciled.
.fig3C_check_text_numbers <- function(cer_data, lor_data) {
  cat("=== Fig 3C: manuscript vs data cross-check ===\n")
  cat("Ceritinib ribosome-biogenesis top-term hits (expected: Neuronal IPC ~1e-4, Retinal pigment ~3e-4, Radial glia ~8e-6):\n")
  cer_hits <- cer_data[!is.na(cer_data$elimFisher), ]
  if (nrow(cer_hits) == 0) {
    cat("  (none)\n")
  } else {
    for (i in seq_len(nrow(cer_hits)))
      cat("  ", format(cer_hits$cell_type[i], width = 30),
          "  p = ", format(cer_hits$elimFisher[i], digits = 2, scientific = TRUE),
          "  |  ", cer_hits$Term[i], "\n", sep = "")
  }
  cat("\nLorlatinib glycolysis top-term hits (expected: ENS glia ~2e-4, Neuronal IPC ~3e-3, Radial glia ~4e-3, Neuroblast ~0.03):\n")
  lor_hits <- lor_data[!is.na(lor_data$elimFisher), ]
  if (nrow(lor_hits) == 0) {
    cat("  (none)\n")
  } else {
    for (i in seq_len(nrow(lor_hits)))
      cat("  ", format(lor_hits$cell_type[i], width = 30),
          "  p = ", format(lor_hits$elimFisher[i], digits = 2, scientific = TRUE),
          "  |  ", lor_hits$Term[i], "\n", sep = "")
  }
  cat("\n")
}

fig3C_motif_go_enrichment <- function() {
  pt    <- load_pattern_table()
  go_df <- .fig3C_load_go()

  # Keyword regexes that pick out the target pathway terms from GO:BP names.
  splic_re <- "splic|spliceosom|mRNA processing"
  ribo_re  <- "ribosom|rRNA|snoRNP|Cajal body|ribonucleoprotein"
  gly_re   <- "glycol|glucose-6-phosphate|glycolytic|glucose metabolic|glucose catabolic"

  ale_p <- .fig3C_subpanel("Ale", splic_re,
                           "Alectinib - splicing / spliceosome",
                           highlight_layer = NULL, pt = pt, go_df = go_df)
  cer_p <- .fig3C_subpanel("Cer", ribo_re,
                           "Ceritinib - ribosome biogenesis",
                           highlight_layer = NULL, pt = pt, go_df = go_df)
  lor_p <- .fig3C_subpanel("Lor", gly_re,
                           "Lorlatinib - glycolysis / glucose metabolism",
                           highlight_layer = "Ectoderm", pt = pt, go_df = go_df)
  sr_p  <- .fig3C_sr_strip(pt)

  # Dump the exact plotted rows to a CSV so reviewers can audit them.
  csv_out <- dplyr::bind_rows(ale_p$data, cer_p$data, lor_p$data)
  write.csv(csv_out, file.path(figures_dir(), "fig3C_values.csv"), row.names = FALSE)
  cat("Saved: fig3C_values.csv  (", nrow(csv_out), "rows)\n\n", sep = "")

  .fig3C_check_text_numbers(cer_p$data, lor_p$data)

  # Assemble via patchwork if present; otherwise emit sub-panels individually.
  if (requireNamespace("patchwork", quietly = TRUE)) {
    top_row  <- patchwork::wrap_plots(ale_p$plot, sr_p$plot,
                                       widths = c(0.55, 0.45), nrow = 1)
    bottom   <- patchwork::wrap_plots(cer_p$plot, lor_p$plot,
                                       widths = c(0.35, 0.65), nrow = 1)
    combined <- patchwork::wrap_plots(top_row, bottom, ncol = 1,
                                       heights = c(0.5, 0.6))
    save_figure(combined, "fig3C_motif_go_enrichment", w = 12, h = 10)
  } else {
    message("patchwork not available - saving sub-panels separately")
    save_figure(ale_p$plot, "fig3C_i_ale_splicing",     w = 6, h = 5)
    save_figure(sr_p$plot,  "fig3C_i_sr_membership",    w = 5, h = 5)
    save_figure(cer_p$plot, "fig3C_ii_cer_ribosome",    w = 5, h = 3)
    save_figure(lor_p$plot, "fig3C_iii_lor_glycolysis", w = 6, h = 6)
  }
}


# Fig 3D: Lor-specific motif overlap upset (main-figure sizing). Same data
# source as the supplementary drug_specific_upset_Lor.png (Supp Fig S8),
# retheme'd for main-figure use and emitted as a single-page vector PDF.

fig3D_lor_motif_upset <- function() {
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    message("Fig 3D skipped: UpSetR not installed."); return(invisible(NULL))
  }
  library(UpSetR)

  pt   <- load_pattern_table()
  rows <- pt[pt$motif_type == "Lor-specific", ]
  if (nrow(rows) < 2) {
    message("Fig 3D: fewer than 2 Lor-specific motifs; skipping.")
    return(invisible(NULL))
  }

  # Read per-motif gene lists, keep only non-empty ones, and shorten cell-type
  # names so the set labels fit inside the plot.
  gene_lists <- lapply(seq_len(nrow(rows)), function(i) {
    ct  <- rows$cell_type[i]; mot <- rows$motif[i]
    f   <- file.path(output_base_dir, ct,
                     paste0(ct, "_motif", mot, "_genes_K2plus.txt"))
    if (!file.exists(f)) { warning("Missing gene list: ", f); return(character(0)) }
    readLines(f)
  })
  names(gene_lists) <- rows$cell_type
  gene_lists <- gene_lists[lengths(gene_lists) > 0]
  names(gene_lists) <- gsub(" cells$| cell$", "", names(gene_lists))

  # Binary membership matrix for UpSetR.
  all_genes <- unique(unlist(gene_lists))
  bin_mat   <- as.data.frame(
    lapply(gene_lists, function(g) as.integer(all_genes %in% g)),
    check.names = FALSE)
  ct_names <- colnames(bin_mat)
  n_sets   <- length(ct_names)
  nintersects <- min(30, 2^n_sets - 1)

  # Cross-check manuscript-cited counts: total motif sizes and the
  # Hepatoblasts x Cardiac mesoderm pairwise overlap. These are set totals,
  # NOT exclusive-intersection bars in the upset (the exclusive bars would
  # be smaller when 3+ way overlaps exist).
  hep_col <- ct_names[ct_names == "Hepatoblasts"]
  cm_col  <- ct_names[ct_names == "Cardiac mesoderm"]
  if (length(hep_col) == 1 && length(cm_col) == 1) {
    hep_total  <- sum(bin_mat[[hep_col]])
    cm_total   <- sum(bin_mat[[cm_col]])
    hep_cm_int <- sum(bin_mat[[hep_col]] == 1 & bin_mat[[cm_col]] == 1)
    cat("=== Fig 3D manuscript cross-checks (set totals + pairwise overlap) ===\n")
    cat("  |Hepatoblasts|      = ", hep_total,  "   (expected 671)\n", sep = "")
    cat("  |Cardiac mesoderm|  = ", cm_total,   "   (expected 2,225)\n", sep = "")
    cat("  |Hep n CM|          = ", hep_cm_int, "   (expected 173)\n\n", sep = "")
  }

  out_pdf <- file.path(figures_dir(), "fig3D_lor_motif_upset.pdf")
  out_png <- file.path(figures_dir(), "fig3D_lor_motif_upset.png")

  render <- function() {
    print(upset(
      bin_mat,
      sets            = ct_names,
      nsets           = n_sets,
      nintersects     = nintersects,
      order.by        = "freq",
      mb.ratio        = c(0.6, 0.4),
      sets.bar.color  = unname(drug_colors["Lor"]),
      point.size      = 3,
      line.size       = 1,
      text.scale      = c(1.4, 1.3, 1.2, 1.1, 1.3, 1.1),
      mainbar.y.label = "Lor-specific genes in intersection",
      sets.x.label    = paste0("Genes in motif (", length(all_genes), " in union)")
    ))
  }

  png(out_png, width = 12, height = 6, units = "in", res = 300)
  render(); dev.off(); cat("Saved:", out_png, "\n")

  # cairo_pdf avoids the spurious blank first page pdf() produces with UpSetR.
  cairo_pdf(out_pdf, width = 12, height = 6)
  render(); dev.off(); cat("Saved:", out_pdf, "\n\n")
}


# Runner: called when the file is Rscript'ed directly (not via
# fig3B_3C_refresh.R, which greps for the "===" banner below and evals
# only the definitions above it).

cat("========================================\n")
cat("CorMotif Fig 3 panels\n")
cat("Output base:", output_base_dir, "\n")
cat("Figures dir:", figures_dir(), "\n")
cat("========================================\n\n")

fig3B_drug_celltype_matrix()
fig3B_drug_celltype_matrix(include_cerlor_column = TRUE)
fig3C_motif_go_enrichment()
fig3D_lor_motif_upset()

cat("========================================\n")
cat("Fig 3B / 3C / 3D complete. Check figures/ for output.\n")
cat("========================================\n")

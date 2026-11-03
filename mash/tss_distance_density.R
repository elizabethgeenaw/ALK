# Supp Fig S12 (3 germ layers): /project/gilad/geena/ALK_draft/V4/Supplementary Figures/S12_tss_distance_density_{ectoderm,endoderm_epithelium_other,mesoderm}.png


#!/usr/bin/env Rscript
#
# Multi-page version of the TSS-distance density plot (supp fig S12).
# Uses build_tss_distance_df() from mash_helpers.R, then splits the
# 21 canonical cell types across 3 pages by germ layer.
#
# Run on gilad-hm (Midway3), R/4.4.1.
#   Rscript plot_tss_distance_multipage.R [lfc]
# lfc defaults to 3 (paper's response-eQTL threshold).

library(tidyverse)

script_dir <- tryCatch(
  dirname(normalizePath(sys.frames()[[1]]$ofile)),
  error = function(e) "/project/gilad/geena/publication_clean/code/mash"
)
source(file.path(script_dir, "mash_helpers.R"))


args           <- commandArgs(trailingOnly = TRUE)
lfc            <- if (length(args) >= 1) as.numeric(args[[1]]) else 3
results_dir    <- "/project/gilad/geena/publication_clean/results/mash/medianUMI1"
gene_loc_file  <- "/project2/gilad/wenhe/EB/combined_3_4/data/ref_set/concise_GRCh38-2020_gtf_strand_considered.csv"
ncol_facets    <- 5
maf            <- 0.05
cisdist        <- "100kb"


germ_groups <- list(
  Ectoderm = c(
    "Radial_glia", "Neuronal_IPC", "Neuroblast", "Neuron", "ENS_glia",
    "Glioblast", "Retinal_neurons", "Retinal_pigment_cells",
    "Retinal_progenitors_and_Muller_glia"
  ),
  Mesoderm = c(
    "Cardiac_mesoderm_cells", "Early_mesoderm_cells",
    "Vascular_endothelial_cells", "Epicardial_fat_cells",
    "Stromal_cells", "Metanephric_cells"
  ),
  `Endoderm/Epithelium/Other` = c(
    "Hepatoblasts", "Bronchiolar_and_alveolar_epithelial_cells",
    "Goblet_cells", "Ciliated_epithelial_cells",
    "Squamous_epithelial_cells", "iPSCs"
  )
)

# Build shared dataframe
message("Building TSS distance df for lfc=", lfc, "...")
lfc_celltypes <- read_tsv(
  file.path(results_dir, paste0("lfc", lfc), "treatment_responsive_eqtls_signfixed.tsv"),
  show_col_types = FALSE
) %>% pull(celltype) %>% unique()

df_tss <- build_tss_distance_df(
  mash_folder   = results_dir,
  lfc           = lfc,
  gene_loc_file = gene_loc_file,
  celltypes     = lfc_celltypes
)

build_page <- function(df_group, group_name) {
  n_facets  <- length(unique(df_group$facet_label))
  nrow_plot <- ceiling(n_facets / ncol_facets)

  ggplot(df_group, aes(x = distance / 1000, fill = source)) +
    geom_histogram(
      aes(fill = source, alpha = source, y = after_stat(density)),
      position = "identity", alpha = 0.3, bins = 50
    ) +
    geom_density(
      aes(color = source, y = after_stat(density)),
      linewidth = 0.8, bounds = c(-100, 100), alpha = 0
    ) +
    scale_fill_manual(values  = c("Response" = "#f53224", "Non-response" = "#619CFF")) +
    scale_color_manual(values = c("Response" = "#f53224", "Non-response" = "#619CFF")) +
    scale_alpha_manual(values = c("Response" = 0.3, "Non-response" = 0.3)) +
    facet_wrap(~facet_label, ncol = ncol_facets, nrow = nrow_plot, scales = "free_y") +
    theme_classic(base_size = 14) +
    theme(legend.position = "bottom", strip.text = element_text(size = 10)) +
    labs(
      x     = "Distance to TSS (Kb)",
      y     = "Density",
      title = paste0(group_name,
                     " — Distance to TSS: fold change ≥ ", lfc,
                     ", MAF ", maf, ", ", cisdist)
    )
}

# Write multi-page pdf
plot_dir <- file.path(results_dir, paste0("lfc", lfc), "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(plot_dir, "tss_distance_density_multipage.pdf")

n_treatments <- length(unique(df_tss$treatment))
max_rows     <- max(vapply(germ_groups, function(cts) {
  ceiling(length(cts) * n_treatments / ncol_facets)
}, numeric(1)))
pdf_width    <- max(12, ncol_facets * 4)
pdf_height   <- max(8,  max_rows * 3)

message("Writing ", out_path, " (", pdf_width, " x ", pdf_height, " in)...")
pdf(out_path, width = pdf_width, height = pdf_height)
for (grp in names(germ_groups)) {
  df_group <- df_tss %>% filter(celltype %in% germ_groups[[grp]])
  if (nrow(df_group) == 0) {
    message("  skipping ", grp, ": no cell types present in df")
    next
  }
  n_cts <- length(unique(df_group$celltype))
  message("  page: ", grp, " (", n_cts, " cell types)")
  print(build_page(df_group, grp))
}
dev.off()

# Also write one PNG per group
# Filename-safe slug: "Endoderm/Epithelium/Other" -> "endoderm_epithelium_other"
slug <- function(x) {
  x %>% tolower() %>% gsub("[^a-z0-9]+", "_", .) %>% gsub("^_|_$", "", .)
}

for (grp in names(germ_groups)) {
  df_group <- df_tss %>% filter(celltype %in% germ_groups[[grp]])
  if (nrow(df_group) == 0) next
  n_cts     <- length(unique(df_group$celltype))
  nrow_grp  <- ceiling(n_cts * n_treatments / ncol_facets)
  png_height <- max(6, nrow_grp * 3)
  png_path  <- file.path(plot_dir,
                         paste0("tss_distance_density_", slug(grp), ".png"))
  message("  png:  ", png_path, " (", pdf_width, " x ", png_height, " in)")
  ggsave(png_path, build_page(df_group, grp),
         width = pdf_width, height = png_height, dpi = 300, limitsize = FALSE)
}

message("Done.")

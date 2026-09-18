# Main Fig 2D: 2D_eye_pathways_tile.pdf

# Eye-development pathway enrichment tile strip: 5 drugs x 30 cell types x
# 2 pathways (lens development, neural retina development). Each cell shows
# the mash-shrunken fgsea NES, with a star if padj < 0.05 in that
# (drug, cell type). Stacked vertically so the two pathways read as one figure.

# The point of this figure is to show the crizotinib eye phenotype: negative
# NES (down-regulation of eye-development pathways) is broad and strong in
# Cri but essentially absent or reversed in the other four ALK inhibitors,
# matching the clinical observation that Cri is the only one with major
# reported visual adverse events.

# Input:
#   results/treatment/dreamlet_mash/gsea/dreamlet_mash_gsea.tsv

# Output:
#   results/treatment/dreamlet_mash/plots/eye_pathways_tile.{pdf,png}

.libPaths(c("/project/gilad/geena/Rlibs", .libPaths()))

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Paths
base_dir <- "/project/gilad/geena/publication_clean"
mash_dir <- file.path(base_dir, "results/treatment/dreamlet_mash")
gsea_tsv <- file.path(mash_dir, "gsea/dreamlet_mash_gsea.tsv")
out_dir <- file.path(mash_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(gsea_tsv))

# Constants
drug_order <- c("Ale", "Bri", "Cer", "Cri", "Lor")

# Germ-layer ordering used for column grouping/annotation in the tile strip.
germ_layer_order <- c("iPSC", "Ectoderm", "Mesoderm", "Endoderm", "Epithelium")
germ_layer_colors <- c(
  iPSC = "#CC6677",
  Ectoderm = "#F16913",
  Mesoderm = "#44AA99",
  Endoderm = "#9970AB",
  Epithelium = "#4292C6"
)
cell_type_germ_layer <- c(
  "Goblet cells" = "Epithelium",
  "Ciliated epithelial cells" = "Epithelium",
  "Intestinal epithelial cells" = "Epithelium",
  "Squamous epithelial cells" = "Epithelium",
  "Bronchiolar and alveolar epithelial cells" = "Epithelium",
  "ENS glia" = "Ectoderm", "Glioblast" = "Ectoderm",
  "Neuroblast" = "Ectoderm", "Neuron" = "Ectoderm",
  "Neuronal IPC" = "Ectoderm", "Radial glia" = "Ectoderm",
  "Retinal neurons" = "Ectoderm", "Retinal pigment cells" = "Ectoderm",
  "Retinal progenitors and Muller glia" = "Ectoderm",
  "Early mesoderm cells" = "Mesoderm", "Cardiac mesoderm cells" = "Mesoderm",
  "Cardiomyocytes" = "Mesoderm", "Epicardial fat cells" = "Mesoderm",
  "Smooth muscle cells" = "Mesoderm", "Vascular endothelial cells" = "Mesoderm",
  "Stromal cells" = "Mesoderm", "Metanephric cells" = "Mesoderm",
  "Skeletal muscle cells" = "Mesoderm", "Lymphoid cells" = "Mesoderm",
  "Erythroblasts" = "Mesoderm", "Megakaryocytes" = "Mesoderm",
  "Myeloid cells" = "Mesoderm",
  "Hepatoblasts" = "Endoderm", "Stellate cells" = "Endoderm",
  "iPSCs" = "iPSC"
)

fdr_thresh <- 0.05

# NES color scale is winsorised at +/- nes_winsor so extreme cells don't blow
# out the range; underlying NES values are not modified, just the display.
nes_winsor <- 4

# The two eye-development pathways shown in Fig 2D.
target_pathways <- c(
  "GOBP_LENS_DEVELOPMENT_IN_CAMERA_TYPE_EYE",
  "GOBP_NEURAL_RETINA_DEVELOPMENT"
)
pathway_labels <- c(
  GOBP_LENS_DEVELOPMENT_IN_CAMERA_TYPE_EYE = "Lens development",
  GOBP_NEURAL_RETINA_DEVELOPMENT = "Neural retina development"
)

ts <- function() format(Sys.time(), "[%H:%M:%S]")

# Load the mash fgsea output (one row per drug x cell type x pathway) and
# restrict to the two eye pathways in the five ALK-inhibitor drugs.
message(ts(), " Reading GSEA: ", gsea_tsv)
gsea_all <- read_tsv(gsea_tsv, show_col_types = FALSE, progress = FALSE)
stopifnot(all(c("pathway", "padj", "NES", "drug", "assay", "collection")
              %in% names(gsea_all)))

eye <- gsea_all %>%
  filter(collection == "GOBP",
         pathway %in% target_pathways,
         drug %in% drug_order) %>%
  mutate(drug = factor(drug, levels = drug_order),
         pathway = factor(pathway, levels = target_pathways))
message(ts(), " Eye-pathway GSEA rows: ", nrow(eye),
        " (drugs x cell types x 2 pathways)")

# Cell-type order: germ layer, then alphabetical within layer, restricted to
# cell types that actually appear in the GSEA table.
ct_meta <- tibble(cell_type = names(cell_type_germ_layer),
                  germ_layer = unname(cell_type_germ_layer)) %>%
  filter(cell_type %in% unique(eye$assay)) %>%
  mutate(germ_layer = factor(germ_layer, levels = germ_layer_order)) %>%
  arrange(germ_layer, cell_type)
ct_order <- ct_meta$cell_type
message(ts(), " Cell types in tile: ", length(ct_order))

# Build one ComplexHeatmap Heatmap object per pathway. We show both a heatmap
# legend and germ-layer annotation legend on the top pathway; hide them on
# the bottom so the two heatmaps stack cleanly.
build_tile <- function(pathway_id, pretty_label, show_column_names,
                       show_heatmap_legend, show_annotation_legend) {
  # Winsorise NES for display (colour) only.
  sub <- eye %>%
    filter(pathway == pathway_id) %>%
    mutate(nes_disp = pmax(pmin(NES, nes_winsor), -nes_winsor))

  # Drug x cell type NES matrix.
  mat <- sub %>%
    select(drug, assay, nes_disp) %>%
    pivot_wider(names_from = assay, values_from = nes_disp) %>%
    tibble::column_to_rownames("drug") %>%
    as.matrix()
  mat <- mat[drug_order, ct_order, drop = FALSE]

  # Boolean sig matrix at the same shape, turned into "*"/"" strings so we
  # can drop them as text in each cell.
  sig_mat <- sub %>%
    mutate(sig = padj < fdr_thresh) %>%
    select(drug, assay, sig) %>%
    pivot_wider(names_from = assay, values_from = sig, values_fill = FALSE) %>%
    tibble::column_to_rownames("drug") %>%
    as.matrix()
  sig_mat <- sig_mat[drug_order, ct_order, drop = FALSE]
  star_mat <- ifelse(sig_mat, "*", "")

  # Diverging NES colors (blue = negative, red = positive).
  col_fun <- colorRamp2(c(-nes_winsor, 0, nes_winsor),
                        c("#1E77B4", "white", "#D62728"))

  # Column top-annotation showing germ layer of each cell type.
  col_germ <- factor(cell_type_germ_layer[ct_order], levels = germ_layer_order)
  col_anno <- HeatmapAnnotation(
    `Germ layer` = col_germ,
    col = list(`Germ layer` = germ_layer_colors),
    show_annotation_name = FALSE,
    show_legend = show_annotation_legend,
    annotation_legend_param = list(
      `Germ layer` = list(title_gp = gpar(fontsize = 8),
                          labels_gp = gpar(fontsize = 7))),
    simple_anno_size = unit(3, "mm")
  )

  Heatmap(
    mat,
    name = "NES",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_split = col_germ,
    column_gap = unit(1.2, "mm"),
    show_row_names = TRUE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 8),
    show_column_names = show_column_names,
    column_names_gp = gpar(fontsize = 6.5),
    column_names_rot = 45,
    top_annotation = col_anno,
    row_title = pretty_label,
    row_title_side = "left",
    row_title_rot = 90,
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    na_col = "grey92",
    border = TRUE,
    show_heatmap_legend = show_heatmap_legend,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8),
      labels_gp = gpar(fontsize = 7),
      direction = "vertical",
      legend_height = unit(3, "cm"),
      at = c(-nes_winsor, 0, nes_winsor),
      labels = c(paste0("<= -", nes_winsor), "0", paste0(">= ", nes_winsor))
    ),
    # Draw the significance stars on top of each cell.
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (nchar(star_mat[i, j]) > 0)
        grid.text(star_mat[i, j], x, y,
                  gp = gpar(fontsize = 8, fontface = "bold", col = "grey15"))
    }
  )
}

# Build the two tile heatmaps (top: lens dev with legends; bottom: retina
# with only column labels, no legends) and stack them vertically.
ht_lens <- build_tile("GOBP_LENS_DEVELOPMENT_IN_CAMERA_TYPE_EYE",
                      pathway_labels[["GOBP_LENS_DEVELOPMENT_IN_CAMERA_TYPE_EYE"]],
                      show_column_names = FALSE,
                      show_heatmap_legend = TRUE,
                      show_annotation_legend = TRUE)
ht_retina <- build_tile("GOBP_NEURAL_RETINA_DEVELOPMENT",
                        pathway_labels[["GOBP_NEURAL_RETINA_DEVELOPMENT"]],
                        show_column_names = TRUE,
                        show_heatmap_legend = FALSE,
                        show_annotation_legend = FALSE)

# Dynamic width so the plot scales with number of cell types.
n_ct <- length(ct_order)
plot_w <- max(11, 0.28 * n_ct + 3)
plot_h <- 6.5

pdf_path <- file.path(out_dir, "eye_pathways_tile.pdf")
png_path <- file.path(out_dir, "eye_pathways_tile.png")

pdf(pdf_path, width = plot_w, height = plot_h)
draw(ht_lens %v% ht_retina,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE,
     column_title = "Eye-development pathway enrichment (drug x cell type)",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

png(png_path, width = plot_w, height = plot_h, units = "in", res = 300)
draw(ht_lens %v% ht_retina,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE,
     column_title = "Eye-development pathway enrichment (drug x cell type)",
     column_title_gp = gpar(fontsize = 12, fontface = "bold"))
dev.off()

message(ts(), " Wrote ", pdf_path)
message(ts(), " Wrote ", png_path)
message(ts(), " fig2D_eye_pathways.R complete.")

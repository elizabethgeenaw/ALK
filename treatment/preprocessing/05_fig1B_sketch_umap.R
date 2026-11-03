# Main Fig 1B: /project/gilad/geena/ALK_draft/V4/Figures/1B_umap.pdf (sketch-based UMAP of 850K cells)

library(BPCells)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(future)

options(future.globals.maxSize = 200e+09)
# Sequential plan: avoids exporting the full BPCells-backed object to parallel
# workers (which blew past the per-plan maxSizeOfObjects cap in the previous run).
plan("sequential")

cellid_rds <- "/project/gilad/geena/publication_clean/results/treatment/celltype_annotation/postQC_lognorm_batch1_2_3_4_all_germlayeradded.rds"
sketch_dir <- "/project/gilad/geena/publication_clean/results/treatment/sketch_umap/"  # BPCells on-disk matrix and sketch RDS files
umap_plot_dir <- "/project/gilad/geena/ALK_paper_code_results/results/treatment/figures/umap_plots/" # UMAP plot output directory

# Cell type color palettes (matching cell type proportion plot)
ectoderm   <- c("Adrenocortical cells", "CNS glia", "CNS neurons", "Neuroendocrine cells", "PNS glia", "PNS neurons", "Photoreceptor cells", "Retinal pigment cells", "Horizontal cells", "Lens fibre cells", "Unipolar brush cells", "Chromaffin cells", "Retinal progenitors and Muller glia", "Ganglion cells", "Amacrine cells", "Bipolar cells", "ENS glia", "Purkinje neurons", "Schwann cells", "Inhibitory neurons", "Astrocytes", "Granule neurons", "Oligodendrocytes", "Sympathoblasts", "Visceral neurons", "Inhibitory interneurons", "Microglia", "Limbic system neurons", "ENS neurons", "Excitatory neurons", "Glioblast", "Neuroblast", "Neuron", "Neuronal IPC", "Oligo", "Radial glia", "Retinal neurons")
mesoderm   <- c("Early mesoderm cells", "Cardiac mesoderm cells", "Cardiac progenitor cells", "Cardiomyocytes", "Epicardial fat cells", "Endocardial cells", "Mesangial cells", "Smooth muscle cells", "Vascular endothelial cells", "Mesothelial cells", "Stromal cells", "Fibroblast", "Metanephric cells", "Ureteric bud cells", "Lymphatic endothelial cells", "Lymphoid cells", "Erythroblasts", "Megakaryocytes", "Hematopoietic stem cells", "Myeloid cells", "Antigen presenting cells", "Thymocytes", "Satellite cells", "Skeletal muscle cells")
endoderm   <- c("Hepatoblasts", "Stellate cells", "Parietal and chief cells", "Acinar cells", "Ductal cells", "Islet endocrine cells")
epithelium <- c("Goblet cells", "Neuroepithelial cells", "Ciliated epithelial cells", "Intestinal epithelial cells", "Squamous epithelial cells", "Bronchiolar and alveolar epithelial cells")

epithelium_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(length(epithelium))
ectoderm_colors   <- colorRampPalette(brewer.pal(9, "Oranges"))(length(ectoderm))
mesoderm_colors   <- colorRampPalette(c("#CCF0E8", "#88CCBB", "#44AA99", "#117777", "#004444"))(length(mesoderm))
endoderm_colors   <- colorRampPalette(brewer.pal(9, "Purples")[3:9])(length(endoderm))

cellid_colors <- c(
  setNames(epithelium_colors, epithelium),
  setNames(ectoderm_colors, ectoderm),
  setNames(mesoderm_colors, mesoderm),
  setNames(endoderm_colors, endoderm),
  setNames("#CC6677", "iPSCs"),
  setNames("grey50", "unassigned"),
  setNames("grey50", "Unassigned")
)

# Germ layer colors (single representative color per layer)
germ_colors <- c(
  "Epithelium" = "#4292C6",
  "Ectoderm"   = "#F16913",
  "Mesoderm"   = "#44AA99",
  "Endoderm"   = "#9970AB",
  "iPSCs"      = "#CC6677",
  "unassigned" = "grey50",
  "Unassigned" = "grey50"
)

# This script generates a sketch-based UMAP of the full dataset (too big to do otherwise)
# Sketching (LeverageScore, n=850,000) allows UMAP to run on a representative subsample, which is then projected back to the full dataset with ProjectData().

# Outputs (in umap_plot_dir):
#  - full_umap_*_dim1_50_res<res>.png (UMAP colored by cluster, cell type, batch, individual, germ layer)
#  - full_umap_*_dim1_50_res<res>_<marker>.png (marker gene FeaturePlots)


################## Functions ##################

# Write BPCells on-disk matrix, normalize, find HVGs, and sketch n=num_cells cells.
make_umap_sketch <- function(samples_object, num_cells) {
  ondisk_path <- file.path(sketch_dir, paste0("filtered_batch1_2_3_4_ncells", num_cells))
  message("Writing on-disk matrix to ", ondisk_path)
  write_matrix_dir(mat = samples_object[["RNA"]]$counts, dir = ondisk_path)
  counts.mat <- open_matrix_dir(dir = ondisk_path)

  ondisk_obj <- CreateSeuratObject(counts = counts.mat, meta.data = samples_object@meta.data)
  ondisk_obj <- NormalizeData(ondisk_obj)
  ondisk_obj <- FindVariableFeatures(ondisk_obj, selection.method = "vst", nfeatures = 3000)

  message("Sketching with n = ", num_cells, " cells (LeverageScore)")
  ondisk_obj <- SketchData(
    object = ondisk_obj,
    ncells = num_cells,
    method = "LeverageScore",
    sketched.assay = "sketch"
  )

  saveRDS(ondisk_obj,
          file.path(sketch_dir, paste0("filtered_batch1_2_3_4_ncells", num_cells, "_sketch.rds")))
  return(ondisk_obj)
}

# Run PCA on the sketched subset (only needs to be done once)
prep_sketch <- function(ondisk_allsamples) {
  DefaultAssay(ondisk_allsamples) <- "sketch"
  ondisk_allsamples <- FindVariableFeatures(ondisk_allsamples, selection.method = "vst", nfeatures = 3000)
  ondisk_allsamples <- ScaleData(ondisk_allsamples)
  ondisk_allsamples <- RunPCA(ondisk_allsamples)
  return(ondisk_allsamples)
}

# Cluster the sketched subset, run UMAP, then project to the full dataset
cluster_sketch <- function(ondisk_allsamples, num_cells, max_dim, res) {
  ondisk_allsamples <- FindNeighbors(ondisk_allsamples, dims = 1:max_dim)
  ondisk_allsamples <- FindClusters(ondisk_allsamples, resolution = res)
  ondisk_allsamples <- RunUMAP(ondisk_allsamples, dims = 1:max_dim, return.model = TRUE)

  # Project clusters and UMAP embedding to full dataset
  ondisk_allsamples <- ProjectData(
    object = ondisk_allsamples,
    assay = "RNA",
    full.reduction = "pca.full",
    sketched.assay = "sketch",
    sketched.reduction = "pca",
    umap.model = "umap",
    dims = 1:max_dim,
    refdata = list(cluster_full = "seurat_clusters")
  )
  return(ondisk_allsamples)
}

# Save UMAP plots (cluster, cell type, batch, individual, germ layer, marker genes).
make_plots <- function(ondisk_allsamples, num_cells, max_dim, res) {
  dir.create(umap_plot_dir, recursive = TRUE, showWarnings = FALSE)
  DefaultAssay(ondisk_allsamples) <- "RNA"
  tag <- paste0("sketchncells_", num_cells, "_dim1_", max_dim, "_res", res)





  # Order germlayer factor to match the cell type proportion plot (panel C) legend
  ondisk_allsamples$germlayer <- factor(
    ondisk_allsamples$germlayer,
    levels = c("Epithelium", "Ectoderm", "Mesoderm", "Endoderm", "iPSCs", "unassigned", "Unassigned")
  )
  germ_legend_order <- c("Epithelium", "Ectoderm", "Mesoderm", "Endoderm", "iPSCs", "unassigned")
  fullumap_germ <- DimPlot(ondisk_allsamples, label = FALSE,
                           reduction = "full.umap", group.by = "germlayer") +
    scale_color_manual(values = germ_colors, breaks = germ_legend_order, drop = TRUE) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    theme(axis.text   = element_blank(),
          axis.ticks  = element_blank(),
          axis.title  = element_text(size = 18),
          legend.text = element_text(size = 16),
          legend.title = element_blank(),
          plot.title  = element_text(size = 20))
  ggsave(file.path(umap_plot_dir, paste0("full_umap_germ_",   tag, ".png")), fullumap_germ, width = 11, height = 11)
  ggsave(file.path(umap_plot_dir, paste0("full_umap_germ_",   tag, ".pdf")), fullumap_germ, width = 11, height = 11)


  return(ondisk_allsamples)
}


################## Run ##################
# batch1_2_3_4 <- readRDS(cellid_rds)  # only needed when re-sketching

# Sketch 850,000 cells for UMAP (only run once, then use readRDS below)
# batch1_2_3_4.sketch <- make_umap_sketch(batch1_2_3_4, 850000)

# Run PCA once
#batch1_2_3_4.sketch <- prep_sketch(batch1_2_3_4.sketch)
#saveRDS(batch1_2_3_4.sketch, file.path(sketch_dir, "filtered_batch1_2_3_4_ncells850000_sketch.rds"))
batch1_2_3_4.sketch <- readRDS(file.path(sketch_dir, paste0("filtered_batch1_2_3_4_ncells850000_sketch.rds")))


# Elbow plot to determine optimal number of dimensions (saved to umap_plot_dir)
#ggsave(file.path(umap_plot_dir, "elbow_plot_sketch850000.png"),
#       ElbowPlot(batch1_2_3_4.sketch, ndims = 50),
#       width = 8, height = 5)

# Cluster and plot at multiple dimension ranges
# Compare to select the final visualization
res <- 0.5
dims <- 30

for (max_dim in dims) {
  message("Running dims 1:", max_dim, " res=", res)
  clustered <- cluster_sketch(batch1_2_3_4.sketch, 850000, max_dim, res)
  make_plots(clustered, 850000, max_dim, res)
}

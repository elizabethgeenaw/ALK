library(BPCells) 
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)
library(scales)
library(SingleCellExperiment)
options(future.globals.maxSize = 50e+09)

library(reticulate)
use_virtualenv("/project/gilad/geena/envs/convert_control", required = TRUE)
library(anndata)
library(CelliD)

# Convert control h5ad subset of controls to Seurat objects
# Define input and output directories
input_dir <- "/project/gilad/geena/ALK_analysis/QC_preprocessing_control"
output_dir <- "/project/gilad/geena/ALK_analysis/QC_preprocessing_control"

# Loop over the subsets
for (i in 1:8) {  # Adjust the range based on the number of splits
  # Construct file paths
  input_file <- file.path(input_dir, paste0("split_by_donor_", i, ".h5ad"))
  bp_dir <- file.path(output_dir, paste0("split_by_donor_", i, "_BP"))
  output_file <- file.path(output_dir, paste0("split_by_donor_", i, ".rds"))
  
  # Read the .h5ad matrix into BPCells format
  mat <- open_matrix_anndata_hdf5(input_file)
  
  # Save BPCells matrix on disk
  write_matrix_dir(
    mat = mat,
    dir = bp_dir
  )
  
  # Reload matrix from BPCells to ensure it's in BPCells format
  bp_mat <- open_matrix_dir(dir = bp_dir)
  
  # Read metadata
  metadata <- reticulate::py_to_r(
    reticulate::import("anndata")$read_h5ad(input_file)$obs
  )
  metadata <- as.data.frame(metadata)
  
  # Create Seurat object using BPCells matrix
  seurat_obj <- CreateSeuratObject(
    counts = bp_mat,
    meta.data = metadata
  )
  
  # Save the Seurat object as .rds
  saveRDS(seurat_obj, output_file)
  
  cat(paste("Processed and saved:", output_file, "\n"))
}

# Import RDS and merge
split1 <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/split_by_donor_1.rds")
split2 <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/split_by_donor_2.rds")
split3 <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/split_by_donor_3.rds")
split4 <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/split_by_donor_4.rds")
split5 <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/split_by_donor_5.rds")
split6 <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/split_by_donor_6.rds")
split7 <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/split_by_donor_7.rds")
split8 <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/split_by_donor_8.rds")

control <- merge(x = split1, y = list(split2, split3, split4, split5, split6, split7, split8)) #36,601 x 1,071,062
control <- JoinLayers(control)
#saveRDS(control, "/project/gilad/geena/ALK_analysis/QC_preprocessing_control/combined_control.rds")

# Preliminary QC
#https://www.cell.com/cell-genomics/fulltext/S2666-979X(24)00330-6 - ebQTL paper
#Filter to keep only those with less than 15% mitochondrial reads and with at least 2500 genes expressed
#Remove cells with very high total counts (keep only those with less than 150,000 total counts)
#Filter to genes expressed in at least 10 cells
# Read combined control data: 36,601 x 1,071,062
control <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/combined_control.rds")

control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^MT-")

# Subset based on filtering criteria -> 36,601 x 914,776
control_filt <- subset(control, subset = percent.mt < 15 & 
                         nFeature_RNA >= 2500 & 
                         nCount_RNA < 150000)
dim(control_filt)

# Subset again based on gene-level filtering criteria
# Get RNA counts matrix
counts_matrix <- GetAssayData(control_filt[["RNA"]], layer = "counts")

# Identify genes expressed in at least 10 cells
genes_to_keep <- rownames(counts_matrix)[rowSums(counts_matrix > 0) >= 10] #35160

# Subset the Seurat object to retain only these genes
control_filt <- subset(control_filt, features = genes_to_keep)
dim(control_filt) # 914,776 cells and 35,160 genes? Paper says 909,536 cells and 35,324 genes. Close enough?
#saveRDS(control_filt, "/project/gilad/geena/ALK_analysis/QC_preprocessing_control/combined_control_filtered.rds")

# Convert to BPCells
#control_filt <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/combined_control_filtered.rds")
# Write the counts layer to a directory
write_matrix_dir(mat = control_filt[["RNA"]]$counts, dir = "/project/gilad/geena/ALK_analysis/QC_preprocessing_all/ebQTL_control/combined_control_filtered")
counts.mat <- open_matrix_dir(dir = "/project/gilad/geena/ALK_analysis/QC_preprocessing_all/ebQTL_control/combined_control_filtered")
control_filt[["RNA"]]$counts <- counts.mat
#saveRDS(control_filt, "/project/gilad/geena/ALK_analysis/QC_preprocessing_control/combined_control_filtered.rds")
control_filt <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/combined_control_filtered.rds")

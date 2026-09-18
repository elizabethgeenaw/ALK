library(BPCells)
library(Seurat)
library(SeuratObject)
library(dplyr)

options(future.globals.maxSize = 50e+09)

preqc_dir <- "/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/intermediate_files" # output of script 2
filtered_dir <- "/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/intermediate_files" # filtered RDS output (can be same as preqc_dir)
ondisk_fixed_dir <- "/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/ondisk_fixedRDS" # normalized/scaled BPCells output


# This script applies QC filters, then log-normalizes and scales data per batch × treatment for downstream CellID annotation and dreamlet DE.
# The fix_object_norm_scale() function re-writes the BPCells matrix to a new on-disk location using the integer type required for RunPCA (workaround for a BPCells PCA compatibility issue).

# Outputs:
# - filtered_<batch>_<condition>.rds (filtered only, for CellID)
# - postQC_lognorm_<batch>_<condition>.rds (filtered + log-normalized, for CellID)
# - ondisk_fixed_dir/<batch>_<condition>/ (BPCells on-disk matrix, normalized/scaled)
# - ondisk_fixed_dir/<batch>_<condition>.rds (Seurat object pointing to above matrix)


################## Functions ##################

# Apply QC filters and save filtered-only Seurat object.
# Filters: nFeature_RNA 1000–8000, percent_mito < 20%, remove Vireo doublets/unassigned.
# Batch 3 removes two low-quality individuals (NA18128, NA19138)
# Input: (1) combined: Seurat object with a sample with lanes from the same sample merged (i.e. batch1_Cri), (2) batch (i.e. "batch1")
# Output: filtered (but not normalized) RDS is saved under "/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/intermediate_files/filtered_<batch>_<condition>.rds"
subset_QC_combined <- function(combined, batch, condition,
                                nFeature_RNA_MIN = 1000,
                                nFeature_RNA_MAX = 8000,
                                percent_mito_MAX = 20) {
  base_filter <- quote(nFeature_RNA > nFeature_RNA_MIN &
                       nFeature_RNA < nFeature_RNA_MAX &
                       percent_mito < percent_mito_MAX &
                       individual != "doublet" &
                       individual != "unassigned")

  if (batch == "batch3") {
    combined_subset <- subset(combined,
      subset = nFeature_RNA > nFeature_RNA_MIN & nFeature_RNA < nFeature_RNA_MAX &
               percent_mito < percent_mito_MAX & individual != "doublet" &
               individual != "unassigned" & individual != "18128" & individual != "19138")
  } else {
    combined_subset <- subset(combined,
      subset = nFeature_RNA > nFeature_RNA_MIN & nFeature_RNA < nFeature_RNA_MAX &
               percent_mito < percent_mito_MAX & individual != "doublet" &
               individual != "unassigned")
  }

  saveRDS(combined_subset,
          file.path(filtered_dir, paste0("filtered_", batch, "_", condition, ".rds")))
  return(combined_subset)
}

# Log-normalize filtered data and save (used as input to CellID annotation).
normalize_combined <- function(combined_subset, batch, condition) {
  normalized <- NormalizeData(combined_subset)
  saveRDS(normalized,
          file.path(filtered_dir, paste0("postQC_lognorm_", batch, "_", condition, ".rds")))
  return(normalized)
}

# Re-write BPCells matrix with correct integer type, then normalize, find variable features, scale, and run PCA. Saves normalized/scaled Seurat RDS to ondisk_fixed_dir.
# This workaround is required because RunPCA fails on certain BPCells matrix types.
fix_object_norm_scale <- function(rds_path, sample) {
  message("Processing: ", sample)
  broken <- readRDS(rds_path)
  broken.counts <- broken[["RNA"]]$counts
  broken.counts <- as(broken.counts, Class = "dgCMatrix")
  broken.counts <- convert_matrix_type(broken.counts, "uint32_t")

  out_dir <- file.path(ondisk_fixed_dir, sample)
  write_matrix_dir(mat = broken.counts, dir = out_dir, overwrite = TRUE)
  fixed.counts <- open_matrix_dir(dir = out_dir)

  obj <- CreateSeuratObject(counts = fixed.counts, meta.data = broken@meta.data)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)

  saveRDS(obj, file.path(ondisk_fixed_dir, paste0(sample, ".rds")))
  return(obj)
}


################## Run ##################
batches <- c("batch1", "batch2", "batch3", "batch4")
conditions <- c("Cri", "Ale", "Cer", "Bri", "Lor")

# Step 1: Load preQC objects, apply QC filters, and save filtered-only RDS
filtered_objects <- list()
for (b in batches) {
  for (cond in conditions) {
    key <- paste0(b, "_", cond)
    message("Filtering: ", key)
    obj <- readRDS(file.path(preqc_dir, paste0("preQC_", key, ".rds")))
    filtered_objects[[key]] <- subset_QC_combined(obj, b, cond)
  }
}

# Step 2: Log-normalize filtered objects (for CellID annotation in script 4)
for (key in names(filtered_objects)) {
  parts <- strsplit(key, "_")[[1]]
  b <- paste(parts[1], parts[2], sep = "_") # e.g., "batch1"
  cond <- parts[3] # e.g., "Cri"

  # Reconstruct batch/condition from key
  b <- sub("_[^_]+$", "", key) # everything before last underscore
  cond <- sub(".*_", "", key) # everything after last underscore
  normalize_combined(filtered_objects[[key]], b, cond)
}

# Step 3: Fix BPCells matrix type, normalize, scale, run PCA per sample
for (key in names(filtered_objects)) {
  b <- sub("_[^_]+$", "", key)
  cond <- sub(".*_", "", key)
  fix_object_norm_scale(
    rds_path = file.path(filtered_dir, paste0("filtered_", key, ".rds")),
    sample   = key
  )
}

# Step 4: Merge all filtered-only objects into a single combined object (used for sketch UMAP in script 5 and dreamlet DE in script 6a)
all_filtered <- lapply(names(filtered_objects), function(key) filtered_objects[[key]])
combined_all <- merge(all_filtered[[1]], y = all_filtered[-1])
combined_all <- JoinLayers(combined_all)
saveRDS(combined_all,
        file.path(filtered_dir, "filtered_batch1_2_3_4.rds"))

normalize_combined(combined_all, "batch1_2_3_4", "all")
fix_object_norm_scale(
  rds_path = file.path(filtered_dir, "filtered_batch1_2_3_4.rds"),
  sample   = "batch1_2_3_4_all"
)

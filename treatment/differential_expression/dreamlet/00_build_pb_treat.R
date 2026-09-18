# Upstream of 01_dreamlet_1v4_DE.R AND both cormotif upstream scripts
# (cormotif/01_setdiff_check.R, cormotif/02_prep_raw_counts.R).
#
# Builds the dreamlet pseudobulk SingleCellExperiment object `pb_treat` from
# the post-QC, log-normalized, cell-type-annotated Seurat object produced by
# treatment/preprocessing/04_celltype_annotation.R.
#
# pb_treat structure:
#   - one assay per cell type (assayNames(pb_treat) = cell type names)
#   - each assay is a genes x samples raw count matrix
#   - sample names are "<drug><individual>" (e.g. "AleNA18486")
#   - cellCounts(pb_treat) returns samples x cell_types cell counts
#
# This aggregation only needs to be run once. Downstream dreamlet and
# cormotif scripts load the cached RDS directly. Re-run if the upstream
# annotation or QC filter changes.
#
# Run on gilad-hm with R/4.4.1 via sbatch 00_build_pb_treat.sh (needs a
# lot of memory — ~200GB — because the input Seurat object is ~850k cells).

library(BPCells)
library(Seurat)
library(SeuratObject)
library(dreamlet)
library(SingleCellExperiment)
options(future.globals.maxSize = 50e+09)

drug_order  <- c("Ale", "Bri", "Cer", "Cri", "Lor")
results_dir <- "/project/gilad/geena/publication_clean/results/treatment/dreamlet_clean"
input_rds   <- "/project/gilad/geena/publication_clean/results/treatment/celltype_annotation/postQC_lognorm_batch1_2_3_4_all_germlayeradded.rds"

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Load the post-QC / log-norm / cell-type-annotated Seurat object.
message("[", format(Sys.time(), "%H:%M:%S"), "] Loading Seurat object...")
seurat_obj <- readRDS(input_rds)

# BPCells-backed count matrices are stored as integer types; dreamlet's
# aggregation expects double. Convert in place before as.SingleCellExperiment.
seurat_obj[['RNA']]$counts <- convert_matrix_type(
  seurat_obj[['RNA']]$counts, type = "double"
)

# Convert to SCE and drop genes with zero total expression across all cells.
sce_treat <- as.SingleCellExperiment(seurat_obj)
sce_treat <- sce_treat[rowSums(counts(sce_treat) > 0) > 0, ]
message("Cells after removing unassigned individuals: ", ncol(sce_treat))

# Set the grouping variables aggregateToPseudoBulk needs:
#   cluster_id = cell type (one assay per cell type in the output)
#   sample_id  = drug + individual (one column per sample per assay)
sce_treat$StimStatus <- factor(sce_treat$treatment, levels = drug_order)
sce_treat$id         <- paste0(sce_treat$StimStatus, sce_treat$individual)

# Aggregate to pseudobulk. Uses raw counts (assay = "counts") so downstream
# steps can apply their own normalization (TMM inside dreamlet, TMM+voom
# inside cormotif's limmafit path).
message("[", format(Sys.time(), "%H:%M:%S"), "] Aggregating to pseudobulk...")
pb_treat <- aggregateToPseudoBulk(
  sce_treat,
  assay      = "counts",
  cluster_id = "cellid_type",
  sample_id  = "id",
  verbose    = TRUE
)

# Cache to disk. Downstream scripts load this RDS directly.
out_path <- file.path(results_dir, "pb_treat")
saveRDS(pb_treat, out_path)
message("[", format(Sys.time(), "%H:%M:%S"), "] Saved: ", out_path)
message("N assays (cell types): ", length(assayNames(pb_treat)))
message("N samples: ",             ncol(pb_treat))

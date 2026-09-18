library(BPCells) 
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(dplyr)
library(ggplot2)
library(harmony) 
library(sctransform)
library(patchwork)
library(scales)
library(qs)
library(biomaRt)
options(future.globals.maxSize = 50e+09)


######### WRITE FUNCTIONS ##########
# Function to add vireo ID to metadata for a Seurat object under column "individual"
## Example to run: add_vireo_id(obj, paste0("/folder/", batch, "/run_count_S", x, "/outs/vireo/donor_ids.tsv"))
add_vireo_id <- function(object, vireo_output) {
  vireo <- read.table(vireo_output, header=TRUE, stringsAsFactors = FALSE) # read in the file
  m <- match(rownames(object@meta.data),vireo$cell)
  vireo <- vireo[m,]
  vireo_id <- vireo$donor_id
  names(vireo_id) <- vireo$cell
  AddMetaData(object, vireo_id, col.name ='individual')
}

# Merge samples from the same batch with the same condition (ex: 2 lanes of batch1_Ale)
# The object is saved under "/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/intermediate_files/preQC_<batch>_<condition>.rds")
# Input: (1) batch (ex: "batch1"), (2) sample names: vector of sample names (each of which must be the name of a directory in ../RDS/<batch>/), (3) treatment/condition name (i.e. Cri)
merge_samples <- function(batch, sample_names, condition){
  eb.obj <- lapply(sample_names, function(x){
    print(x)
    obj <- readRDS(paste0("/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/ondisk_RDS/", batch, "/", x, ".rds")) # RDS from convert_disk_storage
    add_vireo_id(obj, paste0("/project/gilad/geena/ALK_analysis/vireo/", batch, "/", x, "/donor_ids.tsv"))
  })
  #In Seurat v5, merging creates a single object, but keeps the expression information split into different layers for integration. 
  #If not proceeding with integration, rejoin the layers after merging.
  combined <- merge(x = eb.obj[[1]], y = eb.obj[-1], add.cell.ids = paste0(batch, "_", condition, "_", sample_names), project = condition)
  combined <- JoinLayers(combined)
  
  # Add percent_mito, ncell_indiv, batch, treatment to metadata
  combined[["percent_mito"]] <- PercentageFeatureSet(combined, pattern = "^MT-")
  combined[["ncell_indiv"]] <- add_count(combined@meta.data, individual)$n
  combined[["batch"]] <- batch
  combined[["treatment"]] <- condition
  
  # Print QC metrics summaries
  print("Summary of nFeature_RNA:")
  print(summary(combined$nFeature_RNA))
  print("Summary of nCount_RNA:")
  print(summary(combined$nCount_RNA))
  print("Summary of percent_mito:")
  print(summary(combined$percent_mito))
  
  saveRDS(combined, paste0("/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/intermediate_files/preQC_", batch, "_", condition, ".rds"))
  return(combined)}


######### MERGE INDIVIDUAL SAMPLES, 2 LANES PER OBJECT & SAVE PRE-QC OBJECT #########
# Merge two lanes of a sample within a batch (ex: in batch1 YG-KR-10X-ALK10-ALK1-01-1-Cri and 02-2-Cri just become batch1_Cri)
# Save pre-QC object to "/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/intermediate_files/preQC_<batch>_<condition>.rds" 

# Batch1
batch1_Cri <- merge_samples("batch1", c("YG-KR-10X-ALK10-ALK1-01-1-Cri", "YG-KR-10X-ALK10-ALK1-02-2-Cri"), "Cri")
batch1_Ale <- merge_samples("batch1", c("YG-KR-10X-ALK10-ALK1-03-3-Ale", "YG-KR-10X-ALK10-ALK1-04-4-Ale"), "Ale")
batch1_Cer <- merge_samples("batch1", c("YG-KR-10X-ALK10-ALK1-05-5-Cer", "YG-KR-10X-ALK10-ALK1-06-6-Cer"), "Cer")
batch1_Bri <- merge_samples("batch1", c("YG-KR-10X-ALK10-ALK1-07-7-Bri", "YG-KR-10X-ALK10-ALK1-08-8-Bri" ), "Bri")
batch1_Lor <- merge_samples("batch1", c("YG-KR-10X-ALK10-ALK1-09-9-Lor", "YG-KR-10X-ALK10-ALK1-10-10-Lor"), "Lor")

# Batch2 
batch2_Cri <- merge_samples("batch2", c("YG-KR-10X-10s-scR-ALK2-01-Cri", "YG-KR-10X-10s-scR-ALK2-02-Cri"), "Cri")
batch2_Ale <- merge_samples("batch2", c("YG-KR-10X-10s-scR-ALK2-03-Ale", "YG-KR-10X-10s-scR-ALK2-04-Ale"), "Ale")
batch2_Cer <- merge_samples("batch2", c("YG-KR-10X-10s-scR-ALK2-05-Cer", "YG-KR-10X-10s-scR-ALK2-06-Cer"), "Cer")
batch2_Bri <- merge_samples("batch2", c("YG-KR-10X-10s-scR-ALK2-07-Bri", "YG-KR-10X-10s-scR-ALK2-08-Bri"), "Bri")
batch2_Lor <- merge_samples("batch2", c("YG-KR-10X-10s-scR-ALK2-09-Lor", "YG-KR-10X-10s-scR-ALK2-10-Lor"), "Lor")

# Batch3
batch3_Cri <- merge_samples("batch3", c("batch3_01", "YG-KR-10s-ALK3-02"), "Cri")
batch3_Ale <- merge_samples("batch3", c("batch3_03", "batch3_04"), "Ale")
batch3_Cer <- merge_samples("batch3", c("batch3_05", "batch3_06"), "Cer")
batch3_Bri <- merge_samples("batch3", c("batch3_07", "batch3_08"), "Bri")
batch3_Lor <- merge_samples("batch3", c("batch3_09", "batch3_10"), "Lor")

# Batch4
batch4_Cri <- merge_samples("batch4", c("batch4_01", "batch4_02"), "Cri")
batch4_Ale <- merge_samples("batch4", c("batch4_03", "batch4_04"), "Ale")
batch4_Cer <- merge_samples("batch4", c("batch4_05", "batch4_06"), "Cer")
batch4_Bri <- merge_samples("batch4", c("batch4_07", "batch4_08"), "Bri")
batch4_Lor <- merge_samples("batch4", c("batch4_09", "batch4_10"), "Lor")

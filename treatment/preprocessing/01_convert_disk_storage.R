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

# This file is run first to convert output from CellRanger to on-disk matrices using BPCells (necessary since these are large datasets)
# Only need to run this file once (to convert to on-disk)

################## Define functions ###################
# Function to convert ensembl to symbol using gtf (more reliable than BioMart)
# Input: matrix
# Output: matrix with converted ensembl to symbol gene names
ConvertEnsembleToSymbol_gtf <- function(mat, gene_gtf){
  name_df <- data.frame(gene_id = c(rownames(mat)))
  name_df$orig.id <- name_df$gene_id
  name_df$gene_id <- as.character(name_df$gene_id) # Make this a character, otherwise it will throw errors with left_join
  name_df$gene_id <- sub("[.][0-9]*","",name_df$gene_id) # Remove Ensembl suffix
  gene.df <- left_join(name_df, gene_gtf, by = c("gene_id"="ID"))
  
  rownames(gene.df) <- make.unique(gene.df$orig.id)
  gene.df <- gene.df[rownames(mat),]
  gene.df <- gene.df[gene.df[,"gene"] != '',]
  gene.df <- gene.df[ !is.na(gene.df$orig.id),]
  mat.filter <- mat[gene.df$orig.id,]
  rownames(mat.filter) <- make.unique(gene.df[,"gene"])
  return(mat.filter)
}

# Function to convert one sample to on-disk storage
# Input: matrix from CellRanger output
# Output: Seurat object under "/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/ondisk_RDS/<batch>/<sample>.rds"
## Convert ensembl to gene symbol prior to writing on-disk
## BPcells seems to run into problems when renaming dimensions on matrix, so give cells and genes unique names that don't change
## Save this output in "/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/on_disk_RDS/"
convert_disk_storage <- function(batch, sample, gene_gtf) {
  obj.mat <- Read10X(data.dir = paste0("/project/gilad/geena/ALK_analysis/alignment/", batch, "/", sample, "/outs/filtered_feature_bc_matrix"))
  
  # Convert ensembl IDs to gene symbols
  obj.mat <- ConvertEnsembleToSymbol_gtf(mat = obj.mat, gene_gtf)
  
  # Write matrix to a directory
  write_matrix_dir(mat = obj.mat, dir = paste0("/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/ondisk_RDS/", batch, "/", sample, "/outs/counts"))
  
  # Now that we have the matrix on disk, we can load it
  obj.mat <- open_matrix_dir(dir = paste0("/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/ondisk_RDS/", batch, "/", sample, "/outs/counts"))
  obj <- CreateSeuratObject(counts = obj.mat)
  saveRDS(object = obj,  file = paste0("/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/ondisk_RDS/", batch, "/", sample, ".rds")) # save this to separate directory
}


####### Assign sample names based on CellRanger folder and convert to on-disk storage ########
# Import gtf genes from Wenhe's folder (colnames: chr, start, end, ID, gene). NOTE: these are uncorrected for strand but just using them to map Ensembl gene
gene_gtf <- read.csv("/project2/gilad/wenhe/EB/combined_3_4/data/ref_set/concise_GRCh38-2020_gtf.csv") #36,601 genes

batch1_samples <- c("YG-KR-10X-ALK10-ALK1-01-1-Cri",
                    "YG-KR-10X-ALK10-ALK1-02-2-Cri",
                    "YG-KR-10X-ALK10-ALK1-03-3-Ale", 
                    "YG-KR-10X-ALK10-ALK1-04-4-Ale",
                    "YG-KR-10X-ALK10-ALK1-05-5-Cer",
                    "YG-KR-10X-ALK10-ALK1-06-6-Cer",
                    "YG-KR-10X-ALK10-ALK1-07-7-Bri",
                    "YG-KR-10X-ALK10-ALK1-08-8-Bri",
                    "YG-KR-10X-ALK10-ALK1-09-9-Lor",
                    "YG-KR-10X-ALK10-ALK1-10-10-Lor")

print("Batch1")
for (s in batch1_samples) {
  print(s)
  convert_disk_storage("batch1", s, gene_gtf)
}

print("------")

batch2_samples <- c("YG-KR-10X-10s-scR-ALK2-01-Cri", 
                    "YG-KR-10X-10s-scR-ALK2-02-Cri",
                    "YG-KR-10X-10s-scR-ALK2-03-Ale",
                    "YG-KR-10X-10s-scR-ALK2-04-Ale",
                    "YG-KR-10X-10s-scR-ALK2-05-Cer",
                    "YG-KR-10X-10s-scR-ALK2-06-Cer",
                    "YG-KR-10X-10s-scR-ALK2-07-Bri",
                    "YG-KR-10X-10s-scR-ALK2-08-Bri",
                    "YG-KR-10X-10s-scR-ALK2-09-Lor",
                    "YG-KR-10X-10s-scR-ALK2-10-Lor")

print("Batch2")
for (s in batch2_samples) {
  print(s)
  convert_disk_storage("batch2", s, gene_gtf)
}
print("------") 

batch3_samples <- c("batch3_01", #Cri 
                    "YG-KR-10s-ALK3-02", #Cri
                    "batch3_03", #Ale
                    "batch3_04", #Ale
                    "batch3_05", #Cer
                    "batch3_06", #Cer
                    "batch3_07", #Bri
                    "batch3_08", #Bri
                    "batch3_09", #Lor
                    "batch3_10") # Lor

print("Batch3")
for (s in batch3_samples) {
  print(s)
  convert_disk_storage("batch3", s, gene_gtf)
}
print("------")

batch4_samples <- c("batch4_01", #Cri 
                    "batch4_02", #Cri
                    "batch4_03", #Ale
                    "batch4_04", #Ale
                    "batch4_05", #Cer
                    "batch4_06", #Cer
                    "batch4_07", #Bri
                    "batch4_08", #Bri
                    "batch4_09", #Lor
                    "batch4_10") # Lor

print("Batch4")
for (s in batch4_samples) {
  print(s)
  convert_disk_storage("batch4", s, gene_gtf)
}
print("------") 


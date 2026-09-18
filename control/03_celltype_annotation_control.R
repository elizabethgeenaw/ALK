# Split by batch and log-normalize (make input for CellID)
#For my treatment data, I log-normalized each separate batch_treatment combination and processed each of those batches in CellID separately before merging (because CellID crashes when I try to do it for the whole dataset at once)
#Want to split batch by *collection date* not library prep batch or sequencing batch
#Collection date is when we dissociated; the only one that produces an effect; we don't want HVG to be dependent on another factor; probably won't make a huge difference?
# Import filtered object
control_filt <- readRDS("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/combined_control_filtered.rds")

# Get unique collection dates and reorder them chronologically
unique_dates <- unique(control_filt$Collection.Date)
unique_dates_ordered <- sort(as.Date(unique_dates, format = "%m/%d/%Y"))  # Sort dates in chronological order

# Create batch labels based on chronological order
batch_labels <- paste0("batch", seq_along(unique_dates_ordered))

# Split by batch and lognormalize
output_dir <- "/project/gilad/geena/ALK_analysis/QC_preprocessing_control/collectiondatebatch_filt_lognorm/"

# Loop through each ordered batch
# In this case there are 20 batches
for (i in seq_along(unique_dates_ordered)) {
  # Subset the data for the current batch
  batch_date <- unique_dates_ordered[i]
  batch_date_formatted <- format(as.Date(batch_date, format = "%Y-%m-%d"), "%-m/%-d/%Y")
  control_filt_batch <- subset(control_filt, Collection.Date == batch_date_formatted)
  
  # Normalize the data
  control_filt_batch_norm <- NormalizeData(control_filt_batch)
  
  # Define the output file path using batch label
  output_file <- paste0(output_dir, "postQC_lognorm_", batch_labels[i], ".rds")
  print(output_file)
  
  # Save the normalized data
  saveRDS(control_filt_batch_norm, file = output_file)
  cat(paste("Processed and saved normalized data for", batch_labels[i], "\n"))
}

# Print the mapping of batches to dates for reference
batch_date_mapping <- data.frame(
  Batch = batch_labels,
  Date = unique_dates_ordered
)
print(batch_date_mapping)

# Run CellID
# Get rds files in the filtered/normalized rds directory
postQC_lognorm_obj_file <- list.files("/project/gilad/geena/ALK_analysis/QC_preprocessing_control/collectiondatebatch_filt_lognorm", "postQC_lognorm_batch*", full.names = T)

# Load gene signature list generated from FCA
brain_fca <- c("Purkinje neurons", "Inhibitory neurons", "Astrocytes", "Granule neurons", "Oligodendrocytes", "Inhibitory interneurons", "Microglia", "Limbic system neurons", "Excitatory neurons", "Unipolar brush cells")
eye_fca <- c("Ganglion cells", "Amacrine cells", "Retinal progenitors and Muller glia", "Photoreceptor cells", "Bipolar cells", "Horizontal cells", "Retinal pigment cells", "	
Lens fibre cells") 
fca_signatures <- readRDS("/project2/gilad/wenhe/EB/combined_3_4/data/celltype_annotate/train_test/Cao_FCA_1Kcell_train_sce_250_signature_common_20iterate_cut500reads_noCTRiboMT.rds")
fca_signatures[["iPSCs"]] <- c("SOX2", "POU5F1", "NANOG", "NODAL", "ZFP42", "DNMT3B", "UTF1", "LCK", "DPPA2", "IDO1", "DPPA5", "LEFTY1", "ESRG", "CNMD", "SFRP2", "THY1", "USP44", "VSNL1", "SPP1", "KLF4", "GDF3", "L1TD1", "VRTN", "ZSCAN10")

# Add brain atlas (Linnarsson dataset) 
# Linnarson brain atlas has only very few cell types but they have many clusters for each cell type -- Wenhe keeps them in the reference and assign cells to celltype-cluster. Will merge them after cellid assignment
linn_signatures <- readRDS("/project2/gilad/wenhe/EB/combined_3_4/data/celltype_annotate/human_brain/brain_signature_reiterated_1000_Celltype_TopLevelCluster_noRiboCCMT.rds")

# Add CM dataset from our lab
CM_signatures <- readRDS("/project2/gilad/wenhe/EB/combined_3_4/data/celltype_annotate/CM_Reem/CM_signature_reiterated_1000_noRiboCCMT.rds")
CM_signatures <- CM_signatures[2:4]
names(CM_signatures) <- c("Early mesoderm cells", "Cardiac mesoderm cells", "Cardiac progenitor cells")

fca_signatures <- c(c(fca_signatures[!names(fca_signatures)%in%brain_fca], linn_signatures), CM_signatures)
fca_signatures <- lapply(fca_signatures, function(vector) vector[!vector %in% unlist(cc.genes.updated.2019)& vector %in% HgProteinCodingGenes])

all_cellid_orig <- c()
all_cellid_logpval <- c()

for (file in postQC_lognorm_obj_file) {
  #load object
  print(file)
  obj <- readRDS(file) 
  
  # Restrict to protein coding genes
  obj <- obj[rownames(obj) %in% HgProteinCodingGenes,]
  obj_counts <- as(obj[["RNA"]]$counts, "dgCMatrix")
  gene_keep <- rownames(obj_counts)[rowSums(obj_counts>0) > 5]
  obj <- obj[rownames(obj) %in% gene_keep,]
  
  # Get HVG, use top HVG in the obj + genes from the marker list
  hvg <- VariableFeatures(FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000))
  hvg_file <- sub("postQC_lognorm_", "", paste0(sub('\\.rds$', '', file), "_hvg.rds"))
  saveRDS(hvg, hvg_file)
  #hvg <- readRDS(hvg_file)
  
  # Exclude MT genes and ribosomal genes
  hvg <- hvg[!grepl("^RPL|^RPS|^MT-", hvg)&!hvg%in%unlist(cc.genes.updated.2019)]
  
  obj_hvg <- c(unlist(fca_signatures), hvg) %>% unique() %>% intersect(rownames(obj))
  obj <- ScaleData(obj, features = obj_hvg)
  
  # Get an MCA embedding of the EB data
  print("Running MCA")
  
  # Need to include `slot = "RNA`to avoid "`assay` must be one of "RNA", not "data" error
  obj <- RunMCA(obj, nmcs=50, assay = "RNA", slot = "RNA", features=obj_hvg)
  
  # Run hypergeometric test for cell type classification
  # RunCellHGT outputs the -log10 corrected p-value for each cell and each ref cell type
  cellid_enrichments <- RunCellHGT(obj, pathways = fca_signatures, dims = 1:50, n.features=300)
  
  # Output all pvals
  cellid_logpval <- as.data.frame(t(as.matrix(cellid_enrichments)))
  all_cellid_logpval <- rbind(all_cellid_logpval, cellid_logpval)
  
  # Pull labels - highest enrichment w P<=0.01 in each cell
  cellid_labels <- rownames(cellid_enrichments)[apply(cellid_enrichments, 2, which.max)]
  cellid_labels <- ifelse(apply(cellid_enrichments, 2, max)>=2, yes = cellid_labels, "unassigned") 
  all_cellid_orig <- rbind(all_cellid_orig, as_tibble(cellid_labels, rownames="cell_barcode"))
}

# Change the colname of assignment as 'celltype'
all_cellid_orig <- dplyr::select(all_cellid_orig, c(cell_barcode, celltype=value))
unique(all_cellid_orig$celltype)

#write.csv(all_cellid_orig, "/project/gilad/geena/publication_clean/results/control/celltype_annotation/control_all_cellid_orig.csv")
all_cellid_orig <- read.csv("/project/gilad/geena/publication_clean/results/control/celltype_annotation/control_all_cellid_orig.csv") %>% select(-X)


# Merge same neuro cell types from different clusters -> 58 cell types
# Replace numbers and optionally preceded spaces at the end of the strings
all_cellid_orig$celltype <- sub("\\s*\\d+$", "", all_cellid_orig$celltype)

# Now you can check the unique cell types after merging
unique(all_cellid_orig$celltype)

# Merge same neuro cell types from different clusters, discard spurious cell types (too few cells), merge similar cell types -> 42 cell types
all_cellid_orig_merged <- all_cellid_orig %>% 
  mutate(celltype = replace(celltype, celltype == "Visceral neurons" , "ENS neurons")) %>% 
  mutate(celltype = replace(celltype, celltype == "Chromaffin cells" , "ENS neurons"))%>% 
  mutate(celltype = replace(celltype, celltype == "Sympathoblasts" , "ENS neurons"))%>% 
  mutate(celltype = replace(celltype, celltype == "Schwann cells" , "ENS glia")) %>% 
  mutate(celltype = replace(celltype, celltype == "Ureteric bud cells" , "Metanephric cells")) %>% 
  mutate(celltype = replace(celltype, celltype == "Bipolar cells" , "Photoreceptor cells")) %>% 
  mutate(celltype = replace(celltype, celltype == "Thymocytes" , "Lymphoid cells"))%>% 
  mutate(celltype = replace(celltype, celltype == "Mesothelial cells" , "Epicardial fat cells")) %>% 
  mutate(celltype = replace(celltype, celltype == "Amacrine cells" , "Horizontal cells")) %>% 
  mutate(celltype = replace(celltype, celltype == "Ganglion cells" , "Horizontal cells")) %>% 
  mutate(celltype = replace(celltype, celltype == "Lymphatic endothelial cells" , "Vascular endothelial cells")) %>% 
  mutate(celltype = replace(celltype, celltype == "Mesangial cells" , "Stromal cells"))%>% 
  mutate(celltype = replace(celltype, celltype == "Fibroblast" , "Stromal cells"))%>% 
  mutate(celltype = replace(celltype, celltype == "Parietal and chief cells" , "Goblet cells"))%>% 
  mutate(celltype = replace(celltype, celltype == "Cardiac progenitor cells" , "Cardiomyocytes"))%>% 
  mutate(celltype = replace(celltype, celltype == "Horizontal cells" , "Retinal neurons"))%>% 
  mutate(celltype = replace(celltype, celltype == "Photoreceptor cells" , "Retinal neurons"))


# Write labels to a csv file
#write.csv(all_cellid_orig_merged, "/project/gilad/geena/publication_clean/results/control/celltype_annotation/control_allbatches_cellid_label_updated.csv")
all_cellid_orig_merged <- read.csv("/project/gilad/geena/publication_clean/results/control/celltype_annotation/control_allbatches_cellid_label_updated.csv") %>% select(-X)

# write all res of pval to a csv file
#write.csv(all_cellid_logpval, "/project/gilad/geena/publication_clean/results/control/celltype_annotation/control_allbatches_cellid_label_updated_pval.csv", row.names = T)
all_cellid_logpval <- read.csv("/project/gilad/geena/publication_clean/results/control/celltype_annotation/control_allbatches_cellid_label_updated_pval.csv")

#dim(all_cellid_orig)

################ Summary counts ##################
# 42 cell types
# Proportions of cell counts
cell_type_counts <- all_cellid_orig_merged %>%
  group_by(celltype) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  mutate(percent = (n/sum(n))*100)

write.csv(cell_type_counts, "/project/gilad/geena/publication_clean/results/control/celltype_annotation/control_celltype_proportions.csv")

# Add cell type annotations to postQC, filtered, log-normalized object
# Add cellID annotations back to object and save
# Function to add cell ID to sample Seurat object for each specific sample
# cellannotations: cell_id_labels (output from Wenhe's pipeline) (dataframe)
# sample_object: Seurat object from given batch+treatment
add_annot_metadata <- function(cellannotations, sample_object) {
  # Match cell barcodes to the sample's Seurat object
  annotations_subset <- cellannotations[cellannotations$cell_barcode %in% colnames(sample_object), ]
  
  # Create a named vector for annotations
  annotations_vector <- setNames(annotations_subset$celltype, annotations_subset$cell_barcode)
  
  # Add annotations to the Seurat object
  sample_object[["cellid_type"]] <- annotations_vector[colnames(sample_object)]
  
  # Check for missing barcodes (optional for debugging)
  missing_barcodes <- setdiff(colnames(sample_object), annotations_subset$cell_barcode)
  if (length(missing_barcodes) > 0) {
    cat("Missing barcodes in this batch:", length(missing_barcodes), "\n")
  }
  
  return(sample_object)
}

########## Annotate each filtered, log-normalized collection date batch and then combine #########
# Can log-normalize within collection date batch (before combining); order of normalization matters less for log-normalize than SCTransform
input_dir <- "/project/gilad/geena/ALK_analysis/QC_preprocessing_control/collectiondatebatch_filt_lognorm/"

# Annotate and store each batch as an object
for (i in 1:20) {
  # Define file name
  input_file <- file.path(input_dir, paste0("postQC_lognorm_batch", i, ".rds"))
  
  # Load the Seurat object
  batch_obj <- readRDS(input_file)
  
  # Annotate the Seurat object with cellID metadata
  batch_obj_annotated <- add_annot_metadata(all_cellid_orig_merged, batch_obj)
  
  # Dynamically assign the annotated object to a variable
  assign(paste0("batch", i, "_filtnorm_cellid"), batch_obj_annotated)
  
  # Print progress
  cat(paste("Annotated and stored batch", i, "as batch", i, "_filtnorm_cellid\n"))
}

# Create a list of all annotated batch objects
batch_objects <- mget(paste0("batch", 1:20, "_filtnorm_cellid"))

# Merge all batches into one Seurat object
combined_control_filtnorm_cellid <- merge(batch1_filtnorm_cellid, y=list(batch2_filtnorm_cellid, batch3_filtnorm_cellid, batch4_filtnorm_cellid, batch5_filtnorm_cellid,
                                                                         batch6_filtnorm_cellid, batch7_filtnorm_cellid, batch8_filtnorm_cellid, batch9_filtnorm_cellid, batch10_filtnorm_cellid,
                                                                         batch11_filtnorm_cellid, batch12_filtnorm_cellid, batch13_filtnorm_cellid, batch14_filtnorm_cellid, batch15_filtnorm_cellid,
                                                                         batch16_filtnorm_cellid, batch17_filtnorm_cellid, batch18_filtnorm_cellid, batch19_filtnorm_cellid, batch20_filtnorm_cellid))

# Join layers to ensure all data layers are properly combined
combined_control_filtnorm_cellid <- JoinLayers(combined_control_filtnorm_cellid)

# Add germ layer
# Extract metadata
meta_data <- combined_control_filtnorm_cellid@meta.data

meta_data <- meta_data %>%
  mutate(cellid_type = ifelse(is.na(cellid_type), "unassigned", cellid_type))

# Try grouping by germ layer
ectoderm <- c("Adrenocortical cells", "CNS glia", "CNS neurons", "Neuroendocrine cells", "PNS glia", "PNS neurons", "Photoreceptor cells", "Retinal pigment cells",  "Horizontal cells", "Lens fibre cells", "Unipolar brush cells", "Chromaffin cells","Retinal progenitors and Muller glia", "Ganglion cells", "Amacrine cells", "Bipolar cells", "ENS glia", "Purkinje neurons", "Schwann cells", "Inhibitory neurons", "Astrocytes", "Granule neurons", "Oligodendrocytes", "Sympathoblasts", "Visceral neurons", "Inhibitory interneurons", "Microglia", "Limbic system neurons", "ENS neurons", "Excitatory neurons") %>% c("Glioblast", "Neuroblast", "Neuron", "Neuronal IPC", "Oligo", "Radial glia", "Retinal neurons")
mesoderm <- c("Early mesoderm cells", "Cardiac mesoderm cells", "Cardiac progenitor cells", "Cardiomyocytes", "Epicardial fat cells", "Endocardial cells", "Mesangial cells", "Smooth muscle cells", "Vascular endothelial cells", "Mesothelial cells", "Stromal cells", "Fibroblast", "Metanephric cells", "Ureteric bud cells", "Lymphatic endothelial cells", "Lymphoid cells", "Erythroblasts", "Megakaryocytes", "Hematopoietic stem cells", "Myeloid cells", "Antigen presenting cells", "Thymocytes", "Satellite cells","Skeletal muscle cells")
endoderm <- c("Hepatoblasts", "Stellate cells", "Parietal and chief cells", "Acinar cells", "Ductal cells",  "Islet endocrine cells")
epithelium <- c("Goblet cells", "Neuroepithelial cells", "Ciliated epithelial cells", "Intestinal epithelial cells", "Squamous epithelial cells", "Bronchiolar and alveolar epithelial cells")

confirmed_celltypes <- c("iPSCs", "Squamous epithelial cells", "Bronchiolar and alveolar epithelial cells", "Ciliated epithelial cells", "Intestinal epithelial cells", "Hepatoblasts",
                         "Retinal progenitors and Muller glia", "Retinal pigment cells", "Retinal neurons", "Neuroblast","Neuron", "ENS glia", 
                         "Cardiac mesoderm cells", "Cardiomyocytes", "Endocardial cells", "Vascular endothelial cells", "Epicardial fat cells", "Stromal cells", "Erythroblasts", "Megakaryocytes", "Myeloid cells")

meta_data <- meta_data %>% mutate(germlayer=case_when(cellid_type%in% ectoderm ~ "Ectoderm",
                                                      cellid_type%in% endoderm ~ "Endoderm",
                                                      cellid_type%in% mesoderm ~ "Mesoderm",
                                                      cellid_type%in% epithelium ~ "Epithelium",
                                                      .default =cellid_type)) %>%
  mutate(confirmed_celltype=if_else(cellid_type %in% confirmed_celltypes, TRUE, FALSE))

combined_control_filtnorm_cellid@meta.data <- meta_data

# Save RDS
saveRDS(combined_control_filtnorm_cellid, "/project/gilad/geena/publication_clean/results/control/celltype_annotation/combined_control_filtnorm_cellid.rds")

# Print summary of the merged object
cat("Merged control object dimensions:", dim(combined_control_filtnorm_cellid), "\n") # 35,160 x 914,776 (for cell+gene-level filtering)
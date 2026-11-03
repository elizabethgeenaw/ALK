# Main Fig 1C: /project/gilad/geena/ALK_draft/V4/Figures/1C_celltype_proportion_by_treatment.pdf 
# Supp Fig S1: /project/gilad/geena/ALK_draft/V4/Supplementary Figures/S1_celltype_proportion_by_individual.pdf

library(BPCells)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(CelliD)
options(future.globals.maxSize = 50e+09)

source("/project/gilad/geena/ALK_paper_code/shared/figure_defs.R")

filtered_dir <- "/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/"  # output of script 3b (filtered_*.rds and postQC_lognorm_*.rds)
cellid_out_dir <- "/project/gilad/geena/publication_clean/results/treatment/celltype_annotation/" # CellID output directory


# This script runs CellID for automated cell type annotation of all batch × treatment objects.
# It runs MCA embedding + hypergeometric test against three combined reference atlases.
# Cell type labels are merged/renamed, then added back to the Seurat objects.

# Outputs (in cellid_out_dir):
# - all_cellid_orig.csv (raw CellID labels per cell)
# - batch1_2_3_4_cellid_label_updated.csv (merged/renamed labels)
# - batch1_2_3_4_cellid_logpval_updated.csv (-log10 p-values for all reference cell types)
# - batch1_2_3_4_filtered_cellid.rds (filtered + cellid_type + germ layer, for script 5)
# - batch1_2_3_4_filtnorm_cellid.rds (filtered + normalized + cellid_type + germ layer, for scripts 6a+)

# Updated 6/13/24 with Wenhe's newer annotations for unguided HDCs

####### RUN CELLID FOR CELL TYPE ANNOTATION #########
# Get rds files in the processed rds directory
postQC_lognorm_obj_file <- list.files(filtered_dir, "postQC_lognorm_batch*", full.names = T)
# Filter out the specific file to exclude (combined file)
postQC_lognorm_obj_file <- postQC_lognorm_obj_file[!grepl("postQC_lognorm_batch1_2_3_4_all\\.rds$", postQC_lognorm_obj_file)]

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

# Loop through each batch x treatment combination
for (file in postQC_lognorm_obj_file) {
  #load object
  obj <- readRDS(file) 
  
  # Restrict to protein coding genes
  obj <- obj[rownames(obj) %in% HgProteinCodingGenes,]
  obj_counts <- as(obj[["RNA"]]$counts, "dgCMatrix")
  gene_keep <- rownames(obj_counts)[rowSums(obj_counts>0 )>5]
  obj <- obj[rownames(obj) %in% gene_keep,]
  
  # Get HVG, use top HVG in the obj + genes from the marker list
  hvg <- VariableFeatures(FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000))
  hvg_file <- sub("postQC_lognorm_", "", paste0(sub('\\.rds$', '', file), "_hvg.rds"))
  #saveRDS(hvg, hvg_file)
  hvg <- readRDS(hvg_file) # if previously run, can read HVG here
  
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

write.csv(all_cellid_orig, paste0(cellid_out_dir, "all_cellid_orig.csv")) # raw cellID labels per cell
#all_cellid_orig <- read.csv(paste0(cellid_out_dir, "all_cellid_orig.csv"))

# Merge same neuro cell types from different clusters
# Replace numbers and optionally preceded spaces at the end of the strings
all_cellid_orig$celltype <- sub("\\s*\\d+$", "", all_cellid_orig$celltype)

# Now check the unique cell types after merging
unique(all_cellid_orig$celltype)

# Merge same neuro cell types from different clusters, discard spurious cell types (too few cells), merge similar cell types
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


# Write merged/renamed labels to a csv file

# Write all res of pval to a csv file

#dim(all_cellid_orig)

# Proportions of cell counts
cell_type_counts <- all_cellid_orig_merged %>%
  group_by(celltype) %>%
  summarise(n=n()) %>%
  arrange(desc(n)) %>%
  mutate(percent = (n/sum(n))*100)


# Add cellID annotations back to object and save
# Function to add cell ID to sample Seurat object for each specific sample
# cellannotations: cell_id_labels (output from Wenhe's pipeline) (dataframe)
# sample_object: Seurat object from given batch+treatment
# Example to run: add_annot_metadata(all_cellid_orig, batch4_Cri_filtered, "batch4_Cri")
add_annot_metadata <- function(cellannotations, sample_object, sample_name) {
  print(sample_name)
  # Stash cell identity classes in metadata
  sample_object[["old.ident"]] <- Idents(object = sample_object) 
  
  # Get subset of labels for specific sample to merge
  labels_subset <- cellannotations %>% filter(grepl(paste0("^", sample_name), cell_barcode))
  print(head(labels_subset))
  print(dim(labels_subset))
  
  # Check if all your cell barcodes from the annotations are in your Seurat object
  print(all(labels_subset$cell_barcode %in% colnames(sample_object)))
  
  # Create vector with annotations
  labels_subset_vector <- setNames(labels_subset$celltype, labels_subset$cell_barcode)
  missing_cols <- !(colnames(sample_object) %in% names(labels_subset_vector))
  if (any(missing_cols)) {
    print("Missing cell barcodes in the sample object:")
    print(colnames(sample_object)[missing_cols])
  }
  
  # Add the cell type annotations to the Seurat object
  sample_object[['cellid_type']] <- labels_subset_vector[colnames(sample_object)]
  
  print("-----")
  return(sample_object)
}

all_cellid_orig <- read.csv(paste0(cellid_out_dir, "batch1_2_3_4_cellid_label_updated.csv")) %>% dplyr::select(-X)

######## Filtered-only #########
batches <- paste0("batch", 1:4)
drugs <- c("Cri", "Ale", "Bri", "Cer", "Lor")

filtered_list <- list()
for (batch in batches) {
  for (drug in drugs) {
    name <- paste0(batch, "_", drug)
    obj <- readRDS(paste0("/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/intermediate_files/filtered_", name, ".rds"))
    filtered_list[[name]] <- add_annot_metadata(all_cellid_orig, obj, name)
  }
}

batch1_2_3_4_filtered_cellid <- Reduce(function(x, y) merge(x, y), filtered_list)
batch1_2_3_4_filtered_cellid <- JoinLayers(batch1_2_3_4_filtered_cellid)
saveRDS(batch1_2_3_4_filtered_cellid, paste0(cellid_out_dir, "batch1_2_3_4_filtered.rds"))

######## Filtered + normalized ########
filtnorm_list <- list()
for (batch in batches) {
  for (drug in drugs) {
    name <- paste0(batch, "_", drug)
    obj <- readRDS(paste0("/project/gilad/geena/ALK_analysis/QC_preprocessing_treatment/postQC_lognorm_", name, ".rds"))
    filtnorm_list[[name]] <- add_annot_metadata(all_cellid_orig, obj, name)
  }
}

batch1_2_3_4_filtnorm_cellid <- Reduce(function(x, y) merge(x, y), filtnorm_list)
batch1_2_3_4_filtnorm_cellid <- JoinLayers(batch1_2_3_4_filtnorm_cellid)
#saveRDS(batch1_2_3_4_filtnorm_cellid, paste0(cellid_out_dir, "postQC_lognorm_batch1_2_3_4_all.rds"))


######## Add germ layer ##########
# Extract metadata
meta_data <- batch1_2_3_4_filtnorm_cellid@meta.data

meta_data <- meta_data %>%
  mutate(cellid_type = ifelse(is.na(cellid_type), "unassigned", cellid_type))

# Group by germ layer
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

batch1_2_3_4_filtnorm_cellid@meta.data <- meta_data
saveRDS(batch1_2_3_4_filtnorm_cellid, paste0(cellid_out_dir, "postQC_lognorm_batch1_2_3_4_all_germlayeradded.rds"))


######### Make summary tables #########
# Group by cell type
cell_type_counts <- meta_data %>%
  group_by(cellid_type) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n)) %>%
  arrange(desc(n))


# Cells per cell type summary (treatment; exclude unassigned)
cell_type_counts_assigned <- cell_type_counts %>% filter(cellid_type != "unassigned")
print(paste0("Treatment cells per cell type (n_celltypes = ", nrow(cell_type_counts_assigned),
             ") -- median: ", median(cell_type_counts_assigned$n),
             ", min: ",        min(cell_type_counts_assigned$n),
             ", max: ",        max(cell_type_counts_assigned$n),
             ", total cells: ", sum(cell_type_counts_assigned$n)))

# Group by cell type and individual
individual_counts <- meta_data %>%
  group_by(cellid_type, individual) %>%
  summarise(n = n()) %>%
  arrange(desc(n))

# Group by cell type and treatment and individual
ct_tx_indiv_counts <- meta_data %>%
  group_by(cellid_type, treatment, individual) %>%
  summarise(n = n()) %>%
  arrange(desc(n))


# Cell type x treatment total counts (aggregated across individuals)
ct_tx_counts <- meta_data %>%
  group_by(cellid_type, treatment) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(treatment, desc(n))


# Median number of cells per cell type per individual
median_cells_per_ct_per_indiv <- individual_counts %>%
  group_by(cellid_type) %>%
  summarise(
    median_n = median(n),
    mean_n = mean(n),
    min_n = min(n),
    max_n = max(n),
    n_individuals = n()
  ) %>%
  arrange(desc(median_n))


# Total cells per individual across all cell types
total_cells_per_individual <- meta_data %>%
  group_by(individual) %>%
  summarise(total_cells = n())


# Number of cell types detected per individual
ct_per_individual <- meta_data %>%
  filter(cellid_type != "unassigned") %>%
  group_by(individual) %>%
  summarise(n_celltypes = n_distinct(cellid_type)) %>%
  arrange(desc(n_celltypes))

median_ct_per_individual <- median(ct_per_individual$n_celltypes)
min_ct_per_individual    <- min(ct_per_individual$n_celltypes)
max_ct_per_individual    <- max(ct_per_individual$n_celltypes)
print(paste0("Treatment cell types per individual -- median: ", median_ct_per_individual,
             ", min: ", min_ct_per_individual,
             ", max: ", max_ct_per_individual,
             " (n_individuals = ", nrow(ct_per_individual), ")"))

write.csv(ct_per_individual, paste0(cellid_out_dir, "batch1_2_3_4_celltypes_per_individual.csv"))

# Cell type counts by germ layer
germlayer_counts <- meta_data %>%
  group_by(germlayer) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n)) %>%
  arrange(desc(n))


# Cell type counts by treatment
treatment_counts <- meta_data %>%
  group_by(treatment, cellid_type) %>%
  summarise(n = n()) %>%
  mutate(freq = n/sum(n)) %>%
  arrange(treatment, desc(n))




######### Generate cell type proportion plot ############
cellid_out_dir <- "/project/gilad/geena/publication_clean/results/treatment/celltype_annotation/"

# Reload treatment object
batch1_2_3_4_filtnorm_cellid <- readRDS(paste0(cellid_out_dir, "postQC_lognorm_batch1_2_3_4_all_germlayeradded.rds"))
meta_data <- batch1_2_3_4_filtnorm_cellid@meta.data

# Re-define germ layer vectors (needed for ordered_cell_types and color palettes)
ectoderm <- c("Adrenocortical cells", "CNS glia", "CNS neurons", "Neuroendocrine cells", "PNS glia", "PNS neurons", "Photoreceptor cells", "Retinal pigment cells", "Horizontal cells", "Lens fibre cells", "Unipolar brush cells", "Chromaffin cells", "Retinal progenitors and Muller glia", "Ganglion cells", "Amacrine cells", "Bipolar cells", "ENS glia", "Purkinje neurons", "Schwann cells", "Inhibitory neurons", "Astrocytes", "Granule neurons", "Oligodendrocytes", "Sympathoblasts", "Visceral neurons", "Inhibitory interneurons", "Microglia", "Limbic system neurons", "ENS neurons", "Excitatory neurons", "Glioblast", "Neuroblast", "Neuron", "Neuronal IPC", "Oligo", "Radial glia", "Retinal neurons")
mesoderm <- c("Early mesoderm cells", "Cardiac mesoderm cells", "Cardiac progenitor cells", "Cardiomyocytes", "Epicardial fat cells", "Endocardial cells", "Mesangial cells", "Smooth muscle cells", "Vascular endothelial cells", "Mesothelial cells", "Stromal cells", "Fibroblast", "Metanephric cells", "Ureteric bud cells", "Lymphatic endothelial cells", "Lymphoid cells", "Erythroblasts", "Megakaryocytes", "Hematopoietic stem cells", "Myeloid cells", "Antigen presenting cells", "Thymocytes", "Satellite cells", "Skeletal muscle cells")
endoderm <- c("Hepatoblasts", "Stellate cells", "Parietal and chief cells", "Acinar cells", "Ductal cells", "Islet endocrine cells")
epithelium <- c("Goblet cells", "Neuroepithelial cells", "Ciliated epithelial cells", "Intestinal epithelial cells", "Squamous epithelial cells", "Bronchiolar and alveolar epithelial cells")

# Import control cell type counts
control_obj <- readRDS("/project/gilad/geena/publication_clean/results/control/celltype_annotation/combined_control_filtnorm_cellid.rds")
control_meta <- control_obj@meta.data
control_meta <- control_meta %>%
  mutate(
    cellid_type = ifelse(as.character(cellid_type) == "unassigned", "Unassigned", as.character(cellid_type)),
    treatment = "Control"
  )

# Number of cell types detected per individual (control)
control_ct_per_individual <- control_meta %>%
  filter(cellid_type != "Unassigned") %>%
  group_by(donor_id) %>%
  summarise(n_celltypes = n_distinct(cellid_type)) %>%
  arrange(desc(n_celltypes))

median_ct_per_individual_control <- median(control_ct_per_individual$n_celltypes)
min_ct_per_individual_control    <- min(control_ct_per_individual$n_celltypes)
max_ct_per_individual_control    <- max(control_ct_per_individual$n_celltypes)
print(paste0("Control cell types per individual -- median: ", median_ct_per_individual_control,
             ", min: ", min_ct_per_individual_control,
             ", max: ", max_ct_per_individual_control,
             " (n_individuals = ", nrow(control_ct_per_individual), ")"))


# Cells per cell type summary (control; exclude Unassigned)
control_cell_type_counts <- control_meta %>%
  filter(cellid_type != "Unassigned") %>%
  group_by(cellid_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

print(paste0("Control cells per cell type (n_celltypes = ", nrow(control_cell_type_counts),
             ") -- median: ", median(control_cell_type_counts$n),
             ", min: ",        min(control_cell_type_counts$n),
             ", max: ",        max(control_cell_type_counts$n),
             ", total cells: ", sum(control_cell_type_counts$n)))


# Recode treatment meta_data unassigned -> Unassigned
meta_data <- meta_data %>%
  mutate(cellid_type = case_when(
    is.na(cellid_type)                        ~ "Unassigned",
    as.character(cellid_type) == "unassigned" ~ "Unassigned",
    TRUE                                      ~ as.character(cellid_type)
  ))

# Combine treatment and control metadata
meta_data_combined <- bind_rows(
  meta_data %>%
    mutate(cellid_type = as.character(cellid_type)) %>%
    select(individual, treatment, cellid_type, germlayer),
  control_meta %>%
    mutate(cellid_type = as.character(cellid_type)) %>%
    select(individual = donor_id, treatment, cellid_type, germlayer)
)

# Order cellid_type factor levels to group by germ layer
ordered_cell_types <- c(epithelium, ectoderm, mesoderm, endoderm, "iPSCs", "Unassigned")
meta_data_combined$cellid_type <- factor(meta_data_combined$cellid_type, levels = ordered_cell_types)
meta_data_combined$germlayer <- factor(meta_data_combined$germlayer, levels = c("Epithelium", "Ectoderm", "Mesoderm", "Endoderm", "iPSCs", "Unassigned"))

# Generate color palettes for each germ layer
epithelium_colors <- colorRampPalette(brewer.pal(9, "Blues")[3:9])(length(epithelium))
ectoderm_colors   <- colorRampPalette(brewer.pal(9, "Oranges"))(length(ectoderm))
mesoderm_colors   <- colorRampPalette(c("#CCF0E8", "#88CCBB", "#44AA99", "#117777", "#004444"))(length(mesoderm))
endoderm_colors   <- colorRampPalette(brewer.pal(9, "Purples")[3:9])(length(endoderm))
iPSC_color        <- "#CC6677"
unassigned_color  <- "grey50"

cellid_colors <- c(
  setNames(epithelium_colors, epithelium),
  setNames(ectoderm_colors, ectoderm),
  setNames(mesoderm_colors, mesoderm),
  setNames(endoderm_colors, endoderm),
  setNames(iPSC_color, "iPSCs"),
  setNames(unassigned_color, "Unassigned")
)

# Order treatments
treatment_order <- c("Ale", "Bri", "Cer", "Cri", "Lor", "Control")
meta_data_combined$treatment <- factor(meta_data_combined$treatment, levels = treatment_order)

# Summarize data
cell_type_summary <- meta_data_combined %>%
  group_by(individual, treatment, cellid_type) %>%
  summarise(count = n(), .groups = "drop")

# Summarize data without individual (for treatment-level plot)
cell_type_summary_treatment <- meta_data_combined %>%
  group_by(treatment, cellid_type) %>%
  summarise(count = n(), .groups = "drop")

# Publication theme
pub_theme <- theme_minimal(base_size = 14) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y  = element_text(size = 12),
    axis.title   = element_text(size = 14),
    legend.text  = element_text(size = 11),
    legend.title = element_text(size = 13),
    strip.text   = element_text(size = 12),
    plot.title   = element_text(size = 16, face = "bold")
  )

fig_out_dir <- "/project/gilad/geena/publication_clean/results/treatment/figures/celltype_proportion/"
dir.create(fig_out_dir, showWarnings = FALSE, recursive = TRUE)

# Bar plot: proportion by treatment
p_treatment <- ggplot(cell_type_summary_treatment, aes(x = treatment, y = count, fill = cellid_type)) +
  geom_bar(stat = "identity", position = "fill", color = "white", linewidth = 0.08) +
  scale_fill_manual(values = cellid_colors, labels = relabel_celltypes) +
  labs(title = "Cell Type Composition by Treatment",
       x = "Treatment", y = "Proportion of Cells", fill = "Cell Type") +
  pub_theme


ggsave(paste0(fig_out_dir, "celltype_proportion_by_treatment.pdf"), p_treatment, width = 12, height = 7, dpi = 300)
ggsave(paste0(fig_out_dir, "celltype_proportion_by_treatment.png"), p_treatment, width = 12, height = 7, dpi = 300)

# Faceted by individual (shared individuals only)
shared_individuals <- intersect(
  meta_data_combined %>% filter(treatment != "Control") %>% pull(individual) %>% unique(),
  meta_data_combined %>% filter(treatment == "Control") %>% pull(individual) %>% unique()
)

cell_type_summary_shared <- meta_data_combined %>%
  filter(individual %in% shared_individuals) %>%
  group_by(individual, treatment, cellid_type) %>%
  summarise(count = n(), .groups = "drop")


p_individual <- ggplot(cell_type_summary_shared, aes(x = treatment, y = count, fill = cellid_type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = cellid_colors, labels = relabel_celltypes) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(title = "Cell Type Composition by Individual and Treatment",
       x = "Treatment", y = "Proportion of Cells", fill = "Cell Type") +
  pub_theme +
  facet_wrap(~ individual, scales = "free_y") +
  theme(
    legend.text     = element_text(size = 16),
    legend.title    = element_text(size = 18),
    legend.key.size = unit(1, "cm")
  )



ggsave(paste0(fig_out_dir, "celltype_proportion_by_individual.pdf"), p_individual, width = 28, height = 22, dpi = 300)
ggsave(paste0(fig_out_dir, "celltype_proportion_by_individual.png"), p_individual, width = 28, height = 22, dpi = 300)


######### eQTL cell type summary ##########
# 21 cell types used for eQTL mapping (cell types passing >=15 cells per pseudobulk
# and >=20 donors per treatment threshold across both treatment and control)
eqtl_celltypes <- c(
  "Bronchiolar and alveolar epithelial cells", "Cardiac mesoderm cells",
  "Ciliated epithelial cells", "Early mesoderm cells", "ENS glia",
  "Epicardial fat cells", "Glioblast", "Goblet cells", "Hepatoblasts",
  "iPSCs", "Metanephric cells", "Neuroblast", "Neuron", "Neuronal IPC",
  "Radial glia", "Retinal neurons", "Retinal pigment cells",
  "Retinal progenitors and Muller glia", "Squamous epithelial cells",
  "Stromal cells", "Vascular endothelial cells"
)

eqtl_celltypes_df <- tibble(cellid_type = eqtl_celltypes) %>%
  mutate(germlayer = case_when(
    cellid_type %in% ectoderm   ~ "Ectoderm",
    cellid_type %in% endoderm   ~ "Endoderm",
    cellid_type %in% mesoderm   ~ "Mesoderm",
    cellid_type %in% epithelium ~ "Epithelium",
    cellid_type == "iPSCs"      ~ "iPSCs",
    TRUE                        ~ "Other"
  )) %>%
  arrange(germlayer, cellid_type)

print(paste0("Number of eQTL cell types: ", nrow(eqtl_celltypes_df)))
print("Germ layer breakdown of eQTL cell types:")
print(eqtl_celltypes_df %>% count(germlayer))
print(eqtl_celltypes_df)


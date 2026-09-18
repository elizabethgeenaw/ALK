# PAPER OUTPUTS:
# S2_sample_summary.csv
# S3_celltype_counts.csv
# S4_individual_summary_tx_and_control.csv 

library(dplyr)
library(tidyr)

# This script generates supplementary tables from the final annotated Seurat objects.
# Inputs:
# Treatment: postQC_lognorm_batch1_2_3_4_all_germlayeradded.rds (output of script 4)
# Control: combined_control_filtnorm_cellid.rds (output of control script 3)

cellid_out_dir <- "/project/gilad/geena/publication_clean/results/treatment/celltype_annotation/"
control_cellid_dir <- "/project/gilad/geena/publication_clean/results/control/celltype_annotation/"
supp_out_dir <- paste0(cellid_out_dir, "supplementary_tables/")
dir.create(supp_out_dir, recursive = TRUE, showWarnings = FALSE)

# Load the final annotated treatment object
obj <- readRDS(paste0(cellid_out_dir, "postQC_lognorm_batch1_2_3_4_all_germlayeradded.rds"))
meta_data <- obj@meta.data
meta_data <- meta_data %>%
  mutate(cellid_type = ifelse(is.na(cellid_type), "unassigned", cellid_type))

# Load the final annotated control object
control_obj <- readRDS(paste0(control_cellid_dir, "combined_control_filtnorm_cellid.rds"))
control_meta <- control_obj@meta.data
control_meta <- control_meta %>%
  mutate(cellid_type = ifelse(is.na(cellid_type), "unassigned", cellid_type))

######### Sample metadata table (per individual x batch x treatment) #########
sample_metadata <- meta_data %>%
  group_by(individual, batch, treatment) %>%
  summarise(
    n_cells = n(),
    median_genes = median(nFeature_RNA),
    median_UMIs = median(nCount_RNA),
    median_pct_mito = round(median(percent_mito), 2),
    n_celltypes = n_distinct(cellid_type[cellid_type != "unassigned"]),
    .groups = "drop"
  ) %>%
  arrange(batch, treatment, individual)


######### Overall sample summary #########
sample_summary <- sample_metadata %>%
  summarise(
    n_individuals = n_distinct(individual),
    n_samples = n(),
    total_cells = sum(n_cells),
    median_cells_per_sample = median(n_cells),
    median_genes_per_cell = median(median_genes),
    median_UMIs_per_cell = median(median_UMIs),
    median_pct_mito = round(median(median_pct_mito), 2)
  )

print(sample_summary)
write.csv(sample_summary, paste0(cellid_out_dir, "batch1_2_3_4_sample_overall_summary.csv"), row.names = FALSE)

######### Per-individual summary (treatment + control combined) #########
donor_sex <- read.csv("/project/gilad/geena/publication_clean/reference_sets/ALK_sex.csv") %>%
  dplyr::rename("individual" = "line") %>%
  mutate(individual = paste0("NA", as.character(individual)))

individual_summary_treatment <- meta_data %>%
  group_by(individual) %>%
  summarise(
    n_cells = n(),
    median_genes = median(nFeature_RNA),
    median_UMIs = median(nCount_RNA),
    median_pct_mito = round(median(percent_mito), 2),
    n_celltypes = n_distinct(cellid_type[cellid_type != "unassigned"]),
    n_treatments = n_distinct(treatment),
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells)) %>%
  mutate(dataset = "treatment") %>%
  left_join(donor_sex, by = "individual")

individual_summary_control <- control_meta %>%
  group_by(individual = donor_id) %>%
  summarise(
    n_cells = n(),
    median_genes = median(nFeature_RNA),
    median_UMIs = median(nCount_RNA),
    median_pct_mito = round(median(percent.mt), 2),
    n_celltypes = n_distinct(cellid_type[cellid_type != "unassigned"]),
    n_treatments = NA_integer_,
    .groups = "drop"
  ) %>%
  arrange(desc(n_cells)) %>%
  mutate(dataset = "control") %>%
  left_join(donor_sex, by = "individual")

individual_summary <- full_join(
  individual_summary_treatment %>%
    dplyr::select(-dataset) %>%
    dplyr::rename_with(~ paste0("treatment_", .), -c(individual, sex)),
  individual_summary_control %>%
    dplyr::select(-dataset, -sex) %>%
    dplyr::rename_with(~ paste0("control_", .), -individual),
  by = "individual"
) %>%
  dplyr::select(individual, sex, starts_with("treatment_"), starts_with("control_")) %>%
  filter(!is.na(treatment_n_cells)) %>%
  arrange(treatment_n_cells)


# Median number of cell types per individual (excluding unassigned), treatment only
median_celltypes_per_individual <- median(individual_summary_treatment$n_celltypes)
cat("Median cell types per individual:", median_celltypes_per_individual, "\n")

# Median number of cells per cell type (excluding unassigned)
cells_per_ct <- meta_data %>%
  filter(cellid_type != "unassigned") %>%
  group_by(cellid_type) %>%
  summarise(n = n(), .groups = "drop")

median_cells_per_ct <- median(cells_per_ct$n)
cat("Median cells per cell type:", median_cells_per_ct, "\n")

######### Per-batch summary #########
batch_summary <- meta_data %>%
  group_by(batch) %>%
  summarise(
    n_individuals = n_distinct(individual),
    n_cells = n(),
    median_genes = median(nFeature_RNA),
    median_UMIs = median(nCount_RNA),
    median_pct_mito = round(median(percent_mito), 2),
    .groups = "drop"
  )


######### Per-treatment summary #########
treatment_summary <- meta_data %>%
  group_by(treatment) %>%
  summarise(
    n_individuals = n_distinct(individual),
    n_cells = n(),
    median_genes = median(nFeature_RNA),
    median_UMIs = median(nCount_RNA),
    median_pct_mito = round(median(percent_mito), 2),
    n_celltypes = n_distinct(cellid_type[cellid_type != "unassigned"]),
    .groups = "drop"
  )


######### Combined treatment + control cell type counts #########
# Treatment cell type counts
treatment_ct_counts <- meta_data %>%
  group_by(cellid_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::rename(treatment_n = n, treatment_freq = freq)

# Control cell type counts
control_ct_counts <- control_meta %>%
  group_by(cellid_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  mutate(freq = n / sum(n)) %>%
  dplyr::rename(control_n = n, control_freq = freq)

# Join treatment and control counts
combined_ct_counts <- full_join(treatment_ct_counts, control_ct_counts, by = "cellid_type") %>%
  mutate(
    treatment_n = replace_na(treatment_n, 0),
    control_n = replace_na(control_n, 0),
    treatment_freq = replace_na(treatment_freq, 0),
    control_freq = replace_na(control_freq, 0),
    total_n = treatment_n + control_n
  ) %>%
  arrange(desc(total_n))

write.csv(combined_ct_counts, paste0(supp_out_dir, "combined_treatment_control_celltype_counts.csv"), row.names = FALSE)

# Combined summary: total cells in treatment vs control
combined_summary <- data.frame(
  condition = c("Treatment", "Control"),
  total_cells = c(nrow(meta_data), nrow(control_meta)),
  n_individuals = c(n_distinct(meta_data$individual), n_distinct(control_meta$donor_id)),
  n_celltypes = c(n_distinct(meta_data$cellid_type[meta_data$cellid_type != "unassigned"]),
                  n_distinct(control_meta$cellid_type[control_meta$cellid_type != "unassigned"]))
)

print(combined_summary)
write.csv(combined_summary, paste0(supp_out_dir, "combined_treatment_control_summary.csv"), row.names = FALSE)

######### Manuscript statistics (Methods / Results fact-check) #########
# Cell type
ct_counts <- meta_data %>%
  filter(cellid_type != "unassigned") %>%
  group_by(cellid_type) %>%
  summarise(n = n(), .groups = "drop")
cat("Total assigned cell types:", nrow(ct_counts), "\n")
cat("Median cells per cell type:", median(ct_counts$n), "\n")

# Germ layer
germlayer_counts <- meta_data %>%
  filter(cellid_type != "unassigned") %>%
  group_by(germlayer) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))
cat("\nGerm layer counts (treatment):\n")
print(germlayer_counts)

# Treatment cell types per individual
treatment_indiv_ct <- meta_data %>%
  filter(cellid_type != "unassigned") %>%
  group_by(individual) %>%
  summarise(n_celltypes = n_distinct(cellid_type), .groups = "drop")
cat("\nTreatment cell types per individual (pooled across all treatments):\n")
cat("n individuals:", nrow(treatment_indiv_ct), "\n")
cat("Min:", min(treatment_indiv_ct$n_celltypes), "\n")
cat("Max:", max(treatment_indiv_ct$n_celltypes), "\n")
cat("Median:", median(treatment_indiv_ct$n_celltypes), "\n")

# Control cell types per individual
control_indiv_ct <- control_meta %>%
  filter(cellid_type != "unassigned") %>%
  group_by(donor_id) %>%
  summarise(n_celltypes = n_distinct(cellid_type), .groups = "drop")
cat("\nControl cell types per individual:\n")
cat("n individuals:", nrow(control_indiv_ct), "\n")
cat("Min:", min(control_indiv_ct$n_celltypes), "\n")
cat("Max:", max(control_indiv_ct$n_celltypes), "\n")
cat("Median:", median(control_indiv_ct$n_celltypes), "\n")

# Write combined per-individual cell type counts for Table S4 (treatment + control)
treatment_indiv_ct_labeled <- treatment_indiv_ct %>%
  mutate(dataset = "treatment")

control_indiv_ct_labeled <- control_indiv_ct %>%
  dplyr::rename(individual = donor_id) %>%
  mutate(dataset = "control")

combined_indiv_ct <- bind_rows(treatment_indiv_ct_labeled, control_indiv_ct_labeled) %>%
  dplyr::select(dataset, individual, n_celltypes)

# Check overlap between treatment and control cell types
treatment_cts <- unique(meta_data$cellid_type[meta_data$cellid_type != "unassigned"])
control_cts   <- unique(control_meta$cellid_type[control_meta$cellid_type != "unassigned"])
missing_in_control <- setdiff(treatment_cts, control_cts)
cat("\nTreatment cell types absent from control:", length(missing_in_control), "\n")
if (length(missing_in_control) > 0) cat("  Missing:", paste(sort(missing_in_control), collapse = ", "), "\n")

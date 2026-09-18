library(MatrixEQTL)
library(dplyr)
library(tidyverse)

####### MATRIX EQTL ######
# Function to create output directory if it doesn't exist
prepare_output_directory <- function(output_folder_path) {
  if (!dir.exists(output_folder_path)) {
    dir.create(output_folder_path, recursive = TRUE)
    cat("Created output directory:", output_folder_path, "\n")
  }
}

# Function to check if file exists and prompt user to overwrite
check_file_existence <- function(file_path) {
  if (file.exists(file_path)) {
    cat("Warning: The file", file_path, "already exists.\n")
    response <- readline(prompt = "Do you want to overwrite it? (y/n): ")
    return(tolower(response) == "y")
  }
  return(TRUE)
}

# Run MatrixEQTL
single_eqtl_mapping_with_error <- function(celltype, chr, pseudobulk_folder_prefix, output_folder_prefix, maf) {
  print(celltype)
  print(chr)
  
  # Prepare output directory
  output_folder_path <- paste0("/project/gilad/geena/publication_clean/results/control/matrixeqtl/", output_folder_prefix)
  prepare_output_directory(output_folder_path)
  
  tryCatch({
    useModel = modelLINEAR
    SNP_file_name <- paste("/project/gilad/geena/ALK_analysis/matrixeqtl_control_input/geno/pseudoreplicate_control/", chr, "_snp.txt", sep="")
    snps_location_file_name <- paste("/project/gilad/geena/ALK_analysis/matrixeqtl_control_input/geno/pseudoreplicate_control/", chr, "_snpsloc.txt", sep="")
    expression_file_name <- paste("/project/gilad/geena/publication_clean/matrixeqtl_input/control/", pseudobulk_folder_prefix, "/", celltype, "_all_postqc_tmm_logcpm_int.txt", sep="")
    genepos <- read.csv("/project2/gilad/wenhe/EB/combined_3_4/data/ref_set/concise_GRCh38-2020_gtf_strand_considered.csv") %>% dplyr::select(gene, chr, s1 = start, s2 = start)
    snpspos <- read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
    
    pvOutputThreshold_cis <- 1
    cisDist <- 100000 # Ben used 50kb window (for any gene testing, many fewer SNPs; less SNP-gene pairs); strongest eQTL probably in promoter and very proximal
    
    # Load gene expression data
    gene <- SlicedData$new()
    gene$fileDelimiter <- "\t"
    gene$fileOmitCharacters <- "NA"
    gene$fileSkipRows <- 1
    gene$fileSkipColumns <- 1
    gene$fileSliceSize <- 2000
    gene$LoadFile(expression_file_name)
    sample_id <- gene$columnNames
    
    # Ensure sample IDs are formatted consistently
    sample_id <- gsub("_subset[0-9]+", "_Control", sample_id)
    
    # Load genotype data
    snps <- SlicedData$new()
    snps$fileDelimiter <- "\t"
    snps$fileOmitCharacters <- "NA"
    snps$fileSkipRows <- 1
    snps$fileSkipColumns <- 1
    snps$fileSliceSize <- 2000
    snps$LoadFile(SNP_file_name)
    snps$columnNames <- gsub("-", "_", snps$columnNames)
    snps$ColumnSubsample(match(sample_id, snps$columnNames))
    
    covariates_file_name <- list.files(
      paste0("/project/gilad/geena/ALK_analysis_clean/matrixeqtl_input/control/", pseudobulk_folder_prefix, "/cov_elbow/pheno"),
      pattern = paste0("^", celltype, "_\\d+phenoPC\\.txt$"),
      full.names = TRUE)
    
    # Load covariates
    cvrt <- SlicedData$new()
    cvrt$fileDelimiter <- "\t"
    cvrt$fileOmitCharacters <- "NA"
    cvrt$fileSkipRows <- 1
    cvrt$fileSkipColumns <- 1
    if (length(covariates_file_name) > 0) {
      cvrt$LoadFile(covariates_file_name)
    }
    
    # Set result file path and check for existence
    result_file_name_cis <- paste0(output_folder_path, "/Control_", celltype, "_", chr, "_cis.txt")
    if (!check_file_existence(result_file_name_cis)) {
      next  # Skip to the next iteration if the user chooses not to overwrite
    }
    
    # Make temporary directory
    tmpdir_name = "/project/gilad/geena/publication_clean/results/control/matrixeqtl/temp"
    if (!dir.exists(tmpdir_name)) {
      dir.create(tmpdir_name, recursive = TRUE)
      cat("Created output directory:", tmpdir_name, "\n")
    }
    output_file_name_cis <- tempfile(tmpdir = tmpdir_name)
    
    # Prep files for this condition "Control"
    sample_id_cond <- sample_id[grepl("Control", sample_id)]
    gene_cond <- gene$Clone()
    
    gene_cond$columnNames <- gsub("_subset[0-9]+", "_Control", gene_cond$columnNames)
    gene_cond$ColumnSubsample(match(sample_id_cond, gene_cond$columnNames))
    
    cvrt_cond <- cvrt$Clone()
    cvrt_cond$columnNames <- gsub("_subset[0-9]+", "_Control", cvrt_cond$columnNames)
    cvrt_cond$ColumnSubsample(match(sample_id_cond, cvrt_cond$columnNames))
    
    snps_cond <- snps$Clone()
    snps_cond$ColumnSubsample(match(sample_id_cond, snps_cond$columnNames))
    
    # Filter SNPs, do in slices
    cond_snps <- snps_cond$Clone()
    cond_maf.list <- vector('list', length(cond_snps))
    cond_het.list <- vector('list', length(cond_snps))
    cond_homo.list <- vector('list', length(cond_snps))
    
    for(sl in 1:length(cond_snps)) {
      slice <- cond_snps[[sl]]
      cond_maf.list[[sl]] <- rowMeans(slice, na.rm = TRUE) / 2
      cond_maf.list[[sl]] <- pmin(cond_maf.list[[sl]], 1 - cond_maf.list[[sl]]) > maf
      geno_count <- apply(slice, 1, function(x) length(unique(x)))
      cond_het.list[[sl]] <- rowSums(slice == 1, na.rm = TRUE) < (ncol(slice) - 1)
      cond_homo.list[[sl]] <- rowSums(slice == 1, na.rm = TRUE) != 0 | (rowSums(slice == 1, na.rm = TRUE) == 0 & rowSums(slice == 2, na.rm = TRUE) > 1 & rowSums(slice == 2, na.rm = TRUE) < (ncol(slice) - 1))
    }
    
    cond_maf <- unlist(cond_maf.list)
    cond_het <- unlist(cond_het.list)
    cond_homo <- unlist(cond_homo.list)
    
    print(paste("SNPs for Control", celltype, chr, sum(cond_maf & cond_het & cond_homo), "out of", nrow(snps_cond)))
    
    snps_cond$RowReorder(cond_maf & cond_het)
    snpspos_cond <- snpspos[match(snps_cond$GetAllRowNames(), snpspos$snp),]
    
    
    # Run the analysis
    me <- Matrix_eQTL_main(
      snps = snps_cond,
      gene = gene_cond,
      cvrt = cvrt_cond,
      useModel = useModel,
      verbose = TRUE,
      output_file_name.cis = output_file_name_cis,
      pvOutputThreshold.cis = pvOutputThreshold_cis,
      snpspos = snpspos_cond,
      genepos = genepos,
      cisDist = cisDist,
      min.pv.by.genesnp = FALSE,
      noFDRsaveMemory = FALSE
    )
    
    unlink(output_file_name_cis)
    
    # Results
    write.table(me$cis$eqtls, file = result_file_name_cis, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    print(paste("Control", celltype, chr, "done"))
  }, error = function(e) {
    cat("Error in processing Control", "for", celltype, "on chromosome", chr, "\nError Message:", e$message, "\n")
  })
}


# Get p-values for plotting
get_pvalue_distributions <- function(treatment, celltype, files) {
  pattern <- paste(treatment, celltype, sep = "_")
  print(paste("Processing:", pattern))
  
  # Find target files matching the pattern
  target_files <- files[grepl(pattern, files)]
  if (length(target_files) == 0) {
    message("No files found for", pattern)
    next
  }
  
  # Extract p-values from the files
  p_values <- sapply(target_files, function(file) {
    data <- read.table(file, header = TRUE)
    return(data$pvalue)  # ensure this is the correct column name
  }, simplify = FALSE) %>% unlist
  
  # Create a data frame for plotting
  p_values_df <- data.frame(p_value = p_values)
  return(p_values_df)
}

# Get valid cell types
pseudobulk_folder <- "/pseudobulk_15cells_1indiv_summary/"
ct <- read_tsv(paste0("/project/gilad/geena/ALK_analysis_clean/matrixeqtl_input/control/", pseudobulk_folder, "/valid_celltypes_20indiv.tsv"))
colnames(ct) <- "celltypes"
celltypes <- unique(ct$celltypes)


######## RUN IN PARALLEL #######
library(optparse)
library(doParallel)

# Process command-line arguments
parser <- OptionParser()
parser <- add_option(parser, "--chr_range", type = "character", help = "Chromosome range (e.g., 1-3)")
parser <- add_option(parser, "--nc", type = "integer", default = 1, help = "Number of cores")
out <- parse_args(parser)
chr_range <- out$chr_range
nc <- out$nc

# Split chromosome range into individual chromosomes
chr_list <- as.character(strsplit(chr_range, "-")[[1]])
chr_start <- as.numeric(chr_list[1])
chr_end <- as.numeric(chr_list[2])
chr_range_vector <- paste0("chr", chr_start:chr_end)

# Set up parallel execution
registerDoParallel(cores = nc)

# Outer loop over cell types, inner loop over chromosomes
foreach(CT = celltypes, .combine = 'c') %dopar% {
  foreach(chr = chr_range_vector, .combine = 'c') %do% {
    single_eqtl_mapping_with_error(
      celltype = CT,
      chr = chr,
      pseudobulk_folder_prefix = "pseudobulk_15cells_1indiv_summary/normalized_celltypes_medianUMI1/logcpm_int/",
      output_folder_prefix = "medianUMI1_updated4.29.26",
      maf = 0.05
    )
  }
}

# Stop parallel backend
stopImplicitCluster()


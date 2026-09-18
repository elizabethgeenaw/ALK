library(MatrixEQTL)
library(dplyr)
library(readr)

# Function to create output directory if it doesn't exist
prepare_output_directory <- function(output_folder_path) {
  if (!dir.exists(output_folder_path)) {
    dir.create(output_folder_path, recursive = TRUE)
    cat("Created output directory:", output_folder_path, "\n")
  }
}

# Function to run MatrixEQTL
single_eqtl_mapping_with_error <- function(celltype, chr, pseudobulk_folder_prefix, output_folder_prefix, maf) {
  print(celltype)
  print(chr)
  
  # Prepare output directory using helper function
  output_folder_path <- paste0("/project/gilad/geena/publication_clean/results/treatment/matrixeqtl/", output_folder_prefix)
  prepare_output_directory(output_folder_path)
  
  tryCatch({
    useModel = modelLINEAR
    SNP_file_name <- paste("/project/gilad/geena/ALK_analysis/matrixeqtl_treatment_input/geno/pseudoreplicate/", chr, "_snp.txt", sep="")
    snps_location_file_name <- paste("/project/gilad/geena/ALK_analysis/matrixeqtl_treatment_input/geno/pseudoreplicate/", chr, "_snpsloc.txt", sep="") 
    expression_file_name <- paste("/project/gilad/geena/publication_clean/matrixeqtl_input/treatment/", pseudobulk_folder_prefix, "/", celltype, "_all_postqc_tmm_logcpm_int.txt", sep="")
    
    genepos <- read.csv("/project2/gilad/wenhe/EB/combined_3_4/data/ref_set/concise_GRCh38-2020_gtf_strand_considered.csv") %>% 
      dplyr::select(gene, chr, s1 = start, s2 = start)
    snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
    
    pvOutputThreshold_cis = 1
    cisDist = 100000 
    
    ## Load gene expression data
    gene = SlicedData$new()
    gene$fileDelimiter = "\t"
    gene$fileOmitCharacters = "NA"
    gene$fileSkipRows = 1
    gene$fileSkipColumns = 1
    gene$fileSliceSize = 2000
    gene$LoadFile(expression_file_name)
    sample_id <- gene$columnNames
    
    ## Load genotype data
    snps = SlicedData$new()
    snps$fileDelimiter = "\t"
    snps$fileOmitCharacters = "NA"
    snps$fileSkipRows = 1
    snps$fileSkipColumns = 1
    snps$fileSliceSize = 2000
    snps$LoadFile(SNP_file_name)
    column_names = gsub("-", "_", snps$columnNames) # Replace dashes with underscores
    snps$columnNames = column_names
    snps$ColumnSubsample(match(sample_id, snps$columnNames)) # Subset IDs based on exp data
    
    for (condition in c("Ale", "Bri", "Cri", "Cer", "Lor")) {
      print(condition)
      
      covariates_file_name <- list.files(
        paste0("/project/gilad/geena/publication_clean/matrixeqtl_input/treatment/",
               pseudobulk_folder_prefix,
               "/cov_elbow/pheno_by_treatment/"),
        pattern = paste0("^", celltype, "_", condition, "_\\d+phenoPC\\.txt$"),
        full.names = TRUE
      )
      
      ## Load covariates
      cvrt = SlicedData$new()
      cvrt$fileDelimiter = "\t"
      cvrt$fileOmitCharacters = "NA"
      cvrt$fileSkipRows = 1
      cvrt$fileSkipColumns = 1
      if (length(covariates_file_name) > 0) {
        cvrt$LoadFile(covariates_file_name)
      }
      
      result_file_name_cis <- paste0("/project/gilad/geena/publication_clean/results/treatment/matrixeqtl/", output_folder_prefix, "/", condition, "_", celltype, "_", chr, "_cis.txt")
      
      # Make temporary directory
      tmpdir_name = "/project/gilad/geena/publication_clean/results/treatment/matrixeqtl/temp"
      if (!dir.exists(tmpdir_name)) {
        dir.create(tmpdir_name, recursive = TRUE)
        cat("Created output directory:", tmpdir_name, "\n")
      }
      output_file_name_cis <- tempfile(tmpdir = tmpdir_name)
      
      sample_id_cond <- sample_id[grepl(paste0(condition), sample_id)]
      gene_cond <- gene$Clone() # Clone or the original will be changed
      gene_cond$ColumnSubsample(match(sample_id_cond, gene_cond$columnNames)) # Subset IDs to condition
      
      cvrt_cond <- cvrt$Clone() # Clone or the original will be changed
      cvrt_cond$ColumnSubsample(match(sample_id_cond, cvrt_cond$columnNames)) # Subset IDs to condition
      
      snps_cond <- snps$Clone() # Clone or the original will be changed
      snps_cond$ColumnSubsample(match(sample_id_cond, snps_cond$columnNames)) # Subset IDs to condition
      
      ## Filter SNPs
      cond_snps <- snps_cond$Clone()
      cond_maf.list = vector('list', length(cond_snps))
      cond_het.list = vector('list', length(cond_snps))
      cond_homo.list = vector('list', length(cond_snps))
      
      for(sl in 1:length(cond_snps)) {
        slice = cond_snps[[sl]]
        cond_maf.list[[sl]] = rowMeans(slice, na.rm = TRUE)/2 # calculate alt allele freq
        cond_maf.list[[sl]] = pmin(cond_maf.list[[sl]], 1 - cond_maf.list[[sl]]) > maf # MAF > 0.05 (out of 48, 2.4 -> essentially at least 3); TRUE if passes filter
        geno_count = apply(slice, 1, function(x) length(unique(x))) # count number of distinct genotypes for each SNP
        cond_het.list[[sl]] = rowSums(slice == 1, na.rm = TRUE) < (ncol(slice) - 1) # TRUE if not all samples are heterozygotes (makes sure at least one homozygote is present)
        # if heterozygotes -> TRUE; if no heterozygotes, require at least 2 alt homo AND 2 ref homo; ensure enough of both homozygotes when there are no heterozygotes
        cond_homo.list[[sl]] = rowSums(slice == 1, na.rm = TRUE) != 0 | (rowSums(slice == 1, na.rm = TRUE) == 0 & rowSums(slice == 2, na.rm = TRUE) > 1 & rowSums(slice == 2, na.rm = TRUE) < (ncol(slice) - 1))
      }
      cond_maf = unlist(cond_maf.list)
      cond_het = unlist(cond_het.list)
      cond_homo = unlist(cond_homo.list)
      
      print(paste("SNPs for", condition, celltype, chr, sum(cond_maf & cond_het & cond_homo), "out of", nrow(snps_cond)))
      
      snps_cond$RowReorder(cond_maf & cond_het)
      snpspos_cond <- snpspos[match(snps_cond$GetAllRowNames(), snpspos$snp),]
      

      ## Run the analysis
      me = Matrix_eQTL_main(
        snps = snps_cond,
        gene = gene_cond,
        cvrt = cvrt_cond,
        useModel = useModel,
        errorCovariance = numeric(),
        verbose = TRUE,
        output_file_name.cis = output_file_name_cis,
        pvOutputThreshold.cis = pvOutputThreshold_cis,
        snpspos = snpspos_cond,
        genepos = genepos,
        cisDist = cisDist,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE)
      
      unlink(output_file_name_cis)
      
      ## Results:
      write.table(me$cis$eqtls, file = result_file_name_cis, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
      print(paste(condition, celltype, chr, "done"))
    }
  }, error = function(e) {
    # This function will run if there is an error in the try block
    cat("Error in processing", condition, "for", celltype, "on chromosome", chr, "\nError Message:", e$message, "\n")
  })
}


# Get cell types from filtered pseudobulk
ct <- read_tsv("/project/gilad/geena/publication_clean/matrixeqtl_input/treatment/pseudobulk_15cells_20indiv_summary/valid_celltypes_20indiv.tsv")
colnames(ct) <- "celltypes"
celltypes <- unique(ct$celltypes)


# Parallel execution across chromosomes and cell types
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
      pseudobulk_folder_prefix = "pseudobulk_15cells_20indiv_summary/normalized_celltypes_medianUMI1/logcpm_int/",
      output_folder_prefix = "medianUMI1_updated4.29.26",
      maf = 0.05
    )
  }
}

# Stop parallel backend
stopImplicitCluster()



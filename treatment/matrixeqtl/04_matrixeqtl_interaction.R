library(MatrixEQTL)
library(dplyr)
library(readr)
library(optparse)
library(doParallel)
library(foreach)

prepare_output_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

# Genotype × treatment interaction eQTL mapping using modelLINEAR_CROSS.
# Expression files (one per cell type) contain all 5 treatment conditions + control as columns. 
# For each drug, we subset to that drug's samples + control, apply the covariate file (which encodes the treatment indicator as its last row), 
# and run MatrixEQTL with the interaction model.
interaction_eqtl <- function(celltype, chr, umi_label) {
  message(celltype, " | ", chr, " | ", umi_label)

  # Shared genotype files (combined treatment + control individuals)
  # These were built from the intersection of treatment and control YRI lines.
  geno_dir  <- "/project/gilad/geena/ALK_analysis/matrixeqtl_interaction_input/geno/pseudoreplicate/"
  SNP_file  <- paste0(geno_dir, chr, "_snp.txt")
  snpsloc_file <- paste0(geno_dir, chr, "_snpsloc.txt")

  expr_dir  <- paste0("/project/gilad/geena/publication_clean/matrixeqtl_input/interaction/",
                      umi_label, "/")
  cov_dir   <- paste0(expr_dir, "cov_elbow/")
  out_dir   <- paste0("/project/gilad/geena/publication_clean/results/treatment/matrixeqtl/interaction/",
                      umi_label, "/")
  tmp_dir   <- paste0("/project/gilad/geena/publication_clean/results/treatment/matrixeqtl/interaction/tmp/")

  prepare_output_dir(out_dir)
  prepare_output_dir(tmp_dir)

  expr_file <- paste0(expr_dir, celltype, "_all_postqc_tmm_logcpm_int.txt")
  if (!file.exists(expr_file)) {
    message("  Expression file not found, skipping: ", expr_file)
    return(invisible(NULL))
  }

  cov_file <- list.files(cov_dir,
                         pattern = paste0("^", gsub("([()+])", "\\\\\\1", celltype),
                                          "_\\d+phenoPC\\.txt$"),
                         full.names = TRUE)
  if (length(cov_file) == 0) {
    message("  Covariate file not found for: ", celltype)
    return(invisible(NULL))
  }

  snpspos <- read.table(snpsloc_file, header = TRUE, stringsAsFactors = FALSE)
  genepos <- read.csv("/project2/gilad/wenhe/EB/combined_3_4/data/ref_set/concise_GRCh38-2020_gtf_strand_considered.csv") %>%
    dplyr::select(gene, chr, s1 = start, s2 = start)

  # Load expression (shared across all treatments for this cell type)
  gene <- SlicedData$new()
  gene$fileDelimiter     <- "\t"
  gene$fileOmitCharacters <- "NA"
  gene$fileSkipRows      <- 1
  gene$fileSkipColumns   <- 1
  gene$fileSliceSize     <- 2000
  gene$LoadFile(expr_file)
  sample_id <- gene$columnNames

  # Load SNPs and harmonise column names (replace '-' with '_')
  snps <- SlicedData$new()
  snps$fileDelimiter     <- "\t"
  snps$fileOmitCharacters <- "NA"
  snps$fileSkipRows      <- 1
  snps$fileSkipColumns   <- 1
  snps$fileSliceSize     <- 2000
  snps$LoadFile(SNP_file)
  snps$columnNames <- gsub("-", "_", snps$columnNames)
  snps$ColumnSubsample(match(sample_id, snps$columnNames))

  # Load covariates (one file covers all treatments; treatment indicator is last row)
  cvrt <- SlicedData$new()
  cvrt$fileDelimiter     <- "\t"
  cvrt$fileOmitCharacters <- "NA"
  cvrt$fileSkipRows      <- 1
  cvrt$fileSkipColumns   <- 1
  cvrt$LoadFile(cov_file)

  for (drug in c("Ale", "Bri", "Cer", "Cri", "Lor")) {

    out_file <- paste0(out_dir, drug, "_", celltype, "_", chr, "_interact_cis.txt")
    if (file.exists(out_file)) {
      message("  Already done: ", drug, " | ", chr)
      next
    }

    # Subset to this drug + control
    idx_cond <- grepl(paste0("_Control$|_", drug, "$"), sample_id)
    sample_id_cond <- sample_id[idx_cond]

    gene_cond <- gene$Clone()
    gene_cond$ColumnSubsample(match(sample_id_cond, gene_cond$columnNames))

    cvrt_cond <- cvrt$Clone()
    cvrt_cond$ColumnSubsample(match(sample_id_cond, cvrt$columnNames))

    snps_cond <- snps$Clone()
    snps_cond$ColumnSubsample(match(sample_id_cond, snps$columnNames))

    # Identify control / treatment sample indices for MAF filtering
    treat_row  <- cvrt_cond$FindRow("treatment")$row
    ctrl_idx   <- which(treat_row == 0)
    drug_idx   <- which(treat_row == 1)

    # Filter SNPs: MAF > 0.05 in both control and drug subsets
    filter_maf <- function(snps_sd, idx) {
      unlist(lapply(seq_along(snps_sd), function(i) {
        sl  <- snps_sd[[i]]
        maf <- rowMeans(sl[, idx, drop = FALSE], na.rm = TRUE) / 2
        pmin(maf, 1 - maf) > 0.05
      }))
    }
    valid <- filter_maf(snps_cond, ctrl_idx) & filter_maf(snps_cond, drug_idx)
    snps_cond$RowReorder(valid)

    snpspos_cond <- snpspos[match(snps_cond$GetAllRowNames(), snpspos$snp), ]

    tmp_out <- tempfile(tmpdir = tmp_dir)
    tryCatch({
      me <- Matrix_eQTL_main(
        snps              = snps_cond,
        gene              = gene_cond,
        cvrt              = cvrt_cond,
        useModel          = modelLINEAR_CROSS,
        errorCovariance   = numeric(),
        verbose           = FALSE,
        output_file_name.cis = tmp_out,
        pvOutputThreshold.cis = 1,
        snpspos           = snpspos_cond,
        genepos           = genepos,
        cisDist           = 100000,
        min.pv.by.genesnp = TRUE,
        noFDRsaveMemory   = FALSE
      )
      unlink(tmp_out)
      write.table(me$cis$eqtls, file = out_file,
                  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
      message("  Done: ", drug, " | ", chr)
    }, error = function(e) {
      unlink(tmp_out)
      message("  ERROR (", drug, " | ", chr, "): ", e$message)
    })
  }
}

# To run in parallel
parser <- OptionParser()
parser <- add_option(parser, "--chr_range",  type = "character", help = "e.g. 1-3")
parser <- add_option(parser, "--umi_label",  type = "character", default = "medianUMI1",
                     help = "medianUMI1 or medianUMI10")
parser <- add_option(parser, "--nc",         type = "integer",   default = 1,
                     help = "Number of parallel cores")
out <- parse_args(parser)

chr_parts  <- as.integer(strsplit(out$chr_range, "-")[[1]])
chr_vec    <- paste0("chr", seq(chr_parts[1], chr_parts[2]))
umi_label  <- out$umi_label

expr_dir <- paste0("/project/gilad/geena/publication_clean/matrixeqtl_input/interaction/",
                   umi_label, "/")
celltypes <- sub("_all_postqc_tmm_logcpm_int.txt$", "",
                 list.files(expr_dir, pattern = "_all_postqc_tmm_logcpm_int.txt"))

message("Cell types: ", length(celltypes))
message("Chromosomes: ", paste(chr_vec, collapse = ", "))
message("UMI label: ", umi_label)

registerDoParallel(cores = out$nc)

foreach(CT = celltypes, .combine = "c") %dopar% {
  foreach(chr = chr_vec, .combine = "c") %do% {
    interaction_eqtl(CT, chr, umi_label)
  }
}

stopImplicitCluster()

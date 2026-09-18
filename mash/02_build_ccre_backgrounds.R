#!/usr/bin/env Rscript
#
# Pre-computes all cCRE background annotation files needed by post_mash_functional.Rmd. 
#
# Outputs:
# 1) Per-condition backgrounds (one TSV per celltype x treatment): {bg_dir}/{celltype}_{treatment}.tsv
# 2) Global MatrixEQTL background (all SNPs tested in any condition): {meqtl_dir}/all_matrixeqtl_testedsnps_annotations.tsv
# 3) Per-treatment backgrounds (union of SNPs across cell types): {meqtl_dir}/bg_{treatment}_annotated.tsv

library(data.table)
library(stringr)
library(purrr)
library(tibble)
library(dplyr)


# Paths
dir_meqtl  <- "/project/gilad/geena/publication_clean/results/combined/medianUMI1/combine_chr/"
bg_dir     <- "/project/gilad/geena/publication_clean/results/combined/medianUMI1/backgrounds/cttx"
meqtl_dir  <- "/project/gilad/geena/publication_clean/results/combined/medianUMI1/"

dir.create(bg_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(meqtl_dir, recursive = TRUE, showWarnings = FALSE)


# cCRE reference
message("Loading cCRE reference files...")

read_ccre <- function(path, class) {
  fread(path, header = FALSE, select = c(1, 2, 3, 5),
        col.names = c("chromosome", "start", "end", "cCRE_ID")) %>%
    .[, `:=`(chromosome = gsub("^chr", "", chromosome),
             cCRE_class = class)]
}

PLS  <- read_ccre("/project/gilad/geena/publication_clean/reference_sets/encode_ccre_update2026/GRCh38-cCREs.PLS.bed",  "PLS")
pELS <- read_ccre("/project/gilad/geena/publication_clean/reference_sets/encode_ccre_update2026/GRCh38-cCREs.pELS.bed", "pELS")
dELS <- read_ccre("/project/gilad/geena/publication_clean/reference_sets/encode_ccre_update2026/GRCh38-cCREs.dELS.bed", "dELS")

cCREs <- rbindlist(list(PLS, pELS, dELS))
setkey(cCREs, chromosome, start, end)
message("cCRE reference loaded: ", nrow(cCREs), " entries")


# Helpers
normalize_gene_snp <- function(x) {
  gsub("_chr", "_", as.character(x), ignore.case = TRUE)
}

normalize_chr <- function(x) {
  sub("^chr", "", as.character(x), ignore.case = TRUE)
}

parse_gene_snp <- function(x) {
  dt <- data.table(gene_snp = x)
  dt[, snp := sub("^[^_]*_", "", gene_snp)]
  dt[, chr := normalize_chr(sub(":.*", "", snp))]
  dt[, pos := as.numeric(sub(":.*", "", sub("^[^:]*:", "", snp)))]
  dt
}

# Overlap SNPs against cCREs using foverlaps; returns gene_snp + cCRE_class
annotate_ccre_once <- function(gene_snp_vec) {
  gene_snp_vec <- normalize_gene_snp(gene_snp_vec)

  snps <- parse_gene_snp(gene_snp_vec)
  snps <- unique(snps, by = c("gene_snp", "chr", "pos"))
  snps[, pos_end := pos]
  setkey(snps, chr, pos, pos_end)

  ov <- foverlaps(cCREs, snps,
                  by.x = c("chromosome", "start", "end"),
                  by.y = c("chr", "pos", "pos_end"),
                  nomatch = 0L)

  ov[, .(gene_snp, cCRE_class)]
}

# Parse treatment / celltype from filename "<Tx>_<CellType>_combined_cis.txt"
parse_tx_ct <- function(path) {
  b <- basename(path)
  m <- str_match(b, "^([^_]+)_(.*)_combined_cis\\.txt$")
  data.table(path = path, treatment = m[, 2], celltype = m[, 3])
}

# Index MatrixEQTL files
meqtl_files_all <- list.files(dir_meqtl, pattern = "_combined_cis\\.txt$", full.names = TRUE)
meqtl_idx <- rbindlist(lapply(meqtl_files_all, parse_tx_ct))
meqtl_idx[, key := paste(celltype, treatment, sep = "_")]

message("Found ", nrow(meqtl_idx), " MatrixEQTL condition files")


# 1) Per-condition backgrounds
# One TSV per celltype x treatment key, containing all gene_snp pairs tested in that condition with their cCRE annotation (NA if none)
message("1) Per-condition backgrounds")
message("--------------")

n_skipped <- 0
n_built   <- 0

for (i in seq_len(nrow(meqtl_idx))) {
  key      <- meqtl_idx$key[i]
  out_file <- file.path(bg_dir, paste0(key, ".tsv"))

  if (file.exists(out_file)) {
    n_skipped <- n_skipped + 1
    next
  }

  message("[", i, "/", nrow(meqtl_idx), "] Building: ", key)

  snps <- fread(meqtl_idx$path[i],
                select = c("snps", "gene"),
                colClasses = c(snps = "character", gene = "character"))
  snps[, gene_snp := normalize_gene_snp(paste(gene, snps, sep = "_"))]
  snps <- unique(snps[, .(gene_snp)])

  ann <- annotate_ccre_once(snps$gene_snp)

  # Left join so SNPs with no cCRE overlap get NA
  out <- ann[snps, on = "gene_snp"]
  fwrite(out, out_file, sep = "\t")

  n_built <- n_built + 1
  rm(snps, ann, out); gc()
}

message("Per-condition backgrounds: ", n_built, " built, ", n_skipped, " already existed")


# 2) Global MatrixEQTL background
# Union of all gene_snp pairs tested in any condition, annotated with cCRE class
message("2) Global MatrixEQTL background")
message("--------------")

global_ann_path <- file.path(meqtl_dir, "all_matrixeqtl_testedsnps_annotations.tsv")

if (file.exists(global_ann_path)) {
  message("Already exists, skipping: ", global_ann_path)
} else {
  message("Reading all MatrixEQTL files...")

  all_snps <- rbindlist(lapply(meqtl_files_all, function(p) {
    fread(p, select = c("snps", "gene"), colClasses = "character")[
      , .(gene_snp = normalize_gene_snp(paste(gene, snps, sep = "_")))
    ]
  }))
  all_snps <- unique(all_snps)
  message("Total unique gene_snp pairs: ", nrow(all_snps))

  message("Annotating with cCREs...")
  bg_global_annotated <- annotate_ccre_once(all_snps$gene_snp)

  fwrite(bg_global_annotated, global_ann_path, sep = "\t")
  message("Saved: ", global_ann_path)

  rm(all_snps); gc()
}

# 3) Per-treatment backgrounds
# Union of SNPs tested across all cell types for each treatment, annotated via left-join against the global annotation
message("3) Per-treatment backgrounds")
message("--------------")

# Load global annotation (may have just been built above)
message("Loading global annotation for per-treatment joins...")
bg_global_annotated <- fread(global_ann_path)

treatments <- unique(meqtl_idx$treatment)
message("Treatments found: ", paste(treatments, collapse = ", "))

for (tx in treatments) {
  out_path <- file.path(meqtl_dir, paste0("bg_", tolower(tx), "_annotated.tsv"))

  if (file.exists(out_path)) {
    message("Already exists, skipping: ", tx)
    next
  }

  message("Building per-treatment background: ", tx)

  tx_files <- meqtl_idx[treatment == tx, path]

  tx_snps <- rbindlist(lapply(tx_files, function(p) {
    fread(p, select = c("snps", "gene"), colClasses = "character")[
      , .(gene_snp = normalize_gene_snp(paste(gene, snps, sep = "_")))
    ]
  }))
  tx_snps <- unique(tx_snps)
  message("  ", tx, ": ", nrow(tx_snps), " unique gene_snp pairs")

  # Left join against global annotation so non-overlapping SNPs get NA
  out <- bg_global_annotated[tx_snps, on = "gene_snp"]
  out[, gene_snp := sub("_chr", "_", gene_snp)]
  if ("snp" %in% names(out)) out[, snp := gsub("^chr", "", snp)]

  fwrite(out, out_path, sep = "\t")
  message("  Saved: ", out_path)

  rm(tx_snps, out); gc()
}


# Summary
message("--------------")
message("All backgrounds complete")
message("Per-condition TSVs : ", bg_dir)
message("Global annotation  : ", global_ann_path)
message("Per-treatment TSVs : ", meqtl_dir, "bg_*_annotated.tsv")

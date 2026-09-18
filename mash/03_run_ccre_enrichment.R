#!/usr/bin/env Rscript
#
# Computes cCRE enrichment (Fisher tests) across all LFC thresholds in parallel.
# Saves results as RDS for loading in post_mash_functional_v2.Rmd (sections 11/12).
#
# Run via: sbatch run_ccre_enrichment.sbatch

library(tidyverse)
library(data.table)
library(stringr)
library(purrr)
library(parallel)

# Paths
mash_folder <- "/project/gilad/geena/publication_clean/results/mash/medianUMI1/"
dir_meqtl <- "/project/gilad/geena/publication_clean/results/combined/medianUMI1/combine_chr/"
bg_dir <- "/project/gilad/geena/publication_clean/results/combined/medianUMI1/backgrounds/cttx"
alpha_lfsr <- 0.05
ccre_classes <- c("PLS", "pELS", "dELS")
lfcs <- c(2, 3, 3.5, 4, 5, 6, 7, 8, 9, 10)
out_rds_bg <- file.path(mash_folder, "ccre_enrichment_vs_background.rds")
out_rds_rnr <- file.path(mash_folder, "ccre_enrichment_resp_vs_nonresp.rds")


# Helpers 
normalize_gene_snp <- function(x) gsub("_chr", "_", as.character(x), ignore.case = TRUE)

normalize_celltype <- function(x) {
  x %>% as.character() %>% str_replace_all(" ", "_") %>% str_replace_all("__+", "_")
}

split_cttx <- function(x) {
  parts <- strsplit(x, "_", fixed = TRUE)
  tibble(
    treatment = vapply(parts, function(v) v[length(v)], character(1)),
    celltype  = vapply(parts, function(v) paste(v[-length(v)], collapse = "_"), character(1))
  )
}

split_key_cttx <- function(k) {
  parts <- strsplit(k, "_", fixed = TRUE)[[1]]
  list(celltype = paste(parts[-length(parts)], collapse = "_"), treatment = parts[length(parts)])
}

fisher_ccre <- function(target, universe, ccre_set) {
  target <- intersect(target, universe)
  backg  <- setdiff(universe, target)
  YY <- length(intersect(target, ccre_set))
  YN <- length(setdiff(target, ccre_set))
  NY <- length(intersect(backg, ccre_set))
  NN <- length(setdiff(backg, ccre_set))
  if ((YY + YN) == 0 || (NY + NN) == 0 || (YY + NY) == 0 || (YN + NN) == 0) {
    return(tibble(YY = YY, YN = YN, NY = NY, NN = NN,
                  odds_ratio = NA_real_, conf_low = NA_real_,
                  conf_high = NA_real_, p_value = NA_real_))
  }
  ft <- fisher.test(matrix(c(YY, YN, NY, NN), nrow = 2, byrow = TRUE))
  tibble(YY = YY, YN = YN, NY = NY, NN = NN,
         odds_ratio = unname(ft$estimate), conf_low = ft$conf.int[1],
         conf_high = ft$conf.int[2], p_value = ft$p.value)
}

fisher_resp_vs_nonresp <- function(resp, nonresp, universe, ccre_set) {
  resp    <- intersect(resp, universe)
  nonresp <- intersect(nonresp, universe)
  YY <- length(intersect(resp, ccre_set))
  YN <- length(setdiff(resp, ccre_set))
  NY <- length(intersect(nonresp, ccre_set))
  NN <- length(setdiff(nonresp, ccre_set))
  if ((YY + YN) == 0 || (NY + NN) == 0 || (YY + NY) == 0 || (YN + NN) == 0) {
    return(tibble(YY = YY, YN = YN, NY = NY, NN = NN,
                  odds_ratio = NA_real_, conf_low = NA_real_,
                  conf_high = NA_real_, p_value = NA_real_))
  }
  ft <- fisher.test(matrix(c(YY, YN, NY, NN), nrow = 2, byrow = TRUE))
  tibble(YY = YY, YN = YN, NY = NY, NN = NN,
         odds_ratio = unname(ft$estimate), conf_low = ft$conf.int[1],
         conf_high = ft$conf.int[2], p_value = ft$p.value)
}

fisher_from_counts <- function(YY, YN, NY, NN) {
  ft <- fisher.test(matrix(c(YY, YN, NY, NN), 2, byrow = TRUE))
  list(odds_ratio = unname(ft$estimate), conf_low = ft$conf.int[1],
       conf_high = ft$conf.int[2], p_value = ft$p.value)
}

calc_ccre_condition <- function(lfc, bg_cache, sets, ccre_classes = c("PLS", "pELS", "dELS")) {
  map_dfr(names(bg_cache), function(k) {
    ann      <- bg_cache[[k]]
    universe <- unique(ann$gene_snp)
    tmp      <- split_key_cttx(k)
    ct       <- normalize_celltype(tmp$celltype)
    tx       <- tmp$treatment
    resp    <- sets$resp    %>% mutate(celltype = normalize_celltype(celltype)) %>%
               filter(celltype == ct, treatment == tx) %>% pull(gene_snp)
    nonresp <- sets$nonresp %>% mutate(celltype = normalize_celltype(celltype)) %>%
               filter(celltype == ct, treatment == tx) %>% pull(gene_snp)
    map_dfr(ccre_classes, function(cls) {
      ccre_set <- ann[cCRE_class == cls, unique(gene_snp)]
      bind_rows(
        fisher_ccre(resp,    universe, ccre_set) %>% mutate(group = "Response"),
        fisher_ccre(nonresp, universe, ccre_set) %>% mutate(group = "Nonresponse")
      ) %>% mutate(celltype = ct, treatment = tx, facet = paste(ct, tx, sep = " | "),
                   cCRE_class = cls, lfc = lfc, level = "Condition")
    })
  })
}

calc_ccre_condition_resp_vs_nonresp <- function(lfc, bg_cache, sets,
                                                 ccre_classes = c("PLS", "pELS", "dELS")) {
  map_dfr(names(bg_cache), function(k) {
    ann      <- bg_cache[[k]]
    universe <- unique(ann$gene_snp)
    tmp      <- split_key_cttx(k)
    ct       <- normalize_celltype(tmp$celltype)
    tx       <- tmp$treatment
    resp    <- sets$resp    %>% mutate(celltype = normalize_celltype(celltype)) %>%
               filter(celltype == ct, treatment == tx) %>% pull(gene_snp)
    nonresp <- sets$nonresp %>% mutate(celltype = normalize_celltype(celltype)) %>%
               filter(celltype == ct, treatment == tx) %>% pull(gene_snp)
    map_dfr(ccre_classes, function(cls) {
      ccre_set <- ann[cCRE_class == cls, unique(gene_snp)]
      fisher_resp_vs_nonresp(resp, nonresp, universe, ccre_set) %>%
        mutate(celltype = ct, treatment = tx, cCRE_class = cls, lfc = lfc, level = "Condition")
    })
  })
}

calc_global_ccre <- function(cond_res) {
  cond_res %>%
    group_by(group, cCRE_class) %>%
    summarise(YY = sum(YY), YN = sum(YN), NY = sum(NY), NN = sum(NN), .groups = "drop") %>%
    rowwise() %>%
    mutate(ft = list(fisher_from_counts(YY, YN, NY, NN)),
           odds_ratio = ft$odds_ratio, conf_low = ft$conf_low,
           conf_high = ft$conf_high, p_value = ft$p_value,
           facet = "Global", level = "Global") %>%
    dplyr::select(-ft)
}

calc_celltype_ccre <- function(cond_res) {
  cond_res %>%
    group_by(celltype, group, cCRE_class) %>%
    summarise(YY = sum(YY), YN = sum(YN), NY = sum(NY), NN = sum(NN), .groups = "drop") %>%
    rowwise() %>%
    mutate(ft = list(fisher_from_counts(YY, YN, NY, NN)),
           odds_ratio = ft$odds_ratio, conf_low = ft$conf_low,
           conf_high = ft$conf_high, p_value = ft$p_value,
           facet = celltype, level = "Celltype") %>%
    dplyr::select(-ft)
}

calc_treatment_ccre <- function(cond_res) {
  cond_res %>%
    group_by(treatment, group, cCRE_class) %>%
    summarise(YY = sum(YY), YN = sum(YN), NY = sum(NY), NN = sum(NN), .groups = "drop") %>%
    rowwise() %>%
    mutate(ft = list(fisher_from_counts(YY, YN, NY, NN)),
           odds_ratio = ft$odds_ratio, conf_low = ft$conf_low,
           conf_high = ft$conf_high, p_value = ft$p_value,
           facet = treatment, level = "Treatment") %>%
    dplyr::select(-ft)
}

calc_global_resp_vs_nonresp <- function(df) {
  df %>%
    group_by(cCRE_class) %>%
    summarise(A = sum(A), B = sum(B), C = sum(C), D = sum(D), .groups = "drop") %>%
    rowwise() %>%
    mutate(ft = list(fisher_from_counts(A, B, C, D)),
           odds_ratio = ft$odds_ratio, conf_low = ft$conf_low,
           conf_high = ft$conf_high, p_value = ft$p_value,
           level = "Global") %>%
    dplyr::select(-ft)
}

calc_celltype_resp_vs_nonresp <- function(df) {
  df %>%
    group_by(celltype, cCRE_class) %>%
    summarise(A = sum(A), B = sum(B), C = sum(C), D = sum(D), .groups = "drop") %>%
    rowwise() %>%
    mutate(ft = list(fisher_from_counts(A, B, C, D)),
           odds_ratio = ft$odds_ratio, conf_low = ft$conf_low,
           conf_high = ft$conf_high, p_value = ft$p_value,
           level = "Celltype") %>%
    dplyr::select(-ft)
}

calc_treatment_resp_vs_nonresp <- function(df) {
  df %>%
    group_by(treatment, cCRE_class) %>%
    summarise(A = sum(A), B = sum(B), C = sum(C), D = sum(D), .groups = "drop") %>%
    rowwise() %>%
    mutate(ft = list(fisher_from_counts(A, B, C, D)),
           odds_ratio = ft$odds_ratio, conf_low = ft$conf_low,
           conf_high = ft$conf_high, p_value = ft$p_value,
           level = "Treatment") %>%
    dplyr::select(-ft)
}


# Load data
message("Loading MASH results...")
all_mash <- fread(file.path(mash_folder, "combined_mash_results.tsv"))
all_mash[, gene_snp := normalize_gene_snp(gene_snp)]

lfsr_cols <- grep("_lfsr$", names(all_mash), value = TRUE)
tested_annot <- all_mash[, c("gene_snp", lfsr_cols), with = FALSE] %>%
  melt(id.vars = "gene_snp", variable.name = "cttx", value.name = "lfsr") %>%
  .[, cttx := sub("_lfsr$", "", cttx)] %>%
  as_tibble() %>%
  bind_cols(split_cttx(.$cttx)) %>%
  dplyr::select(gene_snp, celltype, treatment, lfsr)
rm(all_mash); gc()

meqtl_files <- tibble(
  path = list.files(dir_meqtl, pattern = "_combined_cis\\.txt$", full.names = TRUE)
) %>%
  mutate(
    file      = basename(path),
    treatment = str_match(file, "^([^_]+)_")[, 2],
    celltype  = str_match(file, "^[^_]+_(.*)_combined_cis")[, 2],
    key       = paste(celltype, treatment, sep = "_")
  )


# Cache backgrounds
message("Caching ", nrow(meqtl_files), " background files...")
bg_cache <- setNames(
  lapply(meqtl_files$key, function(k) {
    dt <- fread(file.path(bg_dir, paste0(k, ".tsv")))
    dt[, gene_snp := normalize_gene_snp(gene_snp)]
    dt
  }),
  meqtl_files$key
)
message("Cached ", length(bg_cache), " conditions.")


# Response / nonresponse helpers
get_resp_tbl <- function(lfc) {
  fread(file.path(mash_folder, paste0("lfc", lfc),
                  "treatment_responsive_eqtls_signfixed.tsv")) %>%
    .[, gene_snp := normalize_gene_snp(gene_snp)] %>%
    as_tibble() %>%
    distinct(gene_snp, celltype, treatment)
}

get_resp_nonresp <- function(lfc) {
  resp <- get_resp_tbl(lfc)
  sig  <- tested_annot %>%
    filter(lfsr < alpha_lfsr) %>%
    distinct(gene_snp, celltype, treatment)
  list(resp = resp, nonresp = anti_join(sig, resp, by = c("gene_snp", "celltype", "treatment")))
}


# Pre-compute response/nonresponse sets for all LFCs (avoids re-reading TSVs)
message("Pre-computing response/nonresponse sets...")
all_sets <- setNames(
  lapply(lfcs, function(lfc) { message("  LFC ", lfc); get_resp_nonresp(lfc) }),
  paste0("lfc", lfcs)
)
rm(tested_annot); gc()


# Parallel enrichment across LFCs
n_cores <- 4L
message("Running enrichment on ", n_cores, " cores...")


# 1) Response / Nonresponse vs Background
if (file.exists(out_rds_bg)) {
  message("Skipping vs-background (already exists): ", out_rds_bg)
} else {
  message("Response/nonresponse vs background")
  message("--------")
  results_all <- setNames(
    mclapply(lfcs, function(lfc) {
      message("  LFC ", lfc)
      sets         <- all_sets[[paste0("lfc", lfc)]]
      cond_res     <- calc_ccre_condition(lfc, bg_cache, sets)
      global_res   <- calc_global_ccre(cond_res)    %>% mutate(lfc = lfc)
      celltype_res <- calc_celltype_ccre(cond_res)  %>% mutate(lfc = lfc)
      treatment_res <- calc_treatment_ccre(cond_res) %>% mutate(lfc = lfc)
      list(condition = cond_res, global = global_res,
           celltype = celltype_res, treatment = treatment_res)
    }, mc.cores = n_cores),
    paste0("lfc", lfcs)
  )
  saveRDS(results_all, out_rds_bg)
  message("Saved: ", out_rds_bg)
}


# 2) Response vs Nonresponse
message("Response vs nonresponse")
message("--------")
results_all_r_nr <- setNames(
  mclapply(lfcs, function(lfc) {
    message("  LFC ", lfc)
    sets     <- all_sets[[paste0("lfc", lfc)]]
    cond_res <- calc_ccre_condition_resp_vs_nonresp(lfc, bg_cache, sets) %>%
      rename(A = YY, B = YN, C = NY, D = NN)
    global_res <- calc_global_resp_vs_nonresp(cond_res) %>% mutate(lfc = lfc)
    celltype_res <- calc_celltype_resp_vs_nonresp(cond_res) %>% mutate(lfc = lfc)
    treatment_res <- calc_treatment_resp_vs_nonresp(cond_res) %>% mutate(lfc = lfc)
    list(condition = cond_res, global = global_res,
         celltype = celltype_res, treatment = treatment_res)
  }, mc.cores = n_cores),
  paste0("lfc", lfcs)
)
saveRDS(results_all_r_nr, out_rds_rnr)
message("Saved: ", out_rds_rnr)

message("Done.")

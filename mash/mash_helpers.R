#!/usr/bin/env Rscript
#
# Helper functions for post_mash_analysis_v2.Rmd
#
# Sections:
# 1) eQTL sharing / specificity pipeline
# 2) TSS distance data prep + plotting

library(tidyverse)
library(dplyr)
library(purrr)

################ eQTL SHARING / SPECIFICITY HELPERS ################
# get_significance_and_effect_sharing_from_df
#
# For a given cell type, adds binary columns indicating:
#   - whether each eQTL is significant (lfsr < thresh) in each treatment
#   - whether each pairwise treatment combination shares effect size (same sign
#     and ratio within `factor`)
#
# Returns a data frame with gene_snp + sig_* + effectsim_* columns.
get_significance_and_effect_sharing_from_df <- function(
  celltype, df,
  lfsr_thresh = 0.05,
  factor = 1.5,
  FUN = identity
) {
  lfsr_cols  <- grep(paste0("^", celltype, "_.*_lfsr$"), colnames(df), value = TRUE)
  beta_cols  <- grep(paste0("^", celltype, "_.*_beta$"), colnames(df), value = TRUE)
  treatments <- gsub(paste0(celltype, "_|_lfsr$"), "", lfsr_cols)

  lfsr_mat <- as.matrix(df[, lfsr_cols])
  beta_mat <- as.matrix(df[, beta_cols])
  rownames(lfsr_mat) <- rownames(beta_mat) <- df$gene_snp
  colnames(lfsr_mat) <- colnames(beta_mat) <- treatments

  # Binary significance per treatment
  sig_df <- as.data.frame(matrix(0, nrow = nrow(lfsr_mat), ncol = length(treatments)))
  colnames(sig_df) <- paste0(celltype, "_sig_", treatments)
  rownames(sig_df) <- rownames(lfsr_mat)
  for (i in seq_along(treatments)) {
    sig_df[, i] <- as.integer(lfsr_mat[, i] < lfsr_thresh)
  }

  # Binary pairwise effect similarity (same sign + magnitude within factor)
  pairwise <- combn(treatments, 2, simplify = FALSE)
  share_df <- as.data.frame(matrix(0, nrow = nrow(lfsr_mat), ncol = length(pairwise)))
  colnames(share_df) <- sapply(pairwise, function(p)
    paste0(celltype, "_effectsim_", p[1], "_", p[2]))
  rownames(share_df) <- rownames(lfsr_mat)

  for (k in seq_along(pairwise)) {
    i <- pairwise[[k]][1]
    j <- pairwise[[k]][2]
    sig_idx <- which(lfsr_mat[, i] < lfsr_thresh | lfsr_mat[, j] < lfsr_thresh)

    same_sign <- sign(beta_mat[sig_idx, i]) == sign(beta_mat[sig_idx, j]) &
                 beta_mat[sig_idx, i] != 0 & beta_mat[sig_idx, j] != 0
    valid_idx <- sig_idx[same_sign]
    ratio <- abs(beta_mat[valid_idx, i] / beta_mat[valid_idx, j])
    shared <- which(ratio > 1 / factor & ratio < factor)
    share_df[valid_idx[shared], k] <- 1
  }

  cbind(gene_snp = rownames(lfsr_mat), sig_df, share_df)
}

# For each treatment pair within a cell type, adds a SHARED column:
# significant in either condition AND effect size is shared.
add_pairwise_shared_columns <- function(sharing_df, celltype) {
  sig_cols   <- grep(paste0("^", celltype, "_sig_"), colnames(sharing_df), value = TRUE)
  treatments <- gsub(paste0(celltype, "_sig_"), "", sig_cols)

  shared_df <- sharing_df
  for (pair in combn(treatments, 2, simplify = FALSE)) {
    t1 <- pair[1]; t2 <- pair[2]
    sig1      <- paste0(celltype, "_sig_", t1)
    sig2      <- paste0(celltype, "_sig_", t2)
    effectsim <- paste0(celltype, "_effectsim_", t1, "_", t2)
    shared_col <- paste0(celltype, "_SHARED_", t1, "_", t2)

    shared_df[[shared_col]] <- as.integer(
      (shared_df[[sig1]] == 1 | shared_df[[sig2]] == 1) &
       shared_df[[effectsim]] == 1
    )
  }
  shared_df
}


# For each treatment within a cell type, adds a specific column:
# significant in this treatment AND no effect sharing with any other treatment.
add_treatment_specific_columns <- function(sharing_df, celltype) {
  sig_cols   <- grep(paste0("^", celltype, "_sig_"), colnames(sharing_df), value = TRUE)
  treatments <- gsub(paste0(celltype, "_sig_"), "", sig_cols)

  output_df <- sharing_df
  for (tx in treatments) {
    sig_this  <- paste0(celltype, "_sig_", tx)
    effect_cols <- intersect(
      c(paste0(celltype, "_effectsim_", tx, "_", setdiff(treatments, tx)),
        paste0(celltype, "_effectsim_", setdiff(treatments, tx), "_", tx)),
      colnames(sharing_df)
    )
    specific_col <- paste0(celltype, "_", tx, "_specific")

    output_df[[specific_col]] <- as.integer(
      output_df[[sig_this]] == 1 &
      rowSums(output_df[, effect_cols, drop = FALSE]) == 0
    )
  }
  output_df
}



# Adds a column indicating eQTLs that are:
#   - significant in at least one of `included` treatments
#   - effect shared among all `included` treatments
#   - effect NOT shared with any `excluded` treatment
add_group_specific_sharing_column <- function(
  df, celltype, included, excluded, prefix = NULL
) {
  if (is.null(prefix)) {
    prefix <- paste0(celltype, "_GROUPSHARED_", paste(included, collapse = "_"))
  }

  sig_included <- rowSums(
    df[, paste0(celltype, "_sig_", included), drop = FALSE] == 1) > 0

  include_pairs <- combn(included, 2, simplify = FALSE)
  expected_include_cols <- paste0(celltype, "_effectsim_",
                                   sapply(include_pairs, paste, collapse = "_"))
  include_cols_exist <- expected_include_cols[expected_include_cols %in% colnames(df)]

  if (length(include_cols_exist) < length(expected_include_cols)) {
    warning(paste("Missing pairwise effectsim columns for", celltype))
    effect_included <- rep(FALSE, nrow(df))
  } else {
    effect_included <- rowSums(
      df[, include_cols_exist, drop = FALSE] == 1) == length(include_cols_exist)
  }

  cross_pairs <- expand.grid(included, excluded, stringsAsFactors = FALSE)
  cross_cols <- intersect(
    c(paste0(celltype, "_effectsim_", cross_pairs$Var1, "_", cross_pairs$Var2),
      paste0(celltype, "_effectsim_", cross_pairs$Var2, "_", cross_pairs$Var1)),
    colnames(df)
  )
  effect_excluded <- if (length(cross_cols) > 0) {
    rowSums(df[, cross_cols, drop = FALSE] == 1) == 0
  } else {
    rep(TRUE, nrow(df))
  }

  df[[prefix]] <- as.integer(sig_included & effect_included & effect_excluded)
  df
}


# For a given cell type × treatment, identifies eQTLs that are treatment-
# responsive: significant (per sig_mode) AND effect size NOT shared with
# Control (sign flip or magnitude outside `factor`).
get_treatment_responsive_eqtls_signfixed <- function(
  celltype, treatment, df,
  lfsr_thresh = 0.05,
  factor = 1.5,
  sig_mode = c("treat_only", "treat_not_control", "treat_or_control"),
  FUN = identity
) {
  sig_mode  <- match.arg(sig_mode)
  lfsr_cols <- c(paste0(celltype, "_", treatment, "_lfsr"),
                 paste0(celltype, "_Control_lfsr"))
  beta_cols <- c(paste0(celltype, "_", treatment, "_beta"),
                 paste0(celltype, "_Control_beta"))

  if (!all(c(lfsr_cols, beta_cols) %in% colnames(df))) {
    warning(paste("Missing columns for", celltype, treatment))
    return(NULL)
  }

  lfsr_mat <- as.matrix(df[, lfsr_cols])
  beta_mat <- as.matrix(df[, beta_cols])
  rownames(lfsr_mat) <- rownames(beta_mat) <- df$gene_snp
  colnames(lfsr_mat) <- colnames(beta_mat) <- c(treatment, "Control")

  sign_flip <- sign(beta_mat[, treatment]) != sign(beta_mat[, "Control"]) &
               beta_mat[, treatment] != 0 & beta_mat[, "Control"] != 0
  ratio     <- abs(beta_mat[, treatment] / beta_mat[, "Control"])
  mag_diff  <- is.finite(ratio) & (ratio < (1 / factor) | ratio > factor)
  no_sharing <- mag_diff | sign_flip

  treat_only_idx <- switch(sig_mode,
    treat_only       = which(lfsr_mat[, treatment] < lfsr_thresh & no_sharing),
    treat_not_control = which(lfsr_mat[, treatment] < lfsr_thresh &
                               lfsr_mat[, "Control"] >= lfsr_thresh & no_sharing),
    treat_or_control  = which((lfsr_mat[, treatment] < lfsr_thresh |
                                lfsr_mat[, "Control"] < lfsr_thresh) & no_sharing)
  )

  if (length(treat_only_idx) == 0) return(NULL)

  out <- df[treat_only_idx, c("gene_snp", "gene")]
  out$celltype  <- celltype
  out$treatment <- treatment
  out$sig_mode  <- sig_mode
  out
}


# Runs the full sharing + treatment-responsiveness pipeline for all cell types
# at a given fold-change threshold. Writes:
#   - combined_celltype_sharing_output.tsv
#   - treatment_specific_eqtls.tsv
#   - treatment_responsive_eqtls_signfixed.tsv
run_sharing_pipeline <- function(
  combined_mash_results,
  lfc_factor,
  output_dir_base,
  sig_mode = c("treat_only", "treat_not_control", "treat_or_control")
) {
  sig_mode <- match.arg(sig_mode)

  lfc_dir <- file.path(output_dir_base, paste0("lfc", lfc_factor))
  dir.create(lfc_dir, recursive = TRUE, showWarnings = FALSE)

  treatments <- c("Ale", "Bri", "Cer", "Cri", "Lor")
  excluded   <- "Control"

  celltypes <- unique(gsub(
    "_(Ale|Bri|Cer|Cri|Lor|Control)_(lfsr|beta|sd)$", "",
    grep("_(lfsr|beta|sd)$", colnames(combined_mash_results), value = TRUE)
  ))

  celltype_sharing_list <- list()
  for (celltype in celltypes) {
    message("  Processing ", celltype)
    tryCatch({
      cols       <- grep(paste0("^", celltype, "_.*_(lfsr|beta|sd)$"),
                         colnames(combined_mash_results), value = TRUE)
      celltype_df <- combined_mash_results[, c(cols, "gene", "gene_snp")]

      sharing_info <- get_significance_and_effect_sharing_from_df(
        celltype, celltype_df, factor = lfc_factor)
      sharing_info <- add_pairwise_shared_columns(sharing_info, celltype)
      sharing_info <- add_treatment_specific_columns(sharing_info, celltype)

      # Ensure all expected pairwise columns exist (fill 0 if missing)
      for (pair in combn(treatments, 2, simplify = FALSE)) {
        col <- paste0(celltype, "_effectsim_", pair[1], "_", pair[2])
        if (!(col %in% colnames(sharing_info))) sharing_info[[col]] <- 0
      }

      sharing_info <- add_group_specific_sharing_column(
        sharing_info, celltype, treatments, excluded)

      celltype_sharing_list[[celltype]] <- left_join(
        sharing_info, celltype_df, by = "gene_snp")
    }, error = function(e) {
      message("  Error in ", celltype, ": ", e$message)
    })
  }

  combined_df <- purrr::reduce(celltype_sharing_list, full_join, by = c("gene_snp", "gene"))
  write_tsv(combined_df, file.path(lfc_dir, "combined_celltype_sharing_output.tsv"))

  # Treatment-specific eQTLs
  treatment_specific_eqtls <- combined_df %>%
    dplyr::select(gene_snp, gene, ends_with("_specific")) %>%
    pivot_longer(ends_with("_specific"), names_to = "celltype_treatment", values_to = "is_specific") %>%
    filter(is_specific == 1) %>%
    mutate(
      celltype  = sub("^(.*)_(.*)_specific$", "\\1", celltype_treatment),
      treatment = sub("^(.*)_(.*)_specific$", "\\2", celltype_treatment)
    ) %>%
    dplyr::select(gene_snp, gene, celltype, treatment) %>%
    arrange(celltype, treatment, gene_snp)
  write_tsv(treatment_specific_eqtls, file.path(lfc_dir, "treatment_specific_eqtls.tsv"))

  # Treatment-responsive eQTLs
  treatment_responsive_list <- list()
  for (celltype in celltypes) {
    for (treatment in treatments) {
      res <- get_treatment_responsive_eqtls_signfixed(
        celltype = celltype, treatment = treatment,
        df = combined_df, lfsr_thresh = 0.05,
        factor = lfc_factor, sig_mode = sig_mode
      )
      if (!is.null(res)) {
        treatment_responsive_list[[paste(celltype, treatment, sep = "_")]] <- res
      }
    }
  }
  write_tsv(
    bind_rows(treatment_responsive_list),
    file.path(lfc_dir, "treatment_responsive_eqtls_signfixed.tsv")
  )
}


# Summarizes how many response eQTLs are seen in 1, 2, ..., N cell types.
summarize_treatment_responsive_by_celltype <- function(mash_folder, lfc) {
  treatment_responsive <- read_tsv(
    file.path(mash_folder, paste0("lfc", lfc), "treatment_responsive_eqtls_signfixed.tsv"),
    show_col_types = FALSE)
  all_mash_eqtls <- read_tsv(
    file.path(mash_folder, "combined_mash_results.tsv"),
    show_col_types = FALSE)

  total_response_eqtls <- all_mash_eqtls %>%
    filter(gene_snp %in% treatment_responsive$gene_snp) %>%
    distinct(gene_snp) %>%
    nrow()

  num_celltypes_per_eqtl <- treatment_responsive %>%
    group_by(gene_snp) %>%
    summarise(n_celltypes = n_distinct(celltype), .groups = "drop") %>%
    arrange(desc(n_celltypes))

  count_summary <- num_celltypes_per_eqtl %>%
    count(n_celltypes, name = "n_eqtls") %>%
    arrange(n_celltypes)

  n_single    <- count_summary %>% filter(n_celltypes == 1) %>% pull(n_eqtls)
  max_ct      <- max(num_celltypes_per_eqtl$n_celltypes)
  n_all       <- count_summary %>% filter(n_celltypes == max_ct) %>% pull(n_eqtls)
  if (!length(n_single)) n_single <- 0
  if (!length(n_all))    n_all    <- 0

  message(
    "LFC ", lfc, ": ", n_single, " eQTLs (of ", total_response_eqtls,
    " unique response eQTLs) response in a single cell type; ",
    n_all, " response in all ", max_ct, " cell types."
  )

  list(summary_table = count_summary,
       eqtl_celltype_counts = num_celltypes_per_eqtl)
}



################## 2) TSS DISTANCE HELPERS ################
# Shared data prep for TSS distance plots. Loads treatment-responsive eQTLs
# and all MASH eQTLs, computes distance to TSS for response and non-response
# eQTLs in each cell type × treatment condition.

build_tss_distance_df <- function(
  mash_folder,
  lfc,
  gene_loc_file,
  celltypes,
  treatments  = c("Ale", "Bri", "Cer", "Cri", "Lor"),
  lfsr_cutoff = 0.05
) {
  lfc_folder <- paste0("lfc", lfc)
  treatment_responsive <- read_tsv(
    file.path(mash_folder, lfc_folder, "treatment_responsive_eqtls_signfixed.tsv"),
    show_col_types = FALSE)
  all_mash_eqtls <- read_tsv(
    file.path(mash_folder, "combined_mash_results.tsv"),
    show_col_types = FALSE)
  gene_loc <- read_csv(gene_loc_file, show_col_types = FALSE) %>%
    distinct(gene, .keep_all = TRUE) %>%
    dplyr::select(gene, chr, start, end) %>%
    mutate(chr = gsub("chr", "", chr))

  make_dist_df <- function(snp_df, gene_loc, ct, tx, label) {
    snp_df %>%
      mutate(
        chr = str_split(gene_snp, "[_:]", simplify = TRUE)[, 2],
        pos = as.numeric(str_split(gene_snp, "[_:]", simplify = TRUE)[, 3])
      ) %>%
      left_join(gene_loc, by = c("gene", "chr")) %>%
      filter(!is.na(start)) %>%
      mutate(distance = pos - start, celltype = ct, treatment = tx, source = label) %>%
      dplyr::select(gene_snp, distance, treatment, celltype, source)
  }

  dt_pos <- list()
  for (ct in celltypes) {
    for (tx in treatments) {
      response_snps <- treatment_responsive %>%
        filter(celltype == ct, treatment == tx) %>%
        pull(gene_snp)

      sig_lfsr_cols <- grep(
        paste0("^", ct, "_(", paste(c("Control", tx), collapse = "|"), ")_lfsr$"),
        colnames(all_mash_eqtls), value = TRUE)

      sig_in_cond <- all_mash_eqtls %>%
        filter(if_any(all_of(sig_lfsr_cols), ~ .x < lfsr_cutoff))

      dt_pos[[paste(ct, tx, sep = "_")]] <- bind_rows(
        make_dist_df(filter(sig_in_cond,  gene_snp %in% response_snps),
                     gene_loc, ct, tx, "Response"),
        make_dist_df(filter(sig_in_cond, !gene_snp %in% response_snps),
                     gene_loc, ct, tx, "Non-response")
      )
    }
  }

  bind_rows(dt_pos) %>%
    mutate(
      source      = factor(source, levels = c("Response", "Non-response")),
      facet_label = paste0(celltype, " — ", treatment)
    )
}


# Density plot of distance to TSS for response vs non-response eQTLs,
# faceted by cell type × treatment. Saves PDF to lfc plots folder.
plot_tss_density <- function(
  df_combined,
  lfc,
  mash_folder,
  ncol       = 5,
  maf        = 0.05,
  cisdist    = "100kb"
) {
  lfc_folder  <- paste0("lfc", lfc)
  n_facets    <- length(unique(df_combined$facet_label))
  nrow_plot   <- ceiling(n_facets / ncol)
  pdf_height  <- max(8,  nrow_plot * 3)
  pdf_width   <- max(12, ncol * 4)

  p <- ggplot(df_combined, aes(x = distance / 1000, fill = source)) +
    geom_histogram(
      aes(fill = source, alpha = source, y = after_stat(density)),
      position = "identity", alpha = 0.3, bins = 50
    ) +
    geom_density(
      aes(color = source, y = after_stat(density)),
      linewidth = 0.8, bounds = c(-100, 100), alpha = 0
    ) +
    scale_fill_manual(values  = c("Response" = "#f53224", "Non-response" = "#619CFF")) +
    scale_color_manual(values = c("Response" = "#f53224", "Non-response" = "#619CFF")) +
    scale_alpha_manual(values = c("Response" = 0.3, "Non-response" = 0.3)) +
    facet_wrap(~facet_label, ncol = ncol, scales = "free_y") +
    theme_classic(base_size = 14) +
    theme(legend.position = "bottom", strip.text = element_text(size = 10)) +
    labs(
      x     = "Distance to TSS (Kb)",
      y     = "Density",
      title = paste0("Distance to TSS: fold change ≥ ", lfc,
                     ", MAF ", maf, ", ", cisdist)
    )

  ggsave(
    file.path(mash_folder, lfc_folder, "plots/tss_distance_density.pdf"),
    p, width = pdf_width, height = pdf_height, limitsize = FALSE
  )

  invisible(p)
}

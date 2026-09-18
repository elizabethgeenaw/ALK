# Supp Fig S14: S14_reQTL_concordance_mashr_vs_interaction_medianUMI1.pdf

# Perform an independent genotype x treatment interaction model.
# We fit a MatrixEQTL modelLINEAR_CROSS on the joint
# treatment+control expression matrix (see 04_matrixeqtl_interaction.R), then
# for every response eQTL (locus x cell type x treatment) we pull the matching
# interaction beta and compare it to the mash-derived delta beta
# (beta_treatment - beta_control). If mash is picking up real treatment
# effects, delta beta and the interaction beta should agree in both magnitude
# and sign; if mash is picking up artifacts we would see poor correlation
# and near-random sign concordance.

# Two stats per drug are printed and shown as scatter panels in Fig S14:
# Pearson r between delta beta and interaction beta
# Sign concordance (% of response eQTLs where the two effects have the same sign; 50% would be chance)

# We also repeat both stats at the stricter lfc5 response-eQTL cutoff as a
# sensitivity check.

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(patchwork)


# Paths and constants
umi_label <- "medianUMI1"
mash_dir <- paste0("/project/gilad/geena/publication_clean/results/mash/", umi_label, "/lfc3/")
mash_root <- paste0("/project/gilad/geena/publication_clean/results/mash/", umi_label, "/")
interact_dir <- paste0("/project/gilad/geena/publication_clean/results/treatment/matrixeqtl/interaction/", umi_label, "/")
out_dir <- file.path(interact_dir, "concordance")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

drugs <- c("Ale", "Bri", "Cer", "Cri", "Lor")

# Load mash response eQTLs (lfc3 = manuscript main threshold).
# delta_beta = beta_treatment - beta_control is the mash effect-size shift we
# want the interaction model to agree with.
message("Loading mash response eQTLs (lfc3)...")
resp <- read_tsv(file.path(mash_dir, "supp_table_response_eqtls.tsv"),
                 show_col_types = FALSE) %>%
  mutate(delta_beta = beta_treat - beta_control,
         snp_id_norm = snp,
         celltype_clean = celltype)
cat("Response eQTLs (lfc3):", nrow(resp), "\n")
print(table(resp$treatment))

# Load the interaction-model results per drug. Each drug has one file per
# chromosome; we read all 22 chr files, tag them with drug and cell type
# (parsed from the filename), then bind into one long table.
load_interaction_results <- function(drug) {
  pattern <- paste0("^", drug, "_.*_chr[0-9]+_interact_cis\\.txt$")
  files <- list.files(interact_dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) {
    message("No interaction result files found for ", drug)
    return(NULL)
  }
  results <- lapply(files, function(f) {
    df <- tryCatch(read_tsv(f, show_col_types = FALSE),
                   error = function(e) { message("Could not read: ", f); NULL })
    if (is.null(df) || nrow(df) == 0) return(NULL)
    base <- basename(f)
    mid <- sub(paste0("^", drug, "_"), "", base)
    mid <- sub("_chr[0-9]+_interact_cis\\.txt$", "", mid)
    df$celltype_clean <- gsub(" ", "_", mid)
    df$drug <- drug
    df
  })
  bind_rows(Filter(Negate(is.null), results))
}

# Load all interaction results from raw files (this takes several minutes).
message("Loading interaction results from raw files (this may take several minutes)...")
interact_all <- bind_rows(lapply(drugs, load_interaction_results))
if (nrow(interact_all) == 0) {
  stop("No interaction results found. Run 04_matrixeqtl_interaction.R first.")
}
# Strip 'chr' prefix so SNP IDs match the mash response eQTL table
interact_all <- interact_all %>% mutate(snp_id_norm = sub("^chr", "", snps))
cat("Interaction results:", nrow(interact_all), "SNP-gene pairs\n")

# Join response eQTLs to interaction betas. Each mash response eQTL has
# exactly one interaction-model beta (matched on gene, SNP, drug, cell type).
message("Matching response eQTLs to interaction betas...")
concordance_df <- resp %>%
  filter(treatment %in% drugs) %>%
  left_join(
    interact_all %>% dplyr::select(drug, celltype_clean, snp_id_norm, gene,
                                   beta_gi = beta, pval_gi = pvalue),
    by = c("treatment" = "drug", "celltype_clean", "snp_id_norm", "gene")
  ) %>%
  filter(!is.na(beta_gi))
cat("Matched response eQTLs:", nrow(concordance_df),
    "of", nrow(resp %>% filter(treatment %in% drugs)), "\n")

rm(interact_all); gc()

# Compute concordance stats per drug: Pearson r, sign concordance, n.
# Response eQTLs where delta_beta and beta_gi have the same sign count as
# concordant (mash and interaction agree on direction).
compute_stats <- function(df, threshold_label) {
  df %>%
    group_by(treatment) %>%
    summarise(
      threshold = threshold_label,
      n = n(),
      pearson_r = cor(delta_beta, beta_gi, use = "complete.obs"),
      sign_concordance = mean(sign(delta_beta) == sign(beta_gi)),
      .groups = "drop"
    )
}

stats_lfc3 <- compute_stats(concordance_df, "lfc3")
cat("\nConcordance stats (lfc3 = manuscript main):\n")
print(stats_lfc3)

# Fig S14: five scatter panels (one per drug) of delta_beta vs beta_gi.
# Each panel shows the Pearson r and sign concordance as an annotation.
plots <- lapply(drugs, function(drug) {
  df <- concordance_df %>% filter(treatment == drug)
  st <- stats_lfc3 %>% filter(treatment == drug)
  label <- paste0("Pearson r = ", round(st$pearson_r, 2),
                  "\nSign concordance = ", round(st$sign_concordance * 100, 1), "%")
  ggplot(df, aes(x = delta_beta, y = beta_gi)) +
    geom_point(size = 0.8, alpha = 0.5) +
    geom_smooth(method = "lm", colour = "red", linetype = "dashed", se = TRUE) +
    annotate("text", x = -Inf, y = Inf, label = label,
             hjust = -0.05, vjust = 1.3, size = 3) +
    labs(title = drug,
         x = expression(Delta * beta ~ "(treat - ctrl, mashr)"),
         y = expression(beta[GI] ~ "(interaction model)")) +
    theme_bw(base_size = 11)
})
fig_s14 <- wrap_plots(plots, nrow = 1) +
  plot_annotation(title = "Effect-size concordance: mashr vs. interaction model (lfc3)")

# Sensitivity analysis: repeat at lfc5 threshold. Response eQTLs at lfc5 are a
# stricter subset of the lfc3 set (larger required effect-size shift). We just
# subset concordance_df on the lfc5 keys and re-compute the same stats.
message("Running lfc5 sensitivity analysis...")
resp5 <- read_tsv(
  file.path(mash_root, "lfc5/treatment_responsive_eqtls_signfixed.tsv"),
  show_col_types = FALSE
)
cat("Response eQTLs (lfc5):", nrow(resp5), "\n")

concordance_lfc5 <- resp5 %>%
  dplyr::select(gene_snp, celltype, treatment) %>%
  inner_join(
    concordance_df %>%
      dplyr::select(gene_snp, celltype_clean, treatment, delta_beta, beta_gi, pval_gi),
    by = c("gene_snp", "celltype" = "celltype_clean", "treatment")
  )
cat("Matched to interaction model (lfc5):", nrow(concordance_lfc5), "\n")

stats_lfc5 <- compute_stats(concordance_lfc5, "lfc5")
cat("\nConcordance stats (lfc5 = sensitivity):\n")
print(stats_lfc5)

# Save outputs: (i) the plot for Fig S14, (ii) a combined stats table with
# both thresholds for the manuscript.
message("Saving outputs...")
ggsave(file.path(out_dir, paste0("concordance_mashr_vs_interaction_", umi_label, ".pdf")),
       plot = fig_s14, width = 12, height = 3.5)

combined_stats <- bind_rows(stats_lfc3, stats_lfc5)
write_tsv(combined_stats,
          file.path(out_dir, paste0("concordance_stats_", umi_label, ".tsv")))

cat("\nDone. Files saved to:", out_dir, "\n")
cat("  concordance_mashr_vs_interaction_", umi_label, ".pdf\n", sep = "")
cat("  concordance_stats_", umi_label, ".tsv\n", sep = "")

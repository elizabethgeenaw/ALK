# Supp Fig S10: S10_std_eqtl_power_scatter.pdf
# Supp Fig S16: S16_response_eqtls_vs_power.pdf 
# Supp Table S20: S20_per_condition_power_and_eqtls.tsv 

# Standalone version of the power_summary chunks from post_mash_analysis_v2.Rmd.
# Loads per-condition counts directly from saved files so we don't need to
# re-run the upstream mash chunks.

library(tidyverse)
library(patchwork)
library(ggrepel)

results_dir <- "/project/gilad/geena/publication_clean/results/mash/medianUMI1"
lfc         <- 3
plot_dir    <- file.path(results_dir, paste0("lfc", lfc), "plots")
pseudobulk_dir <- "/project/gilad/geena/publication_clean/matrixeqtl_input/treatment/pseudobulk_15cells_20indiv_summary"
control_dir    <- "/project/gilad/geena/publication_clean/matrixeqtl_input/control/pseudobulk_15cells_1indiv_summary"


# Load per-condition standard / response eGenes
counts <- read.csv(file.path(plot_dir, "egene_and_response_egene_counts.csv"),
                   stringsAsFactors = FALSE)

control_std <- counts %>%
  filter(treatment == "Control") %>%
  select(celltype, ctrl_std = n_egenes)

power_table <- counts %>%
  filter(treatment != "Control") %>%
  rename(treat_std = n_egenes, response_egenes = n_response_egenes) %>%
  left_join(control_std, by = "celltype")

# Donors per treatment condition
indiv_tx <- read.delim(file.path(pseudobulk_dir, "celltype_individual_counts.tsv")) %>%
  rename(celltype = cellid_type) %>%
  mutate(celltype = gsub(" ", "_", celltype))

power_table <- power_table %>%
  left_join(indiv_tx %>% filter(treatment != "Control") %>%
              rename(treat_donors = num_individuals),
            by = c("celltype", "treatment"))

# Treatment cells per condition
treat_cells <- read.delim(file.path(pseudobulk_dir, "cellcount_filter_summary_combined.tsv")) %>%
  rename(celltype = cellid_type) %>%
  mutate(celltype = gsub(" ", "_", celltype)) %>%
  group_by(celltype, treatment) %>%
  summarise(treat_cells = sum(num_cells_final), .groups = "drop")

power_table <- power_table %>%
  left_join(treat_cells, by = c("celltype", "treatment"))

# Control donors / cells (per cell type, shared across drugs)
ctrl_summary <- read.delim(file.path(control_dir, "cellcount_filter_summary_combined.tsv")) %>%
  rename(celltype = cellid_type) %>%
  mutate(celltype = gsub(" ", "_", celltype)) %>%
  filter(!dropped) %>%
  group_by(celltype) %>%
  summarise(ctrl_donors = n_distinct(individual),
            ctrl_cells  = sum(num_cells_final),
            .groups = "drop")

power_table <- power_table %>%
  left_join(ctrl_summary, by = "celltype")

cat("Power table dimensions: ", nrow(power_table), " rows x ",
    ncol(power_table), " columns\n", sep = "")
cat("Columns: ", paste(colnames(power_table), collapse = ", "), "\n")
cat("Missing values per column:\n")
print(colSums(is.na(power_table)))

# Forest plot + labeled scatter

predictors  <- c("ctrl_donors", "ctrl_std", "ctrl_cells",
                 "treat_cells", "treat_std", "treat_donors")
pred_labels <- c(ctrl_donors  = "Control donors",
                 ctrl_std     = "Control standard eQTLs",
                 ctrl_cells   = "Control cells",
                 treat_cells  = "Treatment cells",
                 treat_std    = "Treatment standard eQTLs",
                 treat_donors = "Treatment donors")

boot_r <- function(x, y, n = 2000, seed = 42) {
  ok <- complete.cases(x, y)
  x  <- x[ok]; y <- y[ok]
  set.seed(seed)
  rs <- replicate(n, {
    i <- sample.int(length(x), replace = TRUE)
    suppressWarnings(cor(x[i], y[i]))
  })
  data.frame(r  = cor(x, y),
             lo = unname(quantile(rs, 0.025, na.rm = TRUE)),
             hi = unname(quantile(rs, 0.975, na.rm = TRUE)))
}

cor_summary <- bind_rows(lapply(predictors, function(p) {
  res <- boot_r(power_table[[p]], power_table$response_egenes)
  res$predictor <- p
  res
})) %>%
  mutate(label = pred_labels[predictor],
         side  = ifelse(grepl("^ctrl", predictor), "Control", "Treatment")) %>%
  arrange(desc(abs(r))) %>%
  mutate(label = factor(label, levels = rev(label)))

cat("\nCorrelations (with 95% bootstrap CI):\n")
print(cor_summary %>% select(label, r, lo, hi))

p_forest <- ggplot(cor_summary, aes(x = r, y = label, color = side)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.7) +
  scale_color_manual(values = c("Control" = "#C0392B", "Treatment" = "#2C7BB6")) +
  coord_cartesian(xlim = c(-0.5, 0.8)) +
  theme_classic(base_size = 13) +
  labs(x = "Pearson r (95% bootstrap CI)",
       y = NULL, color = NULL,
       title = "A. Power proxies vs response eQTL count") +
  theme(legend.position = "top")

scatter_data <- power_table %>%
  group_by(celltype) %>%
  summarise(response_egenes_mean = mean(response_egenes),
            ctrl_std = first(ctrl_std), .groups = "drop") %>%
  mutate(label = gsub("_", " ", celltype))

r_scatter <- cor(scatter_data$ctrl_std, scatter_data$response_egenes_mean)

p_scatter <- ggplot(scatter_data,
                    aes(x = ctrl_std, y = response_egenes_mean)) +
  geom_smooth(method = "lm", se = TRUE, color = "grey40",
              fill = "grey85", linewidth = 0.6) +
  geom_point(size = 2.5, alpha = 0.8, color = "#2C3E50") +
  ggrepel::geom_text_repel(aes(label = label),
                            size = 2.8, max.overlaps = Inf,
                            box.padding = 0.4, segment.color = "grey60") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.8,
           label = sprintf("r = %.2f", r_scatter),
           size = 4, color = "grey20") +
  theme_classic(base_size = 13) +
  labs(x = "Control standard eQTLs",
       y = "Mean response eQTLs per cell type (across 5 drugs)",
       title = "B. Per cell type: response vs control standard eQTL count")

p_power_summary <- p_forest + p_scatter + plot_layout(widths = c(1, 1.3))

ggsave(file.path(plot_dir, "power_check_summary.pdf"),
       p_power_summary, width = 13, height = 5.5)

# ---------------------------------------------------------------------------
# Option 2: variance partitioning + residuals heatmap

marginal_r2 <- sapply(predictors, function(p) {
  cor(power_table[[p]], power_table$response_egenes, use = "complete.obs")^2
})

fmla_full <- as.formula(paste("response_egenes ~", paste(predictors, collapse = " + ")))
m_full    <- lm(fmla_full, data = power_table)
r2_full   <- summary(m_full)$r.squared

unique_r2 <- sapply(predictors, function(p) {
  others <- setdiff(predictors, p)
  fmla   <- as.formula(paste("response_egenes ~", paste(others, collapse = " + ")))
  r2_full - summary(lm(fmla, data = power_table))$r.squared
})

cat(sprintf("\nFull-model R^2 = %.3f\n", r2_full))
cat("Marginal R^2 per predictor:\n"); print(round(marginal_r2, 3))
cat("Unique R^2 per predictor:\n");   print(round(unique_r2, 3))

varexp_long <- tibble(
  predictor = predictors,
  Marginal  = marginal_r2,
  Unique    = unique_r2
) %>%
  mutate(label = pred_labels[predictor],
         side  = ifelse(grepl("^ctrl", predictor), "Control", "Treatment")) %>%
  pivot_longer(c(Marginal, Unique), names_to = "type", values_to = "r2") %>%
  mutate(type = factor(type, levels = c("Marginal", "Unique"),
                       labels = c("Marginal R² (alone)",
                                  "Unique R² (added to others)")))

ord <- varexp_long %>%
  filter(type == "Marginal R² (alone)") %>%
  arrange(desc(r2)) %>% pull(label)
varexp_long$label <- factor(varexp_long$label, levels = rev(ord))

p_varexp <- ggplot(varexp_long, aes(x = r2, y = label, fill = side)) +
  geom_col() +
  scale_fill_manual(values = c("Control" = "#C0392B", "Treatment" = "#2C7BB6")) +
  facet_wrap(~ type, scales = "free_x") +
  theme_classic(base_size = 13) +
  labs(x = expression(R^2), y = NULL, fill = NULL,
       title = sprintf("A. Variance explained (full model R² = %.2f)",
                       r2_full)) +
  theme(legend.position = "top")

power_table_resid <- power_table %>%
  mutate(predicted = predict(m_full, newdata = .),
         residual  = response_egenes - predicted)

ct_order_resid <- power_table_resid %>%
  group_by(celltype) %>%
  summarise(mean_resid = mean(residual), .groups = "drop") %>%
  arrange(mean_resid) %>%
  pull(celltype)

residual_heat <- power_table_resid %>%
  mutate(celltype  = factor(gsub("_", " ", celltype),
                            levels = gsub("_", " ", ct_order_resid)),
         treatment = factor(treatment,
                            levels = c("Ale", "Bri", "Cer", "Cri", "Lor")))

resid_lim <- max(abs(residual_heat$residual))

p_resid <- ggplot(residual_heat,
                  aes(x = treatment, y = celltype, fill = residual)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#C0392B",
                        midpoint = 0, limits = c(-resid_lim, resid_lim),
                        name = "Residual\n(obs - pred)") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 10)) +
  labs(x = "Treatment", y = NULL,
       title = "B. Response eQTLs beyond joint power model")

p_varexp_summary <- p_varexp + p_resid + plot_layout(widths = c(1, 1.1))

ggsave(file.path(plot_dir, "power_variance_partition.pdf"),
       p_varexp_summary, width = 14, height = 6)


# Supplementary table: per-condition power + eQTL counts
supp_per_cond <- power_table %>%
  dplyr::select(
    celltype, treatment,
    treat_donors, treat_cells, treat_std_eqtls = treat_std,
    ctrl_donors,  ctrl_cells,  ctrl_std_eqtls  = ctrl_std,
    response_eqtls = response_egenes
  ) %>%
  arrange(celltype, treatment) %>%
  mutate(celltype = gsub("_", " ", celltype))

supp_table_path <- file.path(results_dir, paste0("lfc", lfc),
                              "supp_table_per_condition_power_and_eqtls.tsv")
write_tsv(supp_per_cond, supp_table_path)

cat("\nSaved:\n")
cat("  ", file.path(plot_dir, "power_check_summary.pdf"), "\n")
cat("  ", file.path(plot_dir, "power_variance_partition.pdf"), "\n")
cat("  ", supp_table_path, "\n")


# Per-condition standard eQTL range and correlations with donor count and cell count 

# Min / max / median across all 126 conditions (105 treatment + 21 control)
ctrl_unique <- power_table %>%
  distinct(celltype, ctrl_std, ctrl_donors, ctrl_cells) %>%
  transmute(celltype, treatment = "Control", n_eqtls = ctrl_std)

all_conds_tbl <- bind_rows(
  power_table %>% transmute(celltype, treatment, n_eqtls = treat_std),
  ctrl_unique
)

cat("\nPer-condition eQTL range (all 126 conditions)\n")
cat(sprintf("Min   : %d — %s / %s\n",
  min(all_conds_tbl$n_eqtls),
  all_conds_tbl$celltype[which.min(all_conds_tbl$n_eqtls)],
  all_conds_tbl$treatment[which.min(all_conds_tbl$n_eqtls)]))
cat(sprintf("Max   : %d — %s / %s\n",
  max(all_conds_tbl$n_eqtls),
  all_conds_tbl$celltype[which.max(all_conds_tbl$n_eqtls)],
  all_conds_tbl$treatment[which.max(all_conds_tbl$n_eqtls)]))
cat(sprintf("Median: %.1f\n", median(all_conds_tbl$n_eqtls)))

# Correlations: treatment standard eQTLs vs donor count and cell count
cat("\nCorrelations: treat standard eQTLs vs donor/cell count\n")
cat(sprintf("Pearson  r (std eQTLs vs donors): %.3f\n",
  cor(power_table$treat_std, power_table$treat_donors)))
cat(sprintf("Pearson  r (std eQTLs vs cells) : %.3f\n",
  cor(power_table$treat_std, power_table$treat_cells)))
cat(sprintf("Spearman r (std eQTLs vs donors): %.3f\n",
  cor(power_table$treat_std, power_table$treat_donors, method = "spearman")))
cat(sprintf("Spearman r (std eQTLs vs cells) : %.3f\n",
  cor(power_table$treat_std, power_table$treat_cells,  method = "spearman")))


# Unique eGenes per condition correlated with donor and cell
# Uses UNIQUE GENES (not gene–SNP pairs) per (celltype, treatment) at lfsr<0.05.

drug_colors_full <- c(Ale = "#E69F00", Bri = "#56B4E9", Cer = "#009E73",
                      Cri = "#CC79A7", Lor = "#0072B2", Control = "#7F7F7F")

sig_all_egene <- read_tsv(
  file.path(results_dir, "supp_table_all_significant_eqtls.tsv"),
  show_col_types = FALSE) %>%
  mutate(celltype = gsub(" ", "_", celltype))

egene_per_cond <- sig_all_egene %>%
  group_by(celltype, treatment) %>%
  summarise(n_egenes = n_distinct(gene), .groups = "drop")

donors_all <- bind_rows(
  indiv_tx %>% filter(treatment != "Control") %>%
    transmute(celltype, treatment, n_donors = num_individuals),
  ctrl_summary %>% transmute(celltype, treatment = "Control", n_donors = ctrl_donors)
)

cells_all <- bind_rows(
  treat_cells %>% transmute(celltype, treatment, n_cells = treat_cells),
  ctrl_summary %>% transmute(celltype, treatment = "Control", n_cells = ctrl_cells)
)

egene_power_all <- egene_per_cond %>%
  left_join(donors_all, by = c("celltype", "treatment")) %>%
  left_join(cells_all,  by = c("celltype", "treatment"))

egene_power_treat <- egene_power_all %>% filter(treatment != "Control")

report_r <- function(tbl, scope_label) {
  r_don  <- cor(tbl$n_egenes, tbl$n_donors, use = "complete.obs")
  r_cell <- cor(tbl$n_egenes, tbl$n_cells,  use = "complete.obs")
  cat(sprintf("\eGenes vs power — %s (n=%d)\n", scope_label, nrow(tbl)))
  cat(sprintf("Pearson r (n_egenes vs n_donors): %.3f\n", r_don))
  cat(sprintf("Pearson r (n_egenes vs n_cells) : %.3f\n", r_cell))
  invisible(list(r_don = r_don, r_cell = r_cell))
}

cors_all   <- report_r(egene_power_all,   "ALL 126 conditions")
cors_treat <- report_r(egene_power_treat, "TREATMENT-ONLY 105 (matches manuscript)")

# Persist tidy table with both scopes flagged
egene_power_out <- egene_power_all %>%
  mutate(celltype = gsub("_", " ", celltype),
         in_treat_only_scope = treatment != "Control") %>%
  arrange(celltype, treatment)
write_tsv(egene_power_out,
          file.path(results_dir, paste0("lfc", lfc),
                    "supp_table_egenes_vs_power.tsv"))

# Scatter panels — one figure per scope
make_egene_scatter <- function(tbl, xvar, xlabel, r_val, log_x = FALSE,
                                panel_title) {
  p <- ggplot(tbl, aes(x = .data[[xvar]], y = n_egenes, color = treatment)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                fill = "grey85", linewidth = 0.5, show.legend = FALSE) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = drug_colors_full, name = "Condition") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = sprintf("r = %.2f", r_val), size = 4, color = "grey20") +
    theme_classic(base_size = 12) +
    labs(x = xlabel, y = "Number of significant eGenes", title = panel_title)
  if (log_x) p <- p + scale_x_log10()
  p
}

pd1 <- make_egene_scatter(egene_power_all, "n_donors", "Number of individuals",
                          cors_all$r_don, log_x = FALSE,
                          "A. All 126 conditions: eGenes vs donors")
pd2 <- make_egene_scatter(egene_power_all, "n_cells", "Cell count (log scale)",
                          cors_all$r_cell, log_x = TRUE,
                          "B. All 126 conditions: eGenes vs cells")
pt1 <- make_egene_scatter(egene_power_treat, "n_donors", "Number of individuals",
                          cors_treat$r_don, log_x = FALSE,
                          "C. Treatment only (n=105): eGenes vs donors")
pt2 <- make_egene_scatter(egene_power_treat, "n_cells", "Cell count (log scale)",
                          cors_treat$r_cell, log_x = TRUE,
                          "D. Treatment only (n=105): eGenes vs cells")

p_egene_power <- (pd1 + pd2) / (pt1 + pt2) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

cat("Saved: ", file.path(plot_dir, "egenes_vs_power.pdf"), "\n", sep = "")
cat("Saved: ",
    file.path(results_dir, paste0("lfc", lfc),
              "supp_table_egenes_vs_power.tsv"),
    "\n", sep = "")


# Per-cell-type: unique eGenes significant in ANY condition (lfsr < 0.05)
mash_results_path <- file.path(results_dir, "combined_mash_results.tsv")
mash_mat <- read_tsv(mash_results_path, show_col_types = FALSE)
lfsr_cols     <- grep("_lfsr$", colnames(mash_mat), value = TRUE)
celltypes_mash <- unique(sub("_[A-Za-z]+_lfsr$", "", lfsr_cols))

per_ct_egenes <- map_int(celltypes_mash, function(ct) {
  cols <- grep(paste0("^", ct, "_"), lfsr_cols, value = TRUE)
  mat  <- mash_mat[, cols, drop = FALSE]
  sum(apply(mat, 1, function(x) any(x < 0.05, na.rm = TRUE)))
})

per_ct_df <- tibble(
  celltype        = gsub("_", " ", celltypes_mash),
  n_unique_egenes = per_ct_egenes
) %>% arrange(n_unique_egenes)

cat("\nPer cell type: unique eGenes (lfsr<0.05 in any condition)\n")
print(per_ct_df, n = Inf)
cat(sprintf("Min   : %d — %s\n", min(per_ct_df$n_unique_egenes),
  per_ct_df$celltype[which.min(per_ct_df$n_unique_egenes)]))
cat(sprintf("Max   : %d — %s\n", max(per_ct_df$n_unique_egenes),
  per_ct_df$celltype[which.max(per_ct_df$n_unique_egenes)]))
cat(sprintf("Median: %.1f\n", median(per_ct_df$n_unique_egenes)))


# Unique response eQTLs per cell type (dedup across treatments) and cell-type rank by treated cell count
resp_path <- file.path(results_dir, paste0("lfc", lfc),
                        "treatment_responsive_eqtls_signfixed.tsv")
resp      <- read_tsv(resp_path, show_col_types = FALSE)

per_ct_resp <- resp %>%
  group_by(celltype) %>%
  summarise(n_unique_response_eqtls = n_distinct(gene_snp), .groups = "drop") %>%
  mutate(celltype = gsub("_", " ", celltype)) %>%
  arrange(desc(n_unique_response_eqtls))

cat("\nUnique response eQTLs per cell type (dedup across treatments)\n")
print(per_ct_resp, n = Inf)
cat(sprintf("Total unique response eQTLs: %d\n", n_distinct(resp$gene_snp)))

# Cell-type rank by mean treated cell count (for power argument)
ct_rank <- power_table %>%
  group_by(celltype) %>%
  summarise(mean_treat_cells = mean(treat_cells), .groups = "drop") %>%
  arrange(desc(mean_treat_cells)) %>%
  mutate(rank = row_number(),
         celltype = gsub("_", " ", celltype))

cat("\nCell type ranks by mean treated cell count\n")
print(ct_rank, n = Inf)

# Control- vs treatment-driven response eQTLs per cell type
# Classifies each unique gene_snp x celltype pair as:
#   control_only  — lfsr<0.05 in control only
#   treatment_only — lfsr<0.05 in treatment only
#   both           — lfsr<0.05 in both
sig_all <- read_tsv(
  file.path(results_dir, "supp_table_all_significant_eqtls.tsv"),
  show_col_types = FALSE)

resp_all <- read_tsv(resp_path, show_col_types = FALSE)

resp_pairs <- resp_all %>% distinct(gene_snp, celltype)

sig_ctrl_pairs  <- sig_all %>% filter(treatment == "Control") %>%
  select(gene_snp, celltype) %>% mutate(sig_ctrl = TRUE)
sig_treat_pairs <- sig_all %>% filter(treatment != "Control") %>%
  distinct(gene_snp, celltype) %>% mutate(sig_treat = TRUE)

resp_classified <- resp_pairs %>%
  left_join(sig_ctrl_pairs,  by = c("gene_snp","celltype")) %>%
  left_join(sig_treat_pairs, by = c("gene_snp","celltype")) %>%
  mutate(
    sig_ctrl  = replace_na(sig_ctrl,  FALSE),
    sig_treat = replace_na(sig_treat, FALSE),
    sig_class = case_when(
      sig_ctrl & sig_treat  ~ "both",
      sig_ctrl & !sig_treat ~ "control_only",
      !sig_ctrl & sig_treat ~ "treatment_only",
      TRUE                  ~ "neither"
    )
  )

tt <- table(resp_classified$sig_class)
cat("\nControl- vs treatment-driven response eQTLs (overall)\n")
print(tt)
cat(sprintf("Control-only : %.1f%%\n", 100*tt["control_only"] / sum(tt)))
cat(sprintf("Treat-only   : %.1f%%\n", 100*tt["treatment_only"] / sum(tt)))
cat(sprintf("Both         : %.1f%%\n", 100*tt["both"] / sum(tt)))

per_ct_class <- resp_classified %>%
  group_by(celltype) %>%
  summarise(
    n_total        = n(),
    n_ctrl_only    = sum(sig_class == "control_only"),
    n_treat_only   = sum(sig_class == "treatment_only"),
    n_both         = sum(sig_class == "both"),
    pct_ctrl_only  = round(100 * n_ctrl_only  / n_total, 1),
    pct_treat_only = round(100 * n_treat_only / n_total, 1),
    pct_both       = round(100 * n_both       / n_total, 1),
    .groups = "drop"
  ) %>%
  mutate(celltype = gsub("_", " ", celltype)) %>%
  arrange(desc(n_total))

cat("\nPer cell type: control- vs treatment-driven response eQTLs\n")
print(per_ct_class, n = Inf)


# Standard eQTL power analysis: correlations and scatter plots
# Treatment conditions (105 rows: 21 cell types x 5 drugs)
# Control conditions  (21 rows:  one per cell type)

drug_colors <- c(Ale = "#E69F00", Bri = "#56B4E9", Cer = "#009E73",
                 Cri = "#CC79A7", Lor = "#0072B2")

ctrl_tbl <- power_table %>%
  distinct(celltype, ctrl_std, ctrl_donors, ctrl_cells) %>%
  mutate(ct_label = gsub("_", " ", celltype))

# Bootstrap correlations
std_cors <- bind_rows(
  boot_r(power_table$treat_donors, power_table$treat_std) %>%
    mutate(predictor = "treat_donors", label = "Treatment donors", side = "Treatment"),
  boot_r(power_table$treat_cells,  power_table$treat_std) %>%
    mutate(predictor = "treat_cells",  label = "Treatment cells",  side = "Treatment"),
  boot_r(ctrl_tbl$ctrl_donors,     ctrl_tbl$ctrl_std) %>%
    mutate(predictor = "ctrl_donors",  label = "Control donors",   side = "Control"),
  boot_r(ctrl_tbl$ctrl_cells,      ctrl_tbl$ctrl_std) %>%
    mutate(predictor = "ctrl_cells",   label = "Control cells",    side = "Control")
) %>%
  arrange(desc(abs(r))) %>%
  mutate(label = factor(label, levels = rev(label)))

cat("\Standard eQTL power correlations (Pearson r, 95% bootstrap CI)\n")
print(std_cors %>% select(label, r, lo, hi))

# Forest plot
p_forest_std <- ggplot(std_cors, aes(x = r, y = label, color = side)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.8) +
  scale_color_manual(values = c(Control = "#C0392B", Treatment = "#2C7BB6")) +
  coord_cartesian(xlim = c(0, 1)) +
  theme_classic(base_size = 13) +
  labs(x = "Pearson r (95% bootstrap CI)", y = NULL, color = NULL,
       title = "Power proxies vs standard eQTL count") +
  theme(legend.position = "top")

# Scatter panels: treatment (colored by drug) and control (labeled by cell type)
make_scatter_treat <- function(xvar, xlabel, log_x = FALSE) {
  r <- cor(power_table[[xvar]], power_table$treat_std, use = "complete.obs")
  p <- ggplot(power_table, aes(x = .data[[xvar]], y = treat_std, color = treatment)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                fill = "grey85", linewidth = 0.5, show.legend = FALSE) +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_color_manual(values = drug_colors, name = "Drug") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = sprintf("r = %.2f", r), size = 4, color = "grey20") +
    theme_classic(base_size = 12) +
    labs(x = xlabel, y = "Standard eQTLs")
  if (log_x) p <- p + scale_x_log10()
  p
}

make_scatter_ctrl <- function(xvar, xlabel, log_x = FALSE) {
  r <- cor(ctrl_tbl[[xvar]], ctrl_tbl$ctrl_std, use = "complete.obs")
  p <- ggplot(ctrl_tbl, aes(x = .data[[xvar]], y = ctrl_std)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                fill = "grey85", linewidth = 0.5) +
    geom_point(size = 3, alpha = 0.8, color = "#C0392B") +
    geom_text_repel(aes(label = ct_label), size = 2.5, max.overlaps = Inf,
                    box.padding = 0.3, segment.color = "grey60") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = sprintf("r = %.2f", r), size = 4, color = "grey20") +
    theme_classic(base_size = 12) +
    labs(x = xlabel, y = "Standard eQTLs")
  if (log_x) p <- p + scale_x_log10()
  p
}

pt1 <- make_scatter_treat("treat_donors", "Treatment donors") +
  labs(title = "A. Treatment: donors vs standard eQTLs")
pt2 <- make_scatter_treat("treat_cells",  "Treatment cells (log scale)", log_x = TRUE) +
  labs(title = "B. Treatment: cells vs standard eQTLs")
pc1 <- make_scatter_ctrl("ctrl_donors", "Control donors") +
  labs(title = "C. Control: donors vs standard eQTLs")
pc2 <- make_scatter_ctrl("ctrl_cells",  "Control cells (log scale)", log_x = TRUE) +
  labs(title = "D. Control: cells vs standard eQTLs")

p_std_power <- (pt1 + pt2) / (pc1 + pc2) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

ggsave(file.path(plot_dir, "std_eqtl_power_scatter.pdf"),
       p_std_power, width = 13, height = 10)

cat("Saved: std_eqtl_power_scatter.pdf\n")
cat("Saved: std_eqtl_power_forest.pdf\n")


# Variance explained by power for standard eQTLs (treatment and control)
cat("\Standard eQTL variance explained — treatment (n=105)\n")
m_std_don  <- lm(treat_std ~ treat_donors,               data = power_table)
m_std_cel  <- lm(treat_std ~ treat_cells,                data = power_table)
m_std_full <- lm(treat_std ~ treat_donors + treat_cells, data = power_table)
cat(sprintf("R² donors only     : %.3f\n", summary(m_std_don)$r.squared))
cat(sprintf("R² cells only      : %.3f\n", summary(m_std_cel)$r.squared))
cat(sprintf("R² donors + cells  : %.3f\n", summary(m_std_full)$r.squared))
cat(sprintf("Unique R² donors   : %.3f\n",
  summary(m_std_full)$r.squared - summary(m_std_cel)$r.squared))
cat(sprintf("Unique R² cells    : %.3f\n",
  summary(m_std_full)$r.squared - summary(m_std_don)$r.squared))

cat("\Standard eQTL variance explained — control (n=21)\n")
m_ctrl_don  <- lm(ctrl_std ~ ctrl_donors,               data = ctrl_tbl)
m_ctrl_cel  <- lm(ctrl_std ~ ctrl_cells,                data = ctrl_tbl)
m_ctrl_full <- lm(ctrl_std ~ ctrl_donors + ctrl_cells,  data = ctrl_tbl)
cat(sprintf("R² donors only     : %.3f\n", summary(m_ctrl_don)$r.squared))
cat(sprintf("R² cells only      : %.3f\n", summary(m_ctrl_cel)$r.squared))
cat(sprintf("R² donors + cells  : %.3f\n", summary(m_ctrl_full)$r.squared))

cat("\Standard eQTL yield: per-drug summary\n")
power_table %>% group_by(treatment) %>%
  summarise(mean = round(mean(treat_std)), median = round(median(treat_std)),
            min = min(treat_std), max = max(treat_std), .groups = "drop") %>%
  arrange(desc(median)) %>% print()

cat("\Per-drug correlation: treat_std vs treat_cells and treat_donors\n")
for (d in c("Ale","Bri","Cer","Cri","Lor")) {
  sub <- filter(power_table, treatment == d)
  cat(sprintf("  %s: r(donors) = %.2f, r(cells) = %.2f  (n = %d)\n", d,
    cor(sub$treat_donors, sub$treat_std),
    cor(sub$treat_cells,  sub$treat_std), nrow(sub)))
}

cat("\Residuals from treatment standard eQTL power model (per cell type)\n")
power_table$std_resid <- residuals(m_std_full)
power_table %>% group_by(celltype) %>%
  summarise(mean_resid = round(mean(std_resid)),
            mean_std   = round(mean(treat_std)), .groups = "drop") %>%
  mutate(celltype = gsub("_", " ", celltype)) %>%
  arrange(desc(mean_resid)) %>% print(n = Inf)

# ---------------------------------------------------------------------------
# Variance explained by power for response eQTLs
cat("\Response eQTL variance explained\n")
m_resp_don   <- lm(response_egenes ~ treat_donors,                               data = power_table)
m_resp_cel   <- lm(response_egenes ~ treat_cells,                                data = power_table)
m_resp_tx    <- lm(response_egenes ~ treat_donors + treat_cells,                 data = power_table)
m_resp_all   <- lm(response_egenes ~ treat_donors + treat_cells +
                     ctrl_donors  + ctrl_cells,                                   data = power_table)
m_resp_stds  <- lm(response_egenes ~ treat_std + ctrl_std,                       data = power_table)
cat(sprintf("R² treat donors only        : %.3f\n", summary(m_resp_don)$r.squared))
cat(sprintf("R² treat cells only         : %.3f\n", summary(m_resp_cel)$r.squared))
cat(sprintf("R² treat donors + cells     : %.3f\n", summary(m_resp_tx)$r.squared))
cat(sprintf("R² all 4 power terms        : %.3f\n", summary(m_resp_all)$r.squared))
cat(sprintf("R² treat_std + ctrl_std     : %.3f\n", summary(m_resp_stds)$r.squared))

cat("\Response eQTL correlations with all power metrics\n")
preds  <- c("treat_donors","treat_cells","treat_std","ctrl_donors","ctrl_cells","ctrl_std")
labs   <- c(treat_donors="Treat donors", treat_cells="Treat cells", treat_std="Treat std eQTLs",
            ctrl_donors ="Ctrl donors",  ctrl_cells ="Ctrl cells",  ctrl_std ="Ctrl std eQTLs")
resp_cors2 <- bind_rows(lapply(preds, function(p) {
  boot_r(power_table[[p]], power_table$response_egenes) %>%
    mutate(predictor = p, label = labs[p])
}))
print(resp_cors2 %>% select(label, r, lo, hi) %>% arrange(desc(abs(r))))

cat("\Response eQTL residuals (per cell type, treatment power model)\n")
power_table$resp_resid <- residuals(m_resp_tx)
power_table %>% group_by(celltype) %>%
  summarise(mean_resid = round(mean(resp_resid)),
            mean_resp  = round(mean(response_egenes)), .groups = "drop") %>%
  mutate(celltype = gsub("_", " ", celltype)) %>%
  arrange(desc(mean_resid)) %>% print(n = Inf)

cat("\Per-drug response eQTL summary\n")
power_table %>% group_by(treatment) %>%
  summarise(median = round(median(response_egenes)), mean = round(mean(response_egenes)),
            min = min(response_egenes), max = max(response_egenes), .groups = "drop") %>%
  arrange(desc(median)) %>% print()

# Response eQTL count per condition correlated with control
# dataset size. Supports the paragraph:
#   "We observed a modest correlation between response eQTL counts and control
#    dataset size (Pearson r = 0.49 for number of control donors, r = 0.25 for
#    total control cells)"
#
# Uses UNIQUE gene_snp (response eQTLs, not genes) per (celltype, treatment).
# Scope: 105 treatment conditions (21 CT × 5 drugs). Control donors/cells are
# per cell type (shared across the 5 drugs).

resp_by_cond <- resp %>%
  group_by(celltype, treatment) %>%
  summarise(n_response_eqtls = n_distinct(gene_snp), .groups = "drop")

# ctrl_summary already loaded above with (celltype, ctrl_donors, ctrl_cells)
reqtl_power_cond <- resp_by_cond %>%
  left_join(ctrl_summary, by = "celltype")

r_ctrl_don_reqtl  <- cor(reqtl_power_cond$n_response_eqtls,
                          reqtl_power_cond$ctrl_donors, use = "complete.obs")
r_ctrl_cell_reqtl <- cor(reqtl_power_cond$n_response_eqtls,
                          reqtl_power_cond$ctrl_cells,  use = "complete.obs")

# Complementary treatment-side correlations
r_tx_don_reqtl <- cor(power_table$treat_donors,
  power_table %>% inner_join(resp_by_cond, by = c("celltype","treatment")) %>%
    pull(n_response_eqtls), use = "complete.obs")
r_tx_cell_reqtl <- cor(power_table$treat_cells,
  power_table %>% inner_join(resp_by_cond, by = c("celltype","treatment")) %>%
    pull(n_response_eqtls), use = "complete.obs")

cat("\Response eQTL vs power (per condition, n = 105)\n")
cat(sprintf("Pearson r (n_response_eqtls vs ctrl_donors): %.3f\n", r_ctrl_don_reqtl))
cat(sprintf("Pearson r (n_response_eqtls vs ctrl_cells) : %.3f\n", r_ctrl_cell_reqtl))
cat(sprintf("Pearson r (n_response_eqtls vs treat_donors): %.3f\n", r_tx_don_reqtl))
cat(sprintf("Pearson r (n_response_eqtls vs treat_cells) : %.3f\n", r_tx_cell_reqtl))

# Supp table: cis-reQTL donor/cell counts per (celltype, treatment)
reqtl_power_out <- reqtl_power_cond %>%
  left_join(power_table %>% select(celltype, treatment, treat_donors, treat_cells),
            by = c("celltype", "treatment")) %>%
  select(celltype, treatment,
         n_response_eqtls,
         treat_donors, treat_cells,
         ctrl_donors,  ctrl_cells) %>%
  mutate(celltype = gsub("_", " ", celltype)) %>%
  arrange(celltype, treatment)

reqtl_supp_path <- file.path(results_dir, paste0("lfc", lfc),
                              "supp_table_cis_reQTL_donor_cell_count.tsv")
write_tsv(reqtl_power_out, reqtl_supp_path)

# Scatter figure — mirrors egenes_vs_power.pdf: control side on top,
# treatment side on bottom
make_reqtl_scatter <- function(tbl, xvar, xlabel, r_val, log_x = FALSE,
                                panel_title) {
  p <- ggplot(tbl, aes(x = .data[[xvar]], y = n_response_eqtls,
                       color = treatment)) +
    geom_smooth(method = "lm", se = TRUE, color = "grey40",
                fill = "grey85", linewidth = 0.5, show.legend = FALSE) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_color_manual(values = drug_colors, name = "Drug") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = sprintf("r = %.2f", r_val), size = 4, color = "grey20") +
    theme_classic(base_size = 12) +
    labs(x = xlabel, y = "Number of response eQTLs", title = panel_title)
  if (log_x) p <- p + scale_x_log10()
  p
}

reqtl_treat_tbl <- reqtl_power_cond %>%
  left_join(power_table %>% select(celltype, treatment, treat_donors, treat_cells),
            by = c("celltype", "treatment"))

pc1 <- make_reqtl_scatter(reqtl_power_cond, "ctrl_donors", "Number of control donors",
                          r_ctrl_don_reqtl, log_x = FALSE,
                          "A. Control donors")
pc2 <- make_reqtl_scatter(reqtl_power_cond, "ctrl_cells", "Control cells (log scale)",
                          r_ctrl_cell_reqtl, log_x = TRUE,
                          "B. Control cells")
pt1 <- make_reqtl_scatter(reqtl_treat_tbl, "treat_donors", "Number of treatment donors",
                          r_tx_don_reqtl, log_x = FALSE,
                          "C. Treatment donors")
pt2 <- make_reqtl_scatter(reqtl_treat_tbl, "treat_cells", "Treatment cells (log scale)",
                          r_tx_cell_reqtl, log_x = TRUE,
                          "D. Treatment cells")

p_reqtl_power <- (pc1 + pc2) / (pt1 + pt2) +
  plot_layout(guides = "collect") & theme(legend.position = "right")

reqtl_pdf_path <- file.path(plot_dir, "response_eqtls_vs_power.pdf")
ggsave(reqtl_pdf_path, p_reqtl_power, width = 13, height = 10)

cat("\nSaved: ", reqtl_pdf_path, "\n", sep = "")
cat("Saved: ", reqtl_supp_path, "\n", sep = "")

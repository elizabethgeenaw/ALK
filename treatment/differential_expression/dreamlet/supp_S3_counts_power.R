# Supp Fig S3: S3_fig_deg_abundance_scatter_preLOO.pdf

# Two things this script does:
#   (1) computes per-(drug, cell type) DE gene counts and the % of tested
#       genes that are DE. Writes deg_counts_by_assay_drug.csv and
#       deg_fraction_by_assay_drug.csv, both consumed by Fig 2A (dotplot).
#   (2) builds the DEG-count vs cell-type-abundance scatter (Fig S3): x =
#       total cells per (drug, cell type), y = number of DEGs (FDR < 0.05),
#       one facet per drug, per-drug Spearman rho annotated on each panel.
#       If rho is near zero the DEG count is not just a function of how many
#       cells were tested, which is what we want to demonstrate.

# Inputs:
#   results/treatment/dreamlet_clean_proteincoding/top_all_1v4.csv
#   results/treatment/dreamlet_clean/pb_treat                       (pseudobulk RDS, for cell counts)

# Outputs (to results/treatment/dreamlet_clean_proteincoding_summaries/):
#   deg_counts_by_assay_drug.csv        (consumed by Fig 2A)
#   deg_fraction_by_assay_drug.csv      (consumed by Fig 2A)
#   plots/deg_abundance_scatter.{pdf,png}  = Fig S3

library(BPCells)
library(Seurat)
library(SeuratObject)
library(dreamlet)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(scales)

# Parameters
fdr_threshold <- 0.05
drug_order <- c("Ale", "Bri", "Cer", "Cri", "Lor")
drug_colors <- c(Ale = "#0072B2", Bri = "#D55E00", Cer = "#009E73",
                 Cri = "#CC79A7", Lor = "#E69F00")

# Paths
base_dir <- "/project/gilad/geena/publication_clean"
in_dir <- file.path(base_dir, "results/treatment/dreamlet_clean_proteincoding")
results_dir <- file.path(base_dir, "results/treatment/dreamlet_clean_proteincoding_summaries")
plots_dir <- file.path(results_dir, "plots")
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

ts <- function() format(Sys.time(), "[%H:%M:%S]")

save_gg <- function(p, name, w, h) {
  ggsave(file.path(plots_dir, paste0(name, ".pdf")), p,
         width = w, height = h, device = cairo_pdf)
  ggsave(file.path(plots_dir, paste0(name, ".png")), p,
         width = w, height = h, dpi = 200)
}

# Load the dreamlet 1v4 DE table (one row per gene x cell type x drug).
message(ts(), " Loading top_all_1v4.csv...")
top_all <- read.csv(file.path(in_dir, "top_all_1v4.csv"),
                    stringsAsFactors = FALSE)
message(ts(), " ", nrow(top_all), " rows, ",
        n_distinct(top_all$assay), " cell types, ",
        n_distinct(top_all$drug), " drugs")

# Load the pseudobulk object to get per-(donor, cell type) cell counts. Sum
# across donors within each drug arm to get the "total cells per (drug, cell
# type)" number that Fig S3's x-axis uses.
message(ts(), " Loading pb_treat (for cell counts)...")
pb_treat <- readRDS(file.path(base_dir, "results/treatment/dreamlet_clean/pb_treat"))
cc <- cellCounts(pb_treat)
arm <- substr(rownames(cc), 1, 3)
rm(pb_treat)

cells_per_cond <- do.call(rbind, lapply(drug_order, function(d) {
  idx <- arm == d
  data.frame(
    drug = d,
    assay = colnames(cc),
    total_cells = as.integer(colSums(cc[idx, , drop = FALSE])),
    stringsAsFactors = FALSE
  )
}))
message(ts(), " Cell counts computed: ", nrow(cells_per_cond),
        " (drug, cell type) combinations")


# DEG counts per (drug, cell type)

# Take the significant subset (FDR < 0.05) and count DEGs per (drug, cell type).
sig <- filter(top_all, adj.P.Val < fdr_threshold)

counts_by_assay <- sig %>%
  group_by(assay, drug) %>%
  summarise(n_sig = n(), .groups = "drop") %>%
  # Fill in zero entries for (drug, cell type) pairs with no DEGs so downstream
  # readers get the full 5 x 30 grid rather than a sparse table.
  complete(assay = unique(top_all$assay), drug = drug_order,
           fill = list(n_sig = 0L)) %>%
  arrange(factor(drug, levels = drug_order), assay)

write.csv(counts_by_assay,
          file.path(results_dir, "deg_counts_by_assay_drug.csv"),
          row.names = FALSE)
message(ts(), " Written: deg_counts_by_assay_drug.csv")


# Fraction of tested genes that are DEGs per (drug, cell type)

# dreamlet applies its gene filter once per cell type, then fits 5 contrasts
# on that same gene set, so every drug within a cell type shares one
# denominator. We take max() defensively in case a rare gene is dropped for
# a particular contrast.
message(ts(), " Computing per-cell-type gene-tested denominator and DEG fractions...")

n_tested_by_ct <- top_all %>%
  group_by(assay, drug) %>%
  summarise(n = n_distinct(ID), .groups = "drop") %>%
  group_by(assay) %>%
  summarise(n_tested = max(n), .groups = "drop")

fraction_by_assay <- counts_by_assay %>%
  left_join(n_tested_by_ct, by = "assay") %>%
  mutate(pct_sig = 100 * n_sig / n_tested)

write.csv(fraction_by_assay,
          file.path(results_dir, "deg_fraction_by_assay_drug.csv"),
          row.names = FALSE)
message(ts(), " Written: deg_fraction_by_assay_drug.csv")


# Fig S3: DEG-count vs abundance scatter

# For each drug, compute Spearman rho between DEG count and total cells across
# the 30 cell types. Near-zero rho means "DEG count doesn't just track how
# many cells you had" -- i.e., biology, not power, drives the differences.
message(ts(), " Computing per-drug abundance correlations...")

abund <- counts_by_assay %>%
  inner_join(cells_per_cond, by = c("assay", "drug"))

rho_per_drug <- abund %>%
  group_by(drug) %>%
  summarise(rho = cor(n_sig, total_cells, method = "spearman",
                      use = "complete.obs"),
            n_pairs = n(),
            .groups = "drop") %>%
  mutate(drug = factor(drug, levels = drug_order)) %>%
  arrange(drug)

message(ts(), " Per-drug Spearman rho:")
print(rho_per_drug)

# Build the scatter: one facet per drug, one point per cell type. Label the
# top-4 cell types per drug by DEG count so readers can see which cell types
# drive the signal.
scatter_dat <- abund %>% mutate(drug = factor(drug, levels = drug_order))

scatter_labels <- scatter_dat %>%
  group_by(drug) %>%
  slice_max(n_sig, n = 4) %>%
  ungroup()

rho_labels <- rho_per_drug %>%
  mutate(label = paste0("rho = ", round(rho, 2)))

p_scatter <- ggplot(scatter_dat,
                    aes(x = total_cells, y = n_sig, color = drug)) +
  geom_point(size = 1.8, alpha = 0.75) +
  # Single dashed regression line (across all drugs) for a visual reference,
  # not a statistical claim; the per-drug rhos above are the real numbers.
  geom_smooth(aes(group = 1), method = "lm", se = FALSE,
              color = "grey40", linewidth = 0.5, linetype = "dashed") +
  geom_text_repel(
    data = scatter_labels,
    aes(label = assay),
    size = 2.0, color = "grey20",
    max.overlaps = Inf, min.segment.length = 0,
    segment.size = 0.25, box.padding = 0.3, force = 2,
    show.legend = FALSE
  ) +
  geom_text(
    data = rho_labels,
    aes(x = -Inf, y = Inf, label = label, color = drug),
    hjust = -0.1, vjust = 1.4, size = 2.8, fontface = "bold",
    inherit.aes = FALSE
  ) +
  facet_wrap(~ drug, nrow = 1, scales = "free_y") +
  scale_color_manual(values = drug_colors, guide = "none") +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "Total cells per cell type (summed across individuals)",
       y = "DEG count (FDR < 0.05)",
       title = "DEG count vs. cell-type abundance") +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 30, hjust = 1))

save_gg(p_scatter, "deg_abundance_scatter", w = 12, h = 3.5)
message(ts(), " Written: plots/deg_abundance_scatter.{pdf,png}")

message(ts(), " supp_S3_counts_power.R complete.")

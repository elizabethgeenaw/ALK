# Supp Fig S4: S4_fig_ctspecific_kequals1_bar.pdf

# Bar chart of cell-type-specific DEGs per drug. A gene is "cell-type-specific"
# for a drug if it hits lfsr < 0.05 in exactly one cell type. Each drug bar is
# split into two segments: the gene was tested in all 30 cell types (bottom,
# dark blue) vs. tested in fewer than 30 (top, pale). The "tested in all 30"
# subset is the more conservative call because those genes had every chance
# to hit multiple cell types and still landed in only one.

# Input:
#   results/treatment/dreamlet_mash/supp_table_ct_specific_with_ntested.csv
#   (or supp_table_ct_specific.csv if the n_ct_tested column is already there)

# Output:
#   results/treatment/dreamlet_mash/plots/fig_ctspecific_kequals1_bar.{pdf,png}
#   results/treatment/dreamlet_mash/plots/fig_ctspecific_kequals1_bar_source.csv

.libPaths(c("/project/gilad/geena/Rlibs/R-4.4.1", .libPaths()))

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# Paths
results_root <- "/project/gilad/geena/publication_clean/results/treatment"
mash_root <- file.path(results_root, "dreamlet_mash")
out_dir <- file.path(mash_root, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Constants
drug_colors <- c(Ale = "#0072B2", Bri = "#D55E00", Cer = "#009E73",
                 Cri = "#CC79A7", Lor = "#E69F00")
drug_full <- c(Ale = "Alectinib", Bri = "Brigatinib", Cer = "Ceritinib",
               Cri = "Crizotinib", Lor = "Lorlatinib")

# Segment colors: deep blue for the conservative subset (all 30 tested),
# pale blue for the "less conservative" subset (fewer than 30 tested).
seg_fills <- c(all30_tested = "#1F4E79",
               fewer_than_30 = "#BDD7EE")

# Load the augmented supp table. Prefer supp_table_ct_specific.csv if it
# already carries n_ct_tested; fall back to the "_with_ntested" variant.
ct_csv <- file.path(mash_root, "supp_table_ct_specific.csv")
if (!"n_ct_tested" %in% names(read_csv(ct_csv, n_max = 0, show_col_types = FALSE))) {
  ct_csv <- file.path(mash_root, "supp_table_ct_specific_with_ntested.csv")
}
ct <- read_csv(ct_csv, show_col_types = FALSE)

# For each drug, split the cell-type-specific DEGs into "tested in all 30
# cell types" vs "tested in fewer". Preserve the total count as any_ntested.
k1_counts <- ct %>%
  group_by(drug) %>%
  summarize(all30_tested = sum(n_ct_tested == 30),
            any_ntested = n(),
            .groups = "drop") %>%
  mutate(fewer_than_30 = any_ntested - all30_tested,
         drug_full = drug_full[drug]) %>%
  arrange(any_ntested) %>%
  mutate(drug_full = factor(drug_full, levels = drug_full))

# Long format for the stacked bar. Factor level order controls stacking:
# level 1 renders at TOP, level 2 at BOTTOM. Put fewer_than_30 first so
# all30_tested becomes the bottom (visual foundation of the bar).
k1_long <- k1_counts %>%
  pivot_longer(c(all30_tested, fewer_than_30),
               names_to = "segment", values_to = "n") %>%
  mutate(segment = factor(segment,
                          levels = c("fewer_than_30", "all30_tested")))

# Combined label above each bar: "total (all30_tested)"
totals <- k1_counts %>%
  mutate(label = paste0(scales::comma(any_ntested),
                        " (", scales::comma(all30_tested), ")"))

p <- ggplot(k1_long, aes(x = drug_full, y = n, fill = segment)) +
  geom_col(width = 0.7, colour = "black", linewidth = 0.3) +
  geom_text(data = totals,
            aes(x = drug_full, y = any_ntested, label = label),
            inherit.aes = FALSE,
            vjust = -0.3, size = 3.5) +
  scale_fill_manual(
    values = seg_fills,
    labels = c(all30_tested = "Tested in all 30 cell types",
               fewer_than_30 = "Tested in fewer than 30 cell types"),
    breaks = c("all30_tested", "fewer_than_30")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)),
                     labels = scales::comma) +
  labs(x = NULL,
       y = "Cell-type-specific DEGs\n(lfsr < 0.05 in exactly one cell type)",
       fill = NULL,
       title = "Cell-type-specific DEGs per ALK inhibitor") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "top",
        axis.text.x = element_text(angle = 0),
        plot.margin = margin(t = 5, r = 10, b = 5, l = 10))

ggsave(file.path(out_dir, "fig_ctspecific_kequals1_bar.pdf"),
       p, width = 8.0, height = 4.6)
ggsave(file.path(out_dir, "fig_ctspecific_kequals1_bar.png"),
       p, width = 8.0, height = 4.6, dpi = 300)

# Also write the source counts
write_csv(k1_counts, file.path(out_dir, "fig_ctspecific_kequals1_bar_source.csv"))

cat("Wrote:\n",
    file.path(out_dir, "fig_ctspecific_kequals1_bar.pdf"), "\n",
    file.path(out_dir, "fig_ctspecific_kequals1_bar.png"), "\n", sep = "")

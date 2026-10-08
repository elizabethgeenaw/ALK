#!/usr/bin/env Rscript
# Figures-only regen of Supp Fig S7 (Cormotif GO BP enrichment pages).
#
# No Cormotif refit, no topGO rerun. Reads the per-cell-type fitted RDS and GO
# BP CSVs already on disk under `output_base_dir`, re-derives the motif pattern
# table (so new motif_type categories like "group-specific" from the updated
# compute_motif_meta() are reflected), and re-assembles the per-cell-type pages
# via assemble_go_supplementary_pages().
#
# Use this after editing compute_motif_meta() or motif_type_display_label() in
# 03_cormotif_run.R without re-running the expensive Cormotif fit or the
# downstream rerun_topgo_from_saved() / replot_all() steps.
#
# Outputs:
#   <output_base_dir>/motif_pattern_table.{csv,xlsx}
#   <output_base_dir>/GO_BP_elim_combined/supp_GO_BP_elim_page01.png, page02.png, ...

# Guard: load 03_cormotif_run.R for its function definitions only (skip the
# main driver that calls rerun_topgo_from_saved(), replot_all(), etc.).
CORMOTIF_HELPERS_ONLY <- TRUE
source("/project/gilad/geena/ALK_paper_code/treatment/differential_expression/cormotif/03_cormotif_run.R")

output_base_dir <- "/project/gilad/geena/publication_clean/results/treatment/cormotif/cormotif_15cells_4arm"

# Refresh the pattern table so new motif_type categories flow through.
generate_pattern_table(output_base_dir)

# Re-assemble per-cell-type GO BP pages (uses motif_type_display_label() so
# group-specific motifs render as "<group> vs others" in panel titles).
assemble_go_supplementary_pages(
  output_base_dir,
  cell_type_order = names(cell_type_germ_layer),
  group_colors    = setNames(germ_layer_colors[cell_type_germ_layer],
                             names(cell_type_germ_layer))
)

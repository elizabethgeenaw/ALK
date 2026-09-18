#!/usr/bin/env Rscript
# Regenerate Fig 3B, 3C, and 3D without running the full cormotif_figures.R
# panel suite. cormotif_figures.R contains function definitions followed by
# a runner block that calls every panel; the runner also invokes an rrvgo
# panel (panel4) that is expensive and not part of the figures we want.
# This wrapper loads only the definitions and then calls the three fig3
# functions directly.
#
# Outputs land under:
#   results/treatment/cormotif/proteincoding_newfilter_15cells_4arm/figures/
#     fig3B_drug_celltype_matrix.{pdf,png}
#     fig3C_motif_go_enrichment.{pdf,png}
#     fig3C_values.csv
#     fig3D_lor_motif_upset.{pdf,png}
#
# Run on gilad-hm with R/4.4.1:
#   module load R/4.4.1
#   Rscript fig3B_3C_refresh.R

fp <- file.path(
  "/project/gilad/geena/publication_clean/code/treatment",
  "differential_expression/cormotif/cormotif_figures.R"
)

# Read the file and split off the trailing runner block. The runner starts
# with a `cat("===...` banner that also serves as a sentinel we can grep for.
# Only the definitions above the sentinel get eval'd so we don't trigger the
# panels we don't want.
lines <- readLines(fp)
main_start <- grep('^cat\\("=+', lines)
if (length(main_start) == 0)
  stop("Could not locate the trailing runner block in cormotif_figures.R")
defs <- lines[1:(main_start[1] - 1)]
eval(parse(text = paste(defs, collapse = "\n")))

cat("Loaded definitions from cormotif_figures.R (runner skipped)\n")
cat("Figures dir:", figures_dir(), "\n\n")

fig3B_drug_celltype_matrix()
fig3C_motif_go_enrichment()
fig3D_lor_motif_upset()

cat("=== Done ===\n")

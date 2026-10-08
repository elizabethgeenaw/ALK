#!/usr/bin/env Rscript
# Utility script: re-derive the motif pattern table from the fitted RDS files
# already written by Cormotif, without refitting Cormotif itself. Useful after
# a bug fix in the classification logic (e.g. drug_pattern, group-specific
# detection) or when regenerating the xlsx with additional columns.
#
# Reads per-cell-type `_motif_fittedK2plus.rds`, re-classifies each motif via
# compute_motif_meta() from 03_cormotif_run.R, then writes
# `motif_pattern_table.{csv,xlsx}` with the current schema (incl. motif_type
# categories like "group-specific").
#
# This script sources 03_cormotif_run.R for its helper definitions only (the
# main driver is skipped via CORMOTIF_HELPERS_ONLY), so the classification
# logic always matches whatever is current in 03_cormotif_run.R — no risk of
# drift from a stale inline duplicate.

CORMOTIF_HELPERS_ONLY <- TRUE
source("/project/gilad/geena/ALK_paper_code/treatment/differential_expression/cormotif/03_cormotif_run.R")

output_base_dir <- "/project/gilad/geena/publication_clean/results/treatment/cormotif/cormotif_15cells_4arm"

generate_pattern_table(output_base_dir)

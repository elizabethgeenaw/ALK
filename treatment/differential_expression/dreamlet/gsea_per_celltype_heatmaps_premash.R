# Pre-mash companion to supp_S5fig_S6_gsea_ct_patterns.R.
#
# Under the current figure numbering, the GO:BP per-drug heatmaps produced
# here are Supp Fig S6 (pre-mash), paired with their post-mash counterpart
# Supp Fig S7 from supp_S5fig_S6_gsea_ct_patterns.R. The Hallmark pre-mash
# heatmap this script also produces is a supplementary companion to Supp
# Fig S5 (Hallmark post-mash) and is not independently numbered.
#
# All heatmaps are built from the pre-mash fgsea table
# (dreamlet_pre_mash_gsea.tsv, ranked by the dreamlet moderated t-statistic).
# Outputs land in a sibling directory with a _premash suffix so the mash
# and pre-mash versions coexist for a direct sensitivity comparison.
#
# The construction logic is identical to the mash version, so this script
# re-uses the mash script by loading it as text, swapping four path tokens,
# and evaluating in a fresh environment. If the mash script changes, this
# companion inherits the change for free.
#
# Inputs:
#   results/treatment/dreamlet_mash/gsea/dreamlet_pre_mash_gsea.tsv
#   reference_sets/celltype_germlayer.csv
#
# Outputs (all to results/treatment/dreamlet_mash/ct_patterns_plots_premash/):
#   fig_supp_full_hallmark_ct_heatmap_all_drugs_stacked_premash.{pdf,png}
#   fig_supp_top_gobp_ct_heatmap_premash_{Ale,Bri,Cer,Cri,Lor}.{pdf,png}

template <- "/project/gilad/geena/ALK_paper_code/treatment/differential_expression/dreamlet/supp_S5fig_S6_gsea_ct_patterns.R"
if (!file.exists(template))
  stop("Expected mash template at: ", template)

src <- readLines(template)
src <- gsub("dreamlet_mash_gsea.tsv",
            "dreamlet_pre_mash_gsea.tsv",
            src, fixed = TRUE)
src <- gsub("ct_patterns_plots",
            "ct_patterns_plots_premash",
            src, fixed = TRUE)
src <- gsub("fig_supp_full_hallmark_ct_heatmap_all_drugs_stacked",
            "fig_supp_full_hallmark_ct_heatmap_all_drugs_stacked_premash",
            src, fixed = TRUE)
src <- gsub("fig_supp_top_gobp_ct_heatmap_",
            "fig_supp_top_gobp_ct_heatmap_premash_",
            src, fixed = TRUE)

eval(parse(text = paste(src, collapse = "\n")),
     envir = new.env(parent = globalenv()))

message("supp_S5fig_S6_premash.R complete.")

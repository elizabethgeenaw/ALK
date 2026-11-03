# Supp Tables S10 (full pre-mash GSEA) and S11 (pre-vs-post sensitivity).
#
# Builds the pre-vs-post sensitivity workbook that supports the manuscript's
# claim that pre-mash GSEA gives consistent results. The main sheet
# "pre_vs_post_cited" is a side-by-side summary for the gene sets the Results
# text names explicitly; "nes_correlations" reports Pearson/Spearman NES
# correlations between pre and post at two grains; "pre_mash_full" and
# "post_mash_full" are the per-(drug, collection, pathway) summaries that back
# the main sheet.
#
# Depends on 04_dreamlet_mash_fgsea.R and 04b_dreamlet_pre_mash_fgsea.R having
# been run.
#
# Inputs:
#   results/treatment/dreamlet_mash/gsea/dreamlet_mash_gsea.tsv       (post-mash)
#   results/treatment/dreamlet_mash/gsea/dreamlet_pre_mash_gsea.tsv   (pre-mash)
#
# Outputs:
#   S11_gsea_pre_vs_post_sensitivity.csv   (pre_vs_post_cited sheet alone)
#   S11_gsea_pre_vs_post_sensitivity.xlsx  (4 sheets: pre_vs_post_cited,
#                                            nes_correlations, pre_mash_full,
#                                            post_mash_full)

library(data.table)
library(openxlsx)

fdr_threshold <- 0.05
drug_order <- c("Ale", "Bri", "Cer", "Cri", "Lor")

base_dir <- "/project/gilad/geena/publication_clean"
gsea_dir <- file.path(base_dir, "results/treatment/dreamlet_mash/gsea")

post_tsv <- file.path(gsea_dir, "dreamlet_mash_gsea.tsv")
pre_tsv  <- file.path(gsea_dir, "dreamlet_pre_mash_gsea.tsv")

stopifnot(file.exists(post_tsv), file.exists(pre_tsv))

ts <- function() format(Sys.time(), "[%H:%M:%S]")

# Cited gene sets (see Results narrative for which pathways are referenced).
cited <- list(
  Hallmark = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
               "HALLMARK_MTORC1_SIGNALING",
               "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
               "HALLMARK_CHOLESTEROL_HOMEOSTASIS"),
  GOBP     = c("GOBP_STEROL_BIOSYNTHETIC_PROCESS",
               "GOBP_REGULATION_OF_PLASMA_LIPOPROTEIN_PARTICLE_LEVELS",
               "GOBP_LENS_DEVELOPMENT_IN_CAMERA_TYPE_EYE",
               "GOBP_NEURAL_RETINA_DEVELOPMENT",
               "GOBP_RNA_SPLICING")
)

message(ts(), " Reading GSEA tables...")
post <- fread(post_tsv)
pre  <- fread(pre_tsv)

# Normalise collection labels across the two files (defensive).
norm_coll <- function(x)
  fifelse(tolower(x) == "hallmark", "Hallmark",
          fifelse(toupper(x) %in% c("GOBP", "GO:BP"), "GOBP", x))
pre[,  collection := norm_coll(collection)]
post[, collection := norm_coll(collection)]

# Per (drug, collection, pathway) summary: cell-type counts (total + sig +
# sig-up + sig-dn) and NES stats (mean + median absolute).
summarize_set <- function(gsea, which_label) {
  gsea <- as.data.table(gsea)
  gsea[, .(
    n_ct_tested    = .N,
    n_ct_sig       = sum(padj < fdr_threshold, na.rm = TRUE),
    n_ct_sig_up    = sum(padj < fdr_threshold & NES > 0, na.rm = TRUE),
    n_ct_sig_dn    = sum(padj < fdr_threshold & NES < 0, na.rm = TRUE),
    mean_NES       = mean(NES, na.rm = TRUE),
    median_abs_NES = median(abs(NES), na.rm = TRUE)
  ), by = .(drug, collection, pathway)][, source := which_label][]
}

message(ts(), " Summarising post-mash...")
post_s <- summarize_set(post, "post")
message(ts(), " Summarising pre-mash...")
pre_s  <- summarize_set(pre,  "pre")

# Side-by-side table for the cited gene sets (long form -> wide-ish merge).
cited_all <- unlist(cited, use.names = FALSE)
post_c <- post_s[pathway %in% cited_all]
pre_c  <- pre_s [pathway %in% cited_all]

reshape_side <- function(dt) {
  dcast(
    melt(dt, measure.vars = c("n_ct_tested", "n_ct_sig",
                               "n_ct_sig_up", "n_ct_sig_dn",
                               "mean_NES", "median_abs_NES"),
         id.vars = c("drug", "collection", "pathway")),
    drug + collection + pathway ~ variable,
    value.var = "value"
  )
}

main <- merge(
  reshape_side(post_c),
  reshape_side(pre_c),
  by = c("drug", "collection", "pathway"),
  suffixes = c("_post", "_pre"), all = TRUE
)

main[, delta_n_sig    := n_ct_sig_post - n_ct_sig_pre]
main[, delta_mean_NES := round(mean_NES_post - mean_NES_pre, 3)]
main[, drug := factor(drug, levels = drug_order)]
setcolorder(main, c("collection", "pathway", "drug",
                    "n_ct_tested_pre",  "n_ct_sig_pre",
                    "n_ct_sig_up_pre",  "n_ct_sig_dn_pre",
                    "mean_NES_pre",     "median_abs_NES_pre",
                    "n_ct_tested_post", "n_ct_sig_post",
                    "n_ct_sig_up_post", "n_ct_sig_dn_post",
                    "mean_NES_post",    "median_abs_NES_post",
                    "delta_n_sig",      "delta_mean_NES"))
main <- main[order(collection, pathway, drug)]

# NES correlations between pre and post:
#   r_dcg  = across (drug, cell type, pathway) within each collection
#   r_mean = across (drug, pathway) using the mean NES across cell types
make_cor_df <- function(pre_raw, post_raw) {
  joint <- merge(
    pre_raw[,  .(drug, collection, pathway, assay, NES_pre = NES)],
    post_raw[, .(drug, collection, pathway, assay, NES_post = NES)],
    by = c("drug", "collection", "pathway", "assay"))
  dcg <- joint[, .(
    n_pairs    = sum(!is.na(NES_pre) & !is.na(NES_post)),
    r_pearson  = round(cor(NES_pre, NES_post, use = "complete.obs"), 3),
    r_spearman = round(cor(NES_pre, NES_post, use = "complete.obs",
                            method = "spearman"), 3)
  ), by = collection][, grain := "drug_x_ct_x_pathway"][]
  mean_dt <- joint[, .(mean_pre  = mean(NES_pre,  na.rm = TRUE),
                       mean_post = mean(NES_post, na.rm = TRUE)),
                   by = .(drug, collection, pathway)]
  mean_cor <- mean_dt[, .(
    n_pairs    = sum(!is.na(mean_pre) & !is.na(mean_post)),
    r_pearson  = round(cor(mean_pre, mean_post, use = "complete.obs"), 3),
    r_spearman = round(cor(mean_pre, mean_post, use = "complete.obs",
                            method = "spearman"), 3)
  ), by = collection][, grain := "drug_x_pathway (mean NES across CTs)"][]
  rbind(dcg, mean_cor)
}
cor_df <- make_cor_df(pre, post)
cat("\n  Pre- vs post-mash NES correlations:\n")
print(cor_df)

# Pre-mash sanity: for each cited pathway, how many CTs are significant for
# at least one drug? (printed, not written.)
cat("\n  Pre-mash sanity (n_ct_sig range across drugs per cited pathway):\n")
pre_check <- pre_c[, .(min_sig = min(n_ct_sig), max_sig = max(n_ct_sig)),
                   by = .(collection, pathway)]
print(pre_check)

# Write outputs. The CSV holds just the pre_vs_post_cited sheet; the XLSX
# holds all four sheets so the full summaries live alongside the comparison.
out_csv  <- file.path(gsea_dir, "S11_gsea_pre_vs_post_sensitivity.csv")
out_xlsx <- file.path(gsea_dir, "S11_gsea_pre_vs_post_sensitivity.xlsx")

fwrite(main, out_csv)
message(ts(), " Wrote: ", out_csv)

wb <- createWorkbook()
addWorksheet(wb, "pre_vs_post_cited")
writeData(wb, "pre_vs_post_cited", main)
addWorksheet(wb, "nes_correlations")
writeData(wb, "nes_correlations", cor_df)
addWorksheet(wb, "pre_mash_full")
writeData(wb, "pre_mash_full",  pre_s)
addWorksheet(wb, "post_mash_full")
writeData(wb, "post_mash_full", post_s)
saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message(ts(), " Wrote: ", out_xlsx)

message(ts(), " supp_S10_S11_pre_mash_sensitivity.R complete.")

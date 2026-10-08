# Supp Table S10 (long form): S10_dreamlet_pre_mash_gsea.tsv

# Pre-mash sensitivity companion to 04_dreamlet_mash_fgsea.R. Runs
# fgseaMultilevel on the raw dreamlet 1v4 DE statistics (no mash shrinkage),
# ranking genes within each (drug, cell type) by t = the dreamlet moderated
# t-statistic from the 1v4 contrast. Same gene-set collections, same fgsea
# arguments as the mash version, so the two outputs are directly comparable
# per (drug, cell type, pathway).

# The output is the pre-mash analog of dreamlet_mash_gsea.tsv and the input
# to supp_S5fig_S6_premash.R and supp_S10_S11_pre_mash_sensitivity.R.

# Inputs:
#   results/treatment/dreamlet_1v4_DE/top_all_1v4.csv

# Outputs (all to results/treatment/dreamlet_mash/gsea/):
#   dreamlet_pre_mash_gsea.tsv           full fgsea results (Supp Table S10)
#   dreamlet_pre_mash_gsea_metadata.txt  msigdbr version and ranking / fgsea params

library(fgsea)
library(msigdbr)
library(dplyr)
library(data.table)

set.seed(42)

# Parameters (match 04_dreamlet_mash_fgsea.R exactly for comparability)
min_set_size <- 15
max_set_size <- 500
n_perm <- 10000
drug_order <- c("Ale", "Bri", "Cer", "Cri", "Lor")

# Paths
base_dir <- "/project/gilad/geena/publication_clean"
de_dir   <- file.path(base_dir, "results/treatment/dreamlet_1v4_DE")
mash_dir <- file.path(base_dir, "results/treatment/dreamlet_mash")
out_dir  <- file.path(mash_dir, "gsea")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ts <- function() format(Sys.time(), "[%H:%M:%S]")

# Fetch MSigDB gene sets.
message(ts(), " Getting gene sets from msigdbr...")

fetch_msig <- function(cat, subcat = NULL) {
  tryCatch(
    msigdbr(species = "Homo sapiens", collection = cat, subcollection = subcat),
    error = function(e)
      msigdbr(species = "Homo sapiens", category = cat, subcategory = subcat)
  )
}

h_df <- fetch_msig("H")
gobp_df <- fetch_msig("C5", "GO:BP")

to_pw_list <- function(df) split(df$gene_symbol, df$gs_name)

collections <- list(
  Hallmark = to_pw_list(h_df),
  GOBP = to_pw_list(gobp_df)
)
message(ts(), " Hallmark: ", length(collections$Hallmark),
        " sets | GO:BP: ", length(collections$GOBP), " sets")

filter_pw <- function(pw_list, universe_genes) {
  pw <- lapply(pw_list, function(g) intersect(g, universe_genes))
  pw[lengths(pw) >= min_set_size & lengths(pw) <= max_set_size]
}

# Run fgsea on one (drug, cell type) ranked t-vector; same per-column helper
# as the mash version so the output schema is identical.
run_fgsea_on_column <- function(t_vec, drug_label, assay_label) {
  t_vec <- t_vec[is.finite(t_vec)]
  if (length(t_vec) == 0) return(NULL)

  ranked <- sort(t_vec, decreasing = TRUE)

  out <- vector("list", length(collections))
  names(out) <- names(collections)

  for (coll_name in names(collections)) {
    valid_pw <- filter_pw(collections[[coll_name]], names(ranked))
    if (length(valid_pw) == 0) next

    res <- tryCatch(
      fgseaMultilevel(pathways = valid_pw,
                      stats = ranked,
                      minSize = min_set_size,
                      maxSize = max_set_size,
                      nPermSimple = n_perm,
                      eps = 0),
      error = function(e) {
        message("  fgsea error for ", drug_label, "/", assay_label,
                "/", coll_name, ": ", e$message)
        NULL
      }
    )
    if (is.null(res)) next

    res <- as.data.frame(res)
    res$drug <- drug_label
    res$assay <- assay_label
    res$collection <- coll_name
    if ("leadingEdge" %in% names(res))
      res$leadingEdge <- sapply(res$leadingEdge,
                                function(x) paste(sort(x), collapse = "|"))
    out[[coll_name]] <- res
  }
  bind_rows(out[!sapply(out, is.null)])
}

# Main loop: read the combined per-drug DE table, rank each (drug x cell type)
# by the dreamlet moderated t-statistic, run fgsea.
message(ts(), " Loading top_all_1v4.csv...")
top_all <- read.csv(file.path(de_dir, "top_all_1v4.csv"),
                    stringsAsFactors = FALSE)

all_results <- vector("list", length(drug_order))
for (i in seq_along(drug_order)) {
  drug <- drug_order[i]
  message(ts(), " [", i, "/", length(drug_order), "] ", drug)
  sub_drug <- top_all[top_all$drug == drug, , drop = FALSE]
  assays   <- sort(unique(sub_drug$assay))

  drug_results <- vector("list", length(assays))
  for (j in seq_along(assays)) {
    assay_label <- assays[j]
    sub <- sub_drug[sub_drug$assay == assay_label, , drop = FALSE]
    # The gene identifier column in top_all is ID (from topTable rownames).
    gene_col <- if ("ID" %in% colnames(sub)) "ID" else "gene"
    t_vec <- setNames(sub$t, sub[[gene_col]])
    drug_results[[j]] <- run_fgsea_on_column(t_vec, drug, assay_label)
  }
  all_results[[i]] <- bind_rows(drug_results[!sapply(drug_results, is.null)])
}

gsea_all <- bind_rows(all_results[!sapply(all_results, is.null)])
message(ts(), " fGSEA complete: ", nrow(gsea_all), " rows")

fwrite(gsea_all,
       file.path(out_dir, "dreamlet_pre_mash_gsea.tsv"),
       sep = "\t")
message(ts(), " Written: dreamlet_pre_mash_gsea.tsv")

# Metadata capture (same style as 04_dreamlet_mash_fgsea.R).
meta_path <- file.path(out_dir, "dreamlet_pre_mash_gsea_metadata.txt")
writeLines(c(
  paste0("# Pre-mash fgsea run metadata (", format(Sys.time()), ")"),
  paste0("ranking: dreamlet moderated t (raw 1v4 contrast, no mash shrinkage)"),
  paste0("min_set_size: ", min_set_size),
  paste0("max_set_size: ", max_set_size),
  paste0("nPermSimple: ", n_perm),
  paste0("eps: 0"),
  paste0("msigdbr: ", as.character(utils::packageVersion("msigdbr"))),
  paste0("Hallmark sets: ", length(collections$Hallmark)),
  paste0("GO:BP sets: ", length(collections$GOBP)),
  paste0("input: ", file.path(de_dir, "top_all_1v4.csv"))
), meta_path)
message(ts(), " Written: dreamlet_pre_mash_gsea_metadata.txt")

message(ts(), " 04b_dreamlet_pre_mash_fgsea.R complete.")

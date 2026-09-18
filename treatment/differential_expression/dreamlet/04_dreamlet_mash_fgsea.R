# Supp Table S9: S9_dreamlet_mash_gsea.tsv

# Runs fgseaMultilevel on the mash-shrunken 1v4 DE posteriors (five per-drug
# mash fits, 30 cell-type columns each). For each (drug, cell type) we rank
# genes by t = posterior_mean / posterior_sd (as in Umans et al.) and test
# every MSigDB Hallmark and GO:BP pathway against that ranking.

# The output is one long fgsea table (one row per drug x cell type x pathway),
# which is Supp Table S9 and the input to every downstream pathway figure
# (Fig 2C, Fig 2D, Supp Fig S5/S6). Also writes a small metadata file
# capturing the msigdbr version, the ranking definition, and fgsea params
# so the run is reproducible, plus a GO ID -> pathway name lookup used by
# downstream scripts that need pretty GO term labels.

# Inputs:
#   results/treatment/dreamlet_mash/by_drug_{Ale,Bri,Cer,Cri,Lor}/{mash_pm.csv, mash_psd.csv}

# Outputs (all to results/treatment/dreamlet_mash/gsea/):
#   dreamlet_mash_gsea.tsv           full fgsea results (Supp Table S9)
#   dreamlet_mash_gsea_metadata.txt  msigdbr version and ranking / fgsea params
#   dreamlet_mash_gsea_go_name_map.csv  GO ID -> GOBP pathway name lookup

library(fgsea)
library(msigdbr)
library(dplyr)
library(data.table)

set.seed(42)

# Parameters
# Gene-set size filter passed to fgsea. Very small sets are noisy; very large
# sets tend to be uninformative (broad "response to X" categories).
min_set_size <- 15
max_set_size <- 500
# fgseaMultilevel's initial permutation count; higher = more precise p-values
# for very small p, at the cost of runtime.
n_perm <- 10000
drug_order <- c("Ale", "Bri", "Cer", "Cri", "Lor")

# Paths
base_dir <- "/project/gilad/geena/publication_clean"
mash_dir <- file.path(base_dir, "results/treatment/dreamlet_mash")
out_dir <- file.path(mash_dir, "gsea")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ts <- function() format(Sys.time(), "[%H:%M:%S]")

# Fetch MSigDB gene sets. The Hallmark collection (H) is a small, curated set
# of 50 well-defined biological programs; GO:BP is a much larger tree of GO
# biological process terms. msigdbr changed its API between versions, so try
# the new argument names first and fall back to the old ones.
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

# Save GO ID -> pretty pathway name lookup for downstream scripts.
go_name_map <- gobp_df %>% distinct(gs_exact_source, gs_name)
write.csv(go_name_map,
          file.path(out_dir, "dreamlet_mash_gsea_go_name_map.csv"),
          row.names = FALSE)


# Restrict a pathway list to genes actually tested, then drop pathways whose
# intersection with the tested-gene universe falls outside the size filter.
filter_pw <- function(pw_list, universe_genes) {
  pw <- lapply(pw_list, function(g) intersect(g, universe_genes))
  pw[lengths(pw) >= min_set_size & lengths(pw) <= max_set_size]
}

# Read a mash CSV (rows = genes, cols = cell types; column names may contain
# spaces, so check.names = FALSE).
read_mash_matrix <- function(path) {
  m <- read.csv(path, row.names = 1, check.names = FALSE,
                stringsAsFactors = FALSE)
  as.matrix(m)
}

# Run fgsea on one (drug, cell type) column's t-vector, then tag results with
# drug + assay + collection so we can bind across everything at the end.
# leadingEdge is flattened to a single pipe-separated string so the output
# writes cleanly to TSV.
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

# Main loop: for each drug, load PM and PSD matrices, compute t = PM/PSD,
# run fgsea on every cell-type column, and collect results.
message(ts(), " Running fgsea on dreamlet_mash/by_drug_* ...")

all_results <- vector("list", length(drug_order))
for (i in seq_along(drug_order)) {
  drug <- drug_order[i]
  dir_path <- file.path(mash_dir, paste0("by_drug_", drug))
  pm_path <- file.path(dir_path, "mash_pm.csv")
  psd_path <- file.path(dir_path, "mash_psd.csv")

  message(ts(), " [", i, "/", length(drug_order), "] ", drug)
  if (!file.exists(pm_path) || !file.exists(psd_path)) {
    message(ts(), " SKIP (missing files): ", dir_path)
    next
  }

  pm <- read_mash_matrix(pm_path)
  psd <- read_mash_matrix(psd_path)

  # Restrict to genes and cell types present in both matrices (defensive).
  common_rows <- intersect(rownames(pm), rownames(psd))
  common_cols <- intersect(colnames(pm), colnames(psd))
  pm <- pm[common_rows, common_cols, drop = FALSE]
  psd <- psd[common_rows, common_cols, drop = FALSE]

  # Ranking statistic: t = posterior mean / posterior sd, one column per CT.
  t_mat <- pm / psd

  drug_results <- vector("list", ncol(t_mat))
  for (j in seq_len(ncol(t_mat))) {
    assay_label <- colnames(t_mat)[j]
    t_vec <- t_mat[, j]
    names(t_vec) <- rownames(t_mat)
    drug_results[[j]] <- run_fgsea_on_column(t_vec, drug, assay_label)
  }
  all_results[[i]] <- bind_rows(drug_results[!sapply(drug_results, is.null)])
}

gsea_all <- bind_rows(all_results[!sapply(all_results, is.null)])
message(ts(), " fGSEA complete: ", nrow(gsea_all), " rows")

fwrite(gsea_all,
       file.path(out_dir, "dreamlet_mash_gsea.tsv"),
       sep = "\t")
message(ts(), " Written: dreamlet_mash_gsea.tsv")

message(ts(), " 04_dreamlet_mash_fgsea.R complete.")

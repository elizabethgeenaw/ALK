# Supp Table S23: /project/gilad/geena/ALK_draft/V4/Supplementary Tables/S23_crisp_fgsea_sig_nominal_p05.tsv

# Run fgseaMultilevel on pre-computed z-scores (from 3_zscore.py)
set.seed(1)
library(fgsea)
library(data.table)
library(qvalue) # https://github.com/StoreyLab/qvalue

pred_dir <- "/project/gilad/geena/chem_perturbation_clean/results/predictions/"
out_dir <- file.path(pred_dir, "zscore_analysis")
#dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Read z-scores calculated in 3_zscore.py
zscore_path <- file.path(out_dir, "gene_zscores_all.tsv")
if (!file.exists(zscore_path)) {
  stop("Z-score file not found: ", zscore_path, "\nRun 3_zscore.py first.")
}
zscores <- fread(zscore_path)
message("Loaded ", nrow(zscores), " z-score rows across ", uniqueN(zscores[, .(cell_type, treatment)]), " conditions")

# Response eQTLs
resp_df <- fread("/project/gilad/geena/publication_clean/results/mash/medianUMI1/lfc3/treatment_responsive_eqtls_signfixed.tsv")

# Get genes tested in mash
mash <- fread("/project/gilad/geena/publication_clean/results/mash/medianUMI1/combined_mash_results.tsv", select = "gene_snp")
mash_universe <- unique(sub("_.*", "", mash$gene_snp))
message("Number of genes tested in mash: ", length(mash_universe))

# Run fGSEA per condition
conditions <- unique(zscores[, .(cell_type, treatment)])
results <- list()

for (i in seq_len(nrow(conditions))) {
  ct <- conditions$cell_type[i]
  tr <- conditions$treatment[i]
  ct_key <- gsub(" ", "_", ct)

  # Get response eGenes for this condition
  resp_genes <- unique(resp_df[celltype == ct_key & treatment == tr, gene])

  # Get z-scores, restrict to genes tested in mash
  sub_z <- zscores[cell_type == ct & treatment == tr]
  sub_z <- sub_z[gene %in% mash_universe]
  sub_z <- sub_z[!duplicated(gene)]

  # Named vector of z-scores
  stats <- setNames(sub_z$z_score, sub_z$gene)
  stats <- sort(stats, decreasing = TRUE)

  # Gene set: response eGenes in mash universe
  resp_in_bg <- intersect(resp_genes, mash_universe)

  pathways <- list(response_eGenes = resp_in_bg)

  res <- fgseaMultilevel(pathways = pathways,
                    stats = stats,
                    minSize = 5,
                    maxSize = 100000,
                    eps = 0)

  if (!is.null(res) && nrow(res) > 0) {
    results[[length(results) + 1]] <- data.table(
      cell_type = ct,
      treatment = tr,
      background = "mash_universe",
      NES = res$NES,
      es = res$ES,
      p_value = res$pval,
      padj_fgsea = res$padj,
      leading_edge_n = res$size
    )
  }
}

fgsea_df <- rbindlist(results)

# BH correction across all conditions, add to fgsea_df
fgsea_df[, adj_p_value := p.adjust(p_value, method = "BH")] 

# Add Storey q-values to fgsea_df
qobj <- qvalue(fgsea_df$p_value)
fgsea_df[, qvalue := qobj$qvalues]
message("Estimated pi0 = ", qobj$pi0) # estimate of the proportion of null p-values
summary(qobj)



out_path <- file.path(out_dir, "fgsea_results.tsv")
fwrite(fgsea_df, out_path, sep = "\t")
message("Saved fgsea -> ", out_path, " (n=", nrow(fgsea_df), ")")

# Print any conditions with significant p-values pre-multiple testing correction
message("\nConditions with nominal p < 0.05:")
sig_nominal <- fgsea_df[p_value < 0.05][order(p_value)]
print(sig_nominal[, .(cell_type, treatment, NES, p_value)])
fwrite(sig_nominal, file.path(out_dir, "fgsea_sig_nominal_p05.tsv"), sep = "\t")

# Print any conditions significant after BH correction
message("\nConditions with BH-adjusted p < 0.05:")
sig_bh <- fgsea_df[adj_p_value < 0.05][order(adj_p_value)]
print(sig_bh[, .(cell_type, treatment, NES, adj_p_value)])

# Print any conditions significant with q-val<0.25
message("\nSignificance summary:")
cutoff = 0.25
n_bh <- sum(fgsea_df$adj_p_value < cutoff)
n_qv <- sum(fgsea_df$qvalue < cutoff)
message("  cutoff <", cutoff, ": BH=", n_bh, " qvalue=", n_qv, " (out of ", nrow(fgsea_df), ")")

sig <- fgsea_df[qvalue < 0.25][order(qvalue)]
print(sig[, .(cell_type, treatment, NES, p_value, adj_p_value, qvalue)])

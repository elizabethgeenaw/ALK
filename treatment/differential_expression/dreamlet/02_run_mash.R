.libPaths(c("/project/gilad/geena/Rlibs/R-4.4.1", .libPaths()))
library(dreamlet)
library(mashr)

# Run mash (multivariate adaptive shrinkage) on the dreamlet 1v4 fit, per
# drug. mash borrows strength across cell types to shrink noisy DE estimates
# toward the patterns supported by the whole coefficient matrix -- so a gene
# that trends the same direction in many cell types gets pulled toward that
# shared signal, while an isolated large estimate in one cell type gets
# shrunk toward zero.

# For each drug we feed its 30-cell-type 1v4 coefficient into dreamlet::run_mash
# and write the outputs to by_drug_<Drug>/: res_mash.rds, mash_pm.csv,
# mash_lfsr.csv, mash_psd.csv, logFC_original.csv, pi_estimated.csv,
# V_estimated.csv. These per-drug outputs feed every downstream mash-based
# figure and supplementary table (Fig 2C, 2D, Supp Figs S5/S6, Supp Tables
# S6/S7/S8/S9).

# Called by 02_run_mash_launcher.sh, which submits the 5 per-drug jobs in
# parallel via 02_run_mash_single.sh.

# Input:  results/treatment/dreamlet_clean_proteincoding/res.dl2_treat (dreamlet fit)
# Output: results/treatment/dreamlet_mash/by_drug_<Drug>/*

drugs <- c("Ale", "Bri", "Cer", "Cri", "Lor")
rds_path <- "/project/gilad/geena/publication_clean/results/treatment/dreamlet_clean_proteincoding/res.dl2_treat"
out_dir <- "/project/gilad/geena/publication_clean/results/treatment/dreamlet_mash"

# Turn a cell-type name into something safe for a directory name (spaces and
# punctuation get collapsed to underscores).
sanitize <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

# Given a dreamlet fit, find the coefficient corresponding to a drug's 1v4
# contrast. The coefficient name looks like "Ale_vs_others = ..." so we grep
# for the drug prefix.
resolve_coef <- function(fit, drug) {
  m <- grep(paste0("^", drug, "_vs_others"), coefNames(fit), value = TRUE)
  if (length(m) != 1)
    stop("Expected exactly one coefficient for ", drug, ", found ", length(m))
  m
}

# mash stores its estimated covariance matrix V in different places depending
# on version. Look in both plausible slots.
extract_V <- function(model) {
  V <- tryCatch(model$fitted_g$V, error = function(e) NULL)
  if (is.null(V)) V <- tryCatch(model$V, error = function(e) NULL)
  V
}

# Save the full set of mash outputs for one per-drug run into sub_dir.
save_outputs <- function(res, sub_dir) {
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
  m <- res$model
  saveRDS(res, file.path(sub_dir, "res_mash.rds"))
  write.csv(get_pm(m), file.path(sub_dir, "mash_pm.csv"))    # posterior means
  write.csv(get_lfsr(m), file.path(sub_dir, "mash_lfsr.csv"))  # local FSR
  write.csv(get_psd(m), file.path(sub_dir, "mash_psd.csv"))   # posterior SDs
  write.csv(res$logFC.original,
            file.path(sub_dir, "logFC_original.csv"))         # pre-mash log2FC
  pi_vec <- get_estimated_pi(m)
  write.csv(data.frame(component = names(pi_vec),
                       weight = as.numeric(pi_vec)),
            file.path(sub_dir, "pi_estimated.csv"), row.names = FALSE)
  V <- extract_V(m)
  if (!is.null(V)) {
    write.csv(V, file.path(sub_dir, "V_estimated.csv"))
    message("V saved. dim = ", paste(dim(V), collapse = " x "))
  } else {
    message("V not found in model object; skipping V csv")
  }
  message("Saved outputs to: ", sub_dir)
}

# Idempotency check: don't recompute if the key outputs are already there.
output_exists <- function(sub_dir) {
  all(file.exists(file.path(sub_dir, c("res_mash.rds", "mash_lfsr.csv"))))
}

# Fit mash on one drug's 30-cell-type coefficient. Writes outputs to
# by_drug_<Drug>/.
run_by_drug <- function(fit, drug) {
  message("=== Starting per-drug mash: ", drug, " ===")
  sub_dir <- file.path(out_dir, paste0("by_drug_", sanitize(drug)))
  if (output_exists(sub_dir)) {
    message("Output already exists -- skipping: ", sub_dir)
    return(invisible(NULL))
  }
  coef <- resolve_coef(fit, drug)
  message("Resolved coef: ", substr(coef, 1, 80))
  message("Calling dreamlet::run_mash ...")
  res <- dreamlet::run_mash(fit, coefList = coef)
  message("run_mash complete. Model class: ",
          paste(class(res$model), collapse = "/"))
  save_outputs(res, sub_dir)
  message("=== Done with drug: ", drug, " ===")
}

# CLI dispatch. With no args, print the drug list (used by the launcher for
# discovery). --drug <D> runs the mash fit for that drug.
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  writeLines(drugs)
  quit(save = "no")
}

if (args[1] != "--drug" || length(args) < 2)
  stop("Usage: Rscript 02_run_mash.R --drug <Drug>")

drug <- args[2]
if (!(drug %in% drugs))
  stop("Unknown drug: ", drug, ". Choose from: ", paste(drugs, collapse = ", "))

message("Loading dreamlet RDS: ", rds_path)
fit <- readRDS(rds_path)
message("Loaded. n cell types: ", length(assayNames(fit)),
        "; n coefs: ", length(coefNames(fit)))

run_by_drug(fit, drug)

message("Done. Outputs in ", out_dir)

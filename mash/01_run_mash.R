library(dplyr)
library(readr)
library(data.table)
library(ashr)
library(mashr)
library(flashier)

set.seed(1)

########### RUN MASH ##########
# Define output directory
output_dir <- paste0("/project/gilad/geena/publication_clean/results/mash/medianUMI1/")


# Create the directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  print("Output directory created")
}

input_mash <- readRDS("/project/gilad/geena/publication_clean/mash_input/medianUMI1_4.30.26/medianUMI1.mash.rds")


# Estimate null correlation using EM method
cat("Estimating null correlation via EM (mash_estimate_corr_em)...\n")
data.temp = mash_set_data(input_mash$random.b, input_mash$random.s)
U.c = cov_canonical(data.temp)
Vhat_em = mash_estimate_corr_em(data.temp, U.c, max_iter = 5, details = TRUE) # setting to max_iter 5 (default is 30)
saveRDS(Vhat_em, file = paste0(output_dir, "Vhat_em.rds")) # Save estimated Vhat_em object
data.random = mash_update_data(data.temp, V = get_fitted_g(Vhat_em)$V)
rm(data.temp)

##### DO WITH EZ (DEFAULT) #####
# Set up mash data
cat("Setting up mash data...\n")
data.strong = mash_set_data(input_mash$strong.b, input_mash$strong.s, V=Vhat_em$V)
saveRDS(data.random, paste0(output_dir, "data_random.rds"))
saveRDS(data.strong, paste0(output_dir, "data_strong.rds"))

Vhat_em <- readRDS(paste0(output_dir, "Vhat_em.rds"))
cat("Vhat_em read...\n")
data.random <- readRDS(paste0(output_dir, "data_random.rds"))
cat("data_random read...\n")
data.strong <- readRDS(paste0(output_dir, "data_strong.rds"))
cat("data_strong read...\n")

# Perform empirical Bayes matrix factorization via flashr
cat("Performing empirical Bayes matrix factorization via flashr...\n")
U.f = cov_flash(data.strong, factors="nonneg", tag="non_neg")
n_conditions <- ncol(data.strong$Bhat)
n_pcs <- min(5, n_conditions)
cat("Performing PCA...\n")
U.pca = cov_pca(data.strong, n_pcs) # make sure number of PCs isn't greater than number of conditions
cat("Learning mixture weights......\n")
U.ed = cov_ed(data.strong, c(U.f, U.pca), maxiter=100) 



cat("Saving U_f...\n")
tryCatch({
  saveRDS(U.f, paste0(output_dir, "U_f.rds"))
}, error = function(e) {
  message("Warning: Failed to save U_f.rds -- ", conditionMessage(e))
})

cat("Saving U_pca...\n")
tryCatch({
  saveRDS(U.pca, paste0(output_dir, "U_pca.rds"))
}, error = function(e) {
  message("Warning: Failed to save U_pca.rds -- ", conditionMessage(e))
})

cat("Saving U_ed...\n")
tryCatch({
  saveRDS(U.ed, paste0(output_dir, "U_ed.rds"))
}, error = function(e) {
  message("Warning: Failed to save U_ed.rds -- ", conditionMessage(e))
})



# Canonical covariances
cat("Canonical covariances...\n")
U.c = cov_canonical(data.random)

# Fit mash model
cat("Fitting mash model...\n")
m = mash(data.random, Ulist = c(U.ed, U.c), outputlevel = 1)
tryCatch({
  saveRDS(m, paste0(output_dir, "m.rds"))
}, error = function(e) {
  message("Warning: Failed to save m.rds -- ", conditionMessage(e))
})

# Use posterior fitted distribution
print("Using posterior fitted distribution...\n")
m2 = mash(data.strong, g=get_fitted_g(m), fixg=TRUE)

# Extract local false sign rate (LFSR)
cat("Writing lfsr results...\n")
results <- get_lfsr(m2)
write_tsv(as.data.frame(results), paste0(output_dir, "lfsr_results.tsv"))

# Identify significant effects
cat("Writing significant effects...\n")
significant_results <- get_significant_results(m2) 
write_tsv(as.data.frame(significant_results), paste0(output_dir, "significant_results.tsv"))

# Compute pairwise sharing (https://github.com/stephenslab/mashr/blob/master/R/get_functions.R)
# default is factor = 0.5, lfsr_thresh = 0.05 (within 2-fold factor difference)
cat("Computing pairwise sharing...\n")
pairwise_sharing <- get_pairwise_sharing(m2, factor = 0.5, lfsr_thresh = 0.05)
write_tsv(as.data.frame(pairwise_sharing), paste0(output_dir, "pairwise_sharing.tsv"))

# Save mash object
saveRDS(m2, paste0(output_dir, "m2.rds"))

# Generate and save meta plot for the first significant effect
if (length(significant_results) > 0) {
  plot_file <- paste0(output_dir, "mash_meta_plot.pdf")
  pdf(plot_file)
  mash_plot_meta(m2, significant_results[1])
  dev.off()
}

# Print summary
cat("Results saved in:", output_dir, "\n")


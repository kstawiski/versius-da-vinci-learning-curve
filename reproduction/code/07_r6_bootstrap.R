# 07_r6_bootstrap.R
# R6: Bootstrap interval for n* of the Versius operative-time curve
# Case-resampling with at least 500 replicates, 1 stochastic imputation per replicate
suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
  library(parallel)
})

Sys.setenv(TMPDIR = "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/tmp")
restricted_dir <- "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro"

imputed_datasets <- readRDS(file.path(restricted_dir, "imputed_datasets.rds"))
M <- length(imputed_datasets)
N_cases_v <- 337
B <- 600 # >= 500 replicates
num_cores <- min(10, parallel::detectCores())

cat("Starting R6 bootstrap with B =", B, "replicates on", num_cores, "cores...\n")

find_n_star <- function(curve, asymp, band) {
  N <- length(curve)
  for (n in 1:N) {
    if (all(abs(curve[n:N] - asymp) <= band)) {
      return(n)
    }
  }
  return(as.integer(N))
}

grid_df <- data.frame(
  platform_n = 1:N_cases_v,
  age = 0, bmi = 0, log_psa = 0, biopsy_isup_cat = "1",
  log_weight = 0, pt_group = "pT2", nerve_sparing = 0, plnd = 0,
  centre = "Bełchatów", first_case = 0
)

# Function to run single replicate
run_rep <- function(b) {
  set.seed(20260925 + b)
  # Pick stochastic imputation
  m_idx <- sample(1:M, 1)
  df_m <- imputed_datasets[[m_idx]]
  v_all <- df_m %>% filter(platform == "Versius")
  
  # Resample cases with replacement
  boot_idx <- sample(1:nrow(v_all), size = nrow(v_all), replace = TRUE)
  boot_all <- v_all[boot_idx, ]
  boot_obs <- boot_all %>% filter(!is.na(or_time))
  
  # Fit GAM
  fit <- tryCatch({
    gam(
      log(or_time) ~ s(platform_n, k = 10) + age + bmi + log_psa + 
        biopsy_isup_cat + log_weight + pt_group + nerve_sparing + plnd + centre + first_case,
      data = boot_obs,
      method = "REML"
    )
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(NULL)
  
  res <- log(boot_obs$or_time) - predict(fit, newdata = boot_obs)
  smear <- mean(exp(res))
  
  terms_all <- tryCatch(predict(fit, newdata = boot_all, type = "terms"), error = function(e) NULL)
  if (is.null(terms_all)) return(NULL)
  
  c_int <- attr(terms_all, "constant")
  cov_terms <- rowSums(terms_all[, colnames(terms_all) != "s(platform_n)", drop = FALSE]) + c_int
  mean_exp_cov <- mean(exp(cov_terms))
  
  s_pred <- tryCatch(predict(fit, newdata = grid_df, type = "terms", terms = "s(platform_n)"), error = function(e) NULL)
  if (is.null(s_pred)) return(NULL)
  
  smooth_vec <- as.vector(s_pred)
  curve_time <- exp(smooth_vec) * mean_exp_cov * smear
  
  last50_idx <- (N_cases_v - 49):N_cases_v
  asymp <- mean(curve_time[last50_idx])
  
  n15 <- find_n_star(curve_time, asymp, 15)
  n10 <- find_n_star(curve_time, asymp, 10)
  
  c(n15 = n15, n10 = n10, asymp = asymp)
}

res_list <- mclapply(1:B, run_rep, mc.cores = num_cores)

# Filter valid results
valid_res <- do.call(rbind, res_list[!sapply(res_list, is.null)])
n_valid <- nrow(valid_res)
cat("Completed replicates:", n_valid, "/", B, "\n")

n15_vec <- valid_res[, "n15"]
n10_vec <- valid_res[, "n10"]

# Percentile intervals
ci95_n15 <- quantile(n15_vec, probs = c(0.025, 0.975))
ci90_n15 <- quantile(n15_vec, probs = c(0.05, 0.95))
ci95_n10 <- quantile(n10_vec, probs = c(0.025, 0.975))
ci90_n10 <- quantile(n10_vec, probs = c(0.05, 0.95))

cat("\n=== R6: Bootstrap Results for Versius n* (Operative Time) ===\n")
cat(sprintf("n* (15-minute band): Median = %3d, Mean = %5.1f, IQR = [%d, %d]\n",
            as.integer(median(n15_vec)), mean(n15_vec), as.integer(quantile(n15_vec, 0.25)), as.integer(quantile(n15_vec, 0.75))))
cat(sprintf("  95%% CI (percentile): [%3d, %3d]\n", as.integer(ci95_n15[1]), as.integer(ci95_n15[2])))
cat(sprintf("  90%% CI (percentile): [%3d, %3d]\n", as.integer(ci90_n15[1]), as.integer(ci90_n15[2])))

cat(sprintf("\nn* (10-minute band): Median = %3d, Mean = %5.1f, IQR = [%d, %d]\n",
            as.integer(median(n10_vec)), mean(n10_vec), as.integer(quantile(n10_vec, 0.25)), as.integer(quantile(n10_vec, 0.75))))
cat(sprintf("  95%% CI (percentile): [%3d, %3d]\n", as.integer(ci95_n10[1]), as.integer(ci95_n10[2])))
cat(sprintf("  90%% CI (percentile): [%3d, %3d]\n", as.integer(ci90_n10[1]), as.integer(ci90_n10[2])))

r6_results <- list(
  n_replicates = n_valid,
  n_star_15 = list(
    median = as.numeric(median(n15_vec)),
    mean = as.numeric(mean(n15_vec)),
    ci_95 = as.numeric(ci95_n15),
    ci_90 = as.numeric(ci90_n15),
    distribution_quantiles = as.list(quantile(n15_vec, probs = seq(0, 1, 0.1)))
  ),
  n_star_10 = list(
    median = as.numeric(median(n10_vec)),
    mean = as.numeric(mean(n10_vec)),
    ci_95 = as.numeric(ci95_n10),
    ci_90 = as.numeric(ci90_n10),
    distribution_quantiles = as.list(quantile(n10_vec, probs = seq(0, 1, 0.1)))
  )
)

saveRDS(r6_results, file.path(restricted_dir, "r6_results.rds"))
Sys.chmod(file.path(restricted_dir, "r6_results.rds"), mode = "0600")
cat("\nR6 analysis completed successfully.\n")

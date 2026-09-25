# 03_r2_versius_curves.R
# R2: Versius learning curve (operative time and pT2 margin)
suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
})

Sys.setenv(TMPDIR = "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/tmp")
restricted_dir <- "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro"

imputed_datasets <- readRDS(file.path(restricted_dir, "imputed_datasets.rds"))
M <- length(imputed_datasets)
N_cases_v <- 337

curves_time <- matrix(0, nrow = M, ncol = N_cases_v)
curves_psm <- matrix(0, nrow = M, ncol = N_cases_v)

cat("Fitting Versius GAMs across", M, "imputations...\n")

grid_df <- data.frame(
  platform_n = 1:N_cases_v,
  age = 0, bmi = 0, log_psa = 0, biopsy_isup_cat = "1",
  log_weight = 0, pt_group = "pT2", nerve_sparing = 0, plnd = 0,
  centre = "Bełchatów", first_case = 0
)

for (m in 1:M) {
  df_m <- imputed_datasets[[m]]
  v_all <- df_m %>% filter(platform == "Versius")
  
  # 1. Operative time (Gaussian on log scale, k=10, REML)
  v_time_obs <- v_all %>% filter(!is.na(or_time))
  
  fit_time <- gam(
    log(or_time) ~ s(platform_n, k = 10) + age + bmi + log_psa + 
      biopsy_isup_cat + log_weight + pt_group + nerve_sparing + plnd + centre + first_case,
    data = v_time_obs,
    method = "REML"
  )
  
  # Smearing factor
  res <- log(v_time_obs$or_time) - predict(fit_time, newdata = v_time_obs)
  smear <- mean(exp(res))
  
  # Terms decomposition for all 337 Versius cases:
  # All non-smooth terms evaluated at each patient's actual covariates
  terms_all <- predict(fit_time, newdata = v_all, type = "terms")
  c_int <- attr(terms_all, "constant")
  cov_terms <- rowSums(terms_all[, colnames(terms_all) != "s(platform_n)", drop = FALSE]) + c_int
  mean_exp_cov <- mean(exp(cov_terms))
  
  # Smooth term evaluated at grid 1..337
  s_pred <- predict(fit_time, newdata = grid_df, type = "terms", terms = "s(platform_n)")
  smooth_vec <- as.vector(s_pred)
  
  # Standardized curve: mean_i(exp(cov_terms[i] + s[n])) * smear = exp(s[n]) * mean_exp_cov * smear
  curves_time[m, ] <- exp(smooth_vec) * mean_exp_cov * smear
  
  # 2. pT2 positive margin (Binomial, k=5, REML)
  v_pt2_all <- v_all %>% filter(pt_group == "pT2")
  v_pt2_obs <- v_pt2_all %>% filter(!is.na(psm))
  
  fit_psm <- gam(
    psm ~ s(platform_n, k = 5) + age + bmi + log_psa + 
      biopsy_isup_cat + log_weight + centre,
    data = v_pt2_obs,
    family = binomial,
    method = "REML"
  )
  
  terms_pt2 <- predict(fit_psm, newdata = v_pt2_all, type = "terms")
  c_int_psm <- attr(terms_pt2, "constant")
  cov_terms_psm <- rowSums(terms_pt2[, colnames(terms_pt2) != "s(platform_n)", drop = FALSE]) + c_int_psm
  
  s_pred_psm <- predict(fit_psm, newdata = grid_df, type = "terms", terms = "s(platform_n)")
  smooth_vec_psm <- as.vector(s_pred_psm)
  
  # Marginal prediction: mean over all 190 pT2 patients of plogis(cov_terms[i] + s[n])
  eta_mat <- outer(cov_terms_psm, smooth_vec_psm, "+")
  prob_mat <- plogis(eta_mat)
  curves_psm[m, ] <- colMeans(prob_mat)
}

# Average curves over imputations
mean_curve_time <- colMeans(curves_time)
mean_curve_psm <- colMeans(curves_psm)

# Asymptote: mean over last 50 cases (288 to 337)
last50_idx <- (N_cases_v - 49):N_cases_v
asymptote_time <- mean(mean_curve_time[last50_idx])
asymptote_psm <- mean(mean_curve_psm[last50_idx])

# Helper for n*
find_n_star <- function(curve, asymp, band) {
  N <- length(curve)
  for (n in 1:N) {
    if (all(abs(curve[n:N] - asymp) <= band)) {
      return(n)
    }
  }
  return(as.integer(N))
}

n_star_time_15 <- find_n_star(mean_curve_time, asymptote_time, 15)
n_star_time_10 <- find_n_star(mean_curve_time, asymptote_time, 10)
n_star_psm_05  <- find_n_star(mean_curve_psm, asymptote_psm, 0.05)

# Slope over final 100 cases
cases_last100 <- 238:337
slope_final100_time <- coef(lm(mean_curve_time[cases_last100] ~ cases_last100))[2]

milestones <- c(1, 25, 50, 100, N_cases_v)

cat("\n=== R2: Versius Operative Time Learning Curve ===\n")
cat("Asymptote (mean over cases 288-337):", round(asymptote_time, 2), "minutes\n")
cat("n* (15-minute band):", n_star_time_15, "\n")
cat("n* (10-minute band):", n_star_time_10, "\n")
cat("Slope over final 100 cases:", round(slope_final100_time, 4), "min/case\n")
cat("Curve values at milestones:\n")
for (ms in milestones) {
  cat(sprintf("  Case %3d: %6.2f min\n", ms, mean_curve_time[ms]))
}

cat("\n=== R2: Versius pT2 Margin Learning Curve ===\n")
cat("Asymptote (mean over cases 288-337):", round(asymptote_psm * 100, 2), "%\n")
cat("n* (5 percentage points band):", n_star_psm_05, "\n")
cat("Curve values at milestones:\n")
for (ms in milestones) {
  cat(sprintf("  Case %3d: %6.2f%%\n", ms, mean_curve_psm[ms] * 100))
}

r2_results <- list(
  mean_curve_time = mean_curve_time,
  asymptote_time = asymptote_time,
  n_star_time_15 = n_star_time_15,
  n_star_time_10 = n_star_time_10,
  slope_final100_time = unname(slope_final100_time),
  time_milestones = setNames(mean_curve_time[milestones], paste0("case_", milestones)),
  mean_curve_psm = mean_curve_psm,
  asymptote_psm = asymptote_psm,
  n_star_psm_05 = n_star_psm_05,
  psm_milestones = setNames(mean_curve_psm[milestones], paste0("case_", milestones))
)

saveRDS(r2_results, file.path(restricted_dir, "r2_results.rds"))
Sys.chmod(file.path(restricted_dir, "r2_results.rds"), mode = "0600")
cat("\nR2 analysis completed successfully.\n")

# 04_r3_davinci_curves.R
# R3: da Vinci restart learning curve (operative time and pT2 margin)
suppressPackageStartupMessages({
  library(dplyr)
  library(mgcv)
})

Sys.setenv(TMPDIR = "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/tmp")
restricted_dir <- "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro"

imputed_datasets <- readRDS(file.path(restricted_dir, "imputed_datasets.rds"))
M <- length(imputed_datasets)
N_cases_dv <- 382

curves_time <- matrix(0, nrow = M, ncol = N_cases_dv)
curves_psm <- matrix(0, nrow = M, ncol = N_cases_dv)

cat("Fitting da Vinci restart GAMs across", M, "imputations...\n")

grid_df <- data.frame(
  dv_restart_n = 1:N_cases_dv,
  age = 0, bmi = 0, log_psa = 0, biopsy_isup_cat = "1",
  log_weight = 0, pt_group = "pT2", nerve_sparing = 0, plnd = 0,
  centre = "Bełchatów", first_case = 0
)

for (m in 1:M) {
  df_m <- imputed_datasets[[m]]
  dv_all <- df_m %>% filter(platform == "da Vinci" & !is.na(dv_restart_n))
  stopifnot(nrow(dv_all) == N_cases_dv)
  
  # 1. Operative time (Gaussian on log scale, k=10, REML)
  dv_time_obs <- dv_all %>% filter(!is.na(or_time))
  
  fit_time <- gam(
    log(or_time) ~ s(dv_restart_n, k = 10) + age + bmi + log_psa + 
      biopsy_isup_cat + log_weight + pt_group + nerve_sparing + plnd + centre + first_case,
    data = dv_time_obs,
    method = "REML"
  )
  
  # Smearing factor
  res <- log(dv_time_obs$or_time) - predict(fit_time, newdata = dv_time_obs)
  smear <- mean(exp(res))
  
  # Terms decomposition for all 382 restart cases:
  terms_all <- predict(fit_time, newdata = dv_all, type = "terms")
  c_int <- attr(terms_all, "constant")
  cov_terms <- rowSums(terms_all[, colnames(terms_all) != "s(dv_restart_n)", drop = FALSE]) + c_int
  mean_exp_cov <- mean(exp(cov_terms))
  
  # Smooth term evaluated at grid 1..382
  s_pred <- predict(fit_time, newdata = grid_df, type = "terms", terms = "s(dv_restart_n)")
  smooth_vec <- as.vector(s_pred)
  
  curves_time[m, ] <- exp(smooth_vec) * mean_exp_cov * smear
  
  # 2. pT2 positive margin (Binomial, k=5, REML)
  dv_pt2_all <- dv_all %>% filter(pt_group == "pT2")
  dv_pt2_obs <- dv_pt2_all %>% filter(!is.na(psm))
  
  fit_psm <- gam(
    psm ~ s(dv_restart_n, k = 5) + age + bmi + log_psa + 
      biopsy_isup_cat + log_weight + centre,
    data = dv_pt2_obs,
    family = binomial,
    method = "REML"
  )
  
  terms_pt2 <- predict(fit_psm, newdata = dv_pt2_all, type = "terms")
  c_int_psm <- attr(terms_pt2, "constant")
  cov_terms_psm <- rowSums(terms_pt2[, colnames(terms_pt2) != "s(dv_restart_n)", drop = FALSE]) + c_int_psm
  
  s_pred_psm <- predict(fit_psm, newdata = grid_df, type = "terms", terms = "s(dv_restart_n)")
  smooth_vec_psm <- as.vector(s_pred_psm)
  
  eta_mat <- outer(cov_terms_psm, smooth_vec_psm, "+")
  prob_mat <- plogis(eta_mat)
  curves_psm[m, ] <- colMeans(prob_mat)
}

# Average curves over imputations
mean_curve_time <- colMeans(curves_time)
mean_curve_psm <- colMeans(curves_psm)

# Asymptote: mean over last 50 cases (333 to 382)
last50_idx <- (N_cases_dv - 49):N_cases_dv
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

# Slope over final 100 cases (283 to 382)
cases_last100 <- (N_cases_dv - 99):N_cases_dv
slope_final100_time <- coef(lm(mean_curve_time[cases_last100] ~ cases_last100))[2]

milestones <- c(1, 25, 50, 100, N_cases_dv)

cat("\n=== R3: da Vinci Restart Operative Time Learning Curve ===\n")
cat("Asymptote (mean over cases 333-382):", round(asymptote_time, 2), "minutes\n")
cat("n* (15-minute band):", n_star_time_15, "\n")
cat("n* (10-minute band):", n_star_time_10, "\n")
cat("Slope over final 100 cases:", round(slope_final100_time, 4), "min/case\n")
cat("Curve values at milestones:\n")
for (ms in milestones) {
  cat(sprintf("  Case %3d: %6.2f min\n", ms, mean_curve_time[ms]))
}

cat("\n=== R3: da Vinci Restart pT2 Margin Learning Curve ===\n")
cat("Asymptote (mean over cases 333-382):", round(asymptote_psm * 100, 2), "%\n")
cat("n* (5 percentage points band):", n_star_psm_05, "\n")
cat("Curve values at milestones:\n")
for (ms in milestones) {
  cat(sprintf("  Case %3d: %6.2f%%\n", ms, mean_curve_psm[ms] * 100))
}

r3_results <- list(
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

saveRDS(r3_results, file.path(restricted_dir, "r3_results.rds"))
Sys.chmod(file.path(restricted_dir, "r3_results.rds"), mode = "0600")
cat("\nR3 analysis completed successfully.\n")

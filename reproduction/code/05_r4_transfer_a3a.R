# 05_r4_transfer_a3a.R
# R4: A3a transfer of proficiency - restart da Vinci 1-100 vs Versius 1-100
suppressPackageStartupMessages({
  library(dplyr)
  library(sandwich)
})

Sys.setenv(TMPDIR = "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/tmp")
restricted_dir <- "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro"

imputed_datasets <- readRDS(file.path(restricted_dir, "imputed_datasets.rds"))
M <- length(imputed_datasets)

# Storage across imputations
# 1. Operative time (linear regression on minutes, HC3)
q_time <- numeric(M)
u_time <- numeric(M)

# 2. pT2 margin risk difference (logistic regression, standardized risk difference)
q_psm <- numeric(M)
u_psm <- numeric(M)

cat("Analyzing A3a (restart da Vinci 1-100 vs Versius 1-100) across", M, "imputations...\n")

for (m in 1:M) {
  df_m <- imputed_datasets[[m]]
  
  # Filter cases 1-100 for each series
  sub_100 <- df_m %>%
    filter(
      (platform == "Versius" & platform_n <= 100) |
      (platform == "da Vinci" & !is.na(dv_restart_n) & dv_restart_n <= 100)
    ) %>%
    mutate(
      dv_group = ifelse(platform == "da Vinci", 1, 0)
    )
  
  # --- 1. Operative time ---
  sub_time <- sub_100 %>% filter(!is.na(or_time))
  
  fit_time <- lm(
    or_time ~ dv_group + age + bmi + log_psa + biopsy_isup_cat + 
      log_weight + pt_group + nerve_sparing + plnd + centre + first_case,
    data = sub_time
  )
  
  vcov_time <- tryCatch(
    sandwich::vcovHC(fit_time, type = "HC3"),
    error = function(e) sandwich::vcovHC(fit_time, type = "HC1")
  )
  if (any(is.nan(vcov_time)) || any(is.na(vcov_time))) {
    vcov_time <- sandwich::vcovHC(fit_time, type = "HC1")
  }
  
  q_time[m] <- coef(fit_time)["dv_group"]
  u_time[m] <- vcov_time["dv_group", "dv_group"]
  
  # --- 2. pT2 Margin risk difference ---
  sub_pt2 <- sub_100 %>% filter(pt_group == "pT2")
  sub_pt2_obs <- sub_pt2 %>% filter(!is.na(psm))
  
  fit_psm <- glm(
    psm ~ dv_group + age + bmi + log_psa + biopsy_isup_cat + 
      log_weight + centre,
    data = sub_pt2_obs,
    family = binomial
  )
  
  vcov_psm <- tryCatch(
    sandwich::vcovHC(fit_psm, type = "HC0"),
    error = function(e) vcov(fit_psm)
  )
  if (any(is.nan(vcov_psm)) || any(is.na(vcov_psm))) {
    vcov_psm <- vcov(fit_psm)
  }
  
  # Marginally standardized risk difference over all pT2 cases in sub_pt2
  # Construct design matrices with dv_group = 1 and dv_group = 0
  nd_dv <- sub_pt2
  nd_dv$dv_group <- 1
  nd_v <- sub_pt2
  nd_v$dv_group <- 0
  
  X_dv <- model.matrix(delete.response(terms(fit_psm)), data = nd_dv)
  X_v  <- model.matrix(delete.response(terms(fit_psm)), data = nd_v)
  
  b <- coef(fit_psm)
  p_dv <- plogis(as.vector(X_dv %*% b))
  p_v  <- plogis(as.vector(X_v %*% b))
  
  rd_m <- mean(p_dv - p_v)
  q_psm[m] <- rd_m
  
  # Delta-method gradient
  # d(p_dv)/db = p_dv * (1 - p_dv) * X_dv
  # d(p_v)/db  = p_v  * (1 - p_v)  * X_v
  grad_dv <- colMeans(p_dv * (1 - p_dv) * X_dv)
  grad_v  <- colMeans(p_v  * (1 - p_v)  * X_v)
  grad <- grad_dv - grad_v
  
  u_psm[m] <- as.numeric(t(grad) %*% vcov_psm %*% grad)
}

# Pool with Rubin's rules
pool_rubin <- function(q, u) {
  M <- length(q)
  q_bar <- mean(q)
  u_bar <- mean(u)
  b <- var(q)
  t_var <- u_bar + (1 + 1/M) * b
  se <- sqrt(t_var)
  
  # Degrees of freedom (Barnard-Rubin)
  df <- if (b > 0) (M - 1) * (1 + u_bar / ((1 + 1/M) * b))^2 else 9999
  
  ci_95 <- c(q_bar - qt(0.975, df) * se, q_bar + qt(0.975, df) * se)
  ci_90 <- c(q_bar - qt(0.95, df) * se, q_bar + qt(0.95, df) * se)
  p_val <- 2 * pt(-abs(q_bar / se), df)
  
  list(
    estimate = q_bar,
    se = se,
    df = df,
    ci_95 = ci_95,
    ci_90 = ci_90,
    p_value = p_val
  )
}

res_time <- pool_rubin(q_time, u_time)
res_psm  <- pool_rubin(q_psm, u_psm)

cat("\n=== R4: A3a Operative Time Difference (Restart da Vinci 1-100 vs Versius 1-100) ===\n")
cat(sprintf("Adjusted difference: %6.2f min (SE: %5.2f)\n", res_time$estimate, res_time$se))
cat(sprintf("95%% CI: [%6.2f, %6.2f]\n", res_time$ci_95[1], res_time$ci_95[2]))
cat(sprintf("90%% CI: [%6.2f, %6.2f]\n", res_time$ci_90[1], res_time$ci_90[2]))
cat(sprintf("p-value: %6.4f\n", res_time$p_value))

cat("\n=== R4: A3a pT2 Margin Risk Difference (Restart da Vinci 1-100 vs Versius 1-100) ===\n")
cat(sprintf("Adjusted risk difference: %6.2f%% (SE: %5.2f%%)\n", res_psm$estimate * 100, res_psm$se * 100))
cat(sprintf("95%% CI: [%6.2f%%, %6.2f%%]\n", res_psm$ci_95[1] * 100, res_psm$ci_95[2] * 100))
cat(sprintf("90%% CI: [%6.2f%%, %6.2f%%]\n", res_psm$ci_90[1] * 100, res_psm$ci_90[2] * 100))
cat(sprintf("p-value: %6.4f\n", res_psm$p_value))

r4_results <- list(
  operative_time = res_time,
  pt2_margin = res_psm
)

saveRDS(r4_results, file.path(restricted_dir, "r4_results.rds"))
Sys.chmod(file.path(restricted_dir, "r4_results.rds"), mode = "0600")
cat("\nR4 analysis completed successfully.\n")

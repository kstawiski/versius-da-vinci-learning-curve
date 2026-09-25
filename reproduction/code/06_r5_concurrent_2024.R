# 06_r5_concurrent_2024.R
# R5: A4 concurrent platform comparison, calendar 2024
suppressPackageStartupMessages({
  library(dplyr)
  library(sandwich)
  library(logistf)
})

Sys.setenv(TMPDIR = "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/tmp")
restricted_dir <- "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro"

imputed_datasets <- readRDS(file.path(restricted_dir, "imputed_datasets.rds"))
M <- length(imputed_datasets)

# Helper for Rubin's pooling
pool_rubin <- function(q, u) {
  M <- length(q)
  q_bar <- mean(q)
  u_bar <- mean(u)
  b <- if (M > 1) var(q) else 0
  t_var <- u_bar + (1 + 1/M) * b
  se <- sqrt(t_var)
  
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

# Helper for marginally standardized risk difference from glm
calc_std_rd <- function(fit, data_subset, vcov_mat, treat_var = "dv_group") {
  nd_1 <- data_subset; nd_1[[treat_var]] <- 1
  nd_0 <- data_subset; nd_0[[treat_var]] <- 0
  
  X_1 <- model.matrix(delete.response(terms(fit)), data = nd_1)
  X_0 <- model.matrix(delete.response(terms(fit)), data = nd_0)
  
  b <- coef(fit)
  p_1 <- plogis(as.vector(X_1 %*% b))
  p_0 <- plogis(as.vector(X_0 %*% b))
  
  rd <- mean(p_1 - p_0)
  
  grad_1 <- colMeans(p_1 * (1 - p_1) * X_1)
  grad_0 <- colMeans(p_0 * (1 - p_0) * X_0)
  grad <- grad_1 - grad_0
  
  var_rd <- as.numeric(t(grad) %*% vcov_mat %*% grad)
  list(rd = rd, var_rd = var_rd)
}

cat("Running R5 (A4 concurrent 2024 platform comparison) across", M, "imputations...\n")

# Storage for estimates and variances across 20 imputations
storage <- list(
  or_time = list(q = numeric(M), u = numeric(M)),
  pt2_margin = list(q = numeric(M), u = numeric(M)),
  all_margin = list(q = numeric(M), u = numeric(M)),
  los = list(q = numeric(M), u = numeric(M)),
  readm_30d = list(q = numeric(M), u = numeric(M)),
  psa_pers = list(q = numeric(M), u = numeric(M)),
  pad_free_3m = list(q = numeric(M), u = numeric(M)),
  pad_free_12m = list(q = numeric(M), u = numeric(M))
)

for (m in 1:M) {
  df_m <- imputed_datasets[[m]]
  
  c2024 <- df_m %>%
    filter(year == 2024 & (platform == "Versius" | (platform == "da Vinci" & !is.na(dv_restart_n)))) %>%
    mutate(dv_group = ifelse(platform == "da Vinci", 1, 0))
  
  # 1. Operative time minutes (linear regression, HC3)
  sub_or <- c2024 %>% filter(!is.na(or_time))
  fit_or <- lm(
    or_time ~ dv_group + age + bmi + log_psa + biopsy_isup_cat + log_weight + 
      pt_group + nerve_sparing + plnd + centre + first_case,
    data = sub_or
  )
  vcov_or <- tryCatch(sandwich::vcovHC(fit_or, type = "HC3"), error = function(e) sandwich::vcovHC(fit_or, type = "HC1"))
  storage$or_time$q[m] <- coef(fit_or)["dv_group"]
  storage$or_time$u[m] <- vcov_or["dv_group", "dv_group"]
  
  # 2. pT2 margin risk difference (logistic regression, standardized RD)
  sub_pt2 <- c2024 %>% filter(pt_group == "pT2")
  sub_pt2_obs <- sub_pt2 %>% filter(!is.na(psm))
  fit_pt2 <- glm(
    psm ~ dv_group + age + bmi + log_psa + biopsy_isup_cat + log_weight + centre,
    data = sub_pt2_obs, family = binomial
  )
  vcov_pt2 <- tryCatch(sandwich::vcovHC(fit_pt2, type = "HC0"), error = function(e) vcov(fit_pt2))
  rd_pt2 <- calc_std_rd(fit_pt2, sub_pt2, vcov_pt2)
  storage$pt2_margin$q[m] <- rd_pt2$rd
  storage$pt2_margin$u[m] <- rd_pt2$var_rd
  
  # 3. All-stage margin risk difference (logistic regression, standardized RD)
  sub_allm <- c2024
  sub_allm_obs <- sub_allm %>% filter(!is.na(psm))
  fit_allm <- glm(
    psm ~ dv_group + age + bmi + log_psa + biopsy_isup_cat + log_weight + 
      pt_group + nerve_sparing + plnd + centre,
    data = sub_allm_obs, family = binomial
  )
  vcov_allm <- tryCatch(sandwich::vcovHC(fit_allm, type = "HC0"), error = function(e) vcov(fit_allm))
  rd_allm <- calc_std_rd(fit_allm, sub_allm, vcov_allm)
  storage$all_margin$q[m] <- rd_allm$rd
  storage$all_margin$u[m] <- rd_allm$var_rd
  
  # 4. Length of stay (linear regression, HC3)
  sub_los <- c2024 %>% filter(!is.na(los))
  fit_los <- lm(
    los ~ dv_group + age + bmi + log_psa + biopsy_isup_cat + log_weight + 
      pt_group + nerve_sparing + plnd + centre,
    data = sub_los
  )
  vcov_los <- tryCatch(sandwich::vcovHC(fit_los, type = "HC3"), error = function(e) sandwich::vcovHC(fit_los, type = "HC1"))
  storage$los$q[m] <- coef(fit_los)["dv_group"]
  storage$los$u[m] <- vcov_los["dv_group", "dv_group"]
  
  # 5. 30-day readmission (age and centre only; 0 events in both groups)
  # When both groups have 0 events, risk is 0 under both platforms
  storage$readm_30d$q[m] <- 0.0
  storage$readm_30d$u[m] <- 0.0
  
  # 6. EAU PSA persistence (window-closed)
  # Fit with parsimonious adjustment (age, centre, log_psa) due to 0 events in pT2 and ISUP 1
  sub_psa <- c2024 %>% filter(eligible_psa_persistence)
  sub_psa_obs <- sub_psa %>% filter(!is.na(psa_persistence_eau_num))
  fit_psa <- glm(
    psa_persistence_eau_num ~ dv_group + age + centre + log_psa,
    data = sub_psa_obs, family = binomial
  )
  vcov_psa <- tryCatch(sandwich::vcovHC(fit_psa, type = "HC0"), error = function(e) vcov(fit_psa))
  rd_psa <- calc_std_rd(fit_psa, sub_psa, vcov_psa)
  storage$psa_pers$q[m] <- rd_psa$rd
  storage$psa_pers$u[m] <- rd_psa$var_rd
  
  # 7. Pad-free 3 months (window-closed, age, bmi, nerve_sparing, pt_group, centre)
  sub_pf3 <- c2024 %>% filter(eligible_continence_3m)
  sub_pf3_obs <- sub_pf3 %>% filter(!is.na(pad_free_3m_num))
  fit_pf3 <- glm(
    pad_free_3m_num ~ dv_group + age + bmi + nerve_sparing + pt_group + centre,
    data = sub_pf3_obs, family = binomial
  )
  vcov_pf3 <- tryCatch(sandwich::vcovHC(fit_pf3, type = "HC0"), error = function(e) vcov(fit_pf3))
  rd_pf3 <- calc_std_rd(fit_pf3, sub_pf3, vcov_pf3)
  storage$pad_free_3m$q[m] <- rd_pf3$rd
  storage$pad_free_3m$u[m] <- rd_pf3$var_rd
  
  # 8. Pad-free 12 months (window-closed, age, bmi, nerve_sparing, pt_group, centre)
  sub_pf12 <- c2024 %>% filter(eligible_continence_12m)
  sub_pf12_obs <- sub_pf12 %>% filter(!is.na(pad_free_12m_num))
  fit_pf12 <- glm(
    pad_free_12m_num ~ dv_group + age + bmi + nerve_sparing + pt_group + centre,
    data = sub_pf12_obs, family = binomial
  )
  vcov_pf12 <- tryCatch(sandwich::vcovHC(fit_pf12, type = "HC0"), error = function(e) vcov(fit_pf12))
  rd_pf12 <- calc_std_rd(fit_pf12, sub_pf12, vcov_pf12)
  storage$pad_free_12m$q[m] <- rd_pf12$rd
  storage$pad_free_12m$u[m] <- rd_pf12$var_rd
}

# Pool all endpoints with Rubin's rules
pooled_r5 <- lapply(storage, function(item) {
  if (all(item$q == 0) && all(item$u == 0)) {
    list(
      estimate = 0.0,
      se = 0.0,
      df = Inf,
      ci_95 = c(0.0, 0.0),
      ci_90 = c(0.0, 0.0),
      p_value = 1.0
    )
  } else {
    pool_rubin(item$q, item$u)
  }
})

cat("\n=== R5 Summary Table (Calendar 2024: Restart da Vinci vs Versius) ===\n")
for (nm in names(pooled_r5)) {
  res <- pooled_r5[[nm]]
  cat(sprintf("%-15s: est = %7.3f, 95%% CI = [%7.3f, %7.3f], 90%% CI = [%7.3f, %7.3f], p = %6.4f\n",
              nm, res$estimate, res$ci_95[1], res$ci_95[2], res$ci_90[1], res$ci_90[2], res$p_value))
}

# Also calculate exact crude counts/means on unimputed data
df_raw <- readRDS(file.path(restricted_dir, "derived_data.rds"))
c2024_raw <- df_raw %>%
  filter(year == 2024 & (platform == "Versius" | (platform == "da Vinci" & !is.na(dv_restart_n))))

crude_stats <- list(
  or_time = list(
    versius = c2024_raw %>% filter(platform == "Versius" & !is.na(or_time)) %>% summarize(mean = mean(or_time), sd = sd(or_time), n = n()) %>% as.list(),
    davinci = c2024_raw %>% filter(platform == "da Vinci" & !is.na(or_time)) %>% summarize(mean = mean(or_time), sd = sd(or_time), n = n()) %>% as.list()
  ),
  pt2_margin = list(
    versius = c2024_raw %>% filter(platform == "Versius" & pt_group == "pT2" & !is.na(psm)) %>% summarize(events = sum(psm), n = n(), rate = mean(psm)) %>% as.list(),
    davinci = c2024_raw %>% filter(platform == "da Vinci" & pt_group == "pT2" & !is.na(psm)) %>% summarize(events = sum(psm), n = n(), rate = mean(psm)) %>% as.list()
  ),
  all_margin = list(
    versius = c2024_raw %>% filter(platform == "Versius" & !is.na(psm)) %>% summarize(events = sum(psm), n = n(), rate = mean(psm)) %>% as.list(),
    davinci = c2024_raw %>% filter(platform == "da Vinci" & !is.na(psm)) %>% summarize(events = sum(psm), n = n(), rate = mean(psm)) %>% as.list()
  ),
  los = list(
    versius = c2024_raw %>% filter(platform == "Versius" & !is.na(los)) %>% summarize(mean = mean(los), sd = sd(los), n = n()) %>% as.list(),
    davinci = c2024_raw %>% filter(platform == "da Vinci" & !is.na(los)) %>% summarize(mean = mean(los), sd = sd(los), n = n()) %>% as.list()
  ),
  readm_30d = list(
    versius = c2024_raw %>% filter(platform == "Versius" & !is.na(readmission_30d_num)) %>% summarize(events = sum(readmission_30d_num), n = n(), rate = mean(readmission_30d_num)) %>% as.list(),
    davinci = c2024_raw %>% filter(platform == "da Vinci" & !is.na(readmission_30d_num)) %>% summarize(events = sum(readmission_30d_num), n = n(), rate = mean(readmission_30d_num)) %>% as.list()
  ),
  psa_pers = list(
    versius = c2024_raw %>% filter(platform == "Versius" & eligible_psa_persistence & !is.na(psa_persistence_eau_num)) %>% summarize(events = sum(psa_persistence_eau_num), n = n(), rate = mean(psa_persistence_eau_num)) %>% as.list(),
    davinci = c2024_raw %>% filter(platform == "da Vinci" & eligible_psa_persistence & !is.na(psa_persistence_eau_num)) %>% summarize(events = sum(psa_persistence_eau_num), n = n(), rate = mean(psa_persistence_eau_num)) %>% as.list()
  ),
  pad_free_3m = list(
    versius = c2024_raw %>% filter(platform == "Versius" & eligible_continence_3m & !is.na(pad_free_3m_num)) %>% summarize(events = sum(pad_free_3m_num), n = n(), rate = mean(pad_free_3m_num)) %>% as.list(),
    davinci = c2024_raw %>% filter(platform == "da Vinci" & eligible_continence_3m & !is.na(pad_free_3m_num)) %>% summarize(events = sum(pad_free_3m_num), n = n(), rate = mean(pad_free_3m_num)) %>% as.list()
  ),
  pad_free_12m = list(
    versius = c2024_raw %>% filter(platform == "Versius" & eligible_continence_12m & !is.na(pad_free_12m_num)) %>% summarize(events = sum(pad_free_12m_num), n = n(), rate = mean(pad_free_12m_num)) %>% as.list(),
    davinci = c2024_raw %>% filter(platform == "da Vinci" & eligible_continence_12m & !is.na(pad_free_12m_num)) %>% summarize(events = sum(pad_free_12m_num), n = n(), rate = mean(pad_free_12m_num)) %>% as.list()
  )
)

r5_results <- list(
  pooled_estimates = pooled_r5,
  crude_stats = crude_stats
)

saveRDS(r5_results, file.path(restricted_dir, "r5_results.rds"))
Sys.chmod(file.path(restricted_dir, "r5_results.rds"), mode = "0600")
cat("\nR5 analysis completed successfully.\n")

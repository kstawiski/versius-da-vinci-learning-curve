# 08_collate_results.R
# Collate all reproduction numbers R1 to R6 into REPRO_RESULTS.json
suppressPackageStartupMessages({
  library(dplyr)
  library(jsonlite)
})

Sys.setenv(TMPDIR = "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/tmp")
restricted_dir <- "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro"
output_json <- "<PROJECT_ROOT>/manuscript6/analysis/repro_results/REPRO_RESULTS.json"

df_raw <- readRDS(file.path(restricted_dir, "derived_data.rds"))
r2_res <- readRDS(file.path(restricted_dir, "r2_results.rds"))
r3_res <- readRDS(file.path(restricted_dir, "r3_results.rds"))
r4_res <- readRDS(file.path(restricted_dir, "r4_results.rds"))
r5_res <- readRDS(file.path(restricted_dir, "r5_results.rds"))
r6_res <- readRDS(file.path(restricted_dir, "r6_results.rds"))

# --- R1: Counts ---
series_defs <- list(
  versius = df_raw$platform == "Versius",
  davinci = df_raw$platform == "da Vinci",
  restart_davinci = df_raw$platform == "da Vinci" & df_raw$platform_n >= 20,
  versius_2024 = df_raw$platform == "Versius" & df_raw$year == 2024,
  restart_davinci_2024 = df_raw$platform == "da Vinci" & df_raw$platform_n >= 20 & df_raw$year == 2024
)

r1_list <- list()
for (nm in names(series_defs)) {
  mask <- series_defs[[nm]]
  n_total <- sum(mask)
  n_pt2 <- sum(mask & df_raw$pt_group == "pT2", na.rm = TRUE)
  n_or_time <- sum(mask & !is.na(df_raw$or_time))
  r1_list[[nm]] <- list(
    total_cases = n_total,
    pt2_cases = n_pt2,
    pt2_denominator = n_total,
    pt2_proportion = n_pt2 / n_total,
    or_time_cases = n_or_time,
    or_time_denominator = n_total,
    or_time_completeness = n_or_time / n_total
  )
}

# --- R2: Versius learning curve ---
r2_list <- list(
  operative_time = list(
    model = "GAM Gaussian on log(or_time), s(platform_n, k=10), REML, smearing retransformation",
    asymptote_minutes = round(r2_res$asymptote_time, 2),
    asymptote_window = "cases 288-337 (final 50 cases)",
    n_star_15min = r2_res$n_star_time_15,
    n_star_10min = r2_res$n_star_time_10,
    slope_final_100_cases_min_per_case = round(r2_res$slope_final100_time, 4),
    standardized_curve_milestones_minutes = lapply(r2_res$time_milestones, function(x) round(x, 2))
  ),
  pt2_positive_margin = list(
    model = "Logistic GAM on psm, s(platform_n, k=5), REML, standardized over pT2 cases",
    asymptote_rate = round(r2_res$asymptote_psm, 4),
    asymptote_percentage = round(r2_res$asymptote_psm * 100, 2),
    asymptote_window = "cases 288-337 (final 50 cases)",
    n_star_5pp = r2_res$n_star_psm_05,
    standardized_curve_milestones_rate = lapply(r2_res$psm_milestones, function(x) round(x, 4)),
    standardized_curve_milestones_percentage = lapply(r2_res$psm_milestones, function(x) round(x * 100, 2))
  )
)

# --- R3: da Vinci restart learning curve ---
r3_list <- list(
  operative_time = list(
    model = "GAM Gaussian on log(or_time), s(dv_restart_n, k=10), REML, smearing retransformation",
    asymptote_minutes = round(r3_res$asymptote_time, 2),
    asymptote_window = "cases 333-382 (final 50 cases)",
    n_star_15min = r3_res$n_star_time_15,
    n_star_10min = r3_res$n_star_time_10,
    slope_final_100_cases_min_per_case = round(r3_res$slope_final100_time, 4),
    standardized_curve_milestones_minutes = lapply(r3_res$time_milestones, function(x) round(x, 2))
  ),
  pt2_positive_margin = list(
    model = "Logistic GAM on psm, s(dv_restart_n, k=5), REML, standardized over pT2 cases",
    asymptote_rate = round(r3_res$asymptote_psm, 4),
    asymptote_percentage = round(r3_res$asymptote_psm * 100, 2),
    asymptote_window = "cases 333-382 (final 50 cases)",
    n_star_5pp = r3_res$n_star_psm_05,
    standardized_curve_milestones_rate = lapply(r3_res$psm_milestones, function(x) round(x, 4)),
    standardized_curve_milestones_percentage = lapply(r3_res$psm_milestones, function(x) round(x * 100, 2))
  )
)

# --- R4: A3a transfer (restart da Vinci 1-100 vs Versius 1-100) ---
r4_list <- list(
  comparison = "da Vinci restart cases 1-100 vs Versius cases 1-100 (reference: Versius)",
  operative_time_minutes = list(
    model = "Linear regression on minutes, HC3 robust errors, pooled across 20 imputations",
    adjusted_difference = round(r4_res$operative_time$estimate, 2),
    se = round(r4_res$operative_time$se, 2),
    ci_95 = round(r4_res$operative_time$ci_95, 2),
    ci_90 = round(r4_res$operative_time$ci_90, 2),
    p_value = signif(r4_res$operative_time$p_value, 4)
  ),
  pt2_margin_risk_difference = list(
    model = "Logistic regression, marginally standardized risk difference, HC0 delta-method, pooled across 20 imputations",
    adjusted_risk_difference = round(r4_res$pt2_margin$estimate, 4),
    adjusted_risk_difference_percentage = round(r4_res$pt2_margin$estimate * 100, 2),
    se = round(r4_res$pt2_margin$se, 4),
    se_percentage = round(r4_res$pt2_margin$se * 100, 2),
    ci_95 = round(r4_res$pt2_margin$ci_95, 4),
    ci_95_percentage = round(r4_res$pt2_margin$ci_95 * 100, 2),
    ci_90 = round(r4_res$pt2_margin$ci_90, 4),
    ci_90_percentage = round(r4_res$pt2_margin$ci_90 * 100, 2),
    p_value = signif(r4_res$pt2_margin$p_value, 4)
  )
)

# --- R5: A4 concurrent 2024 platform comparison ---
endpoints_meta <- list(
  or_time = list(name = "Operative time (minutes)", unit = "minutes"),
  pt2_margin = list(name = "pT2 surgical margin positive", unit = "risk difference"),
  all_margin = list(name = "All-stage surgical margin positive", unit = "risk difference"),
  los = list(name = "Length of stay", unit = "days"),
  readm_30d = list(name = "30-day readmission", unit = "risk difference"),
  psa_pers = list(name = "EAU PSA persistence (window-closed)", unit = "risk difference"),
  pad_free_3m = list(name = "Pad-free continence at 3 months (window-closed)", unit = "risk difference"),
  pad_free_12m = list(name = "Pad-free continence at 12 months (window-closed)", unit = "risk difference")
)

r5_list <- list(
  cohort = "Calendar 2024: Versius (N=99) vs da Vinci restart (N=105)",
  endpoints = list()
)

for (nm in names(endpoints_meta)) {
  pooled <- r5_res$pooled_estimates[[nm]]
  crude <- r5_res$crude_stats[[nm]]
  
  r5_list$endpoints[[nm]] <- list(
    endpoint_name = endpoints_meta[[nm]]$name,
    unit = endpoints_meta[[nm]]$unit,
    crude_versius = crude$versius,
    crude_davinci = crude$davinci,
    adjusted_estimate = round(pooled$estimate, 3),
    se = round(pooled$se, 3),
    ci_95 = round(pooled$ci_95, 3),
    ci_90 = round(pooled$ci_90, 3),
    p_value = signif(pooled$p_value, 4)
  )
}

# --- R6: Bootstrap interval for Versius n* ---
r6_list <- list(
  n_replicates = r6_res$n_replicates,
  method = "Case-resampling bootstrap (B=600), 1 stochastic imputation per replicate, rule applied inside each replicate",
  n_star_15min = list(
    reproduced_point_estimate = r2_res$n_star_time_15,
    bootstrap_median = r6_res$n_star_15$median,
    bootstrap_mean = round(r6_res$n_star_15$mean, 1),
    ci_95_percentile = as.integer(r6_res$n_star_15$ci_95),
    ci_90_percentile = as.integer(r6_res$n_star_15$ci_90),
    distribution_deciles = r6_res$n_star_15$distribution_quantiles
  ),
  n_star_10min = list(
    reproduced_point_estimate = r2_res$n_star_time_10,
    bootstrap_median = r6_res$n_star_10$median,
    bootstrap_mean = round(r6_res$n_star_10$mean, 1),
    ci_95_percentile = as.integer(r6_res$n_star_10$ci_95),
    ci_90_percentile = as.integer(r6_res$n_star_10$ci_90),
    distribution_deciles = r6_res$n_star_10$distribution_quantiles
  )
)

full_repro_results <- list(
  R1 = r1_list,
  R2 = r2_list,
  R3 = r3_list,
  R4 = r4_list,
  R5 = r5_list,
  R6 = r6_list
)

write_json(full_repro_results, output_json, pretty = TRUE, auto_unbox = TRUE)
cat("REPRO_RESULTS.json successfully written to:\n", output_json, "\n")

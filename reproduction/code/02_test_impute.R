# 02_test_impute.R
suppressPackageStartupMessages({
  library(dplyr)
  library(mice)
})

Sys.setenv(TMPDIR = "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/tmp")
restricted_dir <- "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro"
df_pl <- readRDS(file.path(restricted_dir, "derived_data.rds"))

df_rob <- df_pl %>%
  filter(platform %in% c("Versius", "da Vinci")) %>%
  select(
    patient_id = full_polish_series_case_sequence, # internal ID for re-joining
    platform, platform_n, dv_restart_n, year, first_case, centre,
    age, bmi, log_psa, biopsy_isup_cat, log_weight, pt_group,
    nerve_sparing, plnd,
    or_time, psm
  ) %>%
  mutate(
    log_or_time = log(or_time),
    platform = factor(platform),
    centre = factor(centre),
    nerve_sparing = factor(nerve_sparing),
    plnd = factor(plnd),
    first_case = factor(first_case),
    psm_factor = factor(psm)
  )

cat("Robotic rows for imputation:", nrow(df_rob), "\n")
cat("Missingness counts:\n")
sapply(df_rob, function(x) sum(is.na(x)))[sapply(df_rob, function(x) sum(is.na(x))) > 0]

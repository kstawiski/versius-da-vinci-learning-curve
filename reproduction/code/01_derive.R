# 01_derive.R
# Independent cohort derivation and validation of R1 counts
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

Sys.setenv(TMPDIR = "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/tmp")
restricted_dir <- "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro"
dir.create(restricted_dir, recursive = TRUE, showWarnings = FALSE)

parquet_path <- "<PROJECT_ROOT>/data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet/patient_wide.parquet"
df_all <- read_parquet(parquet_path)

# Filter Poland
df_pl <- df_all %>% filter(country == "Poland")
stopifnot(nrow(df_pl) == 838)

# Sorting:
# sort all 838 Polish operations by:
# 1. surgery_date
# 2. surgery_start_minutes_after_midnight (missing last)
# 3. operation_centre_within_day_order (missing last)
# 4. full_polish_series_case_sequence
df_pl <- df_pl %>%
  mutate(
    start_min_sort = ifelse(is.na(surgery_start_minutes_after_midnight), Inf, surgery_start_minutes_after_midnight),
    within_day_sort = ifelse(is.na(operation_centre_within_day_order), Inf, operation_centre_within_day_order)
  ) %>%
  arrange(surgery_date, start_min_sort, within_day_sort, full_polish_series_case_sequence) %>%
  select(-start_min_sort, -within_day_sort)

# platform_n: running count within platform in that order
df_pl <- df_pl %>%
  group_by(platform) %>%
  mutate(platform_n = row_number()) %>%
  ungroup()

# da Vinci restart series: da Vinci cases with platform_n >= 20, dv_restart_n = platform_n - 19
# year = calendar year of surgery_date
df_pl <- df_pl %>%
  mutate(
    dv_restart_n = ifelse(platform == "da Vinci" & platform_n >= 20, platform_n - 19, NA_integer_),
    year = as.integer(format(as.Date(surgery_date), "%Y"))
  )

# day_position = order within the surgeon's operating day (all platforms); first_case = day_position == 1
df_pl <- df_pl %>%
  group_by(surgery_date) %>%
  mutate(
    day_position = row_number(),
    first_case = as.integer(day_position == 1)
  ) %>%
  ungroup()

# months_since_start = days since the first Polish operation / 30.4375
first_surg_date <- min(as.Date(df_pl$surgery_date))
df_pl <- df_pl %>%
  mutate(
    days_since_start = as.numeric(as.Date(surgery_date) - first_surg_date),
    months_since_start = days_since_start / 30.4375
  )

# Data close = the latest (surgery_date + fu_months_last_contact * 30.4375 days) over all 838 Polish rows
contact_dates <- as.Date(df_pl$surgery_date) + (df_pl$fu_months_last_contact * 30.4375)
data_close <- max(contact_dates, na.rm = TRUE)

df_pl <- df_pl %>%
  mutate(
    days_to_close = as.numeric(data_close - as.Date(surgery_date)),
    eligible_continence_3m = days_to_close >= 120,
    eligible_continence_12m = days_to_close >= 425,
    eligible_psa_persistence = days_to_close >= 56
  )

# Outcomes:
# or_time = or_time_min
# psm = margin_status (positive=1, negative=0, blank missing)
# pT group from pT_final: contains T2 -> pT2; T3a -> pT3a; T3b or T4 -> pT3b+
# los = postoperative_length_of_stay_days
# ebl = blood_loss_ml, set to missing where blood_loss_ml_implausible_repeated_value is True
# readmission_30d from readmission_30d (True/False, blank missing)
# psa_persistence_eau from psa_persistence_eau
# pad_free_3m and pad_free_12m from pad_free_3m / pad_free_12m, parsed tolerantly (True/1/1.0 -> 1, False/0/0.0 -> 0)
# Major complication = clavien_dindo_highest_documented grade III or higher (none=0)
df_pl <- df_pl %>%
  mutate(
    or_time = as.numeric(or_time_min),
    psm = case_when(
      margin_status == "positive" ~ 1,
      margin_status == "negative" ~ 0,
      TRUE ~ NA_real_
    ),
    pt_group = case_when(
      grepl("T2", pT_final) ~ "pT2",
      grepl("T3a", pT_final) ~ "pT3a",
      grepl("T3b|T4", pT_final) ~ "pT3b+",
      TRUE ~ NA_character_
    ),
    los = suppressWarnings(as.numeric(postoperative_length_of_stay_days)),
    ebl = ifelse(!is.na(blood_loss_ml_implausible_repeated_value) & blood_loss_ml_implausible_repeated_value == TRUE,
                 NA_real_, as.numeric(blood_loss_ml)),
    readmission_30d_num = case_when(
      readmission_30d == TRUE ~ 1,
      readmission_30d == FALSE ~ 0,
      TRUE ~ NA_real_
    ),
    psa_persistence_eau_num = case_when(
      psa_persistence_eau == TRUE ~ 1,
      psa_persistence_eau == FALSE ~ 0,
      TRUE ~ NA_real_
    ),
    pad_free_3m_num = case_when(
      pad_free_3m %in% c(TRUE, 1, "1", "1.0", "True", "true") ~ 1,
      pad_free_3m %in% c(FALSE, 0, "0", "0.0", "False", "false") ~ 0,
      TRUE ~ NA_real_
    ),
    pad_free_12m_num = case_when(
      pad_free_12m %in% c(TRUE, 1, "1", "1.0", "True", "true") ~ 1,
      pad_free_12m %in% c(FALSE, 0, "0", "0.0", "False", "false") ~ 0,
      TRUE ~ NA_real_
    ),
    major_complication = case_when(
      clavien_dindo_highest_documented %in% c("III", "IIIa", "IIIb", "IVa", "IVb", "V") ~ 1,
      clavien_dindo_highest_documented %in% c("none", "I", "II") ~ 0,
      TRUE ~ NA_real_
    )
  )

# Covariates:
# age; bmi; log(psa_preop); biopsy ISUP category (1, 2, 3, 4 and 5 combined) from biopsy_isup;
# log(prostate_specimen_weight_g); pT group; nerve_sparing_performed;
# pelvic_lymph_node_dissection_state (performed=1, not_performed=0);
# operation_centre; first_case (operative-time models only)
df_pl <- df_pl %>%
  mutate(
    age = as.numeric(age),
    bmi = as.numeric(bmi),
    log_psa = log(as.numeric(psa_preop)),
    biopsy_isup_cat = case_when(
      biopsy_isup %in% c(1) ~ "1",
      biopsy_isup %in% c(2) ~ "2",
      biopsy_isup %in% c(3) ~ "3",
      biopsy_isup %in% c(4, 5) ~ "4+",
      TRUE ~ NA_character_
    ),
    log_weight = log(as.numeric(prostate_specimen_weight_g)),
    pt_group = factor(pt_group, levels = c("pT2", "pT3a", "pT3b+")),
    biopsy_isup_cat = factor(biopsy_isup_cat, levels = c("1", "2", "3", "4+")),
    nerve_sparing = case_when(
      nerve_sparing_performed == TRUE ~ 1,
      nerve_sparing_performed == FALSE ~ 0,
      TRUE ~ NA_real_
    ),
    plnd = case_when(
      pelvic_lymph_node_dissection_state == "performed" ~ 1,
      pelvic_lymph_node_dissection_state == "not_performed" ~ 0,
      TRUE ~ NA_real_
    ),
    centre = factor(operation_centre)
  )

# Save restricted patient-level derived data (mode 0600)
derived_rds <- file.path(restricted_dir, "derived_data.rds")
saveRDS(df_pl, derived_rds)
Sys.chmod(derived_rds, mode = "0600")

cat("Derived data saved successfully to restricted path with mode 0600.\n")
cat("Total Polish rows:", nrow(df_pl), "\n")
cat("Robotic rows (Versius + da Vinci):", sum(df_pl$platform %in% c("Versius", "da Vinci")), "\n")

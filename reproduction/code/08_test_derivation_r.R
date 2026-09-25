# 08_test_derivation_r.R
# Independent derivation of Manuscript 6 cohort in R
suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

Sys.setenv(TMPDIR = "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/tmp")

parquet_path <- "<PROJECT_ROOT>/data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet/patient_wide.parquet"
df <- read_parquet(parquet_path)

# Filter Poland
df_pl <- df %>% filter(country == "Poland")
cat("Total Polish cases:", nrow(df_pl), "\n")

# Sorting rule:
# sort all 838 Polish operations by surgery_date,
# then surgery_start_minutes_after_midnight (missing last),
# then operation_centre_within_day_order,
# then full_polish_series_case_sequence.
# In R: handle NA in start minutes by replacing NA with Inf for sorting
df_pl <- df_pl %>%
  mutate(
    start_min_sort = ifelse(is.na(surgery_start_minutes_after_midnight), Inf, surgery_start_minutes_after_midnight),
    within_day_sort = ifelse(is.na(operation_centre_within_day_order), Inf, operation_centre_within_day_order)
  ) %>%
  arrange(surgery_date, start_min_sort, within_day_sort, full_polish_series_case_sequence) %>%
  select(-start_min_sort, -within_day_sort)

# Running count within platform
df_pl <- df_pl %>%
  group_by(platform) %>%
  mutate(platform_n = row_number()) %>%
  ungroup()

# da Vinci restart series
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

# Data close = latest (surgery_date + fu_months_last_contact * 30.4375 days)
contact_dates <- as.Date(df_pl$surgery_date) + (df_pl$fu_months_last_contact * 30.4375)
data_close <- max(contact_dates, na.rm = TRUE)

df_pl <- df_pl %>%
  mutate(
    days_to_close = as.numeric(data_close - as.Date(surgery_date)),
    eligible_continence_3m = days_to_close >= 120,
    eligible_continence_12m = days_to_close >= 425,
    eligible_psa_persistence = days_to_close >= 56
  )

# Outcomes
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
    readmission_30d = case_when(
      readmission_30d == TRUE ~ 1,
      readmission_30d == FALSE ~ 0,
      TRUE ~ NA_real_
    ),
    psa_persistence_eau = case_when(
      psa_persistence_eau == TRUE ~ 1,
      psa_persistence_eau == FALSE ~ 0,
      TRUE ~ NA_real_
    ),
    pad_free_3m = case_when(
      pad_free_3m %in% c(TRUE, 1, "1", "1.0", "True", "true") ~ 1,
      pad_free_3m %in% c(FALSE, 0, "0", "0.0", "False", "false") ~ 0,
      TRUE ~ NA_real_
    ),
    pad_free_12m = case_when(
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

# Covariates
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

cat("\nCounts validation in R:\n")
cat("Versius count:", sum(df_pl$platform == "Versius"), "\n")
cat("da Vinci count:", sum(df_pl$platform == "da Vinci"), "\n")
cat("da Vinci restart count:", sum(df_pl$platform == "da Vinci" & df_pl$platform_n >= 20), "\n")
cat("2024 Versius count:", sum(df_pl$platform == "Versius" & df_pl$year == 2024), "\n")
cat("2024 restart da Vinci count:", sum(df_pl$platform == "da Vinci" & df_pl$platform_n >= 20 & df_pl$year == 2024), "\n")

cat("\npT2 counts:\n")
cat("Versius pT2:", sum(df_pl$platform == "Versius" & df_pl$pt_group == "pT2", na.rm = TRUE), "\n")
cat("da Vinci pT2:", sum(df_pl$platform == "da Vinci" & df_pl$pt_group == "pT2", na.rm = TRUE), "\n")
cat("Restart da Vinci pT2:", sum(df_pl$platform == "da Vinci" & df_pl$platform_n >= 20 & df_pl$pt_group == "pT2", na.rm = TRUE), "\n")
cat("2024 Versius pT2:", sum(df_pl$platform == "Versius" & df_pl$year == 2024 & df_pl$pt_group == "pT2", na.rm = TRUE), "\n")
cat("2024 restart da Vinci pT2:", sum(df_pl$platform == "da Vinci" & df_pl$platform_n >= 20 & df_pl$year == 2024 & df_pl$pt_group == "pT2", na.rm = TRUE), "\n")

cat("\nor_time counts:\n")
cat("Versius with or_time:", sum(df_pl$platform == "Versius" & !is.na(df_pl$or_time)), "\n")
cat("da Vinci with or_time:", sum(df_pl$platform == "da Vinci" & !is.na(df_pl$or_time)), "\n")
cat("Restart da Vinci with or_time:", sum(df_pl$platform == "da Vinci" & df_pl$platform_n >= 20 & !is.na(df_pl$or_time)), "\n")
cat("2024 Versius with or_time:", sum(df_pl$platform == "Versius" & df_pl$year == 2024 & !is.na(df_pl$or_time)), "\n")
cat("2024 restart da Vinci with or_time:", sum(df_pl$platform == "da Vinci" & df_pl$platform_n >= 20 & df_pl$year == 2024 & !is.na(df_pl$or_time)), "\n")

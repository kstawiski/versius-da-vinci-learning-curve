# 02_impute.R
# Multiple imputation using MICE (20 imputations, PMM)
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
    full_polish_series_case_sequence,
    platform, platform_n, year, first_case, centre,
    age, bmi, log_psa, biopsy_isup_cat, log_weight, pt_group,
    nerve_sparing, plnd,
    or_time, psm
  ) %>%
  mutate(
    log_or_time = log(or_time),
    psm_factor = factor(psm, levels = c(0, 1), labels = c("neg", "pos")),
    platform = factor(platform),
    centre = factor(centre),
    nerve_sparing = factor(nerve_sparing),
    plnd = factor(plnd),
    first_case = factor(first_case),
    biopsy_isup_cat = factor(biopsy_isup_cat, levels = c("1", "2", "3", "4+")),
    pt_group = factor(pt_group, levels = c("pT2", "pT3a", "pT3b+"))
  )

# Prepare mice specs:
# Exclude raw or_time and psm numeric to avoid collinearity with log_or_time and psm_factor
imp_data <- df_rob %>%
  select(
    full_polish_series_case_sequence,
    platform, platform_n, year, first_case, centre,
    age, bmi, log_psa, biopsy_isup_cat, log_weight, pt_group,
    nerve_sparing, plnd,
    log_or_time, psm_factor
  )

ini <- mice(imp_data, maxit = 0)
pred <- ini$predictorMatrix
meth <- ini$method

# Do not predict ID, and do not use ID to predict anything
pred[, "full_polish_series_case_sequence"] <- 0
pred["full_polish_series_case_sequence", ] <- 0

# Set methods: PMM for all incomplete variables
# In mice, pmm works for numeric, binary factor, and polytomous factor
for (v in c("bmi", "log_weight", "biopsy_isup_cat", "pt_group", "log_or_time", "psm_factor")) {
  meth[v] <- "pmm"
}

cat("MICE methods:\n")
print(meth[meth != ""])

cat("Starting MICE (20 imputations, 10 iterations, seed 20260925)...\n")
set.seed(20260925)
imp <- mice(imp_data, m = 20, maxit = 10, method = meth, predictorMatrix = pred, printFlag = FALSE)

cat("Imputation finished.\n")
cat("Logged events in mice:\n")
if (is.null(imp$loggedEvents)) {
  cat("None (clean run).\n")
} else {
  print(imp$loggedEvents)
}

# Create a list of 20 completed full datasets re-joined with the rest of the derived variables
# (such as dv_restart_n, days_to_close, window-closed flags, original outcomes)
# Remember: "Outcomes are never imputed as outcomes."
# So we preserve original or_time and psm as the outcome columns!
imputed_datasets <- vector("list", 20)

for (i in 1:20) {
  comp_i <- complete(imp, i)
  
  # Join imputed covariates back onto df_rob (preserving original un-imputed outcomes)
  joined <- df_rob %>%
    select(full_polish_series_case_sequence, platform, platform_n) %>%
    # join with df_pl to get all original fields
    left_join(
      df_pl %>% select(
        full_polish_series_case_sequence,
        surgery_date, year, day_position, first_case, months_since_start,
        days_to_close, eligible_continence_3m, eligible_continence_12m, eligible_psa_persistence,
        or_time, psm, los, ebl, readmission_30d_num, psa_persistence_eau_num,
        pad_free_3m_num, pad_free_12m_num, major_complication,
        age, log_psa, nerve_sparing, plnd, centre
      ),
      by = "full_polish_series_case_sequence"
    ) %>%
    # replace the missing covariates with the imputed values from comp_i
    mutate(
      bmi = comp_i$bmi,
      log_weight = comp_i$log_weight,
      biopsy_isup_cat = comp_i$biopsy_isup_cat,
      pt_group = comp_i$pt_group,
      dv_restart_n = ifelse(platform == "da Vinci" & platform_n >= 20, platform_n - 19, NA_integer_)
    )
  
  imputed_datasets[[i]] <- joined
}

# Save imputed datasets list to restricted directory with mode 0600
imp_rds <- file.path(restricted_dir, "imputed_datasets.rds")
saveRDS(imputed_datasets, imp_rds)
Sys.chmod(imp_rds, mode = "0600")

cat("Saved 20 completed datasets to restricted path with mode 0600.\n")

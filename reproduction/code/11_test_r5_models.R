# 11_test_r5_models.R
# Test models for R5 (A4 concurrent 2024 platform comparison)
suppressPackageStartupMessages({
  library(dplyr)
  library(sandwich)
})

Sys.setenv(TMPDIR = "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro/tmp")
restricted_dir <- "<PROJECT_ROOT>/manuscript6/analysis/restricted/repro"

imputed_datasets <- readRDS(file.path(restricted_dir, "imputed_datasets.rds"))
df1 <- imputed_datasets[[1]]

c2024 <- df1 %>%
  filter(year == 2024 & (platform == "Versius" | (platform == "da Vinci" & !is.na(dv_restart_n)))) %>%
  mutate(dv_group = ifelse(platform == "da Vinci", 1, 0))

cat("2024 Total cases:", nrow(c2024), "\n")
cat("Versius:", sum(c2024$dv_group == 0), "da Vinci restart:", sum(c2024$dv_group == 1), "\n")

# 1. Operative time
fit_or <- lm(
  or_time ~ dv_group + age + bmi + log_psa + biopsy_isup_cat + log_weight + 
    pt_group + nerve_sparing + plnd + centre + first_case,
  data = c2024 %>% filter(!is.na(or_time))
)
cat("Operative time lm coef on dv_group:", coef(fit_or)["dv_group"], "\n")

# 2. pT2 margin
c2024_pt2 <- c2024 %>% filter(pt_group == "pT2" & !is.na(psm))
fit_pt2 <- glm(
  psm ~ dv_group + age + bmi + log_psa + biopsy_isup_cat + log_weight + centre,
  data = c2024_pt2, family = binomial
)
cat("pT2 psm glm coef on dv_group:", coef(fit_pt2)["dv_group"], "\n")

# 3. All-stage margin
c2024_allm <- c2024 %>% filter(!is.na(psm))
fit_allm <- glm(
  psm ~ dv_group + age + bmi + log_psa + biopsy_isup_cat + log_weight + 
    pt_group + nerve_sparing + plnd + centre,
  data = c2024_allm, family = binomial
)
cat("All-stage psm glm coef on dv_group:", coef(fit_allm)["dv_group"], "\n")

# 4. Length of stay
c2024_los <- c2024 %>% filter(!is.na(los))
fit_los <- lm(
  los ~ dv_group + age + bmi + log_psa + biopsy_isup_cat + log_weight + 
    pt_group + nerve_sparing + plnd + centre,
  data = c2024_los
)
cat("LOS lm coef on dv_group:", coef(fit_los)["dv_group"], "\n")

# 5. Readmission 30d
c2024_readm <- c2024 %>% filter(!is.na(readmission_30d_num))
cat("Readmission cases with data:", nrow(c2024_readm), "events:", sum(c2024_readm$readmission_30d_num), "\n")

# 6. EAU PSA persistence (window-closed)
c2024_psa <- c2024 %>% filter(eligible_psa_persistence & !is.na(psa_persistence_eau_num))
cat("PSA persistence cases:", nrow(c2024_psa), "events:", sum(c2024_psa$psa_persistence_eau_num), "\n")
fit_psa <- glm(
  psa_persistence_eau_num ~ dv_group + age + bmi + log_psa + biopsy_isup_cat + log_weight + 
    pt_group + nerve_sparing + plnd + centre,
  data = c2024_psa, family = binomial
)
cat("PSA persistence glm coef on dv_group:", coef(fit_psa)["dv_group"], "\n")

# 7. Pad-free 3m
c2024_pf3 <- c2024 %>% filter(eligible_continence_3m & !is.na(pad_free_3m_num))
cat("Pad-free 3m cases:", nrow(c2024_pf3), "events:", sum(c2024_pf3$pad_free_3m_num), "\n")
fit_pf3 <- glm(
  pad_free_3m_num ~ dv_group + age + bmi + nerve_sparing + pt_group + centre,
  data = c2024_pf3, family = binomial
)
cat("Pad-free 3m glm coef on dv_group:", coef(fit_pf3)["dv_group"], "\n")

# 8. Pad-free 12m
c2024_pf12 <- c2024 %>% filter(eligible_continence_12m & !is.na(pad_free_12m_num))
cat("Pad-free 12m cases:", nrow(c2024_pf12), "events:", sum(c2024_pf12$pad_free_12m_num), "\n")
fit_pf12 <- glm(
  pad_free_12m_num ~ dv_group + age + bmi + nerve_sparing + pt_group + centre,
  data = c2024_pf12, family = binomial
)
cat("Pad-free 12m glm coef on dv_group:", coef(fit_pf12)["dv_group"], "\n")

#!/usr/bin/env Rscript
# Post hoc sensitivity of the comparisons after the Versius learning phase to the restart operations
# documented as Retzius-sparing, and to all restart operations from October 2025 onwards, the period in which
# all but one of them were done. Uses the same completed covariate datasets as the primary contrasts.
# Writes aggregate results only.
suppressPackageStartupMessages({library(jsonlite); library(arrow)})
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
d <- load_ms6(WS, 20L)$d

# Retzius-sparing surgery is a positive statement in the operative record of the locked release.
# A record that does not mention it is not taken as another approach.
REL <- file.path(dirname(WS), "data/prostatectomy_analysis_database_v1.1_2026-09-04/parquet/patient_wide.parquet")
rel <- as.data.frame(read_parquet(REL, col_select = c("master_patient_id", "retzius_sparing")))
rs <- as.character(d$pid) %in% as.character(rel$master_patient_id[rel$retzius_sparing %in% TRUE])
if (any(rs & is.na(d$dv_restart_n))) stop("A documented Retzius-sparing operation lies outside the restart series.")
pt2 <- function(dd) dd$pT_obs %in% "pT2"
G <- list(
  first100 = function(dd) ifelse(dd$platform == "Versius" & dd$platform_n <= 100, "ref",
                               ifelse(!is.na(dd$dv_restart_n) & dd$dv_restart_n <= 100, "cmp", NA)),
  y2024 = function(dd) ifelse(dd$year == 2024 & dd$platform == "Versius", "ref",
                            ifelse(dd$year == 2024 & !is.na(dd$dv_restart_n), "cmp", NA)),
  after_learning_46 = function(dd) ifelse(dd$platform == "Versius" & dd$platform_n > 46, "ref",
                                         ifelse(!is.na(dd$dv_restart_n), "cmp", NA))
)
restart <- !is.na(d$dv_restart_n)
# The record states Retzius-sparing surgery only where it was done, so the whole period is also left out.
late <- restart & as.Date(d$surgery_date) >= as.Date("2025-10-01")
frac <- function(sel) list(positive = sum(sel & d$psm %in% 1), total = sum(sel & !is.na(d$psm)))
out <- list(label = "Post hoc sensitivity of the comparisons after the Versius learning phase to documented Retzius-sparing surgery",
            contrast = "da Vinci restart minus Versius",
            documented_retzius_sparing = list(
              operations = sum(rs), SalveMedica = sum(rs & d$centre %in% "SM"), Belchatow = sum(rs & d$centre %in% "BE"),
              first_restart_case = min(d$dv_restart_n[rs]), after_restart_case_250 = sum(rs & d$dv_restart_n > 250),
              from_october_2025 = sum(rs & as.Date(d$surgery_date) >= as.Date("2025-10-01")),
              restart_from_october_2025 = sum(restart & as.Date(d$surgery_date) >= as.Date("2025-10-01")),
              in_first100_comparison = sum(rs & G$first100(d) %in% "cmp"), in_2024_comparison = sum(rs & G$y2024(d) %in% "cmp"),
              in_comparison_after_learning = sum(rs & G$after_learning_46(d) %in% "cmp"),
              pT2_margin = frac(rs & pt2(d)), pT2_margin_other_restart = frac(restart & !rs & pt2(d)),
              pT2_margin_by_hospital = lapply(c(SM = "SM", BE = "BE"), function(ct) list(
                documented = frac(rs & pt2(d) & d$centre %in% ct),
                other_restart = frac(restart & !rs & pt2(d) & d$centre %in% ct))),
              from_october_2025_in_first100_comparison = sum(late & !is.na(G$first100(d))),
              from_october_2025_in_2024_comparison = sum(late & !is.na(G$y2024(d))),
              versius_from_october_2025 = sum(d$platform == "Versius" & as.Date(d$surgery_date) >= as.Date("2025-10-01"))))

imps <- impute(d, 20, 11)
excluded <- lapply(imps, function(di) di[!rs, ])
excluded_late <- lapply(imps, function(di) di[!late, ])
margin <- function(completed) {
  x <- pool_contrast(completed, "after_learning_46", G$after_learning_46, "psm", COV_PT2, subset_expr = pt2, eq_margin = 0.10)
  if (isTRUE(x$not_estimable)) stop("pT2 margin contrast not estimable")
  list(estimate_pp = 100 * x$estimate, ci95_pp = 100 * x$ci95,
       crude = list(Versius = list(positive = x$events_ref, total = x$n_ref),
                    da_vinci_restart = list(positive = x$events_cmp, total = x$n_cmp)))
}
as_recorded <- margin(imps)
published <- c(9.8, 1.3, 18.4)
observed <- c(as_recorded$estimate_pp, as_recorded$ci95_pp)
if (any(abs(observed - published) > 0.05))
  stop(sprintf("as-recorded pT2 estimate and 95%% CI (pp) were %s, expected %s within 0.05 pp",
               paste(sprintf("%.4f", observed), collapse = ", "), paste(sprintf("%.1f", published), collapse = ", ")))
out$pT2_margin_after_learning <- list(as_recorded = as_recorded, retzius_sparing_excluded = margin(excluded),
                                      restart_from_october_2025_excluded = margin(excluded_late))

COV_OR_NOCENTRE <- sub(" \\+ centre", "", COV_OR)
published_or <- list(SM = c(11.7, 4.0, 19.3), BE = c(-11.6, -20.0, -3.2))
out$operative_time_after_learning <- list()
for (ct in c("SM", "BE")) {
  gh <- function(dd) ifelse(dd$centre == ct, G$after_learning_46(dd), NA)
  time <- function(completed) {
    x <- pool_contrast(completed, paste("within", ct), gh, "or_time", COV_OR_NOCENTRE)
    if (isTRUE(x$not_estimable)) stop(sprintf("operative-time contrast at %s not estimable", ct))
    list(estimate_min = x$estimate, ci95_min = x$ci95, n_Versius = x$n_ref, n_da_vinci_restart = x$n_cmp)
  }
  a <- time(imps)
  if (any(abs(c(a$estimate_min, a$ci95_min) - published_or[[ct]]) > 0.05))
    stop(sprintf("as-recorded operative-time difference and 95%% CI at %s were %s, expected %s within 0.05 min", ct,
                 paste(sprintf("%.3f", c(a$estimate_min, a$ci95_min)), collapse = ", "), paste(published_or[[ct]], collapse = ", ")))
  out$operative_time_after_learning[[ct]] <- list(as_recorded = a, retzius_sparing_excluded = time(excluded),
                                                  restart_from_october_2025_excluded = time(excluded_late))
}
write_json(out, file.path(WS, "analysis/results/S_retzius_sparing.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("Retzius-sparing sensitivity done\n")

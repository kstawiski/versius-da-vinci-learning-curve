#!/usr/bin/env Rscript
# Post hoc sensitivity of the pT2 margin contrasts to positive margins discordant with pathological R status.
# Uses the same completed covariate datasets for all variants. Writes aggregate results only.
suppressPackageStartupMessages(library(jsonlite))
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
d <- load_ms6(WS, 20L)$d

# The derived flag uses the release concordance field audited in 24_endpoint_verification.py.
if (!"margin_r_discordant" %in% names(d))
  stop("The discordance flag is required. Run scripts/01_derive.py first.")
discordant <- d$psm %in% 1 & d$margin_r_discordant %in% 1
pt2 <- function(dd) dd$pT_obs %in% "pT2"
G <- list(
  first100 = function(dd) ifelse(dd$platform == "Versius" & dd$platform_n <= 100, "ref",
                               ifelse(!is.na(dd$dv_restart_n) & dd$dv_restart_n <= 100, "cmp", NA)),
  y2024 = function(dd) ifelse(dd$year == 2024 & dd$platform == "Versius", "ref",
                            ifelse(dd$year == 2024 & !is.na(dd$dv_restart_n), "cmp", NA)),
  after_learning_46 = function(dd) ifelse(dd$platform == "Versius" & dd$platform_n > 46, "ref",
                                         ifelse(!is.na(dd$dv_restart_n), "cmp", NA))
)
published <- list(first100 = c(20.8, 2.1, 39.6), y2024 = c(8.4, -11.4, 28.2),
                  after_learning_46 = c(9.8, 1.3, 18.4))
imps <- impute(d, 20, 11)

run_variant <- function(comparison, variant, completed, discordant_counts) {
  x <- pool_contrast(completed, comparison, G[[comparison]], "psm", COV_PT2, subset_expr = pt2, eq_margin = 0.10)
  if (isTRUE(x$not_estimable) || !all(is.finite(c(x$estimate, x$ci95))))
    stop(sprintf("pT2 margin contrast not estimable for %s (%s)", comparison, variant))
  list(estimate_pp = 100 * x$estimate, ci95_pp = 100 * x$ci95,
       crude = list(Versius = list(positive = x$events_ref, total = x$n_ref),
                    da_vinci_restart = list(positive = x$events_cmp, total = x$n_cmp)),
       discordant_in_comparison_pT2 = discordant_counts)
}

out <- list(label = "Post hoc sensitivity of pT2 positive margin differences to discordance with pathological R status",
            contrast = "da Vinci restart minus Versius",
            discordant_counts_basis = "Operations with a discordant positive margin before exclusion or recoding")
series <- ifelse(d$platform == "Versius", "Versius", ifelse(!is.na(d$dv_restart_n), "da Vinci restart", NA))
by_hospital_series <- expand.grid(hospital = c("SM", "BE"), series = c("Versius", "da Vinci restart"),
                                  stringsAsFactors = FALSE)
by_hospital_series$all_stages <- mapply(function(h, s) sum(discordant & d$centre %in% h & series %in% s),
                                       by_hospital_series$hospital, by_hospital_series$series)
by_hospital_series$pT2 <- mapply(function(h, s) sum(discordant & pt2(d) & d$centre %in% h & series %in% s),
                                by_hospital_series$hospital, by_hospital_series$series)
out$discordant_by_hospital_and_series <- by_hospital_series

for (comparison in names(G)) {
  grp <- G[[comparison]](d)
  counts <- list(Versius = sum(discordant & pt2(d) & grp %in% "ref"),
                 da_vinci_restart = sum(discordant & pt2(d) & grp %in% "cmp"))
  x <- run_variant(comparison, "as_recorded", imps, counts)
  observed <- c(x$estimate_pp, x$ci95_pp)
  if (any(abs(observed - published[[comparison]]) > 0.05))
    stop(sprintf("%s as_recorded estimate and 95%% CI (pp) were %s, expected %s within 0.05 pp",
                 comparison, paste(sprintf("%.4f", observed), collapse = ", "),
                 paste(sprintf("%.1f", published[[comparison]]), collapse = ", ")))
  out$comparisons[[comparison]] <- list(as_recorded = x)
}

excluded <- lapply(imps, function(di) di[!discordant, ])
recoded <- lapply(imps, function(di) { di$psm[discordant] <- 0; di })
for (comparison in names(G)) {
  counts <- out$comparisons[[comparison]]$as_recorded$discordant_in_comparison_pT2
  out$comparisons[[comparison]]$discordant_excluded <- run_variant(comparison, "discordant_excluded", excluded, counts)
  out$comparisons[[comparison]]$discordant_recoded_negative <- run_variant(comparison, "discordant_recoded_negative", recoded, counts)
}
write_json(out, file.path(WS, "analysis/results/S_margin_r_status.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("margin R-status sensitivity done\n")

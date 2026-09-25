#!/usr/bin/env Rscript
# POST HOC sensitivity of the fixed-cutoff comparison after the Versius learning phase.
# The fixed cutoff (Versius after case 46) was taken from the bootstrap limit of the first segmented fit
# (95% CI 19-46). The corrected fit pooled over 20 imputations gives 19-30, and case 46 equals the unadjusted
# breakpoint. The cutoff is kept, and this script repeats the comparison with the corrected upper limit
# (case 30) and with the case by which 90% of the initial improvement had occurred (case 43). Aggregate output.
suppressPackageStartupMessages({ library(jsonlite) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
dat <- load_ms6(WS, 20L); d <- dat$d
imps <- impute(d, 20, 11)
run <- function(k) {
  g <- function(dd) ifelse(dd$platform == "Versius" & dd$platform_n > k, "ref", ifelse(!is.na(dd$dv_restart_n), "cmp", NA))
  list(label = sprintf("POST HOC: Versius cases %d-337 against da Vinci restart cases 1-382", k + 1),
    or_time = pool_contrast(imps, "fixed or", g, "or_time", COV_OR, eq_margin = 15),
    psm_pt2 = pool_contrast(imps, "fixed pT2", g, "psm", COV_PT2, subset_expr = function(dd) dd$pT_obs %in% "pT2", eq_margin = 0.10),
    psm_all = pool_contrast(imps, "fixed all", g, "psm", COV))
}
out <- list(label = "POST HOC cutoff sensitivity of the comparison after the Versius learning phase",
  after_30 = run(30), after_43 = run(43), after_46 = run(46))
write_json(out, file.path(WS, "analysis/results/B1f_cutoff_sensitivity.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("cutoff sensitivity done\n")

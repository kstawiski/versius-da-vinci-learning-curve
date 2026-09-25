#!/usr/bin/env Rscript
# Pooled tests of the smooth case-number term across the 20 imputations (D2 rule of Li, Meng, Raghunathan and
# Rubin, 1991), replacing the median of per-imputation P values. Each imputation contributes the smooth-term
# test statistic on the chi-square scale (Chi.sq for the binomial margin models, F multiplied by the reference
# degrees of freedom for the operative-time models) with the reference degrees of freedom. Aggregate output.
suppressPackageStartupMessages({ library(jsonlite); library(mgcv) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
dat <- load_ms6(WS, 20L); d <- dat$d
imps <- impute(d, 20, 11)
SERIES <- list(
  Versius   = function(dd) { dd <- dd[dd$platform == "Versius", ]; dd$x <- dd$platform_n; dd },
  dVrestart = function(dd) { dd <- dd[!is.na(dd$dv_restart_n), ]; dd$x <- dd$dv_restart_n; dd })
d2 <- function(stat, k) {
  m <- length(stat); r2 <- (1 + 1 / m) * var(sqrt(stat))
  D <- (mean(stat) / k - (m + 1) / (m - 1) * r2) / (1 + r2)
  df2 <- k^(-3 / m) * (m - 1) * (1 + 1 / r2)^2
  c(D2 = max(D, 0), df1 = k, df2 = df2, p = pf(max(D, 0), k, df2, lower.tail = FALSE))
}
out <- list(label = "Smooth case-number term, pooled across imputations by the D2 rule")
for (sr in names(SERIES)) for (oc in c("or", "psm_pt2")) {
  k <- if (oc == "or") 10 else 5
  st <- t(sapply(imps, function(di) {
    s <- summary(fit_curve(SERIES[[sr]](di), oc, k))$s.table
    stat <- if ("Chi.sq" %in% colnames(s)) s[1, "Chi.sq"] else s[1, "F"] * s[1, "Ref.df"]
    c(stat = stat, ref_df = s[1, "Ref.df"], p = s[1, "p-value"])
  }))
  pooled <- d2(st[, "stat"], mean(st[, "ref_df"]))
  out[[paste(sr, oc)]] <- list(pooled = as.list(pooled), per_imputation_p_median = median(st[, "p"]),
                               per_imputation_p_range = range(st[, "p"]))
}
write_json(out, file.path(WS, "analysis/results/S_pooled_smooth_tests.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("pooled smooth tests done\n")

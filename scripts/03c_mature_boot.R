#!/usr/bin/env Rscript
# B1 mature-phase comparison with cutoff uncertainty (methods-review amendment).
# Each bootstrap replicate resamples the Versius and da Vinci restart series, imputes once, refits both
# operative-time curves, re-selects both plateau cutoffs with the prespecified 15-minute rule, and recomputes
# the operative-time and pT2-margin contrasts. Replicates in which either plateau is not reached are counted,
# not replaced. A post hoc fixed-cutoff comparison (Versius after case 46, the upper bootstrap limit of the
# prespecified segmented breakpoint, against all restart cases) is added with 20 imputations. Aggregate output.
suppressPackageStartupMessages({ library(parallel); library(jsonlite) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
NB <- as.integer(Sys.getenv("MS6_NBOOT_B1", "1000")); NCORE <- as.integer(Sys.getenv("MS6_NCORE", "16"))
dat <- load_ms6(WS, 20L); d <- dat$d
V <- d[d$platform == "Versius", ]; R <- d[!is.na(d$dv_restart_n), ]
one <- function(b) {
  set.seed(30000 + b)
  bs <- rbind(V[sample(nrow(V), replace = TRUE), ], R[sample(nrow(R), replace = TRUE), ])
  di <- tryCatch(impute(bs, 1, 30000 + b)[[1]], error = function(e) NULL); if (is.null(di)) return(c(fail = 1, rep(NA, 6)))
  nstar <- function(sel) { dd <- sel(di); fit <- tryCatch(fit_curve(dd, "or", 10), error = function(e) NULL)
    if (is.null(fit)) return(-1); plateau(std_curve(fit, dd, "or", seq_len(max(dd$x))), 15)$nstar }
  nv <- nstar(function(z) { z <- z[z$platform == "Versius", ]; z$x <- z$platform_n; z })
  nr <- nstar(function(z) { z <- z[!is.na(z$dv_restart_n), ]; z$x <- z$dv_restart_n; z })
  if (identical(nv, -1) || identical(nr, -1)) return(c(fail = 1, rep(NA, 6)))
  if (is.na(nv) || is.na(nr)) return(c(fail = 0, reached = 0, nv = ifelse(is.na(nv), NA, nv), nr = ifelse(is.na(nr), NA, nr), or = NA, psm = NA, n = NA))
  cr <- if (nr <= 1) 0 else nr
  g <- function(dd) ifelse(dd$platform == "Versius" & dd$platform_n > nv, "ref", ifelse(!is.na(dd$dv_restart_n) & dd$dv_restart_n > cr, "cmp", NA))
  o <- tryCatch(fit_contrast(di, g, "or_time", COV_OR)$est, error = function(e) NA)
  p <- tryCatch(fit_contrast(di, g, "psm", COV_PT2, function(dd) dd$pT_obs %in% "pT2")$est, error = function(e) NA)
  c(fail = 0, reached = 1, nv = nv, nr = nr, or = o, psm = p, n = sum(!is.na(g(di))))
}
B <- do.call(rbind, mclapply(seq_len(NB), one, mc.cores = NCORE))
ok <- B[, "fail"] == 0; reached <- ok & B[, "reached"] == 1
pc <- function(v, pr) unname(quantile(v, pr, na.rm = TRUE))
out <- list(label = "B1 mature-phase contrasts with plateau cutoffs re-selected in every bootstrap replicate",
  n_replicates = NB, n_fit_failed = sum(!ok), n_plateau_not_reached = sum(ok & !reached), share_reached = mean(reached[ok]),
  or_time = list(p5 = pc(B[reached, "or"], .05), p95 = pc(B[reached, "or"], .95), p2.5 = pc(B[reached, "or"], .025), p97.5 = pc(B[reached, "or"], .975),
                 share_within_15 = mean(abs(B[reached, "or"]) < 15, na.rm = TRUE)),
  psm_pt2 = list(p5 = pc(B[reached, "psm"], .05), p95 = pc(B[reached, "psm"], .95), p2.5 = pc(B[reached, "psm"], .025), p97.5 = pc(B[reached, "psm"], .975),
                 share_within_10pts = mean(abs(B[reached, "psm"]) < 0.10, na.rm = TRUE)),
  cutoffs_versius = pc(B[reached, "nv"], c(.025, .5, .975)), cutoffs_restart = pc(B[reached, "nr"], c(.025, .5, .975)))
# post hoc fixed cutoff: Versius after case 46 against all restart cases, 20 imputations
imps <- impute(d, 20, 11)
gfix <- function(dd) ifelse(dd$platform == "Versius" & dd$platform_n > 46, "ref", ifelse(!is.na(dd$dv_restart_n), "cmp", NA))
out$posthoc_fixed_cutoff <- list(
  label = "POST HOC: Versius cases 47-337 against da Vinci restart cases 1-382",
  or_time = pool_contrast(imps, "fixed or", gfix, "or_time", COV_OR, eq_margin = 15),
  psm_pt2 = pool_contrast(imps, "fixed pT2", gfix, "psm", COV_PT2, subset_expr = function(dd) dd$pT_obs %in% "pT2", eq_margin = 0.10),
  psm_all = pool_contrast(imps, "fixed all", gfix, "psm", COV))
write_json(out, file.path(WS, "analysis/results/B1_mature_bootstrap.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("mature bootstrap done\n")

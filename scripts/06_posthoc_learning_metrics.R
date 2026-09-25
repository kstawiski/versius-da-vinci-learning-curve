#!/usr/bin/env Rscript
# POST HOC (added after the prespecified plateau rule proved unstable; labelled post hoc in every output).
# Robust descriptions of the steep learning phase from the same standardized GAM curves:
#   (1) first case at which the curve enters the 15-minute band around its asymptote ("first entry");
#   (2) case by which 80% and 90% of the initial improvement (curve at case 1 minus asymptote) is achieved;
#   (3) fitted decline from case 1 to case 50.
# Bootstrap: case resampling within series, one stochastic imputation per replicate.
suppressPackageStartupMessages({ library(parallel); library(jsonlite) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
NB <- as.integer(Sys.getenv("MS6_NBOOT_PH", "1000")); NCORE <- as.integer(Sys.getenv("MS6_NCORE", "16"))
dat <- load_ms6(WS, 20L); d <- dat$d
SERIES <- list(
  Versius   = function(dd) { dd <- dd[dd$platform == "Versius", ]; dd$x <- dd$platform_n; dd },
  dVrestart = function(dd) { dd <- dd[!is.na(dd$dv_restart_n), ]; dd$x <- dd$dv_restart_n; dd })
metrics <- function(cv) {
  asym <- mean(tail(cv, 50)); c1 <- cv[1]; imp <- c1 - asym
  first_entry <- unname(which(abs(cv - asym) <= 15)[1])
  # crossing of a fraction of the initial improvement; undefined (NA) when there is no initial improvement
  f <- function(frac) if (imp <= 0) NA_real_ else unname(which(cv <= c1 - frac * imp)[1])
  c(first_entry = first_entry, frac80 = f(0.8), frac90 = f(0.9), decline_1_50 = unname(cv[1] - cv[50]),
    case1 = unname(cv[1]), case50 = unname(cv[50]), asym = asym, improvement = unname(imp))
}
out <- list(label = "POST HOC learning-phase metrics (not prespecified)")
imps <- impute(d, 20, 11)
for (sr in names(SERIES)) {
  ns <- seq_len(max(SERIES[[sr]](d)$x))
  cm <- rowMeans(sapply(imps, function(di) { dd <- SERIES[[sr]](di); std_curve(fit_curve(dd, "or", 10), dd, "or", ns) }))
  point <- metrics(cm)
  bs <- do.call(rbind, mclapply(seq_len(NB), function(b) {
    set.seed(7000 + b); base <- SERIES[[sr]](d); bb <- base[sample(nrow(base), replace = TRUE), ]
    di <- tryCatch(impute(bb, 1, 7000 + b)[[1]], error = function(e) NULL); if (is.null(di)) return(rep(NA, 8))
    fit <- tryCatch(fit_curve(di, "or", 10), error = function(e) NULL); if (is.null(fit)) return(rep(NA, 8))
    metrics(std_curve(fit, di, "or", ns))
  }, mc.cores = NCORE))
  colnames(bs) <- names(point)
  out[[sr]] <- list(point = as.list(point), n_boot_ok = sum(complete.cases(bs)),
                    ci95 = lapply(as.data.frame(bs), function(v) unname(quantile(v, c(.025, .975), na.rm = TRUE))),
                    median = lapply(as.data.frame(bs), function(v) median(v, na.rm = TRUE)))
}
write_json(out, file.path(WS, "analysis/results/posthoc_learning_metrics.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("posthoc done\n")

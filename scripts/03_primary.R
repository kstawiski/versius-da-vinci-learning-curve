#!/usr/bin/env Rscript
# Manuscript6 primary analyses A1-A4 and B1 (ANALYSIS_PLAN.md v2 and plan/SAP_CHANGELOG.md).
# Reads the restricted analysis dataset; writes aggregate results only.
# Never prints identifiers, dates or row-level values.
suppressPackageStartupMessages({ library(parallel); library(jsonlite) })
set.seed(20260925)
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
RES <- Sys.getenv("MS6_RES", file.path(WS, "analysis/results"))
dir.create(RES, recursive = TRUE, showWarnings = FALSE)
NBOOT <- as.integer(Sys.getenv("MS6_NBOOT", "2000"))
NCORE <- as.integer(Sys.getenv("MS6_NCORE", "12"))
M_IMP <- as.integer(Sys.getenv("MS6_MIMP", "20"))
RESTART_AT <- 20L  # da Vinci restart at platform case 20 (SAP v2)

dat <- load_ms6(WS, RESTART_AT); d <- dat$d
imps <- impute(d, M_IMP, 11)

SERIES <- list(
  Versius   = function(dd) { dd <- dd[dd$platform == "Versius", ]; dd$x <- dd$platform_n; dd },
  dVrestart = function(dd) { dd <- dd[!is.na(dd$dv_restart_n), ]; dd$x <- dd$dv_restart_n; dd })

# ---------- A1/A2 and A3b: standardized curves and plateaus ----------
curves <- list(); plateaus <- list()
for (sr in names(SERIES)) for (oc in c("or", "psm_pt2")) {
  k <- if (oc == "or") 10 else 5
  sel <- SERIES[[sr]]; ns <- seq_len(max(sel(d)$x))
  fits <- lapply(imps, function(di) { dd <- sel(di); list(fit = fit_curve(dd, oc, k), dd = dd) })
  cm <- rowMeans(sapply(fits, function(f) std_curve(f$fit, f$dd, oc, ns)))
  # posterior draws from every imputation (mixture simulation over imputations)
  sims <- do.call(cbind, lapply(fits, function(f) sim_curve(f$fit, f$dd, oc, ns, 100)))
  st <- t(sapply(fits, function(f) { s <- summary(f$fit)$s.table; c(s[1, "edf"], s[1, "p-value"]) }))
  pl <- plateau(cm, if (oc == "or") 15 else 0.05)
  pl10 <- if (oc == "or") plateau(cm, 10)$nstar else NA
  L <- length(cm)
  curves[[paste(sr, oc)]] <- data.frame(series = sr, outcome = oc, n = ns, est = cm,
                                        lo = apply(sims, 1, quantile, 0.025), hi = apply(sims, 1, quantile, 0.975))
  plateaus[[paste(sr, oc)]] <- list(series = sr, outcome = oc, n_model = sum(outcome_rows(fits[[1]]$dd, oc)),
    nstar = pl$nstar, nstar_sens10 = pl10, asymptote = pl$asym,
    value_case1 = unname(cm[1]), value_case25 = unname(cm[min(25, L)]), value_case50 = unname(cm[min(50, L)]),
    value_case100 = unname(cm[min(100, L)]), value_last = unname(cm[L]),
    slope_last100_per10 = if (L >= 100) unname((cm[L] - cm[L - 99]) / 99 * 10) else NA,
    edf_median = median(st[, 1]), smooth_p_median = median(st[, 2]), smooth_p_max = max(st[, 2]))
}
write.csv(do.call(rbind, curves), file.path(RES, "A1_standardized_curves.csv"), row.names = FALSE)

boot_one <- function(b, sr, oc) {
  set.seed(1000 + b)
  base <- SERIES[[sr]](d); bs <- base[sample(nrow(base), replace = TRUE), ]
  di <- tryCatch(impute(bs, 1, 1000 + b)[[1]], error = function(e) NULL)
  if (is.null(di)) return(c(NA, NA, 0))
  fit <- tryCatch(fit_curve(di, oc, if (oc == "or") 10 else 5), error = function(e) NULL)
  if (is.null(fit)) return(c(NA, NA, 0))
  cv <- std_curve(fit, di, oc, seq_len(max(base$x)))
  if (oc == "or") c(plateau(cv, 15)$nstar, plateau(cv, 10)$nstar, 1) else c(plateau(cv, 0.05)$nstar, NA, 1)
}
# Non-attainment (rule never met within the series) is kept as probability mass beyond the last case.
# Unconditional percentiles treat it as +Inf; an upper limit of Inf is reported as "not reached".
bsum <- function(v, ok) {
  vv <- v[ok]; vinf <- ifelse(is.na(vv), Inf, vv)
  q <- unname(quantile(vinf, c(.025, .5, .975), type = 1))
  list(n_ok = sum(ok), n_fit_failed = sum(!ok), n_not_reached = sum(is.na(vv)), share_not_reached = mean(is.na(vv)),
       unconditional = list(p2.5 = q[1], median = q[2], p97.5 = q[3]),
       conditional_on_reached = if (sum(!is.na(vv)) > 10) unname(quantile(vv, c(.025, .5, .975), na.rm = TRUE)) else c(NA, NA, NA),
       share_case1 = mean(vv[!is.na(vv)] == 1))
}
if (NBOOT > 0) for (sr in names(SERIES)) for (oc in c("or", "psm_pt2")) {
  bn <- do.call(rbind, mclapply(seq_len(NBOOT), boot_one, sr = sr, oc = oc, mc.cores = NCORE))
  ok <- bn[, 3] == 1; key <- paste(sr, oc)
  plateaus[[key]]$boot_primary <- bsum(bn[, 1], ok)
  if (oc == "or") plateaus[[key]]$boot_sens10 <- bsum(bn[, 2], ok)
}
write_json(plateaus, file.path(RES, "A2_plateaus.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)

# ---------- A3a initial phase: dV restart 1-100 vs Versius 1-100 ----------
g_init <- function(dd) ifelse(dd$platform == "Versius" & dd$platform_n <= 100, "ref",
                       ifelse(!is.na(dd$dv_restart_n) & dd$dv_restart_n <= 100, "cmp", NA))
A3 <- list(
  a_or = pool_contrast(imps, "A3a dV restart 1-100 vs Versius 1-100, operative time (min)", g_init, "or_time", COV_OR, eq_margin = 15),
  a_psm_pt2 = pool_contrast(imps, "A3a dV restart 1-100 vs Versius 1-100, pT2 PSM (RD)", g_init, "psm", COV_PT2,
                            subset_expr = function(dd) dd$pT_obs %in% "pT2", eq_margin = 0.10),
  a_psm_all = pool_contrast(imps, "A3a dV restart 1-100 vs Versius 1-100, PSM all stages (RD)", g_init, "psm", COV))
a3s <- lapply(imps, function(di) {
  dd <- di[(di$platform == "Versius" & di$platform_n <= 50) | (!is.na(di$dv_restart_n) & di$dv_restart_n <= 50), ]
  dd$n50 <- ifelse(dd$platform == "Versius", dd$platform_n, dd$dv_restart_n); dd$dv <- as.integer(dd$platform == "da Vinci")
  m <- lm(as.formula(paste("or_time ~ dv * n50 +", COV_OR)), data = dd[!is.na(dd$or_time), ])
  V <- robust_vcov(m); b <- coef(m); z <- setNames(rep(0, length(b)), names(b))
  L <- function(w) c(sum(w * b), as.numeric(t(w) %*% V %*% w))
  w1 <- z; w1[c("dv", "dv:n50")] <- 1; w2 <- z; w2["n50"] <- 10; w3 <- z; w3[c("n50", "dv:n50")] <- 10; w4 <- z; w4["dv:n50"] <- 10
  rbind(case1_diff_min = L(w1), slope_versius_per10 = L(w2), slope_dvrestart_per10 = L(w3), slope_diff_per10 = L(w4))
})
A3$slopes_first50 <- lapply(rownames(a3s[[1]]), function(r) {
  pr <- rubin(sapply(a3s, function(x) x[r, 1]), sapply(a3s, function(x) x[r, 2])); q <- qt(.975, pr["df"])
  list(term = r, est = unname(pr["est"]), ci95 = unname(pr["est"] + c(-1, 1) * q * pr["se"]), p = unname(2 * pt(-abs(pr["est"] / pr["se"]), pr["df"])))
})

# ---------- A3c cumulative robotic vs platform-specific experience ----------
a3c <- lapply(imps, function(di) {
  dd <- di[!is.na(di$log_or), ]
  m1 <- gam(as.formula(paste("log_or ~ s(robotic_n, k = 10) + platform +", COV_OR)), data = dd, method = "ML")
  m2 <- gam(as.formula(paste("log_or ~ s(platform_n, by = platform, k = 10) + platform +", COV_OR)), data = dd, method = "ML")
  m3 <- gam(as.formula(paste("log_or ~ s(robotic_n, k = 10) + s(platform_n, by = platform, k = 10) + platform +", COV_OR)), data = dd, method = "ML")
  c(aic_cumulative = AIC(m1), aic_platform = AIC(m2), aic_both = AIC(m3),
    dev_cumulative = summary(m1)$dev.expl, dev_platform = summary(m2)$dev.expl, dev_both = summary(m3)$dev.expl,
    p_platform_given_cumulative = anova(m1, m3, test = "F")$`Pr(>F)`[2],
    p_cumulative_given_platform = anova(m2, m3, test = "F")$`Pr(>F)`[2])
})
A3$c <- as.list(apply(do.call(rbind, a3c), 2, median))

# ---------- A4 concurrent 2024 comparison; B1 all-era mature sensitivity ----------
run_cmp <- function(gdef, tag) {
  pc <- function(lbl, oc, cov, sub = NULL, eq = NA) pool_contrast(imps, paste(tag, lbl), gdef, oc, cov, subset_expr = sub, eq_margin = eq)
  list(or_time = pc("operative time (min)", "or_time", COV_OR, eq = 15),
       psm_pt2 = pc("pT2 PSM (RD)", "psm", COV_PT2, function(dd) dd$pT_obs %in% "pT2", 0.10),
       psm_all = pc("PSM all stages (RD)", "psm", COV),
       los = pc("length of stay (days)", "los", COV),
       ebl = pc("estimated blood loss (mL)", "ebl", COV),
       readmission_30d = pc("30-day readmission (RD)", "readmission_30d", "age + centre"),
       psa_persistence_eau = pc("EAU PSA persistence (RD)", "psa_persistence_eau", COV, function(dd) dd$elig_psa56 == 1),
       pad_free_3m = pc("pad-free 3 months (RD)", "pad_free_3m", "age + bmi + nerve_sparing + pT_group + centre", function(dd) dd$elig_3m == 1),
       pad_free_12m = pc("pad-free 12 months (RD)", "pad_free_12m", "age + bmi + nerve_sparing + pT_group + centre", function(dd) dd$elig_12m == 1),
       rare = rare_summary(d, gdef, c("major_cd_any", "reoperation", "urine_leak", "transfusion_30d")),
       observation = observation_summary(d, gdef))
}
g_2024 <- function(dd) ifelse(dd$year == 2024 & dd$platform == "Versius", "ref", ifelse(dd$year == 2024 & !is.na(dd$dv_restart_n), "cmp", NA))
A4 <- run_cmp(g_2024, "A4 2024 dV restart vs Versius,")
cal <- function(cov) paste(cov, "+ ns(months_since_start, df = 3)")
A4$calendar_adjusted <- list(
  or_time = pool_contrast(imps, "A4cal operative time (min)", g_2024, "or_time", cal(COV_OR), eq_margin = 15),
  psm_pt2 = pool_contrast(imps, "A4cal pT2 PSM (RD)", g_2024, "psm", cal(COV_PT2), subset_expr = function(dd) dd$pT_obs %in% "pT2", eq_margin = 0.10),
  psm_all = pool_contrast(imps, "A4cal PSM all (RD)", g_2024, "psm", cal(COV)))

nv <- plateaus[["Versius or"]]$nstar; ndr <- plateaus[["dVrestart or"]]$nstar
# No fallback: if a plateau was not reached, the mature-phase comparison is not estimable.
cut_v <- if (is.na(nv)) Inf else nv
cut_dr <- if (is.na(ndr)) Inf else if (ndr <= 1) 0 else ndr
g_mature <- function(dd) ifelse(dd$platform == "Versius" & dd$platform_n > cut_v, "ref",
                         ifelse(!is.na(dd$dv_restart_n) & dd$dv_restart_n > cut_dr, "cmp", NA))
B1 <- run_cmp(g_mature, "B1 mature dV restart vs mature Versius,")
B1$definitions <- list(versius_after = cut_v, dv_restart_after = cut_dr)

write_json(list(A3 = A3, A4 = A4, B1 = B1, data_close_days_after_last_surgery = dat$close_gap_days,
                n_imputations = M_IMP, n_bootstrap = NBOOT),
           file.path(RES, "A3_A4_B1_contrasts.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("A1-A4 and B1 complete\n")

#!/usr/bin/env Rscript
# Manuscript6 secondary (B2-B6) and sensitivity (S-*) analyses. Aggregate outputs only,
# except per-case CUSUM series, which are written to analysis/restricted/ for figure rendering.
suppressPackageStartupMessages({ library(parallel); library(jsonlite); library(segmented) })
set.seed(20260926)
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
RES <- Sys.getenv("MS6_RES", file.path(WS, "analysis/results")); dir.create(RES, showWarnings = FALSE, recursive = TRUE)
RESTR <- file.path(WS, "analysis/restricted")
NCORE <- as.integer(Sys.getenv("MS6_NCORE", "12")); M_IMP <- as.integer(Sys.getenv("MS6_MIMP", "20"))
NSEGBOOT <- as.integer(Sys.getenv("MS6_NSEGBOOT", "500")); NPERM <- as.integer(Sys.getenv("MS6_NPERM", "200"))
dat <- load_ms6(WS, 20L); d <- dat$d
imps <- impute(d, M_IMP, 11)
SERIES <- list(
  Versius   = function(dd) { dd <- dd[dd$platform == "Versius", ]; dd$x <- dd$platform_n; dd },
  dVrestart = function(dd) { dd <- dd[!is.na(dd$dv_restart_n), ]; dd$x <- dd$dv_restart_n; dd },
  dVall     = function(dd) { dd <- dd[dd$platform == "da Vinci", ]; dd$x <- dd$platform_n; dd })
out <- list()

curve_nstar <- function(sr, oc, spline = "gam", cov_or = COV_OR, cov_pt2 = COV_PT2, idata = imps) {
  sel <- SERIES[[sr]]; ns <- seq_len(max(sel(d)$x)); k <- if (oc == "or") 10 else 5
  cs <- sapply(idata, function(di) { dd <- sel(di); std_curve(fit_curve(dd, oc, k, cov_or, cov_pt2, spline), dd, oc, ns) })
  cm <- if (is.matrix(cs)) rowMeans(cs) else cs
  m <- if (oc == "or") 15 else 0.05
  list(nstar = plateau(cm, m)$nstar, nstar10 = if (oc == "or") plateau(cm, 10)$nstar else NA, asym = plateau(cm, m)$asym,
       case1 = unname(cm[1]), case50 = unname(cm[min(50, length(cm))]), case100 = unname(cm[min(100, length(cm))]), last = unname(cm[length(cm)]))
}

# ---------- Sensitivity: spline form, adjustment sets, series definitions ----------
out$S_rcs <- list(Versius_or = curve_nstar("Versius", "or", "rcs"), Versius_psm = curve_nstar("Versius", "psm_pt2", "rcs"),
                  dVrestart_or = curve_nstar("dVrestart", "or", "rcs"), dVrestart_psm = curve_nstar("dVrestart", "psm_pt2", "rcs"))
COV_OR_NOLNS <- "age + bmi + log_psa + isup_cat + log_weight + pT_group + centre + first_case"
out$S_noNSLND <- list(Versius_or = curve_nstar("Versius", "or", cov_or = COV_OR_NOLNS), dVrestart_or = curve_nstar("dVrestart", "or", cov_or = COV_OR_NOLNS))
out$S_NS <- list(Versius_psm = curve_nstar("Versius", "psm_pt2", cov_pt2 = paste(COV_PT2, "+ nerve_sparing")),
                 dVrestart_psm = curve_nstar("dVrestart", "psm_pt2", cov_pt2 = paste(COV_PT2, "+ nerve_sparing")))
out$S_dvall <- list(or = curve_nstar("dVall", "or"), psm = curve_nstar("dVall", "psm_pt2"))
# complete case: drop rows with any missing covariate
cc <- d[complete.cases(d[, c("bmi", "log_weight", "isup_cat", "pT_group")]), ]
out$S_cc <- list(n_complete = nrow(cc), Versius_or = curve_nstar("Versius", "or", idata = list(cc)),
                 Versius_psm = curve_nstar("Versius", "psm_pt2", idata = list(cc)),
                 dVrestart_or = curve_nstar("dVrestart", "or", idata = list(cc)), dVrestart_psm = curve_nstar("dVrestart", "psm_pt2", idata = list(cc)))
g_init <- function(dd) ifelse(dd$platform == "Versius" & dd$platform_n <= 100, "ref", ifelse(!is.na(dd$dv_restart_n) & dd$dv_restart_n <= 100, "cmp", NA))
g_2024 <- function(dd) ifelse(dd$year == 2024 & dd$platform == "Versius", "ref", ifelse(dd$year == 2024 & !is.na(dd$dv_restart_n), "cmp", NA))
out$S_cc$A3a_or <- pool_contrast(list(cc), "cc A3a or", g_init, "or_time", COV_OR)
out$S_cc$A3a_psm_pt2 <- pool_contrast(list(cc), "cc A3a pT2", g_init, "psm", COV_PT2, subset_expr = function(dd) dd$pT_obs %in% "pT2")
out$S_cc$A4_or <- pool_contrast(list(cc), "cc A4 or", g_2024, "or_time", COV_OR)
out$S_cc$A4_psm_pt2 <- pool_contrast(list(cc), "cc A4 pT2", g_2024, "psm", COV_PT2, subset_expr = function(dd) dd$pT_obs %in% "pT2")
out$S_cc$A4_psm_all <- pool_contrast(list(cc), "cc A4 all", g_2024, "psm", COV)
# unadjusted contrasts for transparency
out$unadjusted <- list(A3a_or = pool_contrast(imps[1], "unadj A3a or", g_init, "or_time", "centre"),
                       A4_or = pool_contrast(imps[1], "unadj A4 or", g_2024, "or_time", "centre"),
                       A4_psm_pt2 = pool_contrast(imps[1], "unadj A4 pT2", g_2024, "psm", "centre", subset_expr = function(dd) dd$pT_obs %in% "pT2"))

# ---------- S-seg: segmented regression, pooled over imputations, multi-start ----------
seg_fit <- function(dd) {
  dd <- dd[!is.na(dd$log_or), ]
  m0 <- lm(as.formula(paste("log_or ~ x +", COV_OR)), data = dd)
  starts <- unique(round(quantile(dd$x, c(.05, .1, .2, .35, .5, .65))))
  best <- NULL
  for (st in starts) {
    sg <- tryCatch(segmented(m0, seg.Z = ~x, psi = st), error = function(e) NULL)
    if (!is.null(sg) && !is.null(sg$psi) && (is.null(best) || deviance(sg) < deviance(best))) best <- sg
  }
  list(m0 = m0, sg = best)
}
out$S_seg <- list()
for (sr in c("Versius", "dVrestart")) {
  du <- SERIES[[sr]](d); du <- du[!is.na(du$log_or), ]
  mu <- lm(log_or ~ x, data = du); best <- NULL
  for (st in unique(round(quantile(du$x, c(.05, .1, .2, .35, .5, .65))))) {
    sg <- tryCatch(segmented(mu, seg.Z = ~x, psi = st), error = function(e) NULL)
    if (!is.null(sg) && !is.null(sg$psi) && (is.null(best) || deviance(sg) < deviance(best))) best <- sg }
  out$S_seg[[paste0(sr, "_unadjusted")]] <- if (is.null(best)) list(breakpoint = NA) else
    list(breakpoint = unname(best$psi[1, "Est."]), breakpoint_se = unname(best$psi[1, "St.Err"]),
         davies_p = tryCatch(davies.test(mu, seg.Z = ~x)$p.value, error = function(e) NA))
  per <- lapply(imps, function(di) {
    f <- seg_fit(SERIES[[sr]](di))
    dv <- tryCatch(davies.test(f$m0, seg.Z = ~x)$p.value, error = function(e) NA)
    if (is.null(f$sg)) return(c(psi = NA, psi_se = NA, s1 = NA, s1_se = NA, s2 = NA, s2_se = NA, davies = dv))
    sl <- slope(f$sg)$x
    c(psi = f$sg$psi[1, "Est."], psi_se = f$sg$psi[1, "St.Err"], s1 = sl[1, "Est."], s1_se = sl[1, "St.Err."],
      s2 = sl[2, "Est."], s2_se = sl[2, "St.Err."], davies = dv)
  })
  P <- do.call(rbind, per); okm <- !is.na(P[, "psi"])
  pool <- function(e, se) { r <- rubin(P[okm, e], P[okm, se]^2); q <- qt(.975, r["df"]); c(est = unname(r["est"]), lo = unname(r["est"] - q * r["se"]), hi = unname(r["est"] + q * r["se"])) }
  psi <- pool("psi", "psi_se"); s1 <- pool("s1", "s1_se"); s2 <- pool("s2", "s2_se")
  bs <- unlist(mclapply(seq_len(NSEGBOOT), function(b) {
    set.seed(5000 + b); base <- SERIES[[sr]](d); bb <- base[sample(nrow(base), replace = TRUE), ]
    di <- tryCatch(impute(bb, 1, 5000 + b)[[1]], error = function(e) NULL); if (is.null(di)) return(NA)
    g <- seg_fit(di); if (is.null(g$sg)) NA else g$sg$psi[1, "Est."]
  }, mc.cores = NCORE))
  pct <- function(b) (exp(b * 10) - 1) * 100
  out$S_seg[[sr]] <- list(n_imputations_fitted = sum(okm), breakpoint = psi[["est"]], breakpoint_rubin_ci95 = unname(psi[c("lo", "hi")]),
                          slope_before_pct_per10 = pct(s1[["est"]]), slope_before_ci95 = pct(unname(s1[c("lo", "hi")])),
                          slope_after_pct_per10 = pct(s2[["est"]]), slope_after_ci95 = pct(unname(s2[c("lo", "hi")])),
                          davies_p_median = median(P[, "davies"], na.rm = TRUE), davies_p_max = max(P[, "davies"], na.rm = TRUE),
                          boot_ci95 = unname(quantile(bs, c(.025, .975), na.rm = TRUE)), boot_n_ok = sum(!is.na(bs)), boot_n_failed = sum(is.na(bs)))
}

# ---------- S-cusum, S-racusum, S-lccusum ----------
dd_all <- imps[[1]]
risk <- glm(psm ~ pT_group + I(isup_final >= 3) + log_psa + log_weight + centre, data = dd_all[!is.na(dd_all$psm) & !is.na(dd_all$isup_final), ], family = binomial)
rd <- dd_all[!is.na(dd_all$psm) & !is.na(dd_all$isup_final), ]; rd$p <- predict(risk, rd, type = "response")
auc <- function(y, p) { r <- rank(p); n1 <- sum(y == 1); n0 <- sum(y == 0); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
rf <- psm ~ pT_group + I(isup_final >= 3) + log_psa + log_weight + centre
opt <- replicate(200, { bi <- rd[sample(nrow(rd), replace = TRUE), ]; mb <- glm(rf, data = bi, family = binomial)
  auc(bi$psm, predict(mb, bi, type = "response")) - auc(rd$psm, predict(mb, rd, type = "response")) })
rd$q5 <- cut(rd$p, quantile(rd$p, 0:5 / 5), include.lowest = TRUE, labels = FALSE)
rd$one <- 1
oe <- aggregate(cbind(observed = psm, expected = p, n = one) ~ q5, data = rd, FUN = sum)
out$S_racusum_model <- list(n = nrow(rd), events = sum(rd$psm), auc_apparent = auc(rd$psm, rd$p),
                            auc_optimism_corrected = auc(rd$psm, rd$p) - mean(opt),
                            observed_expected_by_risk_quintile = oe,
                            note = "Calibration slope and intercept on the development sample are identities and are not reported.")
steiner_h <- function(p, RA = 2, alpha = 0.05, nsim = 2000) {
  # h such that P(signal within the series | in control) = alpha
  mx <- replicate(nsim, { y <- rbinom(length(p), 1, p); w <- y * log(RA) - log(1 - p + RA * p); s <- 0; m <- 0
    for (wi in w) { s <- max(0, s + wi); m <- max(m, s) }; m })
  unname(quantile(mx, 1 - alpha))
}
lc_h <- function(N, p0 = 0.10, p1 = 0.20, alpha = 0.05, nsim = 5000) {
  s_ok <- log((1 - p0) / (1 - p1)); s_fail <- log(p0 / p1)
  mx <- replicate(nsim, { y <- rbinom(N, 1, p1); s <- 0; m <- 0
    for (yi in y) { s <- max(0, s + ifelse(yi == 1, s_fail, s_ok)); m <- max(m, s) }; m })
  unname(quantile(mx, 1 - alpha))
}
cusum_rows <- list(); out$S_cusum <- list()
for (sr in c("Versius", "dVrestart")) {
  s <- SERIES[[sr]](d); s <- s[order(s$x), ]
  # operative-time CUSUM (deviation from series mean)
  ot <- s[!is.na(s$or_time), ]; cs_ot <- cumsum(ot$or_time - mean(ot$or_time))
  out$S_cusum[[sr]] <- list(peak_case = ot$x[which.max(cs_ot)], n = nrow(ot))
  # RA-CUSUM for PSM (all stages)
  r <- rd[rd$pid %in% s$pid, ]; r <- r[order(match(r$pid, s$pid)), ]
  r$x <- s$x[match(r$pid, s$pid)]
  w <- r$psm * log(2) - log(1 - r$p + 2 * r$p)
  S <- Reduce(function(a, b) max(0, a + b), w, accumulate = TRUE, 0)[-1]
  h <- steiner_h(r$p)
  oe <- cumsum(r$p - r$psm)  # VLAD: expected minus observed (positive = fewer PSM than expected)
  out$S_cusum[[paste0(sr, "_racusum")]] <- list(n = nrow(r), h = h, signals = sum(S >= h), first_signal_case = if (any(S >= h)) r$x[which(S >= h)[1]] else NA,
                                                  observed = sum(r$psm), expected = sum(r$p), vlad_final = tail(oe, 1))
  # LC-CUSUM for pT2 PSM
  p2 <- s[which(s$pT_group == "pT2" & !is.na(s$psm)), ]; p2 <- p2[order(p2$x), ]
  sc <- ifelse(p2$psm == 1, log(0.10 / 0.20), log(0.90 / 0.80))
  L <- Reduce(function(a, b) max(0, a + b), sc, accumulate = TRUE, 0)[-1]; hl <- lc_h(nrow(p2))
  sig_prob <- function(pp, nsim = 20000) mean(replicate(nsim, { y <- rbinom(nrow(p2), 1, pp); s0 <- 0; hit <- FALSE
    for (yi in y) { s0 <- max(0, s0 + ifelse(yi == 1, log(0.10 / 0.20), log(0.90 / 0.80))); if (s0 >= hl) { hit <- TRUE; break } }; hit }))
  achieved_alpha <- sig_prob(0.20); achieved_power <- sig_prob(0.10)
  out$S_cusum[[paste0(sr, "_lccusum_pT2")]] <- list(n_pT2 = nrow(p2), events = sum(p2$psm), h = hl,
                                                     achieved_alpha = achieved_alpha, achieved_beta = 1 - achieved_power,
                                                     competence_signal = any(L >= hl), first_signal_case = if (any(L >= hl)) p2$x[which(L >= hl)[1]] else NA,
                                                     pT2_rank_at_signal = if (any(L >= hl)) which(L >= hl)[1] else NA)
  cusum_rows[[sr]] <- list(ot = data.frame(series = sr, x = ot$x, cusum = cs_ot),
                           ra = data.frame(series = sr, x = r$x, racusum = S, vlad = oe, h = h),
                           lc = data.frame(series = sr, x = p2$x, lccusum = L, h = hl))
}
saveRDS(cusum_rows, file.path(RESTR, "cusum_series.rds")); Sys.chmod(file.path(RESTR, "cusum_series.rds"), "600")

# ---------- B2 surgeon versus team ----------
out$B2 <- list()
for (pf in c("Versius", "da Vinci")) {
  res <- lapply(imps, function(di) {
    dd <- di[di$platform == pf & !is.na(di$log_or), ]
    m1 <- gam(as.formula(paste("log_or ~ s(platform_n, k = 10) +", COV_OR)), data = dd, method = "ML")
    m2 <- gam(as.formula(paste("log_or ~ s(platform_n, k = 10) + s(team_platform_n, k = 10) +", COV_OR)), data = dd, method = "ML")
    m3 <- gam(as.formula(paste("log_or ~ s(team_platform_n, k = 10) +", COV_OR)), data = dd, method = "ML")
    c(cor = cor(dd$platform_n, dd$team_platform_n), aic_surgeon = AIC(m1), aic_both = AIC(m2), aic_team = AIC(m3),
      p_team_given_surgeon = anova(m1, m2, test = "F")$`Pr(>F)`[2], p_surgeon_given_team = anova(m3, m2, test = "F")$`Pr(>F)`[2],
      dev_surgeon = summary(m1)$dev.expl, dev_both = summary(m2)$dev.expl)
  })
  out$B2[[pf]] <- as.list(apply(do.call(rbind, res), 2, median))
}
rs_ <- d[!is.na(d$dv_restart_n), ]
out$B2$note <- "Surgeon and team case numbers correlate above the prespecified 0.8 limit, so joint-model tests are not reported; hospital-stratified team curves (script 07) are the prespecified display."
out$B2$cor_dv_restart_only <- cor(rs_$platform_n, rs_$team_platform_n)
# hospital-stratified Versius curves (descriptive)
out$B2_hospital <- list()
for (ct in c("SM", "BE")) {
  sel <- function(dd) { dd <- dd[dd$platform == "Versius" & dd$centre == ct, ]; dd$x <- dd$platform_n; dd }
  cs <- sapply(imps, function(di) { dd <- sel(di); ns <- sort(unique(dd$x))
    fit <- gam(as.formula(paste("log_or ~ s(x, k = 8) +", sub(" \\+ centre", "", COV_OR))), data = dd[!is.na(dd$log_or), ], method = "REML")
    std_curve(fit, dd, "or", ns) })
  cm <- rowMeans(cs); out$B2_hospital[[ct]] <- list(n = length(cm), first_surgeon_case = min(as.numeric(names(cm))), fitted_first = unname(cm[1]), fitted_last = unname(cm[length(cm)]))
}

# ---------- B3 platform alternation within 2024 ----------
rob <- d[order(d$robotic_n), ]
prev_pf <- c(NA, as.character(rob$platform[-nrow(rob)])); prev_gap <- c(NA, diff(as.numeric(as.Date(rob$surgery_date))))
rob$switch <- as.integer(!is.na(prev_pf) & prev_pf != as.character(rob$platform) & prev_gap <= 7)
sw_map <- setNames(rob$switch, rob$pid)
b3 <- lapply(imps, function(di) {
  dd <- di[di$year == 2024 & !is.na(di$or_time), ]; dd$switch <- sw_map[dd$pid]
  m <- lm(as.formula(paste("or_time ~ switch + platform + platform:platform_n +", COV_OR)), data = dd)
  V <- robust_vcov(m); c(coef(m)["switch"], V["switch", "switch"], sum(dd$switch), nrow(dd))
})
pr <- rubin(sapply(b3, `[`, 1), sapply(b3, `[`, 2)); q <- qt(.975, pr["df"])
out$B3 <- list(n_switch = b3[[1]][3], n_total = b3[[1]][4], est_min = unname(pr["est"]), ci95 = unname(pr["est"] + c(-1, 1) * q * pr["se"]),
               p = unname(2 * pt(-abs(pr["est"] / pr["se"]), pr["df"])))

# ---------- S-ties: permute unresolved order within tie blocks (all 838 operations), imputation held fixed ----------
allops <- as.data.frame(arrow::read_parquet(file.path(WS, "analysis/restricted/ms6_analysis_dataset.parquet")))
allops <- allops[order(allops$series_n), ]
blocks <- split(seq_len(nrow(allops)), allops$tie_block)
blocks <- blocks[names(blocks) != "-1" & lengths(blocks) > 1]
fixed <- imps[[1]]  # covariates fixed so that only the ordering varies
covcols <- c("bmi", "log_weight", "isup_cat", "pT_group")
perm_one <- function(b) {
  set.seed(9000 + b); a <- allops
  for (ix in blocks) {
    slots <- sort(a$series_n[ix]); st <- a$start_min[ix]
    known <- ix[!is.na(st)]; unknown <- ix[is.na(st)]
    # known start times keep their order; identical start times are shuffled among themselves
    known <- known[order(a$start_min[known], runif(length(known)))]
    seqn <- known
    for (u in unknown[sample.int(length(unknown))]) { pos <- sample.int(length(seqn) + 1, 1); seqn <- append(seqn, u, after = pos - 1) }
    a$series_n[seqn] <- slots
  }
  a <- a[order(a$series_n), ]
  a$platform_n <- ave(seq_len(nrow(a)), a$platform, FUN = seq_along)
  a$team_platform_n <- ave(seq_len(nrow(a)), paste(a$centre, a$platform), FUN = seq_along)
  a$day_position <- ave(seq_len(nrow(a)), as.Date(a$surgery_date), FUN = seq_along)
  r <- a[a$robotic == 1, ]; r$robotic_n <- seq_len(nrow(r))
  di <- fixed[match(r$pid, fixed$pid), ]
  for (v in c("platform_n", "team_platform_n", "robotic_n")) di[[v]] <- r[[v]]
  di$first_case <- as.integer(r$day_position == 1)
  di$dv_restart_n <- ifelse(di$platform == "da Vinci" & di$platform_n >= 20, di$platform_n - 19, NA)
  res <- sapply(list(c("Versius", "or"), c("Versius", "psm_pt2"), c("dVrestart", "or")), function(z) {
    sdat <- SERIES[[z[1]]](di); ns <- seq_len(max(sdat$x))
    fit <- tryCatch(fit_curve(sdat, z[2], if (z[2] == "or") 10 else 5), error = function(e) NULL); if (is.null(fit)) return(NA)
    plateau(std_curve(fit, sdat, z[2], ns), if (z[2] == "or") 15 else 0.05)$nstar })
  a3 <- tryCatch(fit_contrast(di, g_init, "or_time", COV_OR)$est, error = function(e) NA)
  c(res, a3)
}
pm <- do.call(rbind, mclapply(seq_len(NPERM), perm_one, mc.cores = NCORE))
rng <- function(v) if (all(is.na(v))) c(NA, NA) else range(v, na.rm = TRUE)
out$S_ties <- list(n_blocks = length(blocks), n_operations_in_blocks = sum(lengths(blocks)),
                   n_robotic_in_blocks = sum(allops$robotic[unlist(blocks)] == 1), n_perm = NPERM,
                   Versius_or_nstar_range = rng(pm[, 1]), Versius_or_nstar_not_reached = sum(is.na(pm[, 1])),
                   Versius_psm_nstar_range = rng(pm[, 2]), dVrestart_or_nstar_range = rng(pm[, 3]),
                   A3a_or_range = rng(pm[, 4]))

# ---------- B4/B5/B6 descriptive tables ----------
d$phase <- with(d, ifelse(platform == "Versius", paste0("V", cut(platform_n, c(0, 50, 100, 200, 337), labels = c("001-050", "051-100", "101-200", "201-337"))),
                  ifelse(platform_n < 20, "dV001-019 early", paste0("dVr", cut(dv_restart_n, c(0, 50, 100, 200, 382), labels = c("001-050", "051-100", "101-200", "201-382"))))))
nm <- function(x) { x <- x[!is.na(x)]; c(n = length(x), mean = mean(x), median = median(x), q1 = unname(quantile(x, .25)), q3 = unname(quantile(x, .75))) }
bn <- function(x, el = NULL) { if (!is.null(el)) x <- x[el == 1]; x <- x[!is.na(x)]; c(events = sum(x), n = length(x), pct = if (length(x)) 100 * mean(x) else NA) }
ph <- split(d, d$phase)
tab <- do.call(rbind, lapply(names(ph), function(p) { g <- ph[[p]]
  data.frame(phase = p, n = nrow(g),
    or_median = nm(g$or_time)["median"], or_q1 = nm(g$or_time)["q1"], or_q3 = nm(g$or_time)["q3"], or_mean = nm(g$or_time)["mean"],
    ebl_median = nm(g$ebl)["median"], los_median = nm(g$los)["median"], los_mean = nm(g$los)["mean"],
    psm_events = bn(g$psm)["events"], psm_n = bn(g$psm)["n"], psm_pct = bn(g$psm)["pct"],
    psm_pt2_events = bn(g$psm_pt2)["events"], psm_pt2_n = bn(g$psm_pt2)["n"], psm_pt2_pct = bn(g$psm_pt2)["pct"],
    psm_pt3_events = bn(g$psm_pt3)["events"], psm_pt3_n = bn(g$psm_pt3)["n"], psm_pt3_pct = bn(g$psm_pt3)["pct"],
    major_cd_events = bn(g$major_cd_any)["events"], major_cd_n = bn(g$major_cd_any)["n"],
    reop_events = bn(g$reoperation)["events"], reop_n = bn(g$reoperation)["n"],
    readm30_events = bn(g$readmission_30d)["events"], readm30_n = bn(g$readmission_30d)["n"],
    transf30_events = bn(g$transfusion_30d)["events"], transf30_n = bn(g$transfusion_30d)["n"],
    leak_events = bn(g$urine_leak)["events"], leak_n = bn(g$urine_leak)["n"],
    ln_yield_median = nm(g$ln_yield)["median"], plnd_pct = 100 * mean(g$plnd, na.rm = TRUE), ns_pct = 100 * mean(g$nerve_sparing, na.rm = TRUE),
    psa_pers_events = bn(g$psa_persistence_eau, g$elig_psa56)["events"], psa_pers_n = bn(g$psa_persistence_eau, g$elig_psa56)["n"], psa_pers_elig = sum(g$elig_psa56),
    pf3_events = bn(g$pad_free_3m, g$elig_3m)["events"], pf3_n = bn(g$pad_free_3m, g$elig_3m)["n"], pf3_elig = sum(g$elig_3m),
    pf12_events = bn(g$pad_free_12m, g$elig_12m)["events"], pf12_n = bn(g$pad_free_12m, g$elig_12m)["n"], pf12_elig = sum(g$elig_12m),
    fu_contact_median = nm(g$fu_months_contact)["median"], row.names = NULL) }))
write.csv(tab, file.path(RES, "B4_outcomes_by_phase.csv"), row.names = FALSE)
# series totals
d$series <- with(d, ifelse(platform == "Versius", "Versius", ifelse(platform_n < 20, "dV early 1-19", "dV restart")))
tot <- do.call(rbind, lapply(split(d, d$series), function(g) data.frame(series = g$series[1], n = nrow(g),
  psm = paste(bn(g$psm)[1:2], collapse = "/"), psm_pt2 = paste(bn(g$psm_pt2)[1:2], collapse = "/"), psm_pt3 = paste(bn(g$psm_pt3)[1:2], collapse = "/"),
  or_median = median(g$or_time, na.rm = TRUE), major_cd = paste(bn(g$major_cd_any)[1:2], collapse = "/"), reop = paste(bn(g$reoperation)[1:2], collapse = "/"),
  readm30 = paste(bn(g$readmission_30d)[1:2], collapse = "/"), transf30 = paste(bn(g$transfusion_30d)[1:2], collapse = "/"),
  pf3 = paste(bn(g$pad_free_3m, g$elig_3m)[1:2], collapse = "/"), pf12 = paste(bn(g$pad_free_12m, g$elig_12m)[1:2], collapse = "/"),
  psa_pers = paste(bn(g$psa_persistence_eau, g$elig_psa56)[1:2], collapse = "/"))))
write.csv(tot, file.path(RES, "B4_outcomes_by_series.csv"), row.names = FALSE)
# margin sites by series (positive cases)
site <- do.call(rbind, lapply(split(d[d$psm %in% 1, ], d$series[d$psm %in% 1]), function(g) data.frame(series = g$series[1], n_psm = nrow(g),
  apex = sum(g$psm_site_apex == 1, na.rm = TRUE), base = sum(g$psm_site_base == 1, na.rm = TRUE), posterior = sum(g$psm_site_posterior == 1, na.rm = TRUE),
  bladder_neck = sum(g$psm_site_bladder_neck == 1, na.rm = TRUE), nvb = sum(g$psm_site_nvb == 1, na.rm = TRUE),
  length_gt3mm = sum(g$psm_length_mm > 3, na.rm = TRUE), length_known = sum(!is.na(g$psm_length_mm)))))
write.csv(site, file.path(RES, "B4_margin_sites.csv"), row.names = FALSE)

# ---------- B5 IPW sensitivity for continence in the 2024 comparison ----------
ipw_rd <- function(oc, el) {
  dd <- imps[[1]]; dd$grp <- g_2024(dd); dd <- dd[!is.na(dd$grp) & dd[[el]] == 1, ]
  dd$obs <- as.integer(!is.na(dd[[oc]]))
  pm <- glm(obs ~ grp + age + centre + pT_group + nerve_sparing + months_since_start, data = dd, family = binomial)
  ps <- predict(pm, type = "response"); w <- 1 / ps; w <- pmin(pmax(w, quantile(w, .01)), quantile(w, .99))
  o <- dd[dd$obs == 1, ]; wo <- w[dd$obs == 1]
  r1 <- weighted.mean(o[[oc]][o$grp == "cmp"], wo[o$grp == "cmp"]); r0 <- weighted.mean(o[[oc]][o$grp == "ref"], wo[o$grp == "ref"])
  list(eligible = nrow(dd), observed = sum(dd$obs), weighted_ref = r0, weighted_cmp = r1, weighted_rd = r1 - r0)
}
out$B5_ipw <- list(pad_free_3m = ipw_rd("pad_free_3m", "elig_3m"), pad_free_12m = ipw_rd("pad_free_12m", "elig_12m"),
                  psa_persistence_eau = ipw_rd("psa_persistence_eau", "elig_psa56"))

# ---------- Table 1 case mix by series ----------
t1 <- function(g) { f <- function(x) sprintf("%.1f (%.1f-%.1f)", median(x, na.rm = TRUE), quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))
  pc <- function(x) sprintf("%d/%d (%.1f%%)", sum(x, na.rm = TRUE), sum(!is.na(x)), 100 * mean(x, na.rm = TRUE))
  data.frame(n = nrow(g), age = f(g$age), bmi = f(g$bmi), psa = f(g$psa), weight = f(g$specimen_weight_g),
    biopsy_isup_ge3 = pc(g$biopsy_isup >= 3), isup_final_ge3 = pc(g$isup_final >= 3), pT3plus = pc(g$pT3plus), pN1 = pc(g$node_positive),
    eau_high_or_la = pc(ifelse(is.na(g$eau_risk), NA, g$eau_risk %in% c("high", "locally_advanced"))), eau_low = pc(ifelse(is.na(g$eau_risk), NA, g$eau_risk == "low")),
    nerve_sparing = pc(g$nerve_sparing), plnd = pc(g$plnd), centre_BE = pc(g$centre == "BE"), row.names = NULL) }
tab1 <- do.call(rbind, lapply(c("Versius", "dV early 1-19", "dV restart"), function(s) cbind(series = s, t1(d[d$series == s, ]))))
tab1 <- rbind(tab1, cbind(series = "2024 Versius", t1(d[d$year == 2024 & d$platform == "Versius", ])), cbind(series = "2024 dV restart", t1(d[d$year == 2024 & !is.na(d$dv_restart_n), ])))
write.csv(tab1, file.path(RES, "table1_by_series.csv"), row.names = FALSE)

write_json(out, file.path(RES, "secondary_sensitivity.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("secondary complete\n")

#!/usr/bin/env Rscript
# POST HOC robustness checks requested in pre-submission review (SAP_CHANGELOG, 27 September 2026).
# 1. Firth penalized logistic regression for the headline pT2 margin contrasts, with the full prespecified
#    pT2 adjustment set, marginally standardized and pooled across the 20 imputations by Rubin's rules.
# 2. Fixed-basis tests of the margin and operative-time curves: natural cubic splines with fixed degrees of
#    freedom in ordinary regression models, pooled across imputations by the D1 rule (multivariate Wald) and
#    by the D2 rule on likelihood-ratio statistics with fixed degrees of freedom. These check the approximate
#    pooled GAM smooth-term tests of 03e.
# 3. The 2024 comparison restricted to hospital-months in which both robots were used, reported descriptively.
# Aggregate output only.
suppressPackageStartupMessages({ library(jsonlite); library(logistf); library(splines) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
dat <- load_ms6(WS, 20L); d <- dat$d
imps <- impute(d, 20, 11)
out <- list(label = "POST HOC robustness checks for the revision (Firth contrasts, fixed-basis curve tests, 2024 hospital-months)")

# ---------- 1. Firth contrasts with the full pT2 adjustment set ----------
G <- list(
  first100 = function(dd) ifelse(dd$platform == "Versius" & dd$platform_n <= 100, "ref", ifelse(!is.na(dd$dv_restart_n) & dd$dv_restart_n <= 100, "cmp", NA)),
  after_learning_46 = function(dd) ifelse(dd$platform == "Versius" & dd$platform_n > 46, "ref", ifelse(!is.na(dd$dv_restart_n), "cmp", NA)),
  y2024 = function(dd) ifelse(dd$year == 2024 & dd$platform == "Versius", "ref", ifelse(dd$year == 2024 & !is.na(dd$dv_restart_n), "cmp", NA)))
firth_rd_full <- function(dd, rhs) {
  dd$grp <- factor(dd$grp, levels = c("ref", "cmp")); dd$centre <- droplevels(dd$centre); dd$isup_cat <- droplevels(dd$isup_cat)
  f <- as.formula(paste("psm ~ grp +", rhs))
  m <- logistf(f, data = dd, plcontrol = logistf.control(maxit = 200), control = logistf.control(maxit = 200))
  tt <- as.formula(paste("~ grp +", rhs))
  X1 <- model.matrix(tt, transform(dd, grp = factor("cmp", levels = c("ref", "cmp"))))
  X0 <- model.matrix(tt, transform(dd, grp = factor("ref", levels = c("ref", "cmp"))))
  b <- coef(m); p1 <- plogis(X1 %*% b); p0 <- plogis(X0 %*% b)
  g <- colMeans(as.vector(p1 * (1 - p1)) * X1) - colMeans(as.vector(p0 * (1 - p0)) * X0)
  list(est = mean(p1) - mean(p0), var = as.numeric(t(g) %*% vcov(m) %*% g), or = unname(exp(b["grpcmp"])))
}
out$firth <- list()
for (cmp in names(G)) {
  base <- d; base$grp <- G[[cmp]](base); base <- base[!is.na(base$grp) & base$pT_obs %in% "pT2" & !is.na(base$psm), ]
  rhs <- if (length(unique(base$centre)) > 1) COV_PT2 else sub(" \\+ centre", "", COV_PT2)
  rs <- lapply(imps, function(di) { dd <- di; dd$grp <- G[[cmp]](dd); dd <- dd[!is.na(dd$grp) & dd$pT_obs %in% "pT2" & !is.na(dd$psm), ]; firth_rd_full(dd, rhs) })
  pr <- rubin(sapply(rs, `[[`, "est"), sapply(rs, `[[`, "var")); q <- qt(.975, pr["df"]); q9 <- qt(.95, pr["df"])
  out$firth[[cmp]] <- list(adjustment = rhs, n_ref = sum(base$grp == "ref"), n_cmp = sum(base$grp == "cmp"),
                           events_ref = sum(base$psm[base$grp == "ref"]), events_cmp = sum(base$psm[base$grp == "cmp"]),
                           rd = unname(pr["est"]), ci95 = unname(pr["est"] + c(-1, 1) * q * pr["se"]), ci90 = unname(pr["est"] + c(-1, 1) * q9 * pr["se"]),
                           p = unname(2 * pt(-abs(pr["est"] / pr["se"]), pr["df"])), or_median = median(sapply(rs, `[[`, "or")))
}

# ---------- 2. Fixed-basis curve tests pooled by D1 and D2 ----------
SERIES <- list(
  Versius   = function(dd) { dd <- dd[dd$platform == "Versius", ]; dd$x <- dd$platform_n; dd },
  dVrestart = function(dd) { dd <- dd[!is.na(dd$dv_restart_n), ]; dd$x <- dd$dv_restart_n; dd })
d1_pool <- function(Q, U) {  # Q: m x k estimates, U: list of k x k covariance matrices
  m <- nrow(Q); k <- ncol(Q); qbar <- colMeans(Q); ubar <- Reduce(`+`, U) / m; B <- cov(Q)
  ui <- solve(ubar); r1 <- (1 + 1 / m) * sum(diag(B %*% ui)) / k
  D1 <- as.numeric(t(qbar) %*% ui %*% qbar) / (k * (1 + r1)); tt <- k * (m - 1)
  v1 <- if (tt > 4) 4 + (tt - 4) * (1 + (1 - 2 / tt) / r1)^2 else tt * (1 + 1 / k) * (1 + 1 / r1)^2 / 2
  c(D1 = D1, df1 = k, df2 = v1, p = pf(D1, k, v1, lower.tail = FALSE))
}
d2_pool <- function(stat, k) {
  m <- length(stat); r2 <- (1 + 1 / m) * var(sqrt(stat))
  D <- (mean(stat) / k - (m + 1) / (m - 1) * r2) / (1 + r2)
  df2 <- k^(-3 / m) * (m - 1) * (1 + 1 / r2)^2
  c(D2 = max(D, 0), df1 = k, df2 = df2, p = pf(max(D, 0), k, df2, lower.tail = FALSE))
}
out$fixed_basis <- list()
for (sr in names(SERIES)) for (oc in c("psm_pt2", "or")) for (dfk in c(3, 4, 5)) {
  per <- lapply(imps, function(di) {
    dd <- SERIES[[sr]](di); dd <- dd[outcome_rows(dd, oc), ]
    cov <- if (oc == "or") COV_OR else COV_PT2
    f1 <- as.formula(paste(if (oc == "or") "log_or" else "psm", "~ ns(x, df =", dfk, ") +", cov))
    f0 <- as.formula(paste(if (oc == "or") "log_or" else "psm", "~", cov))
    if (oc == "or") { m1 <- lm(f1, data = dd); m0 <- lm(f0, data = dd); lr <- as.numeric(2 * (logLik(m1) - logLik(m0))) }
    else { m1 <- glm(f1, data = dd, family = binomial); m0 <- glm(f0, data = dd, family = binomial); lr <- m0$deviance - m1$deviance }
    idx <- grep("^ns\\(x", names(coef(m1)))
    list(q = coef(m1)[idx], u = vcov(m1)[idx, idx], lr = lr)
  })
  Q <- do.call(rbind, lapply(per, `[[`, "q")); U <- lapply(per, `[[`, "u")
  out$fixed_basis[[paste(sr, oc, paste0("df", dfk))]] <- list(D1 = as.list(d1_pool(Q, U)), D2_lr = as.list(d2_pool(sapply(per, `[[`, "lr"), dfk)),
                                                               lr_p_per_imputation_median = median(pchisq(sapply(per, `[[`, "lr"), dfk, lower.tail = FALSE)))
}

# ---------- 3. 2024 restricted to hospital-months with both robots (descriptive) ----------
d$month <- as.integer(format(d$surgery_date, "%m"))
d$robot <- ifelse(d$platform == "Versius", "Versius", ifelse(!is.na(d$dv_restart_n), "restart", NA))
y <- d[d$year == 2024 & !is.na(d$robot), ]
cells <- aggregate(robot ~ centre + month, data = y, FUN = function(r) length(unique(r))); names(cells)[3] <- "n_robots"
y <- merge(y, cells, by = c("centre", "month")); ym <- y[y$n_robots == 2, ]
bn <- function(x) { x <- x[!is.na(x)]; c(events = sum(x), n = length(x)) }
tt <- t.test(ym$or_time[ym$robot == "restart"], ym$or_time[ym$robot == "Versius"])  # restart minus Versius
out$hospital_months_2024 <- list(cells = unique(ym[, c("centre", "month")]), n = nrow(ym), by_robot = as.list(table(ym$robot)),
  or_time_mean = as.list(tapply(ym$or_time, ym$robot, mean, na.rm = TRUE)), or_time_median = as.list(tapply(ym$or_time, ym$robot, median, na.rm = TRUE)),
  or_time_q1 = as.list(tapply(ym$or_time, ym$robot, quantile, .25, na.rm = TRUE)), or_time_q3 = as.list(tapply(ym$or_time, ym$robot, quantile, .75, na.rm = TRUE)),
  unadjusted_mean_difference_restart_minus_versius = unname(tt$estimate[1] - tt$estimate[2]), welch_ci95 = unname(as.numeric(tt$conf.int)),
  psm_pt2 = list(Versius = bn(ym$psm_pt2[ym$robot == "Versius"]), restart = bn(ym$psm_pt2[ym$robot == "restart"])),
  psm_all = list(Versius = bn(ym$psm[ym$robot == "Versius"]), restart = bn(ym$psm[ym$robot == "restart"])))
out$y2024_crude <- list(or_time_mean = as.list(tapply(y$or_time, y$robot, mean, na.rm = TRUE)), or_time_median = as.list(tapply(y$or_time, y$robot, median, na.rm = TRUE)),
  by_hospital_n = as.list(table(paste(y$centre, y$robot))),
  or_time_mean_by_hospital = as.list(tapply(y$or_time, paste(y$centre, y$robot), mean, na.rm = TRUE)))
write_json(out, file.path(WS, "analysis/results/R_revision_checks.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("revision checks done\n")

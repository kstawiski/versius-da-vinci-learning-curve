# Shared functions for manuscript6 learning-curve analyses. Aggregate outputs only.
suppressPackageStartupMessages({ library(arrow); library(mgcv); library(mice); library(sandwich); library(lmtest); library(splines) })

COV    <- "age + bmi + log_psa + isup_cat + log_weight + pT_group + nerve_sparing + plnd + centre"
COV_OR <- paste(COV, "+ first_case")
COV_PT2 <- "age + bmi + log_psa + isup_cat + log_weight + centre"

load_ms6 <- function(WS, restart_at = 20L) {
  raw <- as.data.frame(read_parquet(file.path(WS, "analysis/restricted/ms6_analysis_dataset.parquet")))
  # data close: latest dated postoperative contact in the Polish series
  # data close: the latest dated postoperative contact in the Polish series (source date field)
  close <- max(raw$fu_last_contact_date, na.rm = TRUE)
  d <- raw[raw$robotic == 1, ]
  d$platform <- factor(d$platform, levels = c("Versius", "da Vinci"))
  d$centre <- factor(d$centre, levels = c("SM", "BE"))
  d$log_psa <- log(d$psa)
  d$log_weight <- log(d$specimen_weight_g)
  d$isup_cat <- factor(ifelse(is.na(d$biopsy_isup), NA, pmin(d$biopsy_isup, 4)), levels = 1:4)
  d$pT_group <- factor(d$pT_group, levels = c("pT2", "pT3a", "pT3b+"))
  d$pT_obs <- d$pT_group  # observed stage; defines the pT2 analysis set (never imputed)
  d$first_case <- as.integer(d$day_position == 1)
  d$log_or <- log(d$or_time)
  d$dv_restart_n <- ifelse(d$platform == "da Vinci" & d$platform_n >= restart_at, d$platform_n - restart_at + 1, NA)
  days_open <- as.numeric(difftime(close, d$surgery_date, units = "days"))
  d$elig_3m <- as.integer(days_open >= 120)
  d$elig_12m <- as.integer(days_open >= 425)
  d$elig_psa56 <- as.integer(days_open >= 56)
  list(d = d, close = close, close_gap_days = as.numeric(difftime(close, max(raw$surgery_date), units = "days")))
}

imp_vars <- c("age", "bmi", "log_psa", "isup_cat", "log_weight", "pT_group", "nerve_sparing", "plnd",
              "centre", "platform", "platform_n", "robotic_n", "log_or", "psm", "first_case",
              "los", "ebl", "pad_free_3m", "pad_free_12m", "psa_persistence_eau")
AUX_OUTCOMES <- c("log_or", "psm", "los", "ebl", "pad_free_3m", "pad_free_12m", "psa_persistence_eau")
impute <- function(dd, m, seed) {
  x <- dd[, imp_vars]
  for (v in c("nerve_sparing", "plnd", "psm", "pad_free_3m", "pad_free_12m", "psa_persistence_eau")) x[[v]] <- factor(x[[v]])
  if (length(unique(x$platform)) < 2) x$platform <- NULL
  # drop degenerate auxiliaries (all missing, or a single observed level) so they cannot block imputation
  for (v in intersect(AUX_OUTCOMES, names(x))) {
    obs <- x[[v]][!is.na(x[[v]])]
    if (length(obs) == 0 || length(unique(obs)) < 2) x[[v]] <- NULL else if (is.factor(x[[v]])) x[[v]] <- droplevels(x[[v]])
  }
  meth <- make.method(x)
  for (v in names(meth)) if (!any(is.na(x[[v]]))) meth[v] <- ""
  if ("bmi" %in% names(meth) && any(is.na(x$bmi))) meth["bmi"] <- "pmm"
  if (any(is.na(x$log_weight))) meth["log_weight"] <- "pmm"
  if (any(is.na(x$isup_cat))) meth["isup_cat"] <- "polr"
  if (any(is.na(x$pT_group))) meth["pT_group"] <- "polr"
  # outcomes imputed only as auxiliaries so every covariate row can be completed and outcome-covariate
  # relationships are preserved; imputed outcome values are never copied back into the analysis data
  for (v in intersect(AUX_OUTCOMES, names(x))) if (any(is.na(x[[v]])))
    meth[v] <- if (is.factor(x[[v]])) (if (nlevels(x[[v]]) == 2) "logreg" else "") else "pmm"
  mids <- mice(x, m = m, method = meth, seed = seed, printFlag = FALSE)
  lapply(seq_len(m), function(i) {
    ci <- complete(mids, i); out <- dd
    for (v in c("bmi", "log_weight", "isup_cat", "pT_group")) out[[v]] <- ci[[v]]
    out
  })
}

# HC3 robust variance; HC1 fallback only when a leverage-1 observation makes HC3 undefined
robust_vcov <- function(m) {
  V <- sandwich::vcovHC(m, type = if (inherits(m, "glm")) "HC0" else "HC3")
  if (!all(is.finite(V))) V <- sandwich::vcovHC(m, type = "HC1")
  V
}
rubin <- function(est, var) {
  m <- length(est); qb <- mean(est); ub <- mean(var); b <- if (m > 1) var(est) else 0
  tv <- ub + (1 + 1 / m) * b
  r <- if (ub > 0) (1 + 1 / m) * b / ub else 0
  df <- if (m > 1 && r > 0) (m - 1) * (1 + 1 / r)^2 else 1e6
  c(est = qb, se = sqrt(tv), df = min(df, 1e6))
}

# curve models: s(x) where dd$x holds the series case number
outcome_rows <- function(dd, oc) if (oc == "or") !is.na(dd$log_or) else (dd$pT_obs %in% "pT2" & !is.na(dd$psm))
fit_curve <- function(dd, oc, k, cov_or = COV_OR, cov_pt2 = COV_PT2, spline = "gam") {
  dd <- dd[outcome_rows(dd, oc), ]
  sm <- if (spline == "gam") paste0("s(x, k = ", k, ")") else {
    kn <- if (oc == "or") quantile(dd$x, c(.05, .35, .65, .95)) else quantile(dd$x, c(.10, .50, .90))
    paste0("ns(x, knots = c(", paste(kn[2:(length(kn) - 1)], collapse = ","), "), Boundary.knots = c(", kn[1], ",", kn[length(kn)], "))")
  }
  if (oc == "or") gam(as.formula(paste("log_or ~", sm, "+", cov_or)), data = dd, method = "REML")
  else gam(as.formula(paste("psm ~", sm, "+", cov_pt2)), data = dd, family = binomial, method = "REML")
}
std_newdata <- function(dd, oc, ns) {
  base <- dd[outcome_rows(dd, oc), ]; nb <- nrow(base)
  nd <- base[rep(seq_len(nb), times = length(ns)), ]; nd$x <- rep(ns, each = nb)
  list(nd = nd, idx = rep(ns, each = nb))
}
std_curve <- function(fit, dd, oc, ns) {
  z <- std_newdata(dd, oc, ns)
  p <- if (oc == "or") exp(predict(fit, z$nd, type = "link")) * mean(exp(residuals(fit))) else predict(fit, z$nd, type = "response")
  tapply(p, z$idx, mean)
}
sim_curve <- function(fit, dd, oc, ns, nsim = 100) {
  z <- std_newdata(dd, oc, ns)
  X <- predict(fit, z$nd, type = "lpmatrix")
  B <- rmvn(nsim, coef(fit), vcov(fit, unconditional = TRUE))
  lp <- X %*% t(B)
  p <- if (oc == "or") exp(lp) * mean(exp(residuals(fit))) else plogis(lp)
  apply(p, 2, function(col) tapply(col, z$idx, mean))
}
plateau <- function(curve, margin, tail_n = 50) {
  ns <- as.numeric(names(curve)); asym <- mean(tail(curve, tail_n))
  outside <- abs(curve - asym) > margin
  if (!any(outside)) return(list(nstar = 1, asym = asym))
  last_out <- max(which(outside))
  list(nstar = if (last_out >= length(ns)) NA else ns[last_out + 1], asym = asym)
}

# adjusted contrasts between two groups defined by gdef(dd) -> "ref"/"cmp"/NA
fit_contrast <- function(di, gdef, outcome, cov, subset_expr = NULL) {
  dd <- di; dd$grp <- gdef(dd); dd <- dd[!is.na(dd$grp), ]
  if (!is.null(subset_expr)) dd <- dd[subset_expr(dd), ]
  dd <- dd[!is.na(dd[[outcome]]), ]
  dd$grp <- factor(dd$grp, levels = c("ref", "cmp"))
  if (min(table(dd$grp)) < 5 || (!(outcome %in% c("or_time", "los", "ebl")) && sum(dd[[outcome]]) < 1))
    return(list(est = NA_real_, var = NA_real_, n_ref = sum(dd$grp == "ref"), n_cmp = sum(dd$grp == "cmp"), not_estimable = TRUE))
  # drop factor covariates with a single level in this subset
  cv <- strsplit(cov, " \\+ ")[[1]]
  cv <- cv[sapply(cv, function(v) { v0 <- sub("^ns\\((\\w+).*", "\\1", v); !(v0 %in% names(dd)) || length(unique(dd[[v0]])) > 1 })]
  f <- as.formula(paste(outcome, "~ grp +", paste(cv, collapse = " + ")))
  base <- list(n_ref = sum(dd$grp == "ref"), n_cmp = sum(dd$grp == "cmp"),
               crude_ref = mean(dd[[outcome]][dd$grp == "ref"]), crude_cmp = mean(dd[[outcome]][dd$grp == "cmp"]))
  if (outcome %in% c("or_time", "los", "ebl")) {
    m <- lm(f, data = dd); ct <- coeftest(m, vcov = robust_vcov(m))["grpcmp", ]
    c(list(est = unname(ct[1]), var = unname(ct[2]^2)), base,
      list(median_ref = median(dd[[outcome]][dd$grp == "ref"]), median_cmp = median(dd[[outcome]][dd$grp == "cmp"])))
  } else {
    m <- suppressWarnings(glm(f, data = dd, family = binomial))
    d1 <- dd; d1$grp <- factor("cmp", levels = c("ref", "cmp")); d0 <- dd; d0$grp <- factor("ref", levels = c("ref", "cmp"))
    tt <- delete.response(terms(m))
    X1 <- model.matrix(tt, d1, contrasts.arg = m$contrasts, xlev = m$xlevels); X0 <- model.matrix(tt, d0, contrasts.arg = m$contrasts, xlev = m$xlevels)
    b <- coef(m); keep <- !is.na(b); b <- b[keep]; X1 <- X1[, keep, drop = FALSE]; X0 <- X0[, keep, drop = FALSE]
    p1 <- plogis(X1 %*% b); p0 <- plogis(X0 %*% b)
    g <- colMeans(as.vector(p1 * (1 - p1)) * X1) - colMeans(as.vector(p0 * (1 - p0)) * X0)
    V <- robust_vcov(m)[keep, keep]
    c(list(est = mean(p1) - mean(p0), var = as.numeric(t(g) %*% V %*% g),
           or = unname(exp(coef(m)["grpcmp"])), events_ref = sum(dd[[outcome]][dd$grp == "ref"]),
           events_cmp = sum(dd[[outcome]][dd$grp == "cmp"])), base)
  }
}
pool_contrast <- function(imps, label, gdef, outcome, cov, subset_expr = NULL, eq_margin = NA) {
  rs <- lapply(imps, fit_contrast, gdef = gdef, outcome = outcome, cov = cov, subset_expr = subset_expr)
  if (any(sapply(rs, function(r) isTRUE(r$not_estimable) || !is.finite(r$var))))
    return(list(label = label, outcome = outcome, not_estimable = TRUE, n_ref = rs[[1]]$n_ref, n_cmp = rs[[1]]$n_cmp))
  pr <- rubin(sapply(rs, `[[`, "est"), sapply(rs, `[[`, "var"))
  q95 <- qt(.975, pr["df"]); q90 <- qt(.95, pr["df"])
  out <- list(label = label, outcome = outcome, estimate = unname(pr["est"]), se = unname(pr["se"]),
              ci95 = unname(pr["est"] + c(-1, 1) * q95 * pr["se"]), ci90 = unname(pr["est"] + c(-1, 1) * q90 * pr["se"]),
              p = unname(2 * pt(-abs(pr["est"] / pr["se"]), pr["df"])))
  for (nm in setdiff(names(rs[[1]]), c("est", "var"))) out[[nm]] <- if (nm == "or") median(sapply(rs, `[[`, "or")) else rs[[1]][[nm]]
  if (!is.na(eq_margin)) { out$eq_margin <- eq_margin; out$equivalent <- (out$ci90[1] > -eq_margin) && (out$ci90[2] < eq_margin) }
  out
}
rare_summary <- function(d, gdef, outcomes) {
  res <- list(); dd <- d; dd$grp <- gdef(dd)
  for (oc in outcomes) for (gname in c("ref", "cmp")) {
    x <- dd[[oc]][!is.na(dd$grp) & dd$grp == gname]; nmiss <- sum(is.na(x)); x <- x[!is.na(x)]
    bt <- binom.test(sum(x), max(1, length(x)))
    res[[paste(oc, gname)]] <- list(outcome = oc, group = gname, events = sum(x), n = length(x), n_missing = nmiss, ci95 = unname(bt$conf.int))
  }
  res
}
observation_summary <- function(d, gdef) {
  dd <- d; dd$grp <- gdef(dd); dd <- dd[!is.na(dd$grp), ]
  f <- function(oc, el) sapply(c("ref", "cmp"), function(g) { s <- dd$grp == g & dd[[el]] == 1; c(eligible = sum(s), observed = sum(s & !is.na(dd[[oc]]))) })
  list(pad_free_3m = f("pad_free_3m", "elig_3m"), pad_free_12m = f("pad_free_12m", "elig_12m"), psa_persistence_eau = f("psa_persistence_eau", "elig_psa56"))
}

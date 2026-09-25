#!/usr/bin/env Rscript
# Implementation amendment (SAP_CHANGELOG): sparse binary outcomes.
# When either comparison group has fewer than 10 events (or events equal the group size), the full
# covariate set separates. Such contrasts use Firth penalized logistic regression (logistf) with
# parsimonious adjustment (age, log PSA, hospital) and a marginally standardized risk difference with a
# delta-method variance from the penalized covariance. Zero-event outcomes are reported crude only.
# Applies to A4 (2024 concurrent) and B1 (all-era mature) binary contrasts. Aggregate output only.
suppressPackageStartupMessages({ library(jsonlite); library(logistf) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
dat <- load_ms6(WS, 20L); d <- dat$d
imps <- impute(d, 20, 11)
pl <- fromJSON(file.path(WS, "analysis/results/A2_plateaus.json"))
nv <- pl[["Versius or"]]$nstar; ndr <- pl[["dVrestart or"]]$nstar
cut_v <- if (is.null(nv) || is.na(nv)) Inf else nv; cut_dr <- if (is.null(ndr) || is.na(ndr)) Inf else if (ndr <= 1) 0 else ndr
G <- list(
  A4 = function(dd) ifelse(dd$year == 2024 & dd$platform == "Versius", "ref", ifelse(dd$year == 2024 & !is.na(dd$dv_restart_n), "cmp", NA)),
  B1 = function(dd) ifelse(dd$platform == "Versius" & dd$platform_n > cut_v, "ref", ifelse(!is.na(dd$dv_restart_n) & dd$dv_restart_n > cut_dr, "cmp", NA)))
# rare safety events (grade >=III complications, transfusion) are reported as exact counts only, as prespecified
OUTC <- list(psa_persistence_eau = "elig_psa56", pad_free_3m = "elig_3m", pad_free_12m = "elig_12m", readmission_30d = NULL, psm = NULL)
firth_rd <- function(dd, outcome) {
  dd$grp <- factor(dd$grp, levels = c("ref", "cmp"))
  rhs <- if (length(unique(dd$centre)) > 1) "grp + age + log_psa + centre" else "grp + age + log_psa"
  dd$centre <- droplevels(dd$centre)
  m <- logistf(as.formula(paste(outcome, "~", rhs)), data = dd, plcontrol = logistf.control(maxit = 200))
  X1 <- model.matrix(as.formula(paste("~", rhs)), transform(dd, grp = factor("cmp", levels = c("ref", "cmp"))))
  X0 <- model.matrix(as.formula(paste("~", rhs)), transform(dd, grp = factor("ref", levels = c("ref", "cmp"))))
  b <- coef(m); p1 <- plogis(X1 %*% b); p0 <- plogis(X0 %*% b)
  g <- colMeans(as.vector(p1 * (1 - p1)) * X1) - colMeans(as.vector(p0 * (1 - p0)) * X0)
  list(est = mean(p1) - mean(p0), var = as.numeric(t(g) %*% vcov(m) %*% g))
}
res <- list(label = "Sparse binary contrasts: Firth logistic, adjusted for age, log PSA and hospital; marginal risk difference, da Vinci minus Versius")
for (cmp in names(G)) for (oc in names(OUTC)) {
  el <- OUTC[[oc]]
  base <- d; base$grp <- G[[cmp]](base); base <- base[!is.na(base$grp) & !is.na(base[[oc]]), ]
  if (!is.null(el)) base <- base[base[[el]] == 1, ]
  ev <- tapply(base[[oc]], base$grp, sum); nn <- tapply(base[[oc]], base$grp, length)
  entry <- list(outcome = oc, events_ref = unname(ev["ref"]), n_ref = unname(nn["ref"]), events_cmp = unname(ev["cmp"]), n_cmp = unname(nn["cmp"]))
  if (length(nn) < 2 || any(is.na(nn)) || min(nn) < 5) { entry$note <- "group too small or absent"; res[[paste(cmp, oc)]] <- entry; next }
  if (sum(ev) == 0) { entry$note <- "no events in either group; crude only"; res[[paste(cmp, oc)]] <- entry; next }
  entry$sparse_rule_applies <- min(ev) < 10 || any(ev == nn)
  if (!entry$sparse_rule_applies) { entry$note <- "sparse rule does not apply; the full-model estimate is authoritative"; res[[paste(cmp, oc)]] <- entry; next }
  rs <- lapply(imps, function(di) {
    dd <- di; dd$grp <- G[[cmp]](dd); dd <- dd[!is.na(dd$grp) & !is.na(dd[[oc]]), ]
    if (!is.null(el)) dd <- dd[dd[[el]] == 1, ]
    firth_rd(dd, oc)
  })
  pr <- rubin(sapply(rs, `[[`, "est"), sapply(rs, `[[`, "var")); q <- qt(.975, pr["df"]); q9 <- qt(.95, pr["df"])
  entry$rd <- unname(pr["est"]); entry$ci95 <- unname(pr["est"] + c(-1, 1) * q * pr["se"]); entry$ci90 <- unname(pr["est"] + c(-1, 1) * q9 * pr["se"])
  entry$p <- unname(2 * pt(-abs(pr["est"] / pr["se"]), pr["df"]))
  res[[paste(cmp, oc)]] <- entry
}
write_json(res, file.path(WS, "analysis/results/A4_B1_sparse_binary.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cat("sparse done\n")

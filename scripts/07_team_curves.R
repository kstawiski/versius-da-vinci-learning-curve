#!/usr/bin/env Rscript
# B2 prespecified fallback (surgeon and team experience correlated above 0.8): hospital-stratified,
# case-mix-standardized operative-time curves against each hospital team's own case number on each platform.
# da Vinci curves use restart-series cases only. Aggregate output; Figure 4.
suppressPackageStartupMessages({ library(jsonlite); library(ggplot2); library(patchwork); library(dplyr) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R")); source(file.path(WS, "scripts/fig_style.R"))
dat <- load_ms6(WS, 20L); d <- dat$d
imps <- impute(d, 20, 11)
COV_H <- "age + bmi + log_psa + isup_cat + log_weight + pT_group + nerve_sparing + plnd + first_case"
sel <- list(
  "Versius|SM" = function(dd) dd[dd$platform == "Versius" & dd$centre == "SM", ],
  "Versius|BE" = function(dd) dd[dd$platform == "Versius" & dd$centre == "BE", ],
  "da Vinci|SM" = function(dd) dd[!is.na(dd$dv_restart_n) & dd$centre == "SM", ],
  "da Vinci|BE" = function(dd) dd[!is.na(dd$dv_restart_n) & dd$centre == "BE", ])
curves <- list(); summ <- list()
for (k in names(sel)) {
  base <- sel[[k]](d); base$x <- base$team_platform_n
  ns <- seq(min(base$x), max(base$x))
  kk <- min(8, floor(length(unique(base$x)) / 10))
  fits <- lapply(imps, function(di) { dd <- sel[[k]](di); dd$x <- dd$team_platform_n
    list(fit = gam(as.formula(paste0("log_or ~ s(x, k = ", kk, ") + ", COV_H)), data = dd[!is.na(dd$log_or), ], method = "REML"), dd = dd) })
  cm <- rowMeans(sapply(fits, function(f) std_curve(f$fit, f$dd, "or", ns)))
  sims <- do.call(cbind, lapply(fits, function(f) sim_curve(f$fit, f$dd, "or", ns, 100)))
  pf <- strsplit(k, "\\|")[[1]]
  sc <- sel[[k]](d)
  curves[[k]] <- data.frame(platform = pf[1], hospital = pf[2], x = ns, est = cm, lo = apply(sims, 1, quantile, .025), hi = apply(sims, 1, quantile, .975))
  summ[[k]] <- list(n = nrow(sc), team_case_range = range(sc$team_platform_n), surgeon_platform_case_at_team_start = min(sc$platform_n),
                    surgeon_prior_versius_at_team_start = min(sc$prior_v), first_team_case = min(ns), fitted_first = unname(cm[1]),
                    fitted_team_case_25 = if (25 %in% ns) unname(cm[match(25, ns)]) else NA,
                    fitted_team_case_50 = if (50 %in% ns) unname(cm[match(50, ns)]) else NA,
                    fitted_team_case_100 = if (100 %in% ns) unname(cm[match(100, ns)]) else NA,
                    fitted_24_cases_after_start = unname(cm[min(25, length(cm))]), last_team_case = max(ns), fitted_last = unname(cm[length(cm)]),
                    fitted_min = min(cm), team_case_at_min = ns[which.min(cm)],
                    smooth_p = median(sapply(fits, function(f) summary(f$fit)$s.table[1, "p-value"])))
}
cv <- do.call(rbind, curves)
write.csv(cv, file.path(WS, "analysis/results/B2_team_curves.csv"), row.names = FALSE)
write_json(summ, file.path(WS, "analysis/results/B2_team_curves_summary.json"), auto_unbox = TRUE, digits = 6, pretty = TRUE)
cv$Hospital <- factor(ifelse(cv$hospital == "SM", "SalveMedica", "Bełchatów"), levels = c("SalveMedica", "Bełchatów"))
cv$Panel <- factor(ifelse(cv$platform == "Versius", "A  Versius", "B  da Vinci restart"), levels = c("A  Versius", "B  da Vinci restart"))
lab <- cv %>% group_by(Panel, Hospital) %>% filter(x == min(x))
p <- ggplot(cv, aes(x, est, colour = Hospital, fill = Hospital, linetype = Hospital)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.8) +
  facet_wrap(~Panel, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = c("SalveMedica" = "#1a1a1a", "Bełchatów" = "#B45A1E"), name = NULL) +
  scale_fill_manual(values = c("SalveMedica" = "#1a1a1a", "Bełchatów" = "#B45A1E"), name = NULL) +
  scale_linetype_manual(values = c("SalveMedica" = "solid", "Bełchatów" = "42"), name = NULL) +
  scale_y_continuous(limits = c(100, 320)) +
  labs(x = "Hospital team's case number on the platform", y = "Operative time, minutes") +
  theme_postcddp() + theme(legend.position = "top", strip.text = element_text(size = 9, hjust = 0))
# Team curves are displayed as Figure 2C-D (scripts/05_fig2_curves.R); no standalone figure.
cat("team curves done\n")

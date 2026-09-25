#!/usr/bin/env Rscript
# Figure 2. Case-mix-standardized learning curves for operative time (A) and pT2 positive margins (B):
# Versius from its first case and da Vinci from the 2024 restart, on a shared case-number axis.
# Points are observed means (A) or proportions (B) in consecutive blocks of 25 cases (aggregate).
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(arrow); library(dplyr); library(jsonlite) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/fig_style.R"))
cv <- read.csv(file.path(WS, "analysis/results/A1_standardized_curves.csv"))
cv$Series <- factor(ifelse(cv$series == "Versius", "Versius", "da Vinci (restart)"), levels = c("Versius", "da Vinci (restart)"))
PAL <- c("Versius" = unname(PAL_PLATFORM["Versius"]), "da Vinci (restart)" = unname(PAL_PLATFORM["da Vinci"]))
raw <- as.data.frame(read_parquet(file.path(WS, "analysis/restricted/ms6_analysis_dataset.parquet")))
raw <- raw[raw$robotic == 1, ]
raw$x <- ifelse(raw$platform == "Versius", raw$platform_n, ifelse(raw$platform_n >= 20, raw$platform_n - 19, NA))
raw <- raw[!is.na(raw$x), ]
raw$Series <- factor(ifelse(raw$platform == "Versius", "Versius", "da Vinci (restart)"), levels = levels(cv$Series))
raw$blk <- (raw$x - 1) %/% 25
bo <- raw %>% filter(!is.na(or_time)) %>% group_by(Series, blk) %>%
  summarise(x = mean(x), mean = mean(or_time), n = n(), .groups = "drop") %>% filter(n >= 10)
bp <- raw %>% filter(pT_group == "pT2", !is.na(psm)) %>% group_by(Series, blk) %>%
  summarise(x = mean(x), p = mean(psm), events = sum(psm), n = n(), .groups = "drop") %>% filter(n >= 8)
dir.create(file.path(WS, "figures/source_data"), showWarnings = FALSE)
write.csv(bo, file.path(WS, "figures/source_data/fig2A_block_means.csv"), row.names = FALSE)
write.csv(bp, file.path(WS, "figures/source_data/fig2B_block_proportions.csv"), row.names = FALSE)

co <- cv[cv$outcome == "or", ]
pA <- ggplot(co, aes(n, est, colour = Series, fill = Series)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(data = bo, aes(x, mean, shape = Series), size = 1.4, inherit.aes = FALSE, colour = "#3b3b3b", fill = "white", stroke = 0.4) +
  scale_colour_manual(values = PAL, name = NULL) + scale_fill_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_x_continuous(breaks = c(1, 50, 100, 150, 200, 250, 300, 350), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(100, 350), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Case number on the platform", y = "Operative time, minutes", title = "A") +
  theme_postcddp() + theme(legend.position = "top")
cp <- cv[cv$outcome == "psm_pt2", ]
pB <- ggplot(cp, aes(n, est * 100, colour = Series, fill = Series)) +
  geom_ribbon(aes(ymin = lo * 100, ymax = hi * 100), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(data = bp, aes(x, p * 100, shape = Series), size = 1.4, inherit.aes = FALSE, colour = "#3b3b3b", fill = "white", stroke = 0.4) +
  scale_colour_manual(values = PAL, name = NULL) + scale_fill_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_x_continuous(breaks = c(1, 50, 100, 150, 200, 250, 300, 350), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 80), expand = expansion(mult = c(0.03, 0.02))) +
  labs(x = "Case number on the platform", y = "Positive margin in pT2 disease, %", title = "B") +
  theme_postcddp() + theme(legend.position = "top")
# C and D: hospital-team curves (prespecified display when surgeon and team experience are collinear)
tc <- read.csv(file.path(WS, "analysis/results/B2_team_curves.csv"))
tc$Hospital <- factor(ifelse(tc$hospital == "SM", "SalveMedica", "Bełchatów"), levels = c("SalveMedica", "Bełchatów"))
HPAL <- c("SalveMedica" = "#1a1a1a", "Bełchatów" = "#7a7a7a")
team_panel <- function(pf, title) {
  ggplot(tc[tc$platform == pf, ], aes(x, est, colour = Hospital, fill = Hospital, linetype = Hospital)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.8) +
    scale_colour_manual(values = HPAL, name = NULL) + scale_fill_manual(values = HPAL, name = NULL) +
    scale_linetype_manual(values = c("SalveMedica" = "solid", "Bełchatów" = "42"), name = NULL) +
    scale_x_continuous(breaks = c(1, 50, 100, 150, 200), expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(limits = c(100, 350), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Hospital's case number on the platform", y = "Operative time, minutes", title = title) +
    theme_postcddp() + theme(legend.position = "top")
}
pC <- team_panel("Versius", "C  Versius, by hospital"); pD <- team_panel("da Vinci", "D  da Vinci restart, by hospital")
pA <- pA + labs(title = "A  Operative time") + theme(legend.position = "top")
pB <- pB + labs(title = "B  Positive margin, pT2") + theme(legend.position = "top")
fig <- ((pA | pB) + plot_layout(guides = "collect") & theme(legend.position = "top")) /
       ((pC | pD) + plot_layout(guides = "collect") & theme(legend.position = "top"))
save_fig(fig, "fig2_learning_curves", width = 174, height = 150)
cat("fig2 done\n")

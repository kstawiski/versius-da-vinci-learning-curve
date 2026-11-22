#!/usr/bin/env Rscript
# Figure 1. Quarterly operations (A) and case-mix-standardized learning curves for operative time (B)
# and pT2 positive margins (C), from the first Versius operation and the 2024 da Vinci restart.
# Points are observed means (B) or proportions (C) in consecutive blocks of 25 operations.
# Panels D and E show operative time by hospital.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(arrow); library(dplyr); library(jsonlite) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/fig_style.R"))
cv <- read.csv(file.path(WS, "analysis/results/A1_standardized_curves.csv"))
cv$Series <- factor(ifelse(cv$series == "Versius", "Versius", "da Vinci (restart)"), levels = c("Versius", "da Vinci (restart)"))
PAL <- c("Versius" = unname(PAL_PLATFORM["Versius"]), "da Vinci (restart)" = unname(PAL_PLATFORM["da Vinci"]))
raw <- as.data.frame(read_parquet(file.path(WS, "analysis/restricted/ms6_analysis_dataset.parquet")))
dir.create(file.path(WS, "figures/source_data"), showWarnings = FALSE)
HOSPITAL <- c(SM = "SalveMedica (da Vinci already in use)", BE = "Bełchatów (new to da Vinci)")

# Quarterly counts across all platforms and both hospitals.
quarterly <- raw
quarterly$q <- as.Date(cut(as.Date(quarterly$surgery_date), "quarter"))
quarterly$platform <- factor(quarterly$platform, levels = c("Versius", "da Vinci", "Laparoscopic"))
quarterly$hospital <- factor(ifelse(quarterly$centre == "SM", "SalveMedica, Łódź", "Bełchatów"), levels = c("SalveMedica, Łódź", "Bełchatów"))
qc <- quarterly %>% count(hospital, q, platform, .drop = FALSE)
write.csv(qc %>% mutate(q = paste0(format(as.Date(q), "%Y"), "-Q", (as.integer(format(as.Date(q), "%m")) - 1) %/% 3 + 1)), file.path(WS, "figures/source_data/fig1A_quarterly_counts.csv"), row.names = FALSE)
qc$q <- as.Date(qc$q)
pA <- ggplot(qc, aes(q + 45, n, fill = platform)) +
  geom_col(width = 80, colour = "white", linewidth = 0.15) +
  facet_wrap(~hospital, nrow = 1, labeller = as_labeller(setNames(unname(HOSPITAL), levels(qc$hospital)))) +
  scale_fill_manual(values = PAL_PLATFORM, name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Operations per quarter", title = "A  Operations per quarter") +
  theme_postcddp() + theme(legend.position = "top", panel.grid.major.x = element_blank(), panel.spacing.x = unit(10, "pt"))

raw <- raw[raw$robotic == 1, ]
raw$x <- ifelse(raw$platform == "Versius", raw$platform_n, ifelse(raw$platform_n >= 20, raw$platform_n - 19, NA))
raw <- raw[!is.na(raw$x), ]
raw$Series <- factor(ifelse(raw$platform == "Versius", "Versius", "da Vinci (restart)"), levels = levels(cv$Series))
raw$blk <- (raw$x - 1) %/% 25
bo <- raw %>% filter(!is.na(or_time)) %>% group_by(Series, blk) %>%
  summarise(x = mean(x), mean = mean(or_time), n = n(), .groups = "drop") %>% filter(n >= 10)
bp <- raw %>% filter(pT_group == "pT2", !is.na(psm)) %>% group_by(Series, blk) %>%
  summarise(x = mean(x), p = mean(psm), events = sum(psm), n = n(), .groups = "drop") %>% filter(n >= 8)
write.csv(bo, file.path(WS, "figures/source_data/fig2A_block_means.csv"), row.names = FALSE)
write.csv(bp, file.path(WS, "figures/source_data/fig2B_block_proportions.csv"), row.names = FALSE)

breakpoint <- fromJSON(file.path(WS, "analysis/results/secondary_sensitivity.json"))$S_seg$Versius$breakpoint
if (!is.numeric(breakpoint) || length(breakpoint) != 1L || !is.finite(breakpoint))
  stop("Missing adjusted Versius breakpoint in secondary_sensitivity.json")
co <- cv[cv$outcome == "or", ]
pB <- ggplot(co, aes(n, est, colour = Series, fill = Series)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.18, colour = NA) +
  geom_vline(xintercept = breakpoint, colour = PAL[["Versius"]], linewidth = 0.3, linetype = "dashed") +
  annotate("text", x = breakpoint + 8, y = 330, label = sprintf("Versius breakpoint, case %d", round(breakpoint)),
           hjust = 0, size = 2.3, family = "Open Sans", colour = PAL[["Versius"]]) +
  geom_line(linewidth = 0.8) +
  geom_point(data = bo, aes(x, mean, shape = Series), size = 1.4, inherit.aes = FALSE, colour = "#3b3b3b", fill = "white", stroke = 0.4) +
  scale_colour_manual(values = PAL, name = NULL) + scale_fill_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_x_continuous(breaks = c(1, 50, 100, 150, 200, 250, 300, 350), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(100, 350), expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Case number on the platform", y = "Operative time, minutes", title = "B  Operative time") +
  theme_postcddp() + theme(legend.position = "top")
cp <- cv[cv$outcome == "psm_pt2", ]
pC <- ggplot(cp, aes(n, est * 100, colour = Series, fill = Series)) +
  geom_ribbon(aes(ymin = lo * 100, ymax = hi * 100), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_point(data = bp, aes(x, p * 100, shape = Series), size = 1.4, inherit.aes = FALSE, colour = "#3b3b3b", fill = "white", stroke = 0.4) +
  scale_colour_manual(values = PAL, name = NULL) + scale_fill_manual(values = PAL, name = NULL) +
  scale_shape_manual(values = c(21, 24), name = NULL) +
  scale_x_continuous(breaks = c(1, 50, 100, 150, 200, 250, 300, 350), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(limits = c(0, 80), expand = expansion(mult = c(0.03, 0.02))) +
  labs(x = "Case number on the platform", y = "Positive margin in pT2, %", title = "C  Positive margin, pT2") +
  theme_postcddp() + theme(legend.position = "top")
# D and E: hospital-team curves (prespecified display when surgeon and team experience are collinear)
tc <- read.csv(file.path(WS, "analysis/results/B2_team_curves.csv"))
tc$Hospital <- factor(unname(HOSPITAL[tc$hospital]), levels = unname(HOSPITAL))
HPAL <- setNames(c("#1a1a1a", "#7a7a7a"), HOSPITAL)
HLT <- setNames(c("solid", "42"), HOSPITAL)
team_panel <- function(pf, title) {
  ggplot(tc[tc$platform == pf, ], aes(x, est, colour = Hospital, fill = Hospital, linetype = Hospital)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, colour = NA) + geom_line(linewidth = 0.8) +
    scale_colour_manual(values = HPAL, name = NULL) + scale_fill_manual(values = HPAL, name = NULL) +
    scale_linetype_manual(values = HLT, name = NULL) +
    scale_x_continuous(breaks = c(1, 50, 100, 150, 200), expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_continuous(limits = c(100, 350), expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Hospital's case number on the platform", y = "Operative time, minutes", title = title) +
    theme_postcddp() + theme(legend.position = "top")
}
pD <- team_panel("Versius", "D  Versius, by hospital")
pE <- team_panel("da Vinci", "E  da Vinci restart, by hospital")
row_BC <- (pB | pC) + plot_layout(guides = "collect") & theme(legend.position = "top")
row_DE <- (pD | pE) + plot_layout(guides = "collect") & theme(legend.position = "top")
fig <- pA / row_BC / row_DE + plot_layout(guides = "keep")
save_fig(fig, "fig2_learning_curves", width = 174, height = 205)
cat("fig2 done\n")

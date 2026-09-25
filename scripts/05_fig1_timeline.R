#!/usr/bin/env Rscript
# Figure 1. Structure of the single-surgeon series: quarterly operations by platform and hospital (A)
# and the surgeon's cumulative experience on each platform (B). Aggregate counts only.
suppressPackageStartupMessages({ library(ggplot2); library(ragg); library(patchwork); library(arrow); library(dplyr); library(jsonlite) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/fig_style.R"))
raw <- as.data.frame(read_parquet(file.path(WS, "analysis/restricted/ms6_analysis_dataset.parquet")))
raw$q <- as.Date(cut(as.Date(raw$surgery_date), "quarter"))
raw$platform <- factor(raw$platform, levels = c("Versius", "da Vinci", "Laparoscopic"))
raw$hospital <- factor(ifelse(raw$centre == "SM", "SalveMedica, Łódź", "Bełchatów"), levels = c("SalveMedica, Łódź", "Bełchatów"))
qc <- raw %>% count(hospital, q, platform, .drop = FALSE)
# aggregate table behind panel A
write.csv(qc %>% mutate(q = paste0(format(as.Date(q), "%Y"), "-Q", (as.integer(format(as.Date(q), "%m")) - 1) %/% 3 + 1)), file.path(WS, "figures/source_data/fig1A_quarterly_counts.csv"), row.names = FALSE)
qc$q <- as.Date(qc$q)

pA <- ggplot(qc, aes(q + 45, n, fill = platform)) +
  geom_col(width = 80, colour = "white", linewidth = 0.15) +
  facet_wrap(~hospital, ncol = 1) +
  scale_fill_manual(values = PAL_PLATFORM, name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Operations per quarter", title = "A") +
  theme_postcddp() + theme(legend.position = "top", panel.grid.major.x = element_blank())

# cumulative surgeon experience by platform over calendar time (step function per quarter end)
cum <- raw %>% arrange(series_n) %>% group_by(platform) %>% mutate(cum = row_number()) %>% ungroup()
cq <- cum %>% group_by(platform, q) %>% summarise(cum = max(cum), .groups = "drop")
grid <- expand.grid(platform = levels(raw$platform), q = sort(unique(raw$q)))
cq <- grid %>% left_join(cq, by = c("platform", "q")) %>% group_by(platform) %>% arrange(q) %>%
  mutate(cum = cummax(ifelse(is.na(cum), 0, cum))) %>% ungroup()
cq$platform <- factor(cq$platform, levels = levels(raw$platform))
write.csv(cq %>% mutate(q = paste0(format(as.Date(q), "%Y"), "-Q", (as.integer(format(as.Date(q), "%m")) - 1) %/% 3 + 1)), file.path(WS, "figures/source_data/fig1B_cumulative_by_quarter.csv"), row.names = FALSE)
ends <- cq %>% group_by(platform) %>% filter(q == max(q))
restart_q <- as.Date(cut(as.Date(min(raw$surgery_date[raw$platform == "da Vinci" & raw$platform_n == 20])), "quarter"))
pB <- ggplot(cq, aes(as.Date(q) + 91, cum, colour = platform, linetype = platform)) +
  annotate("rect", xmin = as.Date("2024-01-01"), xmax = as.Date("2025-01-01"), ymin = -Inf, ymax = Inf, fill = "#eeeeee") +
  annotate("text", x = as.Date("2024-07-01"), y = 420, label = "Both robots\nin use (2024)", size = 2.3, colour = "#3b3b3b", lineheight = 0.9) +
  geom_step(linewidth = 0.7, direction = "vh") +
  geom_text(data = ends, aes(x = as.Date("2026-04-29"), label = paste0(platform, " (", cum, ")")), hjust = 0, size = 2.4, show.legend = FALSE, family = "Open Sans") +
  scale_colour_manual(values = PAL_PLATFORM, name = NULL) +
  scale_linetype_manual(values = c("Versius" = "solid", "da Vinci" = "solid", "Laparoscopic" = "22"), name = NULL) +
  scale_x_date(breaks = as.Date(paste0(2021:2026, "-01-01")), date_labels = "%Y",
               limits = c(as.Date("2021-01-01"), as.Date("2026-04-30")), expand = expansion(mult = c(0.01, 0.03))) +
  scale_y_continuous(limits = c(0, 450), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Cumulative operations by the surgeon", title = "B") +
  coord_cartesian(clip = "off") +
  theme_postcddp() + theme(legend.position = "none", plot.margin = margin(2, 25, 2, 2, "mm"))

fig <- pA + pB + plot_layout(widths = c(1.15, 1))
save_fig(fig, "fig1_series_structure", width = 174, height = 88)
cat("fig1 done\n")

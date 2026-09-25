#!/usr/bin/env Rscript
# Figure S2. CUSUM family (prespecified sensitivity): A operative-time CUSUM (deviation from each series mean),
# B risk-adjusted CUSUM (Steiner, odds ratio 2) for positive margin at any stage with its simulated threshold,
# C LC-CUSUM for pT2 margin (acceptable 10%, unacceptable 20%) with its simulated threshold.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/fig_style.R"))
cs <- readRDS(file.path(WS, "analysis/restricted/cusum_series.rds"))
lab <- function(s) factor(ifelse(s == "Versius", "Versius", "da Vinci (restart)"), levels = c("Versius", "da Vinci (restart)"))
PAL <- c("Versius" = unname(PAL_PLATFORM["Versius"]), "da Vinci (restart)" = unname(PAL_PLATFORM["da Vinci"]))
ot <- do.call(rbind, lapply(cs, `[[`, "ot")); ot$Series <- lab(ot$series)
ra <- do.call(rbind, lapply(cs, `[[`, "ra")); ra$Series <- lab(ra$series)
lc <- do.call(rbind, lapply(cs, `[[`, "lc")); lc$Series <- lab(lc$series)
hr <- unique(ra[, c("Series", "h")]); hl <- unique(lc[, c("Series", "h")])
th <- theme_postcddp() + theme(legend.position = "top")
pA <- ggplot(ot, aes(x, cusum / 60, colour = Series)) + geom_hline(yintercept = 0, colour = "#5a5a5a", linewidth = 0.3) +
  geom_line(linewidth = 0.6) + scale_colour_manual(values = PAL, name = NULL) +
  labs(x = "Case number on the platform", y = "Operative-time CUSUM, hours", title = "A") + th
pB <- ggplot(ra, aes(x, racusum, colour = Series)) + geom_line(linewidth = 0.6) +
  geom_hline(data = hr, aes(yintercept = h, colour = Series), linetype = "22", linewidth = 0.4, show.legend = FALSE) +
  scale_colour_manual(values = PAL, name = NULL) +
  labs(x = "Case number on the platform", y = "RA-CUSUM score", title = "B") + th
pC <- ggplot(lc, aes(x, lccusum, colour = Series)) + geom_line(linewidth = 0.6) +
  geom_hline(data = hl, aes(yintercept = h, colour = Series), linetype = "22", linewidth = 0.4, show.legend = FALSE) +
  scale_colour_manual(values = PAL, name = NULL) +
  labs(x = "Case number on the platform", y = "LC-CUSUM score, pT2", title = "C") + th
fig <- (pA | pB | pC) + plot_layout(guides = "collect") & theme(legend.position = "top")
save_fig(fig, "figS2_cusum", width = 174, height = 70)
cat("figS2 done\n")

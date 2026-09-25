#!/usr/bin/env Rscript
# Figure 3. Adjusted differences, da Vinci minus Versius, for three prespecified comparisons:
# initial phase (da Vinci restart cases 1-100 vs Versius cases 1-100), concurrent calendar 2024, and
# all-era mature phases. A: operative time (minutes). B: binary outcomes (percentage points).
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(jsonlite); library(dplyr) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/fig_style.R"))
# every estimate comes from the authoritative contrast table (full model, Firth under the sparse rule, or crude)
ac <- read.csv(file.path(WS, "analysis/results/authoritative_contrasts.csv"), stringsAsFactors = FALSE)
CMP <- c(init = "First 100 cases on each robot", c2024 = "Calendar year 2024", mature = "After the Versius learning phase (post hoc)")
BLK <- c("A3a initial phase" = CMP[["init"]], "A4 concurrent 2024" = CMP[["c2024"]], "B1f post hoc, after the Versius learning phase" = CMP[["mature"]])
OUT <- c(or_time = "Operative time", psm_pt2 = "Positive margin, pT2", psm_all = "Positive margin, all stages",
         psa_persistence_eau = "PSA persistence", pad_free_3m = "Pad-free at 3 months", pad_free_12m = "Pad-free at 12 months")
df <- ac[ac$block %in% names(BLK) & ac$outcome %in% names(OUT) & !is.na(ac$est), ]
df <- data.frame(comparison = unname(BLK[df$block]), outcome = unname(OUT[df$outcome]), est = df$est, lo = df$lo95, hi = df$hi95, lo90 = df$lo90, hi90 = df$hi90,
                 n = sprintf("%d vs %d", df$n_ref, df$n_cmp), source = df$source)
df$comparison <- factor(df$comparison, levels = rev(CMP))
write.csv(df, file.path(WS, "figures/source_data/fig3_contrasts.csv"), row.names = FALSE)
SHP <- c(21, 22, 24); names(SHP) <- CMP
COL <- c("#1a1a1a", "#1E3C78", "#8c8c8c"); names(COL) <- CMP
forest <- function(dd, xlab, title, eq = NA, show_y = TRUE) {
  g <- ggplot(dd, aes(est, comparison))
  if (!is.na(eq)) g <- g + annotate("rect", xmin = -eq, xmax = eq, ymin = -Inf, ymax = Inf, fill = "#eeeeee")
  g <- g + geom_vline(xintercept = 0, colour = "#5a5a5a", linewidth = 0.4) +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = comparison), height = 0.25, linewidth = 0.5) +
    geom_segment(aes(x = lo90, xend = hi90, y = comparison, yend = comparison, colour = comparison), linewidth = 1.6, alpha = 0.55) +
    geom_point(aes(shape = comparison, colour = comparison), fill = "white", size = 2, stroke = 0.7) +
    scale_shape_manual(values = SHP, guide = "none") + scale_colour_manual(values = COL, guide = "none") +
    labs(x = xlab, y = NULL, title = title) + theme_postcddp()
  if (!show_y) g <- g + theme(axis.text.y = element_blank())
  g
}
pA <- forest(df[df$outcome == "Operative time", ], "Minutes", "A  Operative time", 15)
pB <- forest(df[df$outcome == "Positive margin, pT2", ], "Percentage points", "B  Margin, pT2", 10, FALSE)
pC <- forest(df[df$outcome == "Positive margin, all stages", ], "Percentage points", "C  Margin, all stages", NA, FALSE)
fig <- (pA | pB | pC) + plot_layout(widths = c(1, 1, 1))
save_fig(fig, "fig3_contrasts", width = 174, height = 58)
cat("fig3 done\n")

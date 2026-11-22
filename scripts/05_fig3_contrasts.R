#!/usr/bin/env Rscript
# Figure 2. Adjusted differences, da Vinci minus Versius, for the first 100 operations on each robot, calendar 2024
# and the post hoc comparison after the Versius learning phase. A: operative time (minutes), estimated separately
# within each hospital because theatre staff record operative time differently at each hospital (owner decision,
# 27 September 2026; estimates from 23_revision_checks.R). B and C: positive margins (percentage points), pooled
# models adjusted for hospital.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(jsonlite); library(dplyr) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/fig_style.R"))
# every estimate comes from the authoritative contrast table (full model, Firth under the sparse rule, or crude)
ac <- read.csv(file.path(WS, "analysis/results/authoritative_contrasts.csv"), stringsAsFactors = FALSE)
CMP <- c(init = "First 100 cases on each robot", c2024 = "Calendar year 2024", mature = "After the Versius learning phase (post hoc)")
CMP_LABELS <- setNames(c("First 100 Versius and first 100 restart operations", "Calendar year 2024 (primary comparison)",
                        "After the Versius learning phase (post hoc)"), CMP)
BLK <- c("A3a initial phase" = CMP[["init"]], "A4 concurrent 2024" = CMP[["c2024"]], "B1f post hoc, after the Versius learning phase" = CMP[["mature"]])
OUT <- c(or_time = "Operative time", psm_pt2 = "Positive margin, pT2", psm_all = "Positive margin, all stages",
         psa_persistence_eau = "PSA persistence", pad_free_3m = "Pad-free at 3 months", pad_free_12m = "Pad-free at 12 months")
df <- ac[ac$block %in% names(BLK) & ac$outcome %in% names(OUT) & !is.na(ac$est), ]
df <- data.frame(comparison = unname(BLK[df$block]), outcome = unname(OUT[df$outcome]), est = df$est, lo = df$lo95, hi = df$hi95, lo90 = df$lo90, hi90 = df$hi90,
                 n = sprintf("%d vs %d", df$n_ref, df$n_cmp), source = df$source)
df$comparison <- factor(df$comparison, levels = rev(CMP))
write.csv(df, file.path(WS, "figures/source_data/fig3_contrasts.csv"), row.names = FALSE)
# margins panels use one neutral marker so that colour and shape code hospitals only (panel A)
SHP <- c(23, 23, 23); names(SHP) <- CMP
COL <- c("#1a1a1a", "#1a1a1a", "#1a1a1a"); names(COL) <- CMP
# Keep labels below the intervals and the zero line.
direction_cues <- function(left, right) {
  list(
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0.72, fill = "white"),
    annotate("segment", x = 0, xend = 0, y = 0.72, yend = Inf, colour = "#5a5a5a", linewidth = 0.4),
    annotate("text", x = c(-Inf, Inf), y = 0.45, label = c(left, right), hjust = c(-0.1, 1.1),
             size = 2.2, family = "Open Sans", colour = "#5a5a5a"))
}
forest <- function(dd, xlab, title, eq = NA, show_y = TRUE) {
  g <- ggplot(dd, aes(est, comparison))
  if (!is.na(eq)) g <- g + annotate("rect", xmin = -eq, xmax = eq, ymin = -Inf, ymax = Inf, fill = "#eeeeee")
  g <- g + direction_cues("fewer on da Vinci", "more on da Vinci") +
    geom_errorbarh(aes(xmin = lo, xmax = hi, colour = comparison), height = 0.25, linewidth = 0.5) +
    geom_segment(aes(x = lo90, xend = hi90, y = comparison, yend = comparison, colour = comparison), linewidth = 1.6, alpha = 0.55) +
    geom_point(aes(shape = comparison, colour = comparison), fill = "white", size = 2, stroke = 0.7) +
    scale_shape_manual(values = SHP, guide = "none") + scale_colour_manual(values = COL, guide = "none") +
    scale_y_discrete(labels = CMP_LABELS, expand = expansion(add = c(0.75, 0.4))) +
    labs(x = xlab, y = NULL, title = title) + theme_postcddp()
  if (!show_y) g <- g + theme(axis.text.y = element_blank())
  g
}
rv <- fromJSON(file.path(WS, "analysis/results/R_revision_checks.json"), simplifyVector = FALSE)$within_hospital_or_time
KEY <- c(first100 = CMP[["init"]], y2024 = CMP[["c2024"]], after_learning_46 = CMP[["mature"]])
HOS <- c(SM = "SalveMedica", BE = "Bełchatów")
wh <- do.call(rbind, lapply(names(KEY), function(k) do.call(rbind, lapply(names(HOS), function(h) {
  x <- rv[[paste(k, h)]]
  data.frame(comparison = KEY[[k]], hospital = HOS[[h]], est = x$estimate, lo = x$ci95[[1]], hi = x$ci95[[2]], lo90 = x$ci90[[1]], hi90 = x$ci90[[2]],
             n = sprintf("%d vs %d", x$n_ref, x$n_cmp))
}))))
wh$y <- match(wh$comparison, rev(CMP)) + ifelse(wh$hospital == "SalveMedica", 0.17, -0.17)
write.csv(wh, file.path(WS, "figures/source_data/fig3A_within_hospital_operative_time.csv"), row.names = FALSE)
HCOL <- c("SalveMedica" = "#1E3C78", "Bełchatów" = "#8c8c8c"); HSHP <- c("SalveMedica" = 21, "Bełchatów" = 24)
HOS_LABELS <- setNames(c("SalveMedica (da Vinci already in use)", "Bełchatów (new to da Vinci)"), HOS)
pA <- ggplot(wh, aes(est, y)) + annotate("rect", xmin = -15, xmax = 15, ymin = -Inf, ymax = Inf, fill = "#eeeeee") +
  direction_cues("da Vinci shorter", "da Vinci longer") +
  geom_errorbarh(aes(xmin = lo, xmax = hi, colour = hospital), height = 0.14, linewidth = 0.5) +
  geom_segment(aes(x = lo90, xend = hi90, y = y, yend = y, colour = hospital), linewidth = 1.6, alpha = 0.55) +
  geom_point(aes(shape = hospital, colour = hospital), fill = "white", size = 2, stroke = 0.7) +
  scale_colour_manual(values = HCOL, labels = HOS_LABELS, name = NULL) +
  scale_shape_manual(values = HSHP, labels = HOS_LABELS, name = NULL) +
  scale_y_continuous(breaks = seq_along(CMP), labels = unname(CMP_LABELS[rev(CMP)]),
                     limits = c(0.25, length(CMP) + 0.4), expand = expansion()) +
  labs(x = "Difference in operative time, minutes (da Vinci minus Versius)", y = NULL,
       title = "A  Operative time, within each hospital") + theme_postcddp() +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0))
MARGIN_X <- "Percentage points\n(da Vinci minus Versius)"
pB <- forest(df[df$outcome == "Positive margin, pT2", ], MARGIN_X, "B  Margin, pT2", 10)
pC <- forest(df[df$outcome == "Positive margin, all stages", ], MARGIN_X, "C  Margin, all stages", show_y = FALSE)
fig <- pA / (pB | pC) + plot_layout(heights = c(1.25, 1))
save_fig(fig, "fig3_contrasts", width = 174, height = 118)
cat("fig3 done\n")

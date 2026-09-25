#!/usr/bin/env Rscript
# Davies test for a breakpoint in the adjusted log operative-time model, run at top level in each of the
# 20 imputations (the in-function call in 04_secondary.R failed to resolve the model data). Aggregate output.
suppressPackageStartupMessages({ library(jsonlite); library(segmented) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
dat <- load_ms6(WS, 20L); d <- dat$d; imps <- impute(d, 20, 11)
res <- list()
for (sr in c("Versius", "dVrestart")) {
  p <- c()
  for (i in seq_along(imps)) {
    dd <- imps[[i]]
    dd <- if (sr == "Versius") dd[dd$platform == "Versius", ] else dd[!is.na(dd$dv_restart_n), ]
    dd$x <- if (sr == "Versius") dd$platform_n else dd$dv_restart_n
    dd <- dd[!is.na(dd$log_or), ]
    m0 <- lm(as.formula(paste("log_or ~ x +", COV_OR)), data = dd)
    p <- c(p, davies.test(m0, seg.Z = ~x)$p.value)
  }
  res[[sr]] <- list(davies_p_median = median(p), davies_p_max = max(p), n = length(p))
}
write_json(res, file.path(WS, "analysis/results/S_seg_davies_adjusted.json"), auto_unbox = TRUE, digits = 8, pretty = TRUE)
print(res)

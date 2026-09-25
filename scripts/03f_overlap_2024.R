#!/usr/bin/env Rscript
# POST HOC. The 2024 comparison by hospital and month: platform availability, the share of operations in
# hospital-quarters where both robots were used, and a comparison restricted to those hospital-quarters.
# Also observed and eligible denominators for the functional outcomes, and complication-grade availability.
# Aggregate output only.
suppressPackageStartupMessages({ library(jsonlite) })
args_file <- sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))
WS <- normalizePath(file.path(dirname(args_file), ".."))
source(file.path(WS, "scripts/lc_common.R"))
dat <- load_ms6(WS, 20L); d <- dat$d
d$year <- as.integer(format(d$surgery_date, "%Y")); d$month <- as.integer(format(d$surgery_date, "%m"))
d$q <- (d$month - 1) %/% 3 + 1
d$robot <- ifelse(d$platform == "Versius", "Versius", ifelse(!is.na(d$dv_restart_n), "restart", NA))
y <- d[d$year == 2024 & !is.na(d$robot), ]
monthly <- as.data.frame(table(hospital = y$centre, month = y$month, robot = y$robot))
cells <- aggregate(robot ~ centre + q, data = y, FUN = function(r) length(unique(r)))
names(cells)[3] <- "n_robots"
y <- merge(y, cells, by = c("centre", "q"))
both <- y$n_robots == 2
out <- list(label = "POST HOC 2024 overlap and ascertainment",
  monthly_counts = monthly[monthly$Freq > 0, ],
  n_2024 = nrow(y), n_in_overlap_cells = sum(both),
  overlap_cells = unique(y[both, c("centre", "q")]),
  overlap_by_robot = as.list(table(y$robot[both])))
# comparison restricted to hospital-quarters with both robots
imps <- impute(d, 20, 11)
keys <- paste(out$overlap_cells$centre, out$overlap_cells$q)
g_ov <- function(dd) {
  yr <- as.integer(format(dd$surgery_date, "%Y")); qq <- (as.integer(format(dd$surgery_date, "%m")) - 1) %/% 3 + 1
  inov <- yr == 2024 & paste(dd$centre, qq) %in% keys
  ifelse(inov & dd$platform == "Versius", "ref", ifelse(inov & !is.na(dd$dv_restart_n), "cmp", NA))
}
COV_SMALL <- "age + bmi + log_psa + isup_cat + log_weight + pT_group"
out$overlap_or_time <- pool_contrast(imps, "overlap operative time", g_ov, "or_time", COV_SMALL, eq_margin = 15)
ov <- d[!is.na(g_ov(d)), ]; ov$grp <- g_ov(ov)
cnt <- function(v, g) { x <- ov[[v]][ov$grp == g]; x <- x[!is.na(x)]; c(events = sum(x), n = length(x)) }
out$overlap_counts <- list(or_median = tapply(ov$or_time, ov$grp, median, na.rm = TRUE),
  psm_pt2 = list(ref = cnt("psm_pt2", "ref"), cmp = cnt("psm_pt2", "cmp")),
  psm = list(ref = cnt("psm", "ref"), cmp = cnt("psm", "cmp")))
# functional outcome observation among eligible men, 2024 and whole series
obs <- function(dd, v, el) c(observed = sum(!is.na(dd[[v]]) & dd[[el]] == 1), eligible = sum(dd[[el]] == 1))
for (nm in c("Versius", "restart")) {
  y2 <- y[y$robot == nm, ]; s2 <- d[!is.na(d$robot) & d$robot == nm, ]
  out$functional_2024[[nm]] <- list(pad_free_3m = obs(y2, "pad_free_3m", "elig_3m"), pad_free_12m = obs(y2, "pad_free_12m", "elig_12m"),
                                    psa_persistence = obs(y2, "psa_persistence_eau", "elig_psa56"))
  out$functional_series[[nm]] <- list(pad_free_3m = obs(s2, "pad_free_3m", "elig_3m"), pad_free_12m = obs(s2, "pad_free_12m", "elig_12m"),
                                      psa_persistence = obs(s2, "psa_persistence_eau", "elig_psa56"))
  out$grade_available[[nm]] <- c(graded = sum(!is.na(s2$major_cd_any)), n = nrow(s2),
                                 reoperations = sum(s2$reoperation == 1, na.rm = TRUE),
                                 reop_graded_below3 = sum(s2$reoperation == 1 & s2$major_cd_any == 0, na.rm = TRUE),
                                 readmission_ascertained = sum(!is.na(s2$readmission_30d)))
}
j <- toJSON(out, auto_unbox = TRUE, digits = 6, pretty = TRUE)
if (grepl("\\b(19|20)[0-9]{2}-[0-9]{2}-[0-9]{2}\\b", j)) stop("date pattern in output")
writeLines(j, file.path(WS, "analysis/results/S_overlap_2024.json"))
cat("overlap 2024 done\n")

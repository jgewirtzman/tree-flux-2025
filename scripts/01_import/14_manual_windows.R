# ============================================================
# 14_manual_windows.R  --  INTERACTIVE; run in an R console, not with Rscript
#
# Review the closures on the trace-QC shortlist one at a time. For each you
# see CO2 (top) and CH4 (bottom) with 2 min of context and the current
# window, then choose:
#   k = keep the current window
#   c = click a new window (START then END on the CO2 panel)
#   x = exclude the closure (e.g. leak, breath contamination, no usable closure)
#   s = skip for now (decide later)      q = quit (progress is saved)
# Decisions are appended to data/input/manual_windows.csv (closure_id, decision,
# start, end, note). 09_goflux_reprocess.R applies them on the next run:
# clicked windows replace the automatic window (no deadband, no end trim) and
# excluded closures are kept in the data table but flagged manual_exclude = TRUE
# and dropped from the analysis dataset by 10_quality_flags.R.
#
# Usage (from the project root, in R or RStudio with a pop-up graphics device):
#   source("scripts/01_import/14_manual_windows.R")
#   review_windows()                         # the shortlist
#   review_windows(ids = c("LGR_20240626_431_155638"))   # any closures
# ============================================================
suppressPackageStartupMessages({ library(fluxqc) })

review_windows <- function(ids = NULL,
                           shortlist = "outputs/tables/trace_qc_clickpeak_shortlist.csv",
                           out = "data/input/manual_windows.csv") {
  tr <- readRDS("data/processed/goflux_traces.rds")
  traces <- tr$traces; clo <- tr$closures
  if (is.null(ids)) ids <- read.csv(shortlist, stringsAsFactors = FALSE)$closure_id
  done <- if (file.exists(out)) read.csv(out, stringsAsFactors = FALSE) else
    data.frame(closure_id = character(), decision = character(), start = character(), end = character(), note = character())
  todo <- setdiff(ids, done$closure_id[done$decision %in% c("keep", "click", "exclude")])
  message(length(todo), " closures to review (", length(ids) - length(todo), " already decided)")
  tz <- attr(traces$POSIX.time, "tzone")
  for (id in todo) {
    d <- traces[traces$UniqueID == id, ]
    if (!nrow(d)) { message(id, ": no trace"); next }
    d$start.time <- min(d$POSIX.time[d$flag == 1])       # blue line = current window start
    inw <- d$flag == 1
    op <- par(mfrow = c(2, 1), mar = c(2, 4, 2, 1))
    for (g in c("CO2dry_ppm", "CH4dry_ppb")) {
      plot(d$POSIX.time, d[[g]], pch = 16, cex = .4, col = ifelse(inw, "red", "grey50"), ylab = g,
           main = if (g == "CO2dry_ppm") paste(id, "- current window in red") else "")
    }
    par(op)
    a <- tolower(trimws(readline(sprintf("%s  [k]eep / [c]lick / e[x]clude / [s]kip / [q]uit: ", id))))
    if (a == "q") break
    if (a == "s" || !a %in% c("k", "c", "x")) next
    note <- if (a == "x") readline("  reason (optional): ") else ""
    row <- data.frame(closure_id = id, decision = c(k = "keep", c = "click", x = "exclude")[[a]],
                      start = NA_character_, end = NA_character_, note = note)
    if (a == "c") {
      m <- click_peak2_stacked(list(d), gases = c("CO2dry_ppm", "CH4dry_ppb"), sleep = 2)
      if (!nrow(m)) { message("  no window clicked; skipped"); next }
      row$start <- format(min(m$POSIX.time[m$flag == 1]), "%Y-%m-%d %H:%M:%S")
      row$end <- format(max(m$POSIX.time[m$flag == 1]), "%Y-%m-%d %H:%M:%S")
    }
    done <- rbind(done[done$closure_id != id, ], row)
    write.csv(done, out, row.names = FALSE)
    message("  saved: ", row$decision)
  }
  message("Decisions in ", out, ". Rerun the pipeline: bash scripts/run_pipeline.sh")
  invisible(done)
}

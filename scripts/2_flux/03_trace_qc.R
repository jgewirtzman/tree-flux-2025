# ============================================================
# 03_trace_qc.R  (flux step 3 of 6)
#
# Automated red-flag checks on every fitted closure, plus CO2 and CH4 trace
# plots for all closures, so manual window selection (click-peak) is needed
# only for a short list.
#
# Uses the traces saved by 02_fit_fluxes.R (fitting window flag = 1,
# with 120 s of context on either side) and its closure table.
#
# Checks (per closure):
#   window problems (candidates for manual re-windowing)
#     co2_not_rising    CO2 slope <= 0 or rise < 3 sigma over the window (goFlux co2.tracer)
#     co2_drop_in_window CO2 peaks and then falls back toward ambient before the window ends
#                       (chamber lifted / leak before the logged end time)
#     flat_start        CO2 flat in the first quarter of the window, rising afterwards
#                       (window starts before the chamber was sealed)
#     rise_mismatch     automatic rise detection (goFlux::find.rise on CO2 incl. context)
#                       covers < 70 % of the fitting window, or finds no rise at all
#     overlap           window overlaps another closure's window on the same analyzer
#   data problems (inspect; usually keep)
#     ch4_step          a single CH4 step > 10 x the day's sigma and > 20 % of the CH4 change
#                       (ebullition or analyzer glitch)
#     gap               missing data > 5 s inside the window
#     noisy / c0 / convex  from goFlux qc.flags (02_fit_fluxes.R)
#     lm_hm_disagree    LM and HM fluxes differ in sign, or g-factor > 2
#   other
#     short_window / long_window  < 90 s or > 600 s
#
# Outputs
#   data/interim/flux_trace_qc.csv             one row per closure, all checks (joined to the
#                                              final dataset by 06_flux_dataset.R)
#   outputs/tables/flux_processing/trace_qc_clickpeak_shortlist.csv  closures for inspection
#   outputs/figures/flux_processing/trace_qc/traces_all_<analyzer>_<year>.pdf   every closure, date order
#   outputs/figures/trace_qc/traces_flagged.pdf                 flagged closures, worst first
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(goFlux) })
# rough noise level of a trace for the screening thresholds below: MAD of first differences / sqrt(2)
precision_mad <- function(x) { d <- diff(x[!is.na(x)]); if (length(d) < 2) NA_real_ else stats::mad(d, constant = 1.4826) / sqrt(2) }
SCRIPT <- "03_trace_qc"
source("scripts/2_flux/flux_settings.R")

tr <- readRDS(PATH_TRACES)
traces <- tr$traces
g <- read.csv(PATH_FITS, stringsAsFactors = FALSE)
g <- g[g$fitted, ]
FIG <- file.path(FIG_DIR, "trace_qc"); dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

sp <- split(traces, traces$UniqueID)
win <- g %>% select(closure_id, instrument, window_start, window_end) %>%
  mutate(ws = as.POSIXct(window_start, tz = "America/New_York"), we = as.POSIXct(window_end, tz = "America/New_York"))

# overlaps between windows on the same analyzer
win <- win %>% arrange(instrument, ws) %>% group_by(instrument) %>%
  mutate(overlap = (!is.na(lag(we)) & ws < lag(we)) | (!is.na(lead(ws)) & we > lead(ws))) %>% ungroup()

check <- function(d, id) {
  w <- d[d$flag == 1, ]; w <- w[order(w$POSIX.time), ]
  n <- nrow(w); tt <- as.numeric(w$POSIX.time - w$POSIX.time[1], units = "secs")
  out <- list(closure_id = id, n_win = n)
  if (n < 10) return(as.data.frame(c(out, list(note = "too few points"))))
  s_co2 <- suppressWarnings(precision_mad(d$CO2dry_ppm)); s_ch4 <- suppressWarnings(precision_mad(d$CH4dry_ppb))
  co2 <- w$CO2dry_ppm; ch4 <- w$CH4dry_ppb
  sl <- coef(lm(co2 ~ tt))[2]; rise <- sl * diff(range(tt))
  out$co2_slope <- sl; out$co2_rise_ppm <- rise
  # CO2 falling through the whole window: the window sits on the flush-down after another closure
  out$co2_falling <- is.finite(rise) && rise < -max(10, 5 * s_co2)
  # drop back toward ambient inside the window
  k <- which.max(co2); tail10 <- mean(co2[tt >= max(tt) - 10])
  gain <- co2[k] - mean(co2[tt <= 10])
  out$co2_drop_in_window <- k < 0.85 * n & gain > 0 & (co2[k] - tail10) > max(5, 0.2 * gain)
  # flat start (window begins before sealing)
  q1 <- tt <= 0.25 * max(tt); mid <- tt > 0.25 * max(tt) & tt <= 0.75 * max(tt)
  s1 <- if (sum(q1) > 5) coef(lm(co2[q1] ~ tt[q1]))[2] else NA
  s2 <- if (sum(mid) > 5) coef(lm(co2[mid] ~ tt[mid]))[2] else NA
  out$flat_start <- !is.na(s1) && !is.na(s2) && s2 > 0 && s1 < 0.25 * s2 && (s2 * diff(range(tt[mid]))) > 5 * s_co2
  # start transient: a sharp jump (> 10 x sigma of differences) in the first 30 s of the window
  first <- which(tt[-1] <= 30)
  dco2 <- diff(co2); s_d <- function(v) mad(v[tt[-1] > 30], constant = 1.4826)
  out$start_transient <- length(first) > 0 &&
    (any(abs(diff(ch4)[first] - median(diff(ch4))) > 10 * s_d(diff(ch4))) ||
     any(abs(dco2[first] - median(dco2)) > 10 * s_d(dco2)))
  # CH4 declining from an elevated start: the chamber was not flushed to ambient (e.g. after the
  # previous tree) and CH4 relaxes toward ambient -> spurious "uptake"
  amb <- quantile(d$CH4dry_ppb, 0.05, na.rm = TRUE)
  ch4_sl <- coef(lm(ch4 ~ tt))[2]
  out$ch4_start_excess <- mean(ch4[tt <= 10]) - amb
  out$elevated_decline <- is.finite(ch4_sl) && ch4_sl < 0 && out$ch4_start_excess > max(20 * s_ch4, 20) &&
    (-ch4_sl * max(tt)) > 5 * s_ch4
  # CH4 step
  dch4 <- abs(diff(ch4)); rng <- diff(range(ch4))
  out$ch4_max_step <- max(dch4)
  out$ch4_step <- is.finite(s_ch4) && max(dch4) > 10 * s_ch4 && max(dch4) > 0.2 * rng
  # gaps
  out$max_gap_s <- max(diff(as.numeric(w$POSIX.time)))
  out$gap <- out$max_gap_s > 5
  # automatic rise detection on CO2 including the 120-s context
  r <- tryCatch(find.rise(d$POSIX.time, d$CO2dry_ppm, rise = 6, rise.secs = 30, drop = 8,
                          min.dur = 60, gap.secs = 12, min.n = 30), error = function(e) NULL)
  if (is.null(r)) { out$rise_found <- FALSE; out$rise_cover <- NA_real_; out$rise_start_offset_s <- NA_real_; out$rise_end_offset_s <- NA_real_ }
  else {
    ws <- min(w$POSIX.time); we <- max(w$POSIX.time)
    ov <- as.numeric(min(we, r$end) - max(ws, r$start), units = "secs")
    out$rise_found <- TRUE
    out$rise_cover <- max(0, ov) / as.numeric(we - ws, units = "secs")
    out$rise_start_offset_s <- as.numeric(ws - r$start, units = "secs")   # >0: window starts after rise onset
    out$rise_end_offset_s <- as.numeric(r$end - we, units = "secs")       # >0: rise continues after window
  }
  as.data.frame(out)
}
message("Checking ", length(sp), " closures ...")
chk <- bind_rows(lapply(names(sp), function(id) check(sp[[id]], id)))

f <- g %>% select(closure_id, date, Tree, SPECIES, location, instrument, in_legacy_dataset, t_sec,
                  CH4_flux_goflux, CH4_LM_flux, CH4_HM_flux, CH4_model, CH4_g_factor, CH4_MDF, CH4_below_MDF,
                  CO2_flux_goflux, CO2_LM_r2, qc_c0, qc_co2_tracer, qc_convex, qc_noisy, qc_noisy_ratio, fl_check, fl_notes) %>%
  left_join(chk, by = "closure_id") %>% left_join(win %>% select(closure_id, overlap), by = "closure_id") %>%
  mutate(
    month = as.integer(substr(date, 6, 7)), growing = month %in% 5:10,
    # a flat CO2 trace indicates a poor seal only in the growing season (guidelines paper);
    # dormant stems can have tiny CO2 efflux while CH4 accumulates cleanly
    co2_not_rising_any = coalesce(qc_co2_tracer, FALSE) | coalesce(co2_rise_ppm, 0) <= 0,
    co2_not_rising = co2_not_rising_any & growing,
    # the rise detector is only trusted when it clearly disagrees (covers < 30 % of the window)
    rise_mismatch = rise_found %in% TRUE & coalesce(rise_cover, 1) < 0.3,
    lm_hm_disagree = (sign(CH4_LM_flux) != sign(CH4_HM_flux) & abs(CH4_LM_flux) > coalesce(CH4_MDF, 0)) | coalesce(CH4_g_factor, 0) > 2,
    short_window = t_sec < 90, long_window = t_sec > 600,
    # the automatic rise detector proved unreliable on inspection; it is reported but not used
    legacy_disagree = in_legacy_dataset & !is.na(CH4_flux_goflux) &
      abs(CH4_flux_goflux - g$CH4_flux_legacy[match(closure_id, g$closure_id)]) >
        pmax(0.5, 0.5 * abs(g$CH4_flux_legacy[match(closure_id, g$closure_id)])),
    window_problem = co2_not_rising | coalesce(co2_falling, FALSE) | coalesce(co2_drop_in_window, FALSE) |
      coalesce(elevated_decline, FALSE) | coalesce(legacy_disagree, FALSE) |
      coalesce(flat_start, FALSE) | coalesce(overlap, FALSE),
    data_problem = coalesce(ch4_step, FALSE) | coalesce(gap, FALSE) | coalesce(qc_noisy, FALSE) |
      coalesce(start_transient, FALSE),
    reasons = trimws(paste(
      ifelse(co2_not_rising, "CO2 not rising (growing season);", ifelse(co2_not_rising_any, "CO2 flat (dormant);", "")), ifelse(coalesce(co2_drop_in_window, FALSE), "CO2 drops before window end;", ""),
      ifelse(coalesce(co2_falling, FALSE), "CO2 falling throughout;", ""),
      ifelse(coalesce(elevated_decline, FALSE), "CH4 declining from elevated start (not flushed);", ""),
      ifelse(coalesce(legacy_disagree, FALSE), "differs from legacy flux by >50%;", ""),
      ifelse(coalesce(flat_start, FALSE), "flat start;", ""), ifelse(coalesce(rise_mismatch, FALSE), "auto-rise mismatch;", ""),
      ifelse(coalesce(overlap, FALSE), "overlaps another window;", ""), ifelse(coalesce(ch4_step, FALSE), "CH4 step;", ""),
      ifelse(coalesce(gap, FALSE), "data gap;", ""), ifelse(coalesce(start_transient, FALSE), "start transient;", ""), ifelse(coalesce(qc_noisy, FALSE), "noisy;", ""),
      ifelse(coalesce(qc_c0, FALSE), "C0 high;", ""), ifelse(coalesce(lm_hm_disagree, FALSE), "LM/HM disagree;", ""),
      ifelse(short_window, "short window;", ""), ifelse(fl_check %in% "bad", "field log 'bad';", ""))),
    reasons = gsub("\\s+", " ", reasons),
    # effect on the result: does the flag touch a flux that matters? (above detection, in the analysis set)
    priority = 3 * window_problem + 2 * data_problem + coalesce(lm_hm_disagree, FALSE) +
      2 * (!coalesce(CH4_below_MDF, TRUE)) + in_legacy_dataset)
write.csv(f, PATH_TRACE_QC, row.names = FALSE)
log_step("8 trace QC", "fitted closures with a window problem flagged for inspection", sum(f$window_problem & f$in_legacy_dataset))
log_step("8 trace QC", "fitted closures with a data problem flagged (CH4 step, gap, noise, start transient)", sum(f$data_problem & f$in_legacy_dataset))

# click-peak shortlist: window problems where a better window plausibly exists
# (a CO2 rise is visible in the trace or context) and the flux is in the analysis set
short <- f %>% filter(window_problem, in_legacy_dataset,
                      co2_drop_in_window %in% TRUE | flat_start %in% TRUE | overlap %in% TRUE | co2_falling %in% TRUE |
                        elevated_decline %in% TRUE | legacy_disagree %in% TRUE |
                        (co2_not_rising & rise_found %in% TRUE)) %>%
  arrange(desc(priority), date)
# days where most closures have window problems point to a clock/timing error for the whole day:
# fix with one offset for the day (goFlux::find.clock.offset, checked by eye), not by clicking
days <- f %>% group_by(date, instrument) %>%
  summarise(n = n(), n_window_problem = sum(window_problem), n_co2_falling = sum(co2_falling %in% TRUE),
            frac = n_window_problem / n, .groups = "drop") %>%
  mutate(suspect_day = n >= 5 & frac >= 0.4) %>% arrange(desc(frac))
write.csv(days, file.path(TAB_DIR, "trace_qc_days.csv"), row.names = FALSE)
message("Suspect days (>= 40% of closures with window problems): ",
        paste(days$date[days$suspect_day], collapse = ", "))
short <- short %>% mutate(on_suspect_day = paste(date, instrument) %in% paste(days$date, days$instrument)[days$suspect_day])
write.csv(short, file.path(TAB_DIR, "trace_qc_clickpeak_shortlist.csv"), row.names = FALSE)
log_step("8 trace QC", "closures on the shortlist for inspection", nrow(short))
write_log()

summ <- f %>% summarise(closures = n(), in_analysis = sum(in_legacy_dataset),
  across(c(co2_not_rising, co2_falling, co2_drop_in_window, elevated_decline, legacy_disagree, flat_start, rise_mismatch, overlap, start_transient, ch4_step, gap, qc_noisy, qc_c0,
           lm_hm_disagree, short_window, window_problem, data_problem), ~ sum(.x %in% TRUE)))
print(t(summ))
message("Window problems among analysis-set closures: ", sum(f$window_problem & f$in_legacy_dataset),
        " (above MDF: ", sum(f$window_problem & f$in_legacy_dataset & !f$CH4_below_MDF %in% TRUE), ")")
message("Click-peak shortlist: ", nrow(short), " closures -> outputs/tables/trace_qc_clickpeak_shortlist.csv")

# ---------------- plots ----------------
plot_closure <- function(id, row) {
  d <- sp[[id]]; if (is.null(d)) return(invisible())
  t0 <- min(d$POSIX.time[d$flag == 1]); x <- as.numeric(d$POSIX.time - t0, units = "secs")
  inw <- d$flag == 1; xr <- range(x[inw])
  for (gas in c("CO2dry_ppm", "CH4dry_ppb")) {
    y <- d[[gas]]; ok <- is.finite(y)
    yl <- range(c(quantile(y[ok], c(0.01, 0.99)), y[ok & inw]))
    par(mar = c(if (gas == "CH4dry_ppb") 2 else 0.4, 3.2, if (gas == "CO2dry_ppm") 2.2 else 0.3, 0.4), mgp = c(1.9, 0.5, 0))
    plot(x[ok], y[ok], type = "n", xlab = "", ylim = yl, ylab = if (gas == "CO2dry_ppm") "CO2 (ppm)" else "CH4 (ppb)",
         xaxt = if (gas == "CH4dry_ppb") "s" else "n", cex.axis = 0.7, cex.lab = 0.75)
    rect(xr[1], par("usr")[3], xr[2], par("usr")[4], col = "#E8F1FA", border = NA)
    points(x[ok & !inw], y[ok & !inw], pch = 16, cex = 0.25, col = "grey60")
    points(x[ok & inw], y[ok & inw], pch = 16, cex = 0.3, col = if (gas == "CO2dry_ppm") "#B03A2E" else "#1B4F72")
    if (sum(ok & inw) > 3) abline(lm(y[ok & inw] ~ x[ok & inw]), col = "black", lwd = 0.8)
    if (gas == "CO2dry_ppm") {
      title(main = sprintf("%s | %s %s | CH4 %.3g (%s) | CO2 R2 %.2f", id, row$SPECIES, row$location,
                           row$CH4_flux_goflux, row$CH4_model, row$CO2_LM_r2), cex.main = 0.62, font.main = 1, line = 1.1)
      if (nzchar(row$reasons)) mtext(row$reasons, side = 3, line = 0.2, cex = 0.5, col = "#B03A2E")
    }
  }
}
pdf_pages <- function(rows, file) {
  pdf(file, width = 11, height = 8.5)
  layout(matrix(1:24, nrow = 8, ncol = 3, byrow = FALSE), heights = rep(c(1, 1), 4))
  # 4 closures per column x 3 columns = 12 closures per page, 2 panels each
  lay <- matrix(0, 8, 3); k <- 1
  for (cc in 1:3) for (rr in seq(1, 8, 2)) { lay[rr, cc] <- k; lay[rr + 1, cc] <- k + 1; k <- k + 2 }
  layout(lay)
  for (i in seq_len(nrow(rows))) plot_closure(rows$closure_id[i], rows[i, ])
  dev.off()
}
f <- f %>% mutate(year = substr(date, 1, 4), analyzer = ifelse(instrument == "LI-7810", "LI7810", "LGR"))
for (k in split(f, paste(f$analyzer, f$year))) {
  k <- k[order(k$date, k$closure_id), ]
  pdf_pages(k, file.path(FIG, sprintf("traces_all_%s_%s.pdf", k$analyzer[1], k$year[1])))
}
fl <- f %>% filter(window_problem | data_problem) %>% arrange(desc(priority), date)
pdf_pages(fl, file.path(FIG, "traces_flagged.pdf"))
pdf_pages(short, file.path(FIG, "traces_clickpeak_shortlist.pdf"))
message("Trace plots: ", FIG)

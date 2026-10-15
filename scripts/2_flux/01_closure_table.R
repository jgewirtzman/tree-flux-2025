# ============================================================
# 01_closure_table.R  (flux step 1 of 5)
#
# One row per measurement from data/package/measurements.csv (tree, analyzer, chamber
# closure start and end on the analyzer clock, deadband, sampling time, air temperature and
# pressure), with the chamber geometry and the trace segment cut from the raw 1-Hz record.
#
# Window: closure start + deadband (20 s; 0 where the start was set on inspection of the
# record) to the closure end, ended earlier where an abrupt, sustained fall in CO2 shows
# the chamber was lifted. Analyzer precision per analyzer x field day from the windows
# (goFlux::empirical.prec, method "hadamard").
#
# Inputs : data/package/ (measurements.csv, analyzer_raw/, chamber_volumes.csv, trees.csv)
# Output : data/interim/flux_closures.rds
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(goFlux)
})
stopifnot(packageVersion("goFlux") >= "0.5.0.9002")   # jgewirtzman/goFlux fork (doi:10.5281/zenodo.23256675)

SCRIPT <- "01_closure_table"
source("scripts/2_flux/flux_settings.R")

# ============================================================
# PART 1: MEASUREMENTS
# ============================================================
meas <- read.csv(PATH_MEASUREMENTS, stringsAsFactors = FALSE)
meas$date <- as.Date(meas$date)
meas$closure_start <- as.POSIXct(meas$closure_start, format = "%Y-%m-%d %H:%M:%S", tz = TZ)
meas$closure_end   <- as.POSIXct(meas$closure_end,   format = "%Y-%m-%d %H:%M:%S", tz = TZ)
stopifnot(!anyNA(meas$closure_start), !anyNA(meas$closure_end), !anyDuplicated(meas$measurement_id))
message("Measurements: ", nrow(meas), " (", paste(names(table(meas$analyzer)), table(meas$analyzer), sep = " ", collapse = ", "), ")")
log_step("1 measurements", "measurements (data/package/measurements.csv)", nrow(meas))

# ============================================================
# PART 2: RAW RECORDS
# ============================================================
# LGR/UGGA: analyzer_raw/lgr/LGRn/YYYY-MM-DD/micro_YYYY-MM-DD_fNNNN.txt(.zip)
lgr_files <- list.files(PATH_LGR, pattern = "^micro_\\d{4}-\\d{2}-\\d{2}_f\\d+\\.txt(\\.zip)?$",
                        recursive = TRUE, full.names = TRUE)
lgr_catalog <- data.frame(
  file = lgr_files,
  machine = sub("^([^/]+)/.*$", "\\1", sub(paste0("^", PATH_LGR, "/"), "", lgr_files)),
  date = as.Date(sub("^micro_(\\d{4}-\\d{2}-\\d{2})_.*$", "\\1", basename(lgr_files))),
  stringsAsFactors = FALSE)
message("LGR data files: ", nrow(lgr_catalog), " on ", length(unique(lgr_catalog$date)), " dates")

parse_lgr_file <- function(fp) {
  lines <- tryCatch({
    if (grepl("\\.zip$", fp)) {
      inner <- utils::unzip(fp, list = TRUE)$Name
      inner <- inner[grepl("\\.txt$", inner) & !grepl("^__MACOSX|/\\._", inner)][1]
      con <- unz(fp, inner); on.exit(close(con)); readLines(con, warn = FALSE)
    } else readLines(fp, warn = FALSE)
  }, error = function(e) character(0))
  if (length(lines) < 3) return(NULL)
  dat <- tryCatch(read.csv(text = paste(lines[-1], collapse = "\n"), stringsAsFactors = FALSE,
                           strip.white = TRUE, check.names = FALSE), error = function(e) NULL)
  if (is.null(dat) || !nrow(dat)) return(NULL)
  names(dat) <- trimws(names(dat))
  ch4 <- grep("^\\[CH4\\]d_ppm$", names(dat)); co2 <- grep("^\\[CO2\\]d_ppm$", names(dat))
  h2o <- grep("^\\[H2O\\]_ppm$", names(dat))
  if (!length(ch4) || !length(co2)) return(NULL)
  ts <- as.POSIXct(trimws(dat[[1]]), format = "%m/%d/%Y %H:%M:%OS", tz = TZ)
  out <- data.frame(POSIX.time = ts, CH4dry_ppb = as.numeric(dat[[ch4[1]]]) * 1000,
                    CO2dry_ppm = as.numeric(dat[[co2[1]]]),
                    H2O_ppm = if (length(h2o)) as.numeric(dat[[h2o[1]]]) else NA_real_)
  out[!is.na(out$POSIX.time), ]
}
load_lgr_day <- function(date, machine) {
  f <- lgr_catalog$file[lgr_catalog$date == date & lgr_catalog$machine == machine]
  if (!length(f)) return(NULL)
  d <- do.call(rbind, lapply(f, parse_lgr_file))
  if (is.null(d) || !nrow(d)) return(NULL)
  d <- d[order(d$POSIX.time), ]
  d[!duplicated(d$POSIX.time), ]
}

# LI-7810: one .data file per export
li_files <- list.files(PATH_LI7810, pattern = "\\.data$", full.names = TRUE)
li_raw <- do.call(rbind, lapply(li_files, function(f) {
  x <- import.LI7810(f, timezone = TZ, keep_all = TRUE)
  as.data.frame(x)[, c("POSIX.time", "CH4dry_ppb", "CO2dry_ppm", "H2O_ppm")]
}))
li_raw <- li_raw[order(li_raw$POSIX.time), ]
li_raw <- li_raw[!duplicated(li_raw$POSIX.time), ]
message("LI-7810 data files: ", length(li_files))

# ============================================================
# PART 3: WINDOWS AND TRACE SEGMENTS
# ============================================================
# End the window at chamber removal: the first abrupt fall of the smoothed CO2 (7-point running
# median) by more than max(10 ppm, 8 sigma, 25 % of the rise) within 10 s, after >= 60 s, that
# never recovers.
window_trims <- list()
trim_at_opening <- function(seg) {
  w <- which(seg$flag == 1 & is.finite(seg$CO2dry_ppm)); if (length(w) < 70) return(seg)
  co2 <- seg$CO2dry_ppm[w]; tt <- as.numeric(seg$POSIX.time[w] - seg$POSIX.time[w[1]], units = "secs")
  sm <- runmed(co2, 7, endrule = "median"); sd_d <- mad(diff(co2), constant = 1.4826) / sqrt(2)
  base <- median(sm[tt <= 10]); thr <- max(10, 8 * sd_d, 0.25 * max(max(sm) - base, 0))
  for (i in which(tt >= 60)) {
    j <- which(tt > tt[i] & tt <= tt[i] + 10); if (!length(j)) next
    if (sm[i] - min(sm[j]) > thr) {
      rest <- which(tt > tt[i] + 10)
      if (!length(rest) || max(sm[rest]) < sm[i] - thr / 2) {
        new_end <- seg$POSIX.time[w[max(which(tt <= tt[i] - 2))]]
        window_trims[[seg$UniqueID[1]]] <<- as.numeric(max(seg$POSIX.time[w]) - new_end, units = "secs")
        seg$flag <- as.numeric(seg$flag == 1 & seg$POSIX.time <= new_end)
        seg$end.time <- new_end; seg$obs.length <- as.numeric(new_end - seg$start.time[1], units = "secs")
        return(seg)
      }
    }
  }
  seg
}
segment_trace <- function(day, start, end, uid) {
  seg <- day[day$POSIX.time >= start - SHOULDER_S & day$POSIX.time <= end + SHOULDER_S, ]
  if (!nrow(seg)) return(NULL)
  seg$UniqueID <- uid
  seg$flag <- as.numeric(seg$POSIX.time >= start & seg$POSIX.time <= end)
  if (sum(seg$flag) < 5) return(NULL)
  seg$start.time <- start; seg$end.time <- end
  seg$Etime <- as.numeric(seg$POSIX.time - start, units = "secs")
  seg$obs.length <- as.numeric(end - start, units = "secs")
  trim_at_opening(seg)
}

meas$window_start <- meas$closure_start + meas$deadband_s
segments <- list()
lgr_days <- unique(meas[meas$analyzer != "LI-7810", c("date", "analyzer")])
for (k in seq_len(nrow(lgr_days))) {
  day <- load_lgr_day(lgr_days$date[k], lgr_days$analyzer[k])
  if (is.null(day)) next
  for (i in which(meas$date == lgr_days$date[k] & meas$analyzer == lgr_days$analyzer[k])) {
    seg <- segment_trace(day, meas$window_start[i], meas$closure_end[i], meas$measurement_id[i])
    if (!is.null(seg)) { seg$instrument <- meas$analyzer[i]; segments[[meas$measurement_id[i]]] <- seg }
  }
}
for (i in which(meas$analyzer == "LI-7810")) {
  day <- li_raw[li_raw$POSIX.time >= meas$closure_start[i] - SHOULDER_S & li_raw$POSIX.time <= meas$closure_end[i] + SHOULDER_S, ]
  seg <- segment_trace(day, meas$window_start[i], meas$closure_end[i], meas$measurement_id[i])
  if (!is.null(seg)) { seg$instrument <- "LI-7810"; segments[[meas$measurement_id[i]]] <- seg }
}
message("Measurements with a trace in the window: ", length(segments), " / ", nrow(meas))
log_step("2 windows", "measurements with a raw trace in the window", length(segments))

meas$end_trim_s <- unlist(window_trims)[meas$measurement_id]
meas$end_trim_s[is.na(meas$end_trim_s)] <- 0
meas$window_end <- meas$closure_end - meas$end_trim_s
log_step("2 windows", "windows ended early at a detected chamber opening (abrupt, sustained CO2 fall)", sum(meas$end_trim_s > 0),
         sprintf("median %.0f s removed", median(meas$end_trim_s[meas$end_trim_s > 0])))

# ============================================================
# PART 4: CHAMBER GEOMETRY AND TREE METADATA
# ============================================================
tv <- read.csv(PATH_VOLUMES, stringsAsFactors = FALSE, check.names = FALSE)
tv$tree <- as.integer(round(tv$Tag))
trees <- read.csv(PATH_TREES, stringsAsFactors = FALSE)
closures <- meas %>%
  transmute(closure_id = measurement_id, date, tree = Tree, PLOT, instrument = analyzer, closure_start, window_start,
            window_end, closure_end, deadband_s, end_trim_s, sample_time_local, Tcham = air_temp_C, Pcham = pressure_kPa) %>%
  left_join(tv %>% transmute(tree, Vcollar = Volume), by = "tree") %>%
  left_join(trees %>% transmute(tree = Tree, SPECIES = species, DBH = DBH_cm, location = site), by = "tree") %>%
  mutate(Vextra = ifelse(instrument == "LI-7810", EXTRA_VOL_7810_L, EXTRA_VOL_LGR_L), Vtot = Vcollar + Vextra,
         Area = AREA_CM2, has_trace = closure_id %in% names(segments))
stopifnot(!anyNA(closures$Vtot), !anyNA(closures$SPECIES))

# ============================================================
# PART 5: ANALYZER PRECISION PER ANALYZER x FIELD DAY
# ============================================================
# goFlux::empirical.prec(method = "hadamard"): for each closure, the MAD-based SD of the second
# differences of the concentration inside the fitting window, / sqrt(6); the median over the closures
# of that analyzer and day (at the day's dominant logging interval). Second differences remove the
# flux trend, so the estimate is not inflated by the rise itself.
pd <- do.call(rbind, lapply(segments, function(s) data.frame(
  UniqueID = s$UniqueID, flag = s$flag, POSIX.time = s$POSIX.time, CH4dry_ppb = s$CH4dry_ppb,
  CO2dry_ppm = s$CO2dry_ppm, instrument = s$instrument,
  date = as.Date(min(s$POSIX.time[s$flag == 1]), tz = TZ))))
pd$inst_day <- paste(pd$instrument, pd$date, sep = "|")
prec_day <- function(gas) {
  r <- suppressWarnings(empirical.prec(pd, gas, method = "hadamard", by = "inst_day", warn = FALSE))
  r %>% group_by(inst_day) %>% mutate(n_interval_runs = n()) %>%
    slice_max(n.closures, n = 1, with_ties = FALSE) %>% ungroup()
}
p_ch4 <- prec_day("CH4dry_ppb"); p_co2 <- prec_day("CO2dry_ppm")
day_sigma <- p_ch4 %>% transmute(inst_day, n_closures = n.closures, CH4_sigma_day = prec, dt_day = dt_s, n_interval_runs) %>%
  left_join(p_co2 %>% transmute(inst_day, CO2_sigma_day = prec), by = "inst_day") %>%
  mutate(instrument = sub("\\|.*", "", inst_day), date = as.Date(sub(".*\\|", "", inst_day))) %>%
  select(instrument, date, n_closures, CH4_sigma_day, CO2_sigma_day, dt_day, n_interval_runs) %>% as.data.frame()
message("Precision (goFlux empirical.prec, hadamard), median CH4 ppb by analyzer: ",
        paste(sprintf("%s %.3f", names(tapply(day_sigma$CH4_sigma_day, day_sigma$instrument, median)),
                      tapply(day_sigma$CH4_sigma_day, day_sigma$instrument, median)), collapse = "; "))

saveRDS(list(closures = closures, segments = segments, day_sigma = day_sigma), PATH_CLOSURES)
write_log()
message("Wrote ", PATH_CLOSURES, ": ", nrow(closures), " measurements (with trace: ", sum(closures$has_trace), ")")

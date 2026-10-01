# ============================================================
# 01_closure_table.R  (flux step 1 of 6)
#
# Builds one row per chamber closure from the field logs and the raw analyzer
# records, with the fitting window, chamber geometry and met, and cuts the trace
# segment for each closure. Every cleaning rule logs how many records it touched
# (data/interim/flux_log_01_closure_table.csv).
#
#   1 earlier dataset   curated tree IDs, met and the fluxes of measurements with no raw record
#   2 field logs        read, apply the team's timing corrections, map old tags, drop unusable
#                       or duplicate entries (LGR/UGGA, Jun 2023 - Mar 2025)
#   3 matching          field-log closures <-> measurements of the earlier dataset
#   4 LGR traces        cut each window from the 1-Hz day files; end it at chamber opening
#   5 LI-7810 remarks   closures from the tagged remarks (2025); aborted starts and superseded
#                       redos removed
#   6 windows           manual window decisions, curated tags, geometry and met
#
# Windows: field-log start + 20 s deadband to the logged end (LGR/UGGA), or remark start
# + 20 s to remark end (LI-7810). Nothing is clicked by hand; the QC screens of
# 03_trace_qc.R flag closures for inspection, and the decisions are recorded in
# data/package/qc_decisions/manual_windows.csv (04_review_windows.R).
#
# Inputs : data/package/ (field_logs, analyzer_raw, chamber_volumes.csv,
#          previous_processing, qc_decisions), data/interim/wtd_met.csv
# Output : data/interim/flux_closures.rds
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(goFlux)
  library(fluxqc)
})
stopifnot(packageVersion("fluxqc") >= "0.2.3")

SCRIPT <- "01_closure_table"
source("scripts/2_flux/flux_settings.R")
INCLUDE_NEW_CLOSURES <- FALSE   # closures with a trace but no row in the earlier dataset (Jan-Mar 2025)

# ============================================================
# PART 1: LEGACY DATASET
# ============================================================

message("=== Part 1: legacy dataset ===")

legacy <- read.csv(PATH_V1,
                   stringsAsFactors = FALSE)
legacy$legacy_row <- seq_len(nrow(legacy))
# datetime_posx is local wall time written with a spurious "Z"
legacy$legacy_time <- as.POSIXct(sub("Z$", "", legacy$datetime_posx),
                                 format = "%Y-%m-%dT%H:%M:%S", tz = TZ)
legacy$date <- as.Date(legacy$date)
legacy$UniqueID_legacy <- legacy$UniqueID
message("Legacy rows: ", nrow(legacy), " (",
        paste(names(table(legacy$year)), table(legacy$year), sep = ": ", collapse = ", "), ")")
log_step("1 earlier dataset", "measurements in the earlier processed dataset (curated tree IDs, met, fluxes of measurements without a raw record)", nrow(legacy))

# ============================================================
# PART 2: FIELD LOGS -> CLOSURE WINDOWS (LGR/UGGA, Jun 2023 - Mar 2025)
# ============================================================

message("\n=== Part 2: field logs ===")

fl_dir   <- PATH_FIELD_LOGS
xlsx_dir <- PATH_FIELD_LOGS
pick <- function(a, b) ifelse(nz(a), as.character(a), as.character(b))
hms_of <- function(x) {            # readxl gives 1899-12-31 hh:mm:ss POSIXct; csv gives "hh:mm:ss";
                                   # two xlsx files hold Excel fractions of a day as text ("0.51310")
  if (inherits(x, "POSIXt")) return(format(x, "%H:%M:%S"))
  x <- trimws(as.character(x)); x <- sub("^.*\\s+", "", x)
  frac <- nz(x) & grepl("^0?\\.\\d+$", x)
  if (any(frac)) {
    secs <- round(as.numeric(x[frac]) * 86400)
    x[frac] <- sprintf("%02d:%02d:%02d", secs %/% 3600, (secs %% 3600) %/% 60, secs %% 60)
  }
  x[nz(x) & grepl("^\\d{1,2}:\\d{2}$", x)] <- paste0(x[nz(x) & grepl("^\\d{1,2}:\\d{2}$", x)], ":00")
  x
}

fl1 <- read.csv(file.path(fl_dir, "Field_Data_Monthly_Summer2023.csv"), stringsAsFactors = FALSE,
                check.names = FALSE)
fl1_std <- data.frame(
  UniqueID_fl = fl1$UniqueID, date_raw = fl1$Date, tree_raw = as.character(fl1$Tree_Tag),
  comp_start = hms_of(fl1$comp_start_time), comp_end = hms_of(fl1$comp_end_time),
  real_start = hms_of(fl1$`Real start`), machine = trimws(fl1$Analyzer),
  fl_check = NA_character_, fl_notes = fl1$Notes, source = "summer2023",
  stringsAsFactors = FALSE)

fl2 <- read.csv(file.path(fl_dir, "Field_Data_Monthly_updated2.csv"), stringsAsFactors = FALSE,
                check.names = FALSE)
fl2_std <- data.frame(
  UniqueID_fl = fl2$UniqueID, date_raw = fl2$format_Date, tree_raw = as.character(fl2$`Tree Tag`),
  comp_start = hms_of(pick(fl2$`Updated Start Time`, fl2$comp_start_time)),
  comp_end   = hms_of(pick(fl2$`Updated End Time`,   fl2$comp_end_time)),
  real_start = hms_of(fl2$`Real start`), machine = trimws(fl2$Machine),
  fl_check = fl2$check, fl_notes = paste(fl2$check_notes, fl2$Notes), source = "updated2",
  stringsAsFactors = FALSE)

fl3 <- read.csv(file.path(fl_dir, "summer_2024.csv"), stringsAsFactors = FALSE, check.names = FALSE)
fl3_std <- data.frame(
  UniqueID_fl = fl3$UniqueID, date_raw = fl3$format_Date, tree_raw = as.character(fl3$Tree),
  comp_start = hms_of(fl3$comp_start_time), comp_end = hms_of(fl3$comp_end_time),
  real_start = hms_of(fl3$`Real start`), machine = "LGR1",
  fl_check = NA_character_, fl_notes = NA_character_, source = "summer2024",
  stringsAsFactors = FALSE)

# Team's manual window review for 22 May - 17 Jul 2024 (Timing Updates.xlsx): hand-adjusted
# start/end times and good/bad calls. These are the equivalent of click-peak windows and override
# the summer_2024.csv times; "bad" closures were already left out of the published dataset.
tu_path <- file.path(PATH_FIELD_LOGS, "Timing Updates.xlsx")
if (file.exists(tu_path)) {
  tu <- as.data.frame(readxl::read_excel(tu_path, sheet = 1))
  tu <- tu[!is.na(tu$`Unique ID`), ]
  tu_start <- hms_of(tu$Updated_start); tu_end <- hms_of(tu$Updated_end)
  k <- match(fl3_std$UniqueID_fl, tu$`Unique ID`)
  has <- !is.na(k)
  fl3_std$start_adjusted <- FALSE
  fl3_std$start_adjusted[has] <- nz(tu_start[k[has]]) & tu_start[k[has]] != fl3_std$comp_start[has]
  fl3_std$comp_start[has] <- ifelse(nz(tu_start[k[has]]), tu_start[k[has]], fl3_std$comp_start[has])
  fl3_std$comp_end[has]   <- ifelse(nz(tu_end[k[has]]),   tu_end[k[has]],   fl3_std$comp_end[has])
  fl3_std$fl_check[has]   <- tu$quality[k[has]]
  fl3_std$fl_notes[has]   <- tu$Notes[k[has]]
  message("Timing Updates.xlsx applied to ", sum(has), " summer-2024 closures (bad: ",
          sum(tu$quality[k[has]] %in% "bad"), ")")
  log_step("2 field logs", "summer-2024 windows replaced by the team's hand-adjusted start/end (Timing Updates.xlsx)", sum(has),
           sprintf("start changed for %d; marked bad %d", sum(fl3_std$start_adjusted), sum(tu$quality[k[has]] %in% "bad")))
}

# Sep 2024 - Mar 2025 xlsx logs (the LGR was used until the LI-7810 arrived in Apr 2025;
# the Jan-Mar 2025 files must NOT be skipped)
xlsx_files <- list.files(xlsx_dir, pattern = "^Stem_flux.*[.]xlsx$", full.names = TRUE)
fl4_std <- do.call(rbind, lapply(xlsx_files, function(xf) {
  xd <- as.data.frame(readxl::read_excel(xf))
  data.frame(
    UniqueID_fl = as.character(xd$UniqueID), date_raw = as.character(xd$format_Date),
    tree_raw = as.character(xd$`Tree Tag`),
    comp_start = ifelse(nz(xd$`Updated Start Time`), hms_of(xd$`Updated Start Time`), hms_of(xd$comp_start_time)),
    comp_end   = ifelse(nz(xd$`Updated End Time`),   hms_of(xd$`Updated End Time`),   hms_of(xd$comp_end_time)),
    real_start = hms_of(xd$`Real start`),
    machine = if ("Machine" %in% names(xd)) as.character(xd$Machine) else NA_character_,
    fl_check = if ("check" %in% names(xd)) as.character(xd$check) else NA_character_,
    fl_notes = paste(if ("check_notes" %in% names(xd)) xd$check_notes else "",
                     if ("Notes" %in% names(xd)) xd$Notes else ""),
    source = paste0("xlsx:", basename(xf)), stringsAsFactors = FALSE)
}))
message("xlsx field logs: ", length(xlsx_files), " files, ", nrow(fl4_std), " entries (",
        paste(range(fl4_std$date_raw), collapse = " to "), ")")

field_logs <- bind_rows(fl1_std, fl2_std, fl3_std, fl4_std)
log_step("2 field logs", "entries read from the field logs", nrow(field_logs),
         paste(sprintf("%s %d", c("summer 2023", "2023-24 updated", "summer 2024", "Sep 2024-Mar 2025 xlsx"),
                       c(nrow(fl1_std), nrow(fl2_std), nrow(fl3_std), nrow(fl4_std))), collapse = "; "))
# hand-adjusted starts (Timing Updates.xlsx) already skip the transient: no extra deadband
field_logs$start_adjusted <- field_logs$start_adjusted %in% TRUE
field_logs$deadband_s <- ifelse(field_logs$start_adjusted, 0, DEADBAND_LGR)
field_logs$machine <- gsub("LGR #", "LGR", field_logs$machine)
field_logs$machine <- gsub("\\s+", "", field_logs$machine)
log_step("2 field logs", "analyzer not recorded; set to LGR1 (the logs' default)", sum(!nz(field_logs$machine) | field_logs$machine == "NA"))
field_logs$machine[!nz(field_logs$machine) | field_logs$machine == "NA"] <- "LGR1"

parse_date_multi <- function(x) {
  x <- trimws(as.character(x)); out <- as.Date(rep(NA, length(x)))
  for (f in c("%Y-%m-%d", "%m/%d/%Y", "%m/%d/%y", "%m-%d-%Y", "%m-%d-%y", "%y-%m-%d")) {
    i <- is.na(out) & nz(x)
    if (!any(i)) break
    d <- suppressWarnings(as.Date(x[i], format = f))
    d[!is.na(d) & (as.integer(format(d, "%Y")) < 2000 | as.integer(format(d, "%Y")) > 2030)] <- NA
    out[i] <- d
  }
  out
}
field_logs$date <- parse_date_multi(field_logs$date_raw)
field_logs$tree <- suppressWarnings(as.integer(gsub("[^0-9]", "", field_logs$tree_raw)))
# The 2023 and Jan-2024 logs used the wetland trees' old tags (1-32); map them to the ForestGEO tags
tv0 <- read.csv(PATH_VOLUMES, stringsAsFactors = FALSE, check.names = FALSE)
old_map <- setNames(as.integer(round(tv0$Tag)), as.integer(tv0$`Old tag`))
old_map <- old_map[!is.na(names(old_map))]
is_old <- !is.na(field_logs$tree) & field_logs$tree %in% as.integer(names(old_map)) &
  !field_logs$tree %in% as.integer(round(tv0$Tag))
field_logs$tree_logged <- field_logs$tree
field_logs$tree[is_old] <- unname(old_map[as.character(field_logs$tree[is_old])])
message("Field-log old wetland tags mapped to ForestGEO tags: ", sum(is_old))
log_step("2 field logs", "old wetland tags (1-32) mapped to ForestGEO tags", sum(is_old))
log_step("2 field logs", "entries without a readable date dropped", sum(is.na(field_logs$date)))
field_logs <- field_logs %>%
  filter(!is.na(date)) %>%
  mutate(
    comp_start_posix = as.POSIXct(paste(date, comp_start), format = "%Y-%m-%d %H:%M:%S", tz = TZ),
    comp_end_posix   = as.POSIXct(paste(date, comp_end),   format = "%Y-%m-%d %H:%M:%S", tz = TZ),
    real_start_posix = as.POSIXct(paste(date, real_start), format = "%Y-%m-%d %H:%M:%S", tz = TZ),
    t_sec_fl = as.numeric(difftime(comp_end_posix, comp_start_posix, units = "secs"))
  )
bad_win <- is.na(field_logs$comp_start_posix) | is.na(field_logs$comp_end_posix) |
  field_logs$t_sec_fl <= MIN_WINDOW_S | field_logs$t_sec_fl >= 1800
message("Field-log entries: ", nrow(field_logs), "; unusable windows (missing/short/long): ", sum(bad_win))
log_step("2 field logs", sprintf("entries without a usable window dropped (start or end missing, <= %d s or >= 30 min)", MIN_WINDOW_S), sum(bad_win))
field_logs_bad <- field_logs[bad_win, ]
field_logs <- field_logs[!bad_win, ]

# De-duplicate entries logged in more than one file (prefer the later, QC'd log)
n_before_dedup <- nrow(field_logs)
field_logs <- field_logs %>%
  mutate(priority = case_when(grepl("^xlsx", source) ~ 4, source == "summer2024" ~ 3,
                              source == "updated2" ~ 2, TRUE ~ 1),
         dup_key = ifelse(nz(UniqueID_fl), UniqueID_fl, paste(date, tree, comp_start))) %>%
  group_by(dup_key) %>% slice_max(priority, n = 1, with_ties = FALSE) %>% ungroup() %>%
  # also collapse the same closure logged under different IDs
  group_by(date, tree, comp_start) %>% slice_max(priority, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(-priority, -dup_key) %>% arrange(comp_start_posix)
log_step("2 field logs", "duplicate entries (same closure in more than one log) collapsed, keeping the later QC'd log", n_before_dedup - nrow(field_logs))
log_step("2 field logs", "unique field-log closures (LGR/UGGA)", nrow(field_logs))
message("Unique field-log closures: ", nrow(field_logs), " (", min(field_logs$date), " to ",
        max(field_logs$date), "); machines: ",
        paste(names(table(field_logs$machine)), table(field_logs$machine), sep = "=", collapse = ", "))

# ============================================================
# PART 3: MATCH LEGACY 2023-24 ROWS TO FIELD-LOG WINDOWS
# ============================================================

message("\n=== Part 3: legacy <-> field log matching ===")

leg_lgr <- legacy %>% filter(year < 2025)
field_logs$legacy_row <- NA_integer_

# 3a. exact UniqueID (2024 style "45296_895"); the 2023 IDs were mangled by Excel ("1-Jan")
m1 <- match(leg_lgr$UniqueID_legacy, field_logs$UniqueID_fl)
ok1 <- !is.na(m1) & is.na(field_logs$legacy_row[m1]) & !duplicated(m1)
field_logs$legacy_row[m1[ok1]] <- leg_lgr$legacy_row[ok1]
# 3b. date + tree + nearest Real-start time (legacy datetime_posx = field-log "Real start")
for (i in which(!ok1)) {
  cand <- which(field_logs$date == leg_lgr$date[i] & field_logs$tree == leg_lgr$Tree[i] &
                  is.na(field_logs$legacy_row))
  if (!length(cand)) next
  if (length(cand) > 1 && !is.na(leg_lgr$legacy_time[i])) {
    dt <- abs(as.numeric(difftime(field_logs$real_start_posix[cand], leg_lgr$legacy_time[i], units = "mins")))
    dt[is.na(dt)] <- Inf
    cand <- cand[which.min(dt)]
    if (!is.finite(min(dt)) || min(dt) > 30) next
  } else cand <- cand[1]
  field_logs$legacy_row[cand] <- leg_lgr$legacy_row[i]
}
n_leg_matched <- sum(!is.na(field_logs$legacy_row))
log_step("3 matching", "2023-24 measurements matched to a field-log window", n_leg_matched,
         sprintf("by ID %d, by date + tree + time %d; of %d", sum(ok1), n_leg_matched - sum(ok1), nrow(leg_lgr)))
message("Legacy 2023-24 rows matched to a field-log window: ", n_leg_matched, " / ", nrow(leg_lgr))
unmatched_leg <- leg_lgr[!leg_lgr$legacy_row %in% field_logs$legacy_row, ]
message("  unmatched legacy rows by month: ",
        paste(names(table(substr(unmatched_leg$date, 1, 7))), table(substr(unmatched_leg$date, 1, 7)),
              sep = "=", collapse = ", "))
message("  field-log closures without a legacy row: ", sum(is.na(field_logs$legacy_row)), " (",
        paste(names(table(substr(field_logs$date[is.na(field_logs$legacy_row)], 1, 7))),
              table(substr(field_logs$date[is.na(field_logs$legacy_row)], 1, 7)), sep = "=", collapse = ", "), ")")

# ============================================================
# PART 4: LGR RAW TRACES
# ============================================================

message("\n=== Part 4: LGR raw traces ===")

LGR_BASE <- PATH_LGR
lgr_files <- list.files(LGR_BASE, pattern = "^micro_\\d{4}-\\d{2}-\\d{2}_f\\d+\\.txt(\\.zip)?$",
                        recursive = TRUE, full.names = TRUE)
lgr_files <- lgr_files[!grepl("/\\._", lgr_files)]
lgr_catalog <- data.frame(
  file = lgr_files,
  machine = sub("^([^/]+)/.*$", "\\1", sub(paste0("^", LGR_BASE, "/"), "", lgr_files)),
  date = as.Date(sub("^micro_(\\d{4}-\\d{2}-\\d{2})_.*$", "\\1", basename(lgr_files))),
  stringsAsFactors = FALSE)
message("LGR data files: ", nrow(lgr_catalog), " on ", length(unique(lgr_catalog$date)), " dates (",
        paste(names(table(lgr_catalog$machine)), table(lgr_catalog$machine), sep = "=", collapse = ", "), ")")

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
  d <- d[!duplicated(d$POSIX.time), ]           # "(1)" copies of the same file
  d
}

# Which machine's record to use for each field-log closure: the logged machine if it has a
# file that day, otherwise whichever LGR does (the logs default to LGR1)
field_logs$trace_machine <- NA_character_
for (i in seq_len(nrow(field_logs))) {
  have <- unique(lgr_catalog$machine[lgr_catalog$date == field_logs$date[i]])
  if (!length(have)) next
  field_logs$trace_machine[i] <- if (field_logs$machine[i] %in% have) field_logs$machine[i] else have[1]
}
message("Field-log closures with an LGR file that day: ", sum(!is.na(field_logs$trace_machine)),
        " / ", nrow(field_logs))
log_step("4 LGR traces", "closures recorded on a different LGR than logged (the logged analyzer had no file that day)",
         sum(!is.na(field_logs$trace_machine) & field_logs$trace_machine != field_logs$machine))

# End the window at chamber removal: the first abrupt CO2 drop from the running maximum
# (> max(10 ppm, 20 % of the rise so far) within 10 s). The LI-7810 remark often continues
# after the chamber is lifted; the same rule is applied to the UGGA field-log windows.
window_trims <- list()
trim_at_opening <- function(seg) {
  # opening = abrupt fall of the smoothed CO2 (7-point running median) by more than
  # max(10 ppm, 8 sigma, 25 % of the rise) within 10 s, after >= 60 s, that never recovers
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
# Manual review decisions (2_flux/04_review_windows.R): clicked windows replace the automatic ones
manual <- if (file.exists(PATH_MANUAL))
  read.csv(PATH_MANUAL, stringsAsFactors = FALSE) else
  data.frame(closure_id = character(), decision = character(), start = character(), end = character())
manual_click <- manual[manual$decision == "click", ]
apply_manual <- function(uid, start, end) {
  k <- match(uid, manual_click$closure_id)
  if (is.na(k)) return(list(start = start, end = end, manual = FALSE))
  list(start = as.POSIXct(manual_click$start[k], tz = TZ), end = as.POSIXct(manual_click$end[k], tz = TZ), manual = TRUE)
}
segment_trace <- function(day, start, end, uid) {
  mw <- apply_manual(uid, start, end); start <- mw$start; end <- mw$end
  seg <- day[day$POSIX.time >= start - SHOULDER_S & day$POSIX.time <= end + SHOULDER_S, ]
  if (!nrow(seg)) return(NULL)
  seg$UniqueID <- uid
  seg$flag <- as.numeric(seg$POSIX.time >= start & seg$POSIX.time <= end)
  if (sum(seg$flag) < 5) return(NULL)
  seg$start.time <- start; seg$end.time <- end
  seg$Etime <- as.numeric(seg$POSIX.time - start, units = "secs")
  seg$obs.length <- as.numeric(end - start, units = "secs")
  if (mw$manual) seg else trim_at_opening(seg)
}

# Clock check: the field logs record the LGR's own clock, so no shift is applied. A
# fluxqc::find_clock_offset() pass was tried and rejected: comp_start marks the start of the
# fitting segment (after the closure onset), so the onset score is not diagnostic here; the
# CO2 traces show the rises inside the logged windows. Misaligned or disturbed windows are
# left to the fluxqc/goFlux CO2 screens (co2_tracer), which flag them for review.
field_logs$clock_offset_s <- 0
lgr_segments <- list(); lgr_meta <- list()

# Analyzer precision per analyzer x field day, from the WHOLE day record (not just the closures):
# MAD of first differences / sqrt(2) per constant-interval run (fluxqc::precision_mad_runs); the
# dominant run is used. This follows the filtering paper (sigma once per analyzer x campaign) and
# absorbs drift, analyzer swaps and logging-interval changes.
day_sigma <- list()
sigma_day <- function(x, instrument, date) {
  one <- function(v) {
    r <- suppressWarnings(fluxqc::precision_mad_runs(x[[v]], x$POSIX.time))
    r <- r[which.max(r$n), ]
    c(sigma = r$sigma, dt = r$dt_s, runs = nrow(suppressWarnings(fluxqc::precision_mad_runs(x[[v]], x$POSIX.time))))
  }
  a <- one("CH4dry_ppb"); b <- one("CO2dry_ppm")
  data.frame(instrument = instrument, date = as.Date(date), n_rows = nrow(x),
             CH4_sigma_day = a[["sigma"]], CO2_sigma_day = b[["sigma"]], dt_day = a[["dt"]], n_interval_runs = a[["runs"]])
}
field_logs$closure_id <- sprintf("LGR_%s_%s_%s", format(field_logs$date, "%Y%m%d"), field_logs$tree,
                                 gsub(":", "", field_logs$comp_start))
day_keys <- unique(field_logs[!is.na(field_logs$trace_machine), c("date", "trace_machine")])
for (k in seq_len(nrow(day_keys))) {
  day <- load_lgr_day(day_keys$date[k], day_keys$trace_machine[k])
  if (is.null(day)) next
  day_sigma[[length(day_sigma) + 1]] <- sigma_day(day, day_keys$trace_machine[k], day_keys$date[k])
  rows <- which(field_logs$date == day_keys$date[k] & field_logs$trace_machine == day_keys$trace_machine[k])
  for (i in rows) {
    seg <- segment_trace(day, field_logs$comp_start_posix[i] + field_logs$deadband_s[i], field_logs$comp_end_posix[i],
                         field_logs$closure_id[i])
    if (is.null(seg)) next
    seg$instrument <- day_keys$trace_machine[k]
    lgr_segments[[field_logs$closure_id[i]]] <- seg
  }
  if (k %% 20 == 0) message("  ", k, "/", nrow(day_keys), " LGR days read; ", length(lgr_segments), " closures with traces")
}
message("LGR closures with a trace in the window: ", length(lgr_segments), " / ", nrow(field_logs))

# ============================================================
# PART 5: LI-7810 (2025) RAW TRACES -> REMARK WINDOWS
# ============================================================

message("\n=== Part 5: LI-7810 traces ===")

li_files <- list.files(PATH_LI7810, pattern = "\\.data$", full.names = TRUE)
message("LI-7810 day files: ", length(li_files))
li_raw <- do.call(rbind, lapply(li_files, function(f) {
  x <- tryCatch(import.LI7810(f, timezone = TZ, keep_all = TRUE), error = function(e) NULL)
  if (is.null(x)) { message("  failed: ", basename(f)); return(NULL) }
  x <- as.data.frame(x)[, c("POSIX.time", "DATE", "REMARK", "CH4dry_ppb", "CO2dry_ppm", "H2O_ppm")]
  x$source_file <- basename(f); x
}))
li_raw <- li_raw[order(li_raw$POSIX.time), ]
li_raw <- li_raw[!duplicated(li_raw$POSIX.time), ]
li_raw$REMARK <- trimws(as.character(li_raw$REMARK))
li_raw$REMARK[is.na(li_raw$REMARK)] <- ""

li_raw$instrument <- "LI-7810"
li_raw$day <- as.Date(li_raw$POSIX.time, tz = TZ)
for (dd in split(li_raw, li_raw$day)) {
  dd <- dd[is.finite(dd$CH4dry_ppb) & is.finite(dd$CO2dry_ppm), ]
  if (nrow(dd) > 100) day_sigma[[length(day_sigma) + 1]] <- sigma_day(dd, "LI-7810", dd$day[1])
}

# Remark runs: consecutive rows with the same non-empty REMARK and no gap > 5 s
r <- li_raw$REMARK; nzr <- nz(r)
brk <- c(TRUE, r[-1] != r[-length(r)] | diff(as.numeric(li_raw$POSIX.time)) > 5)
li_raw$run <- cumsum(brk)
runs <- li_raw %>% filter(nz(REMARK)) %>% group_by(run) %>%
  summarise(REMARK = REMARK[1], date = as.Date(POSIX.time[1], tz = TZ), start = min(POSIX.time),
            end = max(POSIX.time), n = n(), .groups = "drop") %>%
  mutate(dur = as.numeric(end - start, units = "secs"),
         tree = suppressWarnings(as.integer(sub("R$", "", REMARK))),
         redo = grepl("R$", REMARK), short = dur < MIN_REMARK_S)
message("Remark runs: ", nrow(runs), "; tree-tag remarks: ", sum(!is.na(runs$tree)),
        "; redo (…R): ", sum(runs$redo), "; < ", MIN_REMARK_S, " s: ", sum(runs$short))
# The team's rule: a "412R" redo supersedes the earlier "412" that day; short runs are aborted starts
log_step("5 LI-7810 remarks", "remark runs that are not a tree tag dropped", sum(is.na(runs$tree)))
log_step("5 LI-7810 remarks", sprintf("aborted starts (remark < %d s) dropped", MIN_REMARK_S), sum(!is.na(runs$tree) & runs$short))
runs <- runs %>% filter(!is.na(tree), !short)
superseded <- runs %>% filter(redo) %>% select(date, tree) %>% mutate(sup = TRUE)
runs <- runs %>% left_join(superseded, by = c("date", "tree")) %>%
  mutate(superseded = !redo & !is.na(sup)) %>% select(-sup)
message("  superseded by a redo: ", sum(runs$superseded))
log_step("5 LI-7810 remarks", "measurement superseded by a redo the same day (\"412R\" replaces \"412\")", sum(runs$superseded))
runs <- runs %>% filter(!superseded) %>% arrange(start) %>%
  group_by(date, tree) %>% mutate(k = row_number()) %>% ungroup() %>%
  mutate(closure_id = sprintf("LI7810_%s_%s_%d", format(date, "%Y%m%d"), tree, k),
         win_start = start + DEADBAND_7810, win_end = end)

li_segments <- list()
for (i in seq_len(nrow(runs))) {
  day <- li_raw[li_raw$POSIX.time >= runs$start[i] - SHOULDER_S & li_raw$POSIX.time <= runs$end[i] + SHOULDER_S,
                c("POSIX.time", "CH4dry_ppb", "CO2dry_ppm", "H2O_ppm", "instrument")]
  seg <- segment_trace(day, runs$win_start[i], runs$win_end[i], runs$closure_id[i])
  if (!is.null(seg)) li_segments[[runs$closure_id[i]]] <- seg
}
message("LI-7810 closures with windows: ", length(li_segments))

# Match to legacy 2025 rows: same date and tree, nearest legacy time (datetime_posx = remark start
# floored to the minute), within 10 min
leg_li <- legacy %>% filter(year == 2025)
runs$legacy_row <- NA_integer_
for (i in seq_len(nrow(runs))) {
  cand <- which(leg_li$date == runs$date[i] & leg_li$Tree == runs$tree[i] &
                  !leg_li$legacy_row %in% runs$legacy_row)
  if (!length(cand)) next
  dt <- abs(as.numeric(difftime(leg_li$legacy_time[cand], runs$start[i], units = "mins")))
  j <- which.min(dt); if (dt[j] <= 10) runs$legacy_row[i] <- leg_li$legacy_row[cand[j]]
}
log_step("5 LI-7810 remarks", "LI-7810 closures (one per tree and remark)", nrow(runs))
log_step("5 LI-7810 remarks", "2025 measurements matched to a remark (same date and tree, within 10 min)", sum(!is.na(runs$legacy_row)),
         sprintf("of %d", nrow(leg_li)))
message("Legacy 2025 rows matched to a remark: ", sum(!is.na(runs$legacy_row)), " / ", nrow(leg_li),
        "; unmatched legacy 2025 dates: ",
        paste(names(table(leg_li$date[!leg_li$legacy_row %in% runs$legacy_row])),
              table(leg_li$date[!leg_li$legacy_row %in% runs$legacy_row]), sep = "=", collapse = ", "))

# ============================================================
# PART 6: CLOSURE TABLE + CHAMBER GEOMETRY
# ============================================================

message("\n=== Part 6: closure table and geometry ===")

closures <- bind_rows(
  field_logs %>% transmute(closure_id, legacy_row, date, tree, instrument = trace_machine, clock_offset_s,
                           machine_logged = machine,
                           window_src = ifelse(start_adjusted, "field log, start hand-adjusted by team (Timing Updates.xlsx)",
                                               sprintf("field log start/end + %d s deadband", DEADBAND_LGR)),
                           closure_start = comp_start_posix,
                           window_start = comp_start_posix + deadband_s, window_end = comp_end_posix,
                           real_start = real_start_posix, fl_source = source, fl_check, fl_notes,
                           UniqueID_fl),
  runs %>% transmute(closure_id, legacy_row, date, tree, instrument = "LI-7810", clock_offset_s = 0,
                     machine_logged = "LI-7810", closure_start = start,
                     window_src = sprintf("LI-7810 remark + %d s deadband", DEADBAND_7810),
                     window_start = win_start, window_end = win_end, real_start = start,
                     fl_source = "LI-7810 REMARK", fl_check = NA_character_,
                     fl_notes = ifelse(redo, "redo remark (…R)", NA_character_), UniqueID_fl = REMARK)
)
closures$has_trace <- closures$closure_id %in% c(names(lgr_segments), names(li_segments))
closures$manual_window <- closures$closure_id %in% manual_click$closure_id
closures$manual_exclude <- closures$closure_id %in% manual$closure_id[manual$decision == "exclude"]
if (nrow(manual)) message("Manual decisions applied: ", sum(closures$manual_window), " clicked windows, ",
                          sum(closures$manual_exclude), " exclusions")
log_step("6 windows", "windows set by hand after inspection (qc_decisions/manual_windows.csv)", sum(closures$manual_window))
log_step("6 windows", "measurements marked for exclusion after inspection (applied in 06_flux_dataset.R)", sum(closures$manual_exclude))
closures$end_trim_s <- unlist(window_trims)[closures$closure_id]
closures$end_trim_s[is.na(closures$end_trim_s)] <- 0
closures$window_end <- closures$window_end - closures$end_trim_s
message("Windows ended early at a detected chamber opening: ", sum(closures$end_trim_s > 0),
        " (", paste(names(table(closures$instrument[closures$end_trim_s > 0])),
                    table(closures$instrument[closures$end_trim_s > 0]), sep = "=", collapse = ", "),
        "); median trim ", median(closures$end_trim_s[closures$end_trim_s > 0]), " s")
log_step("6 windows", "windows ended early at a detected chamber opening (abrupt, sustained CO2 fall)", sum(closures$end_trim_s > 0),
         sprintf("median %.0f s removed", median(closures$end_trim_s[closures$end_trim_s > 0])))
# the legacy Tree column is curated (typos fixed); prefer it where a legacy row is matched
lt <- legacy$Tree[match(closures$legacy_row, legacy$legacy_row)]
closures$tree_logged <- closures$tree
closures$tree <- ifelse(!is.na(lt), lt, closures$tree)
log_step("6 windows", "tree tag corrected to the curated tag of the matched measurement", sum(!is.na(lt) & lt != closures$tree_logged, na.rm = TRUE))
# legacy rows with no closure at all (no field-log entry / no remark) are carried as their own rows
orphan <- legacy %>% filter(!legacy_row %in% closures$legacy_row) %>%
  transmute(closure_id = sprintf("LEGACY_%04d", legacy_row), legacy_row, date, tree = Tree,
            instrument = ifelse(year == 2025, "LI-7810", NA_character_), clock_offset_s = NA_real_,
            machine_logged = NA_character_, window_src = "none (legacy row without field-log window)",
            window_start = as.POSIXct(NA, tz = TZ), window_end = as.POSIXct(NA, tz = TZ),
            real_start = legacy_time, fl_source = NA_character_, fl_check = NA_character_,
            fl_notes = NA_character_, UniqueID_fl = NA_character_, has_trace = FALSE)
closures <- bind_rows(closures, orphan)
log_step("6 windows", "measurements without any recorded window or raw record (earlier flux carried)", nrow(orphan))
message("Closures: ", nrow(closures), " (with trace: ", sum(closures$has_trace),
        "; legacy rows without window: ", nrow(orphan), ")")

# geometry: per-tree measured collar volume (+ tubing/analyzer), fixed collar area
tv <- read.csv(PATH_VOLUMES, stringsAsFactors = FALSE, check.names = FALSE)
tv$tree <- as.integer(round(tv$Tag))
tv$Vcollar <- tv$Volume
# met: from the legacy row (Fisher met matched upstream) else the hourly Fisher record
met <- read.csv(PATH_MET, stringsAsFactors = FALSE)
met$t <- as.POSIXct(sub("Z$", "", met$datetime), format = "%Y-%m-%dT%H:%M:%S", tz = TZ)
met <- met[!is.na(met$t) & !is.na(met$tair_C) & !is.na(met$p_kPa), ]
nearest_met <- function(t) {
  j <- findInterval(as.numeric(t), as.numeric(met$t))
  j <- pmax(1, pmin(j, nrow(met)))
  j2 <- pmin(j + 1, nrow(met))
  use2 <- abs(as.numeric(met$t[j2]) - as.numeric(t)) < abs(as.numeric(met$t[j]) - as.numeric(t))
  j[use2] <- j2[use2]
  data.frame(tair = met$tair_C[j], p_kPa = met$p_kPa[j])
}
closures <- closures %>%
  left_join(tv[, c("tree", "Vcollar")], by = "tree") %>%
  left_join(legacy %>% select(legacy_row, airt, bar, PLOT, SPECIES, DBH, location, year), by = "legacy_row") %>%
  mutate(is_7810 = instrument %in% "LI-7810" | (is.na(instrument) & !is.na(year) & year == 2025),
         Vextra = ifelse(is_7810, EXTRA_VOL_7810_L, EXTRA_VOL_LGR_L),
         Vtot = Vcollar + Vextra,
         Area = AREA_CM2,
         Pcham = bar / 10, Tcham = airt)      # legacy bar = hPa (millibar) -> kPa
need_met <- is.na(closures$Pcham) | is.na(closures$Tcham)
if (any(need_met)) {
  nm <- nearest_met(coalesce(closures$window_start, closures$real_start)[need_met])
  closures$Pcham[need_met] <- coalesce(closures$Pcham[need_met], nm$p_kPa)
  closures$Tcham[need_met] <- coalesce(closures$Tcham[need_met], nm$tair)
  closures$met_src <- ifelse(need_met, "Fisher met hourly (wtd_met.csv)", "legacy row (Fisher met)")
} else closures$met_src <- "legacy row (Fisher met)"
# tree metadata for closures not in the legacy set
tree_info <- legacy %>% filter(!is.na(PLOT)) %>% group_by(Tree) %>%
  summarise(PLOT_t = first(na.omit(PLOT)), SPECIES_t = first(na.omit(SPECIES)),
            DBH_t = first(na.omit(DBH)), location_t = first(na.omit(location)), .groups = "drop")
closures <- closures %>% left_join(tree_info, by = c("tree" = "Tree")) %>%
  mutate(PLOT = coalesce(PLOT, PLOT_t), SPECIES = coalesce(SPECIES, SPECIES_t),
         DBH = coalesce(DBH, DBH_t), location = coalesce(location, location_t)) %>%
  select(-PLOT_t, -SPECIES_t, -DBH_t, -location_t, -is_7810, -year)
message("Geometry: Vtot missing for ", sum(is.na(closures$Vtot) & closures$has_trace),
        " closures with traces; Pcham/Tcham from hourly met for ", sum(need_met), " closures")
log_step("6 windows", "chamber air temperature and pressure from the hourly Fisher record (no met in the earlier dataset)", sum(need_met))
log_step("6 windows", "closures with a raw trace but no chamber volume (not fitted)", sum(is.na(closures$Vtot) & closures$has_trace))

# day_sigma: analyzer precision per analyzer x day, from the whole day record
saveRDS(list(closures = closures, field_logs = field_logs, runs = runs, lgr_segments = lgr_segments,
             li_segments = li_segments, legacy = legacy, day_sigma = do.call(rbind, day_sigma)), PATH_CLOSURES)
write_log()
message("Wrote ", PATH_CLOSURES, ": ", nrow(closures), " closures (with trace: ", sum(closures$has_trace), ")")


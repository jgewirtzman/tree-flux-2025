# ============================================================
# 09_goflux_reprocess.R
#
# Reprocess every Harvard Forest stem-chamber closure that has a raw
# analyzer trace through goFlux (LM / HM, best.flux with the Hüppi et al.
# 2018 criteria) and fluxqc (empirical precision, minimum detectable flux,
# physical QC screens), so that the published dataset and the papers that
# depend on it share one flux fit.
#
# Why this exists (found 25 Sep 2026, see ch4-data-filtering/WORKLOG_2026-09.md):
#   1. The legacy 2023-24 fluxes (Matthes_Lab/stem-CH4-flux, linear model on
#      field-log windows) wrote CH4_SE as SE_slope * nmol in ppm with the
#      "/surfarea" commented out; the old "x1000" correction was insufficient.
#   2. The old MDF code used t = nb.obs (only seconds at 1 Hz), a per-closure
#      x3.t "Christiansen" MDF (misattributed) and datasheet precision.
#   3. The 2025 legacy fluxes came from a different pipeline (best 60-s
#      window by CO2 R2), so the two analyzers were not fitted the same way.
#
# Conventions (shared with the guidelines paper / fluxqc >= 0.2.3):
#   - Windows come from the field-log start/end times (2023-24, LGR/UGGA) or
#     from the analyzer REMARK span with a fixed deadband (2025, LI-7810).
#     Nothing is clicked by hand; QC screens flag closures for review.
#   - Flux = goFlux::best.flux (LM or HM, Hüppi criteria, g.limit = 2).
#   - MDF = 1.96 * sigma / t * flux.term (95 %, two-sided), sigma = MAD of
#     first differences over each analyzer's whole record per constant-
#     interval run (fluxqc::flag_detection, precision = "mad"), t = closure
#     length in seconds (fluxqc::closure_seconds). Retain-and-flag.
#
# Inputs
#   data/input/HF_2023-2025_tree_flux_corrected.csv      legacy dataset (n = 1,640)
#   data/raw/upland_wetland/processing_csvs/data/*.csv    field logs 2023 - Sep 2024
#   data/raw/upland_wetland/Sept2024_onwards/*.xlsx       field logs Sep 2024 - Mar 2025
#   data/raw/upland_wetland/processing_csvs/tree_volumes.csv
#   data/raw/7810_Processed/Tree_Fluxes/raw-7810-tree_flux_data/*.data   LI-7810 (1 Hz)
#   Matthes_Lab/stem-CH4-flux/Raw LGR Data/LGR{1,2,3}/   LGR day files (1 Hz; read-only)
#   Matthes_Lab/stem-CH4-flux/raw-7810-data/             Aug 2025 LI-7810 files (read-only)
#   data/processed/wtd_met.csv                            Fisher met for closures not in legacy
#
# Outputs
#   data/input/HF_2023-2025_tree_flux_goflux.csv          one row per closure; legacy and
#                                                         goFlux fluxes side by side
#   data/processed/goflux_traces.rds                      the trace segments used (for plots)
#   data/processed/goflux_settings.json                   fluxqc settings record
#   outputs/tables/goflux_vs_legacy_by_analyzer_season.csv
#   outputs/tables/goflux_qc_review.csv                   closures the QC screens caught
#   outputs/tables/goflux_match_report.csv                matching / coverage summary
#   outputs/figures/goflux/flux_plots_*.pdf               goFlux trace plots (all closures)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(goFlux)
  library(fluxqc)
})
stopifnot(packageVersion("fluxqc") >= "0.2.3")

# ============================================================
# CONSTANTS
# ============================================================

TZ            <- "America/New_York"   # wall clock of both analyzers and the field logs
SURFAREA_M2   <- pi * 0.0508^2        # m2, collar inner radius 5.08 cm (all collars)
AREA_CM2      <- SURFAREA_M2 * 1e4    # goFlux wants cm2
# Analyzer internal volume added to the measured collar volume (tree_volumes.csv already holds
# collar + cap + 28.9 cm3 of tubing). One value for both analyzers (user decision, 30 Sep 2026,
# shared with ch4-data-filtering; see its WORKLOG_2026-09.md "Resolved 2026-09-30 - analyzer
# volume"): 0.028 L, the LI-COR LI-7810 total sample volume; goFlux lists 25 cm3 for the LGR
# microportable GLA131 cavity, within 10 %.
# The legacy LGR pipeline (Matthes_Lab .../diurnal_flux_processing.R) used lgr_volume <- .2 with
# no documented source, so 2023-24 legacy fluxes are ~1/0.70 of the goFlux values by volume alone.
EXTRA_VOL_LGR_L  <- 0.028
EXTRA_VOL_7810_L <- 0.028
SHOULDER_S    <- 120                  # s of trace kept outside the window (flag = 0), for plots/QC context
DEADBAND_7810 <- 20                   # s dropped after the LI-7810 REMARK starts (team convention)
DEADBAND_LGR  <- 20                   # s dropped after the field-log start (chamber closure) for the UGGA:
                                      # same deadband for both analyzers; removes the placement transient
                                      # (scripts/01_import/13_deadband_sensitivity.R: 0-45 s tested)
MIN_REMARK_S  <- 90                   # LI-7810 remarks shorter than this are aborted starts
MIN_WINDOW_S  <- 30                   # windows shorter than this are not fitted
PREC_LGR      <- c(CO2 = 0.35, CH4 = 0.9, H2O = 100)   # datasheet 1-s precision, GLA131 (ppm, ppb, ppm)
PREC_7810     <- c(CO2 = 3.5,  CH4 = 0.6, H2O = 45)    # LI-7810 (goFlux default)
INCLUDE_NEW_CLOSURES <- FALSE         # closures with a trace but no legacy row (Jan-Mar 2025 LGR
                                      # days) are written to the output but, by default, not
                                      # promoted into the analysis dataset (keeps n comparable)

LGR_BASE <- "/Users/jongewirtzman/My Drive/Matthes_Lab/stem-CH4-flux/Raw LGR Data"
LI7810_DIRS <- c(
  "data/raw/7810_Processed/Tree_Fluxes/raw-7810-tree_flux_data",
  "/Users/jongewirtzman/My Drive/Matthes_Lab/stem-CH4-flux/raw-7810-data"   # Aug 2025 only lives here
)
OUT_CSV      <- file.path("data", "input", "HF_2023-2025_tree_flux_goflux.csv")
OUT_TRACES   <- file.path("data", "processed", "goflux_traces.rds")
OUT_SETTINGS <- file.path("data", "processed", "goflux_settings.json")
TAB_DIR      <- file.path("outputs", "tables")
FIG_DIR      <- file.path("outputs", "figures", "goflux")
for (d in c(dirname(OUT_CSV), dirname(OUT_TRACES), TAB_DIR, FIG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# ============================================================
# PART 1: LEGACY DATASET
# ============================================================

message("=== Part 1: legacy dataset ===")

legacy <- read.csv("data/input/HF_2023-2025_tree_flux_corrected.csv",
                   stringsAsFactors = FALSE)
legacy$legacy_row <- seq_len(nrow(legacy))
# datetime_posx is local wall time written with a spurious "Z"
legacy$legacy_time <- as.POSIXct(sub("Z$", "", legacy$datetime_posx),
                                 format = "%Y-%m-%dT%H:%M:%S", tz = TZ)
legacy$date <- as.Date(legacy$date)
legacy$UniqueID_legacy <- legacy$UniqueID
message("Legacy rows: ", nrow(legacy), " (",
        paste(names(table(legacy$year)), table(legacy$year), sep = ": ", collapse = ", "), ")")

# ============================================================
# PART 2: FIELD LOGS -> CLOSURE WINDOWS (LGR/UGGA, Jun 2023 - Mar 2025)
# ============================================================

message("\n=== Part 2: field logs ===")

fl_dir   <- file.path("data", "raw", "upland_wetland", "processing_csvs", "data")
xlsx_dir <- file.path("data", "raw", "upland_wetland", "Sept2024_onwards")

nz <- function(x) !is.na(x) & nchar(trimws(as.character(x))) > 0
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
tu_path <- file.path("data", "raw", "upland_wetland", "May23_Sept24", "Timing Updates.xlsx")
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
}

# Sep 2024 - Mar 2025 xlsx logs (the LGR was used until the LI-7810 arrived in Apr 2025;
# the Jan-Mar 2025 files must NOT be skipped)
xlsx_files <- list.files(xlsx_dir, pattern = "[.]xlsx$", full.names = TRUE)
xlsx_files <- xlsx_files[!grepl("Template", xlsx_files)]
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
# hand-adjusted starts (Timing Updates.xlsx) already skip the transient: no extra deadband
field_logs$start_adjusted <- field_logs$start_adjusted %in% TRUE
field_logs$deadband_s <- ifelse(field_logs$start_adjusted, 0, DEADBAND_LGR)
field_logs$machine <- gsub("LGR #", "LGR", field_logs$machine)
field_logs$machine <- gsub("\\s+", "", field_logs$machine)
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
tv0 <- read.csv(file.path("data", "raw", "upland_wetland", "processing_csvs", "tree_volumes.csv"),
                stringsAsFactors = FALSE, check.names = FALSE)
old_map <- setNames(as.integer(round(tv0$Tag)), as.integer(tv0$`Old tag`))
old_map <- old_map[!is.na(names(old_map))]
is_old <- !is.na(field_logs$tree) & field_logs$tree %in% as.integer(names(old_map)) &
  !field_logs$tree %in% as.integer(round(tv0$Tag))
field_logs$tree_logged <- field_logs$tree
field_logs$tree[is_old] <- unname(old_map[as.character(field_logs$tree[is_old])])
message("Field-log old wetland tags mapped to ForestGEO tags: ", sum(is_old))
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
field_logs_bad <- field_logs[bad_win, ]
field_logs <- field_logs[!bad_win, ]

# De-duplicate entries logged in more than one file (prefer the later, QC'd log)
field_logs <- field_logs %>%
  mutate(priority = case_when(grepl("^xlsx", source) ~ 4, source == "summer2024" ~ 3,
                              source == "updated2" ~ 2, TRUE ~ 1),
         dup_key = ifelse(nz(UniqueID_fl), UniqueID_fl, paste(date, tree, comp_start))) %>%
  group_by(dup_key) %>% slice_max(priority, n = 1, with_ties = FALSE) %>% ungroup() %>%
  # also collapse the same closure logged under different IDs
  group_by(date, tree, comp_start) %>% slice_max(priority, n = 1, with_ties = FALSE) %>% ungroup() %>%
  select(-priority, -dup_key) %>% arrange(comp_start_posix)
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
# Manual review decisions (14_manual_windows.R): clicked windows replace the automatic ones
manual <- if (file.exists(file.path("data", "input", "manual_windows.csv")))
  read.csv(file.path("data", "input", "manual_windows.csv"), stringsAsFactors = FALSE) else
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

li_files <- unlist(lapply(LI7810_DIRS, list.files, pattern = "\\.data$", full.names = TRUE))
li_files <- li_files[!duplicated(basename(li_files))]
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
runs <- runs %>% filter(!is.na(tree), !short)
superseded <- runs %>% filter(redo) %>% select(date, tree) %>% mutate(sup = TRUE)
runs <- runs %>% left_join(superseded, by = c("date", "tree")) %>%
  mutate(superseded = !redo & !is.na(sup)) %>% select(-sup)
message("  superseded by a redo: ", sum(runs$superseded))
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
closures$end_trim_s <- unlist(window_trims)[closures$closure_id]
closures$end_trim_s[is.na(closures$end_trim_s)] <- 0
closures$window_end <- closures$window_end - closures$end_trim_s
message("Windows ended early at a detected chamber opening: ", sum(closures$end_trim_s > 0),
        " (", paste(names(table(closures$instrument[closures$end_trim_s > 0])),
                    table(closures$instrument[closures$end_trim_s > 0]), sep = "=", collapse = ", "),
        "); median trim ", median(closures$end_trim_s[closures$end_trim_s > 0]), " s")
# the legacy Tree column is curated (typos fixed); prefer it where a legacy row is matched
lt <- legacy$Tree[match(closures$legacy_row, legacy$legacy_row)]
closures$tree_logged <- closures$tree
closures$tree <- ifelse(!is.na(lt), lt, closures$tree)
# legacy rows with no closure at all (no field-log entry / no remark) are carried as their own rows
orphan <- legacy %>% filter(!legacy_row %in% closures$legacy_row) %>%
  transmute(closure_id = sprintf("LEGACY_%04d", legacy_row), legacy_row, date, tree = Tree,
            instrument = ifelse(year == 2025, "LI-7810", NA_character_), clock_offset_s = NA_real_,
            machine_logged = NA_character_, window_src = "none (legacy row without field-log window)",
            window_start = as.POSIXct(NA, tz = TZ), window_end = as.POSIXct(NA, tz = TZ),
            real_start = legacy_time, fl_source = NA_character_, fl_check = NA_character_,
            fl_notes = NA_character_, UniqueID_fl = NA_character_, has_trace = FALSE)
closures <- bind_rows(closures, orphan)
message("Closures: ", nrow(closures), " (with trace: ", sum(closures$has_trace),
        "; legacy rows without window: ", nrow(orphan), ")")

# geometry: per-tree measured collar volume (+ tubing/analyzer), fixed collar area
tv <- read.csv(file.path("data", "raw", "upland_wetland", "processing_csvs", "tree_volumes.csv"),
               stringsAsFactors = FALSE, check.names = FALSE)
tv$tree <- as.integer(round(tv$Tag))
tv$Vcollar <- tv$Volume
# met: from the legacy row (Fisher met matched upstream) else the hourly Fisher record
met <- read.csv(file.path("data", "processed", "wtd_met.csv"), stringsAsFactors = FALSE)
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
if (Sys.getenv("GOFLUX_DRYRUN") == "1") { saveRDS(list(closures = closures, field_logs = field_logs, runs = runs,
  lgr_segments = lgr_segments, li_segments = li_segments, legacy = legacy), "/private/tmp/claude-501/-Users-jongewirtzman-My-Drive-Research-tree-flux-2025/e55182b0-acd5-4758-82f7-33e91b24a757/scratchpad/dryrun.rds"); message("DRY RUN: stopping after Part 6"); quit(save = "no") }

# ============================================================
# PART 7: goFlux + fluxqc
# ============================================================

message("\n=== Part 7: goFlux / fluxqc ===")

segs <- c(lgr_segments, li_segments)
fit_ids <- closures$closure_id[closures$has_trace & !is.na(closures$Vtot)]
segs <- segs[fit_ids]
manID <- do.call(rbind, lapply(segs, function(s) {
  uid <- s$UniqueID[1]; c1 <- closures[closures$closure_id == uid, ]
  s$Vtot <- c1$Vtot; s$Area <- c1$Area; s$Pcham <- c1$Pcham; s$Tcham <- c1$Tcham
  s$DATE <- format(s$POSIX.time[1], "%Y-%m-%d")
  pr <- if (s$instrument[1] == "LI-7810") PREC_7810 else PREC_LGR
  s$CO2_prec <- pr[["CO2"]]; s$CH4_prec <- pr[["CH4"]]; s$H2O_prec <- pr[["H2O"]]
  s$H2O_ppm[is.na(s$H2O_ppm)] <- 0
  s$obs.length_corr <- s$obs.length; s$start.time_corr <- s$start.time; s$end.time_corr <- s$end.time
  s
}))
rownames(manID) <- NULL
manID$UniqueID <- as.character(manID$UniqueID)
message("Trace rows: ", nrow(manID), " for ", length(unique(manID$UniqueID)), " closures")
aux <- closures %>% filter(closure_id %in% fit_ids) %>% transmute(UniqueID = closure_id, instrument)

t0 <- Sys.time()
res_co2 <- process_fluxes(manID, aux = aux, gastype = "CO2dry_ppm", precision = "mad",
                          mdf = "wassmann", conf = 0.95, group = "instrument", qc = FALSE,
                          extra = list(precision = c("mad", "allan", "datasheet"), mdf = c("wassmann", "goflux")))
message("CO2 done in ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
t0 <- Sys.time()
res_ch4 <- process_fluxes(manID, aux = aux, gastype = "CH4dry_ppb", precision = "mad",
                          mdf = "wassmann", conf = 0.95, group = "instrument",
                          extra = "all", co2 = res_co2$fluxes,
                          # ambient_start is off: the field-log windows begin after the closure onset,
                          # so the first seconds of a window are never at ambient by design
                          qc = list(c0 = TRUE, co2_tracer = TRUE, convex = TRUE, min_window = list(secs = 60),
                                    ambient_start = FALSE, noisy = TRUE))
message("CH4 done in ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
print(res_ch4)
message("Campaign sigma (ppb) by analyzer: ",
        paste(capture.output(print(res_ch4$fluxes %>% group_by(instrument) %>%
          summarise(n = n(), sigma_mad_ppb = round(median(sigma_emp), 3),
                    sigma_allan_med = round(median(sigma_allan, na.rm = TRUE), 3), .groups = "drop"))),
          collapse = "\n"))

# closure length and logging interval actually used
cs <- closure_seconds(manID)

# ============================================================
# PART 8: ASSEMBLE OUTPUT
# ============================================================

message("\n=== Part 8: output table ===")

pick_se <- function(f) ifelse(f$model == "HM" & !is.na(f$HM.SE), f$HM.SE, f$LM.SE)
gf_ch4 <- res_ch4$fluxes %>% transmute(
  closure_id = UniqueID,
  CH4_flux_goflux = best.flux, CH4_model = model, CH4_LM_flux = LM.flux, CH4_HM_flux = HM.flux,
  CH4_SE_goflux = pick_se(res_ch4$fluxes), CH4_LM_SE = LM.SE, CH4_LM_r2 = LM.r2, CH4_LM_p = LM.p.val,
  CH4_HM_r2 = HM.r2, CH4_g_factor = g.fact, CH4_HM_k = HM.k, CH4_k_max = k.max,
  CH4_C0 = C0, CH4_Ct = Ct, CH4_MAE_LM = LM.MAE, CH4_MAE_HM = HM.MAE,
  CH4_quality_check = quality.check, CH4_LM_diagnose = LM.diagnose, CH4_HM_diagnose = HM.diagnose,
  nb_obs = nb.obs, flux_term = flux.term,
  CH4_sigma_mad_record = sigma_emp, CH4_MDF_record = MDF_emp,
  CH4_sigma_allan = sigma_allan, CH4_sigma_datasheet = sigma_datasheet, CH4_sigma_rolling = sigma_rolling,
  CH4_MDF_wass90 = MDF_wassmann_mad_90, CH4_MDF_wass95 = MDF_wassmann_mad_95, CH4_MDF_wass99 = MDF_wassmann_mad_99,
  CH4_MDF_chr90 = MDF_christiansen_allan_90, CH4_MDF_chr95 = MDF_christiansen_allan_95,
  CH4_MDF_chr99 = MDF_christiansen_allan_99,
  CH4_MDF_datasheet = MDF_goflux_datasheet, CH4_MDF_2sigma_allan = MDF_two_sigma_allan,
  qc_c0, qc_c0_ratio, qc_co2_tracer, qc_co2_slope, qc_convex, qc_min_window,
  qc_noisy, qc_noisy_ratio, qc_any, qc_note)
gf_co2 <- res_co2$fluxes %>% transmute(
  closure_id = UniqueID,
  CO2_flux_goflux = best.flux, CO2_model = model, CO2_LM_flux = LM.flux, CO2_HM_flux = HM.flux,
  CO2_SE_goflux = pick_se(res_co2$fluxes), CO2_LM_r2 = LM.r2, CO2_LM_p = LM.p.val, CO2_g_factor = g.fact,
  CO2_C0 = C0, CO2_quality_check = quality.check,
  CO2_sigma_mad_record = sigma_emp, CO2_sigma_allan = sigma_allan, CO2_MDF_record = MDF_emp)

day_sigma <- do.call(rbind, day_sigma)
write.csv(day_sigma, file.path(TAB_DIR, "goflux_precision_by_analyzer_day.csv"), row.names = FALSE)
message("Per-day precision (median CH4 ppb): ",
        paste(tapply(round(day_sigma$CH4_sigma_day, 3), day_sigma$instrument, median), collapse = " / "),
        " for ", paste(names(table(day_sigma$instrument)), collapse = " / "),
        "; days with >1 logging interval: ", sum(day_sigma$n_interval_runs > 1))
out <- closures %>%
  left_join(day_sigma %>% select(instrument, date, CH4_sigma_day, CO2_sigma_day), by = c("instrument", "date")) %>%
  left_join(cs %>% transmute(closure_id = UniqueID, t_sec = t, dt_s = dt), by = "closure_id") %>%
  left_join(gf_ch4, by = "closure_id") %>%
  left_join(gf_co2, by = "closure_id") %>%
  left_join(legacy %>% transmute(
    legacy_row, UniqueID_legacy, datetime_posx, year, jday, rh, parr, sampling_round, new_round,
    days_since_last, species_label, qa_check_legacy = qa_check,
    CH4_flux_legacy = CH4_flux_nmolpm2ps, CH4_r2_legacy = CH4_r2, CH4_SE_legacy_raw = CH4_SE,
    # the 2023-24 SE was SE_slope * n(mol) in ppm without /surfarea -> x1000/surfarea puts it in
    # nmol m-2 s-1 like the 2025 SE (checked against the white-noise expectation in
    # ch4-data-filtering/12_se_analysis.R)
    CH4_SE_legacy = ifelse(year < 2025, CH4_SE * 1000 / SURFAREA_M2, CH4_SE),
    CO2_flux_legacy = CO2_flux_umolpm2ps, CO2_r2_legacy = CO2_r2, CO2_SE_legacy_raw = CO2_SE,
    CO2_SE_legacy = ifelse(year < 2025, CO2_SE / SURFAREA_M2, CO2_SE),
    legacy_nmol = nmol, legacy_vol_system = vol_system, legacy_REMARK_LENGTH = REMARK_LENGTH),
    by = "legacy_row") %>%
  mutate(
    # reference detection limit: 1.96 x per-day analyzer sigma / closure seconds x flux term
    CH4_sigma_mad = coalesce(CH4_sigma_day, CH4_sigma_mad_record),
    CO2_sigma_mad = coalesce(CO2_sigma_day, CO2_sigma_mad_record),
    CH4_MDF = abs(qnorm(0.975) * CH4_sigma_mad / t_sec * flux_term),
    CO2_MDF = abs(qnorm(0.975) * CO2_sigma_mad / t_sec * flux_term),
    CH4_below_MDF = ifelse(is.na(CH4_MDF), NA, abs(CH4_flux_goflux) < CH4_MDF),
    CO2_below_MDF = ifelse(is.na(CO2_MDF), NA, abs(CO2_flux_goflux) < CO2_MDF),
    CH4_det_class = ifelse(is.na(CH4_MDF), NA_character_, ifelse(CH4_flux_goflux > CH4_MDF, "emission",
                           ifelse(CH4_flux_goflux < -CH4_MDF, "uptake", "below detection"))),
    CH4_MDF_method = ifelse(is.na(CH4_MDF), NA_character_, "1.96 x MAD sigma (analyzer x day) / t x flux term"),
    CH4_MDF_wass90 = abs(qnorm(0.95) * CH4_sigma_mad / t_sec * flux_term),
    CH4_MDF_wass95 = CH4_MDF,
    CH4_MDF_wass99 = abs(qnorm(0.995) * CH4_sigma_mad / t_sec * flux_term),
    year = coalesce(year, as.integer(format(date, "%Y"))),
    jday = coalesce(jday, as.integer(format(date, "%j"))),
    in_legacy_dataset = !is.na(legacy_row),
    fitted = !is.na(CH4_flux_goflux),
    # implied legacy chamber moles: legacy flux / (my LM slope on the same window in ppb/s / area);
    # ~1 means the legacy pipeline used the same volume/moles as this run
    legacy_vol_ratio = ifelse(year < 2025 & fitted & !is.na(CH4_flux_legacy) & abs(CH4_LM_flux) > 0.05,
                              CH4_flux_legacy / CH4_LM_flux, NA_real_),
    t_src = case_when(fitted ~ "trace (closure_seconds)", !is.na(window_start) ~ "field log only",
                      TRUE ~ "none"),
    # canonical flux for the analysis dataset: goFlux where a trace exists, legacy otherwise
    flux_source = ifelse(fitted, "goFlux", "legacy"),
    CH4_flux_nmolpm2ps = ifelse(fitted, CH4_flux_goflux, CH4_flux_legacy),
    CH4_SE = ifelse(fitted, CH4_SE_goflux, CH4_SE_legacy),
    CH4_r2 = ifelse(fitted, CH4_LM_r2, CH4_r2_legacy),
    CO2_flux_umolpm2ps = ifelse(fitted, CO2_flux_goflux, CO2_flux_legacy),
    CO2_SE = ifelse(fitted, CO2_SE_goflux, CO2_SE_legacy),
    CO2_r2 = ifelse(fitted, CO2_LM_r2, CO2_r2_legacy),
    inst_label = case_when(!is.na(instrument) & instrument == "LI-7810" ~ "LI-7810",
                           !is.na(instrument) ~ "LGR/UGGA",
                           year == 2025 & !is.na(window_start) ~ "LGR/UGGA",   # Jan-Mar 2025 field-log rows
                           year == 2025 ~ "LI-7810",
                           TRUE ~ "LGR/UGGA"),
    season = case_when(month(date) %in% c(12, 1, 2) ~ "winter", month(date) %in% 3:5 ~ "spring",
                       month(date) %in% 6:8 ~ "summer", TRUE ~ "fall")
  ) %>%
  rename(Tree = tree, Vtot_L = Vtot, Vcollar_L = Vcollar, Vextra_L = Vextra, Area_cm2 = Area,
         Pcham_kPa = Pcham, Tcham_C = Tcham) %>%
  arrange(coalesce(window_start, real_start), Tree)

fmt_t <- function(x) format(x, "%Y-%m-%d %H:%M:%S")
out_csv <- out %>% mutate(across(where(function(v) inherits(v, "POSIXct")), fmt_t),
                          date = format(date, "%Y-%m-%d"))
write.csv(out_csv, OUT_CSV, row.names = FALSE, na = "NA")
saveRDS(list(traces = manID, closures = closures), OUT_TRACES)
writeLines(fluxqc:::to_json(list(CH4 = res_ch4$settings, CO2 = res_co2$settings,
                                 constants = list(SURFAREA_M2 = SURFAREA_M2, EXTRA_VOL_LGR_L = EXTRA_VOL_LGR_L, EXTRA_VOL_7810_L = EXTRA_VOL_7810_L,
                                                  SHOULDER_S = SHOULDER_S, DEADBAND_7810 = DEADBAND_7810, DEADBAND_LGR = DEADBAND_LGR,
                                                  MIN_REMARK_S = MIN_REMARK_S, TZ = TZ))), OUT_SETTINGS)
message("Wrote ", OUT_CSV, ": ", nrow(out), " rows (fitted: ", sum(out$fitted),
        "; legacy rows: ", sum(out$in_legacy_dataset), "; new closures: ", sum(!out$in_legacy_dataset), ")")

# ============================================================
# PART 9: REPORTS
# ============================================================

message("\n=== Part 9: reports ===")

# 9a. matching / coverage
cov <- out %>% group_by(inst_label, year) %>%
  summarise(closures = n(), in_legacy = sum(in_legacy_dataset), with_trace = sum(has_trace),
            n_fitted = sum(fitted), legacy_without_fit = sum(in_legacy_dataset & !fitted),
            new_closures = sum(!in_legacy_dataset), new_fitted = sum(!in_legacy_dataset & fitted),
            .groups = "drop")
print(as.data.frame(cov))
write.csv(cov, file.path(TAB_DIR, "goflux_match_report.csv"), row.names = FALSE)
# legacy rows without a trace, by date
miss <- out %>% filter(in_legacy_dataset, !fitted) %>% count(inst_label, date, window_src, name = "n_rows")
write.csv(miss, file.path(TAB_DIR, "goflux_legacy_rows_without_trace.csv"), row.names = FALSE)
message("Legacy rows without a goFlux fit: ", sum(miss$n_rows), " on ", nrow(miss), " date/source combos")

# 9b. legacy vs goFlux by analyzer and season
cmp <- out %>% filter(fitted, in_legacy_dataset, !is.na(CH4_flux_legacy)) %>%
  mutate(diff = CH4_flux_goflux - CH4_flux_legacy,
         ratio = ifelse(abs(CH4_flux_legacy) > 0.05, CH4_flux_goflux / CH4_flux_legacy, NA),
         same_sign = sign(CH4_flux_goflux) == sign(CH4_flux_legacy),
         within20 = abs(diff) <= 0.2 * abs(CH4_flux_legacy),
         se_ratio = CH4_SE_goflux / CH4_SE_legacy)
cmp_tab <- bind_rows(
  cmp %>% group_by(inst_label, season) %>% summarise(
    n = n(), r_pearson = cor(CH4_flux_goflux, CH4_flux_legacy), rho = cor(CH4_flux_goflux, CH4_flux_legacy, method = "spearman"),
    median_legacy = median(CH4_flux_legacy), median_goflux = median(CH4_flux_goflux),
    mean_legacy = mean(CH4_flux_legacy), mean_goflux = mean(CH4_flux_goflux),
    median_diff = median(diff), median_abs_diff = median(abs(diff)), median_ratio = median(ratio, na.rm = TRUE),
    pct_same_sign = 100 * mean(same_sign), pct_within_20pct = 100 * mean(within20, na.rm = TRUE),
    pct_HM = 100 * mean(CH4_model == "HM"), median_SE_ratio = median(se_ratio, na.rm = TRUE),
    pct_below_MDF = 100 * mean(CH4_below_MDF, na.rm = TRUE), pct_qc_any = 100 * mean(qc_any, na.rm = TRUE),
    .groups = "drop"),
  cmp %>% group_by(inst_label) %>% summarise(
    season = "all", n = n(), r_pearson = cor(CH4_flux_goflux, CH4_flux_legacy), rho = cor(CH4_flux_goflux, CH4_flux_legacy, method = "spearman"),
    median_legacy = median(CH4_flux_legacy), median_goflux = median(CH4_flux_goflux),
    mean_legacy = mean(CH4_flux_legacy), mean_goflux = mean(CH4_flux_goflux),
    median_diff = median(diff), median_abs_diff = median(abs(diff)), median_ratio = median(ratio, na.rm = TRUE),
    pct_same_sign = 100 * mean(same_sign), pct_within_20pct = 100 * mean(within20, na.rm = TRUE),
    pct_HM = 100 * mean(CH4_model == "HM"), median_SE_ratio = median(se_ratio, na.rm = TRUE),
    pct_below_MDF = 100 * mean(CH4_below_MDF, na.rm = TRUE), pct_qc_any = 100 * mean(qc_any, na.rm = TRUE),
    .groups = "drop")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(cmp_tab))
write.csv(cmp_tab, file.path(TAB_DIR, "goflux_vs_legacy_by_analyzer_season.csv"), row.names = FALSE)
message("Implied legacy/goFlux LM volume ratio (2023-24, |LM flux| > 0.05): median ",
        round(median(out$legacy_vol_ratio, na.rm = TRUE), 3), ", IQR ",
        paste(round(quantile(out$legacy_vol_ratio, c(.25, .75), na.rm = TRUE), 3), collapse = "-"))

# 9c. QC review list: everything a screen or goFlux diagnostic caught, plus unmatched/short windows
review <- out %>%
  filter(fitted) %>%
  mutate(reasons = paste0(
    ifelse(!is.na(qc_any) & qc_any, paste0("fluxqc: ", qc_note, "; "), ""),
    ifelse(!is.na(CO2_flux_goflux) & CO2_flux_goflux <= 0, "CO2 flux <= 0; ", ""),
    ifelse(!is.na(fl_check) & fl_check == "bad", paste0("field log check = bad (", trimws(fl_notes), "); "), ""),
    ifelse(!is.na(legacy_vol_ratio) & (legacy_vol_ratio < 0.75 | legacy_vol_ratio > 1.33),
           sprintf("legacy/goFlux LM ratio %.2f (window or volume differs); ", legacy_vol_ratio), ""),
    ifelse(fitted & in_legacy_dataset & !is.na(CH4_flux_legacy) &
             sign(CH4_flux_goflux) != sign(CH4_flux_legacy) & abs(CH4_flux_legacy) > 0.05,
           "sign differs from legacy; ", ""))) %>%
  filter(nz(reasons)) %>%
  mutate(n_reasons = lengths(regmatches(reasons, gregexpr("; ", reasons)))) %>%
  arrange(desc(n_reasons), date) %>%
  transmute(closure_id, date, Tree, PLOT, SPECIES, instrument, in_legacy_dataset, window_start, window_end,
            t_sec, CH4_flux_legacy, CH4_flux_goflux, CH4_model, CH4_g_factor, CH4_MDF, CH4_det_class,
            CO2_flux_goflux, qc_any, n_reasons, reasons)
write.csv(review %>% mutate(across(where(function(v) inherits(v, "POSIXct")), fmt_t)),
          file.path(TAB_DIR, "goflux_qc_review.csv"), row.names = FALSE)
message("QC review list: ", nrow(review), " closures -> outputs/tables/goflux_qc_review.csv")
message("  screens fired: ", paste(
  vapply(c("qc_c0", "qc_co2_tracer", "qc_convex", "qc_min_window", "qc_noisy"),
         function(s) sprintf("%s=%d", sub("qc_", "", s), sum(out[[s]], na.rm = TRUE)), character(1)), collapse = ", "))

# 9d. goFlux trace plots for every fitted closure (one pdf per analyzer x year)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  grp <- aux %>% left_join(out %>% select(closure_id, year), by = c("UniqueID" = "closure_id")) %>%
    mutate(g = paste0(gsub("[^A-Za-z0-9]", "", instrument), "_", year))
  for (gg in unique(grp$g)) {
    ids <- grp$UniqueID[grp$g == gg]
    pl <- tryCatch(flux.plot(res_ch4$fluxes[res_ch4$fluxes$UniqueID %in% ids, ],
                             manID[manID$UniqueID %in% ids, ], "CH4dry_ppb", shoulder = SHOULDER_S),
                   error = function(e) { message("flux.plot failed for ", gg, ": ", conditionMessage(e)); NULL })
    if (!is.null(pl)) invisible(capture.output(flux2pdf(pl, outfile = file.path(FIG_DIR, paste0("flux_plots_CH4_", gg, ".pdf")))))
  }
  message("Trace plots written to ", FIG_DIR)
}

message("\nDone. New dataset: ", normalizePath(OUT_CSV))

# ============================================================
# 07_filter_sensitivity_ridges.R
#
# CH4 flux distribution under quality filter criteria including
# MDF-based thresholds (manufacturer, Wassmann, Christiansen),
# empirical SNR (Allan deviation), SE-based SNR, and R^2 filters.
#
# Loads pre-computed quality flags from 09_quality_flags.R
# (scripts/01_import/09_quality_flags.R must be run first).
#
# Two instruments:
#   - 2025 LI-7810
#   - 2023-24 LGR/UGGA
# ============================================================

library(dplyr)
library(ggplot2)
library(ggridges)
library(tidyr)
library(scales)

# ============================================================
# CONSTANTS (for reference lines in plots)
# ============================================================

PREC_CH4_LGR   <- 0.9    # ppb, GLA131 (UGGA) manufacturer spec
PREC_CO2_LGR   <- 0.35   # ppm, GLA131
PREC_CH4_7810  <- 0.6    # ppb, LI-7810 manufacturer spec
PREC_CO2_7810  <- 3.5    # ppm, LI-7810

# ============================================================
# LOAD PRE-COMPUTED QUALITY FLAGS
# (from scripts/01_import/09_quality_flags.R)
# ============================================================

flagged_path <- file.path("data", "processed", "flux_with_quality_flags.csv")
if (!file.exists(flagged_path)) {
  stop("flux_with_quality_flags.csv not found. Run scripts/01_import/09_quality_flags.R first.")
}

df <- read.csv(flagged_path, stringsAsFactors = FALSE)
n_total <- nrow(df)
message("Loaded pre-computed quality flags: ", n_total, " measurements")
message("  2023-24 (LGR/UGGA): ", sum(df$year < 2025),
        " | 2025 (LI-7810): ", sum(df$year == 2025))

# Reconstruct per-instrument Allan deviation summaries for plotting
allan_df_lgr <- df %>%
  filter(year < 2025, !is.na(allan_sd_CH4)) %>%
  select(allan_sd_CO2, allan_sd_CH4)
allan_df_7810 <- df %>%
  filter(year == 2025, !is.na(allan_sd_CH4)) %>%
  select(allan_sd_CO2, allan_sd_CH4)

# NOTE: Parts 1-7 (Allan deviation, chamber geometry, field log parsing, MDF
# computation) have been moved to scripts/01_import/09_quality_flags.R.
# The code below was formerly Part 8+.

# DEAD CODE MARKER — everything between here and Part 8 is skipped
if (FALSE) {
# Preserved for git history reference only
fluxes_UNUSED <- data.frame()
df_UNUSED <- fluxes_UNUSED %>%
  mutate(
    CH4_SE_corr = if_else(year < 2025, CH4_SE * 1000, CH4_SE),
    # SE-based SNR
    CO2_snr_se = abs(CO2_flux_umolpm2ps) / CO2_SE,
    CH4_snr_se = abs(CH4_flux_nmolpm2ps) / CH4_SE_corr
  )

n_total <- nrow(df)
message("After NA removal: ", n_total, " measurements")
message("  2023-24 (LGR/UGGA): ", sum(df$year < 2025),
        " | 2025 (LI-7810): ", sum(df$year == 2025))

# ============================================================
# PART 2: CHAMBER GEOMETRY FROM tree_volumes.csv
# ============================================================

message("\n=== Loading chamber geometry ===")

tv <- read.csv(file.path("data", "raw", "upland_wetland",
                          "processing_csvs", "tree_volumes.csv"),
               stringsAsFactors = FALSE)
tv$Tag_int <- as.integer(round(tv$Tag))
# vol_system = Volume (L) + 0.028 L extra connecting tubing
tv$vol_system_tv <- tv$Volume + EXTRA_TUBE_VOL

message("tree_volumes.csv: ", nrow(tv), " trees")

# Join tree_volumes to corrected flux by Tree = Tag_int
df <- df %>%
  left_join(tv[, c("Tag_int", "vol_system_tv")],
            by = c("Tree" = "Tag_int"))

# Fill vol_system for ALL rows using tree_volumes.csv
# (2025 data already has vol_system; 2023-24 is NA)
df <- df %>%
  mutate(
    vol_system = if_else(is.na(vol_system), vol_system_tv, vol_system),
    surfarea   = if_else(is.na(surfarea), SURFAREA_M2, surfarea),
    nmol       = if_else(is.na(nmol) & !is.na(bar) & !is.na(airt) & !is.na(vol_system),
                          bar * vol_system / (R_hPa_L * (airt + 273.15)),
                          nmol)
  ) %>%
  select(-vol_system_tv)

message("Chamber geometry filled:")
message("  vol_system NAs: ", sum(is.na(df$vol_system)), " / ", n_total)
message("  nmol NAs:       ", sum(is.na(df$nmol)), " / ", n_total)
message("  surfarea NAs:   ", sum(is.na(df$surfarea)), " / ", n_total)

# ============================================================
# PART 3: LOAD FIELD LOGS (for LGR measurement timestamps)
# ============================================================

message("\n=== Loading field logs for LGR measurement timestamps ===")

field_data_dir <- file.path("data", "raw", "upland_wetland",
                             "processing_csvs", "data")

# --- Field log 1: Summer 2023 ---
fl1 <- read.csv(file.path(field_data_dir, "Field_Data_Monthly_Summer2023.csv"),
                stringsAsFactors = FALSE)
fl1_std <- fl1 %>%
  transmute(
    UniqueID       = UniqueID,
    date_raw       = Date,
    tree_raw       = as.character(Tree_Tag),
    comp_start     = comp_start_time,
    comp_end       = comp_end_time,
    real_start     = Real.start,
    machine        = trimws(Analyzer),
    source         = "summer2023"
  )

# --- Field log 2: Updated2 (Sep 2023 - Sep 2024) ---
fl2 <- read.csv(file.path(field_data_dir, "Field_Data_Monthly_updated2.csv"),
                stringsAsFactors = FALSE)
# Use the "Updated Start Time" / "Updated End Time" columns if available
# (these incorporate QC corrections), else fall back to comp_start_time
fl2_std <- fl2 %>%
  transmute(
    UniqueID       = UniqueID,
    date_raw       = format_Date,
    tree_raw       = as.character(Tree.Tag),
    comp_start     = if_else(nchar(trimws(Updated.Start.Time)) > 0,
                              Updated.Start.Time, comp_start_time),
    comp_end       = if_else(nchar(trimws(Updated.End.Time)) > 0,
                              Updated.End.Time, comp_end_time),
    real_start     = Real.start,
    machine        = trimws(Machine),
    source         = "updated2"
  )

# --- Field log 3: Summer 2024 ---
fl3 <- read.csv(file.path(field_data_dir, "summer_2024.csv"),
                stringsAsFactors = FALSE)
fl3_std <- fl3 %>%
  transmute(
    UniqueID       = UniqueID,
    date_raw       = format_Date,
    tree_raw       = as.character(Tree),
    comp_start     = comp_start_time,
    comp_end       = comp_end_time,
    real_start     = Real.start,
    machine        = "LGR1",   # no Machine column; all 2024 used LGR1
    source         = "summer2024"
  )

# --- Field log 4: Sept-Dec 2024 Excel files ---
xlsx_dir <- file.path("data", "raw", "upland_wetland", "Sept2024_onwards")
fl4_std <- data.frame()
if (requireNamespace("readxl", quietly = TRUE) && dir.exists(xlsx_dir)) {
  xlsx_files <- list.files(xlsx_dir, pattern = "[.]xlsx$", full.names = TRUE)
  xlsx_files <- xlsx_files[!grepl("Template", xlsx_files)]
  xlsx_files <- xlsx_files[!grepl("2025", xlsx_files)]  # skip 2025 (LI-7810)

  for (xf in xlsx_files) {
    tryCatch({
      xd <- as.data.frame(readxl::read_excel(xf))
      # Normalize column names (spaces → dots)
      names(xd) <- make.names(names(xd))
      xd_std <- xd %>%
        transmute(
          UniqueID   = as.character(UniqueID),
          date_raw   = as.character(format_Date),
          tree_raw   = as.character(Tree.Tag),
          comp_start = as.character(comp_start_time),
          comp_end   = as.character(comp_end_time),
          real_start = as.character(Real.start),
          machine    = if ("Machine" %in% names(xd))
                         as.character(Machine) else "LGR1",
          source     = "xlsx_late2024"
        )
      fl4_std <- bind_rows(fl4_std, xd_std)
    }, error = function(e) {
      message("  Warning: could not read ", basename(xf), ": ", e$message)
    })
  }
  # Clean up comp_start/comp_end that may have Excel datetime artifacts
  # (e.g., "1899-12-31 12:18:52" → extract just "12:18:52")
  fl4_std$comp_start <- sub("^.*\\s+", "", fl4_std$comp_start)
  fl4_std$comp_end   <- sub("^.*\\s+", "", fl4_std$comp_end)
  fl4_std$real_start <- sub("^.*\\s+", "", fl4_std$real_start)

  message("Excel field logs (Sept-Dec 2024): ", nrow(fl4_std), " entries")
} else {
  message("readxl not available or xlsx dir missing; skipping Excel field logs")
}

# Combine all field logs
field_logs <- bind_rows(fl1_std, fl2_std, fl3_std, fl4_std)
message("Combined field logs: ", nrow(field_logs), " entries")

# Standardize machine names
field_logs$machine <- gsub("LGR #", "LGR", field_logs$machine)
field_logs$machine <- gsub("\\s+", "", field_logs$machine)
field_logs$machine[is.na(field_logs$machine) |
                    field_logs$machine == ""] <- "LGR1"  # default

# Parse dates — handle multiple formats:
#   "9/21/2023" or "01/05/2024" (%m/%d/%Y)
#   "2/27/24" or "3/6/24"       (%m/%d/%y — 2-digit year!)
#   "6-29-2023"                  (%m-%d-%Y)
#   "24-02-27"                   (%y-%m-%d)
field_logs$date_parsed <- as.Date(field_logs$date_raw,
                                   format = "%m/%d/%Y")
# Fix 2-digit years that %Y parsed as year 24 AD instead of 2024
bad_year <- !is.na(field_logs$date_parsed) &
            as.integer(format(field_logs$date_parsed, "%Y")) < 100
if (any(bad_year)) {
  field_logs$date_parsed[bad_year] <- as.Date(
    field_logs$date_raw[bad_year], format = "%m/%d/%y")
}
# Try alternate format for dates like "6-29-2023"
na_dates <- is.na(field_logs$date_parsed)
if (any(na_dates)) {
  field_logs$date_parsed[na_dates] <- as.Date(
    field_logs$date_raw[na_dates], format = "%m-%d-%Y")
}
# Handle 2-digit years like "24-02-27"
na_dates <- is.na(field_logs$date_parsed)
if (any(na_dates)) {
  field_logs$date_parsed[na_dates] <- as.Date(
    field_logs$date_raw[na_dates], format = "%y-%m-%d")
}
# Handle ISO format "2024-09-17" (from readxl Date objects converted to char)
na_dates <- is.na(field_logs$date_parsed)
if (any(na_dates)) {
  field_logs$date_parsed[na_dates] <- as.Date(
    field_logs$date_raw[na_dates], format = "%Y-%m-%d")
}

message("Field log date range: ",
        min(field_logs$date_parsed, na.rm = TRUE), " to ",
        max(field_logs$date_parsed, na.rm = TRUE))
message("Dates with NAs: ", sum(is.na(field_logs$date_parsed)))

# Build comp_start/comp_end as POSIXct (date + time)
field_logs <- field_logs %>%
  filter(!is.na(date_parsed)) %>%
  mutate(
    comp_start_posix = as.POSIXct(paste(date_parsed, comp_start),
                                    format = "%Y-%m-%d %H:%M:%S"),
    comp_end_posix   = as.POSIXct(paste(date_parsed, comp_end),
                                    format = "%Y-%m-%d %H:%M:%S"),
    t_sec_fl         = as.numeric(difftime(comp_end_posix,
                                            comp_start_posix,
                                            units = "secs")),
    date_str         = format(date_parsed, "%Y-%m-%d")
  )

# Remove entries with unparseable times or unreasonable durations
field_logs <- field_logs %>%
  filter(!is.na(comp_start_posix), !is.na(comp_end_posix),
         t_sec_fl > 30, t_sec_fl < 1800)

message("After cleaning: ", nrow(field_logs), " field log entries with valid timestamps")
message("  Measurement duration range: ",
        round(min(field_logs$t_sec_fl)), " - ",
        round(max(field_logs$t_sec_fl)), " sec")

# De-duplicate: prefer updated2 over summer2023 (has QC corrections), and
# prefer summer2024 over updated2 (more complete for those dates)
field_logs <- field_logs %>%
  mutate(priority = case_when(
    source == "summer2024" ~ 3,
    source == "updated2"   ~ 2,
    source == "summer2023" ~ 1
  )) %>%
  group_by(UniqueID) %>%
  slice_max(priority, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(-priority)

message("After de-duplication: ", nrow(field_logs), " unique field log entries")

# ============================================================
# PART 4: ALLAN DEVIATION - 2025 LI-7810
# ============================================================

message("\n=== Allan deviation: 2025 LI-7810 ===")

raw_dir_7810 <- file.path("data", "raw", "7810_Processed", "Tree_Fluxes")
raw_files_7810 <- list.files(raw_dir_7810, pattern = "\\.csv$",
                              recursive = TRUE, full.names = TRUE)
raw_files_7810 <- raw_files_7810[grepl("/tables/", raw_files_7810)]
message("Found ", length(raw_files_7810), " raw 1-Hz LI-7810 files")

allan_7810 <- lapply(raw_files_7810, function(f) {
  tryCatch({
    raw <- read.csv(f, stringsAsFactors = FALSE)
    if (nrow(raw) < 5) return(NULL)
    data.frame(
      REMARK       = raw$REMARK[1],
      DATE         = raw$DATE[1],
      allan_sd_CO2 = allan_sd(raw$CO2),
      allan_sd_CH4 = allan_sd(raw$CH4),
      n_pts        = nrow(raw),
      t_sec        = nrow(raw),
      instrument   = "LI-7810",
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
})

allan_df_7810 <- do.call(rbind, Filter(Negate(is.null), allan_7810))

if (is.null(allan_df_7810) || nrow(allan_df_7810) == 0) {
  warning("No LI-7810 Allan results. Continuing with LGR only.")
  allan_df_7810 <- data.frame(
    REMARK = character(), DATE = character(),
    allan_sd_CO2 = numeric(), allan_sd_CH4 = numeric(),
    n_pts = integer(), t_sec = integer(), instrument = character(),
    stringsAsFactors = FALSE
  )
} else {
  message("LI-7810 Allan deviation: ", nrow(allan_df_7810), " measurements")
  message("  CH4 Allan SD: median ", round(median(allan_df_7810$allan_sd_CH4, na.rm = TRUE), 4),
          " ppb, range ", round(min(allan_df_7810$allan_sd_CH4, na.rm = TRUE), 4),
          " - ", round(max(allan_df_7810$allan_sd_CH4, na.rm = TRUE), 4))
}

# Build match key for 7810: REMARK_DATE
allan_df_7810$match_key <- paste(allan_df_7810$REMARK, allan_df_7810$DATE, sep = "_")

# ============================================================
# PART 5: ALLAN DEVIATION - 2023-24 LGR/UGGA
# ============================================================

message("\n=== Allan deviation: 2023-24 LGR/UGGA ===")

lgr_base <- normalizePath(
  file.path("..", "..", "..", "Matthes_Lab", "stem-CH4-flux", "Raw LGR Data"),
  mustWork = FALSE
)

# Check if the path exists, try alternate
if (!dir.exists(lgr_base)) {
  # Try absolute path
  lgr_base <- "/Users/jongewirtzman/My Drive/Matthes_Lab/stem-CH4-flux/Raw LGR Data"
}
message("LGR raw data directory: ", lgr_base)
message("  Exists: ", dir.exists(lgr_base))

# Catalog all available LGR date folders
lgr_catalog <- data.frame(machine = character(), date_str = character(),
                           folder = character(), stringsAsFactors = FALSE)
for (lgr_name in c("LGR1", "LGR2", "LGR3")) {
  lgr_dir <- file.path(lgr_base, lgr_name)
  if (!dir.exists(lgr_dir)) next
  date_folders <- list.dirs(lgr_dir, recursive = FALSE, full.names = TRUE)
  for (folder in date_folders) {
    folder_name <- basename(folder)
    # Folder names like "2023-06-29" or "2023-06-06 (1)"
    date_part <- substr(folder_name, 1, 10)
    lgr_catalog <- rbind(lgr_catalog, data.frame(
      machine = lgr_name, date_str = date_part, folder = folder,
      stringsAsFactors = FALSE
    ))
  }
}
message("LGR catalog: ", nrow(lgr_catalog), " date-folders across ",
        length(unique(lgr_catalog$machine)), " instruments")

# Parse one LGR raw file (txt or txt.zip) and return 1-Hz data frame
parse_lgr_file <- function(filepath) {
  tryCatch({
    if (grepl("\\.zip$", filepath)) {
      # Read from zip without extracting
      txt_name <- sub("\\.zip$", "", basename(filepath))
      con <- unz(filepath, txt_name)
      lines <- readLines(con, warn = FALSE)
      close(con)
    } else {
      lines <- readLines(filepath, warn = FALSE)
    }

    if (length(lines) < 3) return(NULL)

    # Skip line 1 (serial/metadata), read from line 2 (headers) onwards
    dat <- tryCatch(
      read.csv(text = paste(lines[-1], collapse = "\n"),
               stringsAsFactors = FALSE, strip.white = TRUE),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0) return(NULL)

    # Parse timestamp from first column (SysTime)
    # Format: "06/29/2023 10:06:11.087"
    ts_col <- names(dat)[1]  # usually "SysTime" or similar
    dat$timestamp <- as.POSIXct(dat[[ts_col]], format = "%m/%d/%Y %H:%M:%OS")

    # Get CH4 and CO2 dry mole fractions
    ch4_col <- grep("CH4.*d_ppm", names(dat), value = TRUE)[1]
    co2_col <- grep("CO2.*d_ppm", names(dat), value = TRUE)[1]
    if (is.na(ch4_col) || is.na(co2_col)) return(NULL)

    data.frame(
      timestamp = dat$timestamp,
      CH4_ppm   = as.numeric(dat[[ch4_col]]),
      CO2_ppm   = as.numeric(dat[[co2_col]]),
      stringsAsFactors = FALSE
    ) %>% filter(!is.na(timestamp))
  }, error = function(e) {
    # message("  Error parsing: ", filepath, " - ", e$message)
    NULL
  })
}

# Load and concatenate all LGR files for a given folder
load_lgr_day <- function(folder_path) {
  # Find all txt and txt.zip files with _f prefix (data files, not _b, _l, _p)
  all_files <- list.files(folder_path, full.names = TRUE)
  data_files <- all_files[grepl("_f\\d+\\.txt(\\.zip)?$", all_files)]

  if (length(data_files) == 0) return(NULL)

  # Parse each file and concatenate
  day_data <- lapply(data_files, parse_lgr_file)
  day_data <- do.call(rbind, Filter(Negate(is.null), day_data))

  if (is.null(day_data) || nrow(day_data) == 0) return(NULL)

  # Sort by timestamp
  day_data <- day_data[order(day_data$timestamp), ]
  day_data
}

# Process each field log entry: find LGR data, segment, compute Allan deviation
message("\nProcessing LGR measurements...")

# Group field log entries by date + machine for efficiency
# (load LGR file once per date, segment all measurements)
field_logs_lgr <- field_logs %>%
  filter(date_parsed < as.Date("2025-01-01"))

# Track unique dates to load
unique_day_machine <- field_logs_lgr %>%
  distinct(date_str, machine) %>%
  arrange(date_str)

message("  Unique LGR date-machine combos: ", nrow(unique_day_machine))

allan_lgr_all <- list()
n_matched <- 0
n_missing <- 0

for (i in seq_len(nrow(unique_day_machine))) {
  d <- unique_day_machine$date_str[i]
  m <- unique_day_machine$machine[i]

  # Find matching folder in catalog
  folders <- lgr_catalog$folder[lgr_catalog$machine == m &
                                  lgr_catalog$date_str == d]

  if (length(folders) == 0) {
    # Try other LGR instruments as fallback
    folders <- lgr_catalog$folder[lgr_catalog$date_str == d]
  }

  if (length(folders) == 0) {
    n_missing <- n_missing + 1
    next
  }

  # Load LGR data for this day (concatenate all folders if multiple)
  day_data <- NULL
  for (folder in folders) {
    dd <- load_lgr_day(folder)
    if (!is.null(dd)) {
      day_data <- rbind(day_data, dd)
    }
  }

  if (is.null(day_data) || nrow(day_data) < 10) {
    n_missing <- n_missing + 1
    next
  }

  # Get all field log entries for this date+machine
  entries <- field_logs_lgr %>%
    filter(date_str == d, machine == m)

  # Segment each measurement and compute Allan deviation
  for (j in seq_len(nrow(entries))) {
    start_t <- entries$comp_start_posix[j]
    end_t   <- entries$comp_end_posix[j]

    segment <- day_data[day_data$timestamp >= start_t &
                          day_data$timestamp <= end_t, ]

    if (nrow(segment) < 5) next

    n_matched <- n_matched + 1
    allan_lgr_all[[n_matched]] <- data.frame(
      UniqueID     = entries$UniqueID[j],
      date_str     = d,
      allan_sd_CO2 = allan_sd(segment$CO2_ppm),          # ppm (same as LI-7810)
      allan_sd_CH4 = allan_sd(segment$CH4_ppm) * 1000,  # convert ppm -> ppb to match LI-7810
      n_pts        = nrow(segment),
      t_sec        = nrow(segment),
      instrument   = m,
      stringsAsFactors = FALSE
    )
  }

  # Progress message every 10 dates
  if (i %% 10 == 0 || i == nrow(unique_day_machine)) {
    message("  Processed ", i, "/", nrow(unique_day_machine),
            " dates (", n_matched, " measurements matched so far)")
  }
}

allan_df_lgr <- do.call(rbind, allan_lgr_all)

if (is.null(allan_df_lgr) || nrow(allan_df_lgr) == 0) {
  warning("No LGR Allan results computed.")
  allan_df_lgr <- data.frame(
    UniqueID = character(), date_str = character(),
    allan_sd_CO2 = numeric(), allan_sd_CH4 = numeric(),
    n_pts = integer(), t_sec = integer(), instrument = character(),
    stringsAsFactors = FALSE
  )
} else {
  message("\nLGR Allan deviation: ", nrow(allan_df_lgr), " measurements")
  message("  CH4 Allan SD: median ", round(median(allan_df_lgr$allan_sd_CH4, na.rm = TRUE), 4),
          " ppb, range ", round(min(allan_df_lgr$allan_sd_CH4, na.rm = TRUE), 4),
          " - ", round(max(allan_df_lgr$allan_sd_CH4, na.rm = TRUE), 4))
  message("  Missing date-folders: ", n_missing)
}

# ============================================================
# PART 6: MERGE ALLAN RESULTS INTO MAIN DATA
# ============================================================

message("\n=== Merging Allan deviation into corrected flux data ===")

# --- Match 2025 LI-7810 by Tree_Date key ---
df$match_key_7810 <- paste(df$Tree, df$date, sep = "_")

allan_7810_unique <- allan_df_7810 %>%
  group_by(match_key) %>%
  slice(1) %>%
  ungroup() %>%
  select(match_key, allan_sd_CO2, allan_sd_CH4, n_pts, t_sec, instrument)

df <- merge(df, allan_7810_unique,
            by.x = "match_key_7810", by.y = "match_key",
            all.x = TRUE, suffixes = c("", "_7810"))

# --- Match 2023-24 LGR by UniqueID ---
if (nrow(allan_df_lgr) > 0) {
  allan_lgr_unique <- allan_df_lgr %>%
    group_by(UniqueID) %>%
    slice(1) %>%
    ungroup() %>%
    select(UniqueID, allan_sd_CO2, allan_sd_CH4, n_pts, t_sec, instrument)

  df <- merge(df, allan_lgr_unique,
              by = "UniqueID", all.x = TRUE,
              suffixes = c("", "_lgr"))

  # Combine: prefer existing (7810) values, fill in LGR where missing
  df <- df %>%
    mutate(
      allan_sd_CO2 = if_else(is.na(allan_sd_CO2), allan_sd_CO2_lgr, allan_sd_CO2),
      allan_sd_CH4 = if_else(is.na(allan_sd_CH4), allan_sd_CH4_lgr, allan_sd_CH4),
      n_pts        = if_else(is.na(n_pts), n_pts_lgr, n_pts),
      t_sec        = if_else(is.na(t_sec), t_sec_lgr, t_sec),
      instrument   = if_else(is.na(instrument) | instrument == "",
                              if_else(!is.na(instrument_lgr), instrument_lgr, NA_character_),
                              instrument)
    ) %>%
    select(-ends_with("_lgr"))
}

# --- For 2023 garbled UniqueIDs, try matching by datetime ---
# The corrected flux file datetime_posx = "Real start" from field logs
# Match 2023-24 unmatched rows by (date + tree + datetime proximity)
unmatched_idx <- which(is.na(df$allan_sd_CH4) & df$year < 2025)
message("Unmatched 2023-24 rows before datetime fallback: ", length(unmatched_idx))

if (length(unmatched_idx) > 0 && nrow(field_logs_lgr) > 0 && nrow(allan_df_lgr) > 0) {
  # Build lookup from field logs: real_start datetime + tree -> UniqueID
  fl_lookup <- field_logs %>%
    filter(!is.na(date_parsed)) %>%
    mutate(
      real_start_posix = as.POSIXct(paste(date_parsed, real_start),
                                     format = "%Y-%m-%d %H:%M:%S")
    ) %>%
    filter(!is.na(real_start_posix))

  # For each unmatched row, try to find a field log entry with matching
  # datetime_posx (within 2 minutes tolerance)
  for (idx in unmatched_idx) {
    dt_str <- df$datetime_posx[idx]
    if (is.na(dt_str) || dt_str == "") next

    # Parse datetime from corrected flux file
    dt <- as.POSIXct(dt_str, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")
    if (is.na(dt)) dt <- as.POSIXct(dt_str, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
    if (is.na(dt)) next

    # Compare with field log real_start times on same date
    same_date <- fl_lookup$date_str == df$date[idx]
    if (!any(same_date, na.rm = TRUE)) next

    fl_sub <- fl_lookup[same_date & !is.na(same_date), ]
    if (nrow(fl_sub) == 0) next

    # Find closest match within 2 minutes
    time_diffs <- abs(as.numeric(difftime(fl_sub$real_start_posix, dt, units = "secs")))
    best <- which.min(time_diffs)
    if (time_diffs[best] > 120) next  # no match within 2 min

    # Get UniqueID from field log and look up Allan result
    matched_uid <- fl_sub$UniqueID[best]
    allan_row <- allan_df_lgr[allan_df_lgr$UniqueID == matched_uid, ]
    if (nrow(allan_row) == 0) next

    df$allan_sd_CO2[idx] <- allan_row$allan_sd_CO2[1]
    df$allan_sd_CH4[idx] <- allan_row$allan_sd_CH4[1]
    df$n_pts[idx]        <- allan_row$n_pts[1]
    df$t_sec[idx]        <- allan_row$t_sec[1]
    df$instrument[idx]   <- allan_row$instrument[1]
  }
}

unmatched_final <- sum(is.na(df$allan_sd_CH4) & df$year < 2025)
message("Unmatched 2023-24 after datetime fallback: ", unmatched_final)

# Summary
message("\n=== Allan deviation coverage ===")
message("  2025 LI-7810: ", sum(!is.na(df$allan_sd_CH4) & df$year == 2025),
        " / ", sum(df$year == 2025))
message("  2023-24 LGR:  ", sum(!is.na(df$allan_sd_CH4) & df$year < 2025),
        " / ", sum(df$year < 2025))
message("  Total:        ", sum(!is.na(df$allan_sd_CH4)), " / ", n_total)

# ============================================================
# PART 7: COMPUTE FLUX TERM, NOISE, MDF THRESHOLDS
# ============================================================

message("\n=== Computing MDF thresholds ===")

df <- df %>%
  mutate(
    # flux_term = nmol / surfarea
    flux_term = nmol / surfarea,

    # Actual measurement duration; fall back to 300 s if unknown
    t_sec_est = if_else(!is.na(t_sec), as.numeric(t_sec), 300),

    # Empirical noise floor in flux units (Allan SD / t * flux_term)
    CH4_noise_floor = ifelse(!is.na(allan_sd_CH4) & !is.na(flux_term),
                              allan_sd_CH4 / t_sec_est * flux_term,
                              NA_real_),
    CO2_noise_floor = ifelse(!is.na(allan_sd_CO2) & !is.na(flux_term),
                              allan_sd_CO2 / t_sec_est * flux_term,
                              NA_real_),

    # Empirical SNR (Allan deviation based)
    CH4_snr_allan = ifelse(!is.na(CH4_noise_floor) & CH4_noise_floor > 0,
                            abs(CH4_flux_nmolpm2ps) / CH4_noise_floor,
                            NA_real_),

    # Manufacturer precision (instrument-specific)
    prec_ch4 = if_else(year == 2025, PREC_CH4_7810, PREC_CH4_LGR),
    prec_co2 = if_else(year == 2025, PREC_CO2_7810, PREC_CO2_LGR),

    # Manufacturer MDF = precision / t * flux_term
    CH4_MDF_manufacturer = ifelse(!is.na(flux_term),
                                   prec_ch4 / t_sec_est * flux_term,
                                   NA_real_),

    # Christiansen MDF = allan_sd * 3 * t_crit / t * flux_term
    df_meas = pmax(t_sec_est - 2, 1),
    t99 = qt(0.995, df = df_meas),
    t95 = qt(0.975, df = df_meas),
    t90 = qt(0.95,  df = df_meas),

    CH4_MDF_chr99 = ifelse(!is.na(allan_sd_CH4) & !is.na(flux_term),
                            (allan_sd_CH4 * 3 * t99) / t_sec_est * flux_term,
                            NA_real_),
    CH4_MDF_chr95 = ifelse(!is.na(allan_sd_CH4) & !is.na(flux_term),
                            (allan_sd_CH4 * 3 * t95) / t_sec_est * flux_term,
                            NA_real_),
    CH4_MDF_chr90 = ifelse(!is.na(allan_sd_CH4) & !is.na(flux_term),
                            (allan_sd_CH4 * 3 * t90) / t_sec_est * flux_term,
                            NA_real_),

    # Below-MDF flags
    CH4_below_MDF_manuf = ifelse(!is.na(CH4_MDF_manufacturer),
                                   abs(CH4_flux_nmolpm2ps) < CH4_MDF_manufacturer,
                                   NA),
    CH4_below_MDF_chr99 = ifelse(!is.na(CH4_MDF_chr99),
                                   abs(CH4_flux_nmolpm2ps) < CH4_MDF_chr99,
                                   NA),
    CH4_below_MDF_chr95 = ifelse(!is.na(CH4_MDF_chr95),
                                   abs(CH4_flux_nmolpm2ps) < CH4_MDF_chr95,
                                   NA),
    CH4_below_MDF_chr90 = ifelse(!is.na(CH4_MDF_chr90),
                                   abs(CH4_flux_nmolpm2ps) < CH4_MDF_chr90,
                                   NA)
  )

# --- Wassmann MDF: use INSTRUMENT-SPECIFIC global Allan SD ---
global_sd_ch4_lgr  <- if (nrow(allan_df_lgr) > 0)
  median(allan_df_lgr$allan_sd_CH4, na.rm = TRUE) else NA_real_
global_sd_ch4_7810 <- if (nrow(allan_df_7810) > 0)
  median(allan_df_7810$allan_sd_CH4, na.rm = TRUE) else NA_real_

message("Global CH4 precision (median Allan SD):")
message("  LGR:    ", ifelse(is.na(global_sd_ch4_lgr), "N/A",
                              round(global_sd_ch4_lgr, 4)), " ppb")
message("  LI-7810: ", ifelse(is.na(global_sd_ch4_7810), "N/A",
                               round(global_sd_ch4_7810, 4)), " ppb")

z99 <- qnorm(0.995); z95 <- qnorm(0.975); z90 <- qnorm(0.95)

df <- df %>%
  mutate(
    # Instrument-specific global precision for Wassmann
    global_sd_ch4 = if_else(year == 2025, global_sd_ch4_7810, global_sd_ch4_lgr),

    CH4_MDF_wass99 = ifelse(!is.na(global_sd_ch4) & !is.na(flux_term),
                             (z99 * global_sd_ch4) / t_sec_est * flux_term,
                             NA_real_),
    CH4_MDF_wass95 = ifelse(!is.na(global_sd_ch4) & !is.na(flux_term),
                             (z95 * global_sd_ch4) / t_sec_est * flux_term,
                             NA_real_),
    CH4_MDF_wass90 = ifelse(!is.na(global_sd_ch4) & !is.na(flux_term),
                             (z90 * global_sd_ch4) / t_sec_est * flux_term,
                             NA_real_),

    CH4_below_MDF_wass99 = ifelse(!is.na(CH4_MDF_wass99),
                                    abs(CH4_flux_nmolpm2ps) < CH4_MDF_wass99,
                                    NA),
    CH4_below_MDF_wass95 = ifelse(!is.na(CH4_MDF_wass95),
                                    abs(CH4_flux_nmolpm2ps) < CH4_MDF_wass95,
                                    NA),
    CH4_below_MDF_wass90 = ifelse(!is.na(CH4_MDF_wass90),
                                    abs(CH4_flux_nmolpm2ps) < CH4_MDF_wass90,
                                    NA)
  )

# ============================================================
# MDF SUMMARY
# ============================================================

message("\n=== MDF Summary ===")
# Helper to report filter stats accounting for NAs
report_mdf <- function(flag_col, label) {
  vals <- df[[flag_col]]
  n_eval <- sum(!is.na(vals))
  n_below <- sum(vals, na.rm = TRUE)
  pct <- if (n_eval > 0) round(100 * n_below / n_eval, 1) else NA
  message(sprintf("  %-25s %d/%d evaluated, %d below (%.1f%%)",
                  label, n_eval, n_total, n_below,
                  ifelse(is.na(pct), 0, pct)))
}

report_mdf("CH4_below_MDF_manuf", "Manufacturer MDF")
report_mdf("CH4_below_MDF_wass90", "Wassmann 90%")
report_mdf("CH4_below_MDF_wass95", "Wassmann 95%")
report_mdf("CH4_below_MDF_wass99", "Wassmann 99%")
report_mdf("CH4_below_MDF_chr90", "Christiansen 90%")
report_mdf("CH4_below_MDF_chr95", "Christiansen 95%")
report_mdf("CH4_below_MDF_chr99", "Christiansen 99%")
} # end if (FALSE) — dead code from old Parts 1-7

# ============================================================
# DEFINE ALL QUALITY FILTERS
# ============================================================

r2_sym <- "\u00B2"

quality_filters <- list()
quality_filters[["No filter"]] <- function(d) rep(FALSE, nrow(d))

# Manufacturer MDF (instrument-specific)
quality_filters[["Manufacturer MDF"]] <- function(d) {
  ifelse(!is.na(d$CH4_below_MDF_manuf), d$CH4_below_MDF_manuf, FALSE)
}

# CO2-based filters
quality_filters[["CO2 flux > 0"]]                    <- function(d) d$CO2_flux_umolpm2ps <= 0
quality_filters[["CO2 SNR (SE) > 2"]]                <- function(d) d$CO2_snr_se <= 2
quality_filters[[paste0("CO2 R", r2_sym, " > 0.7")]] <- function(d) d$CO2_r2 <= 0.7
quality_filters[[paste0("CO2 R", r2_sym, " > 0.8")]] <- function(d) d$CO2_r2 <= 0.8

# CH4 SE-based SNR
quality_filters[["CH4 SNR (SE) > 2"]] <- function(d) d$CH4_snr_se <= 2
quality_filters[["CH4 SNR (SE) > 3"]] <- function(d) d$CH4_snr_se <= 3

# CH4 Allan deviation SNR
# For measurements without Allan deviation, they pass by default (not excluded)
quality_filters[["CH4 SNR (Allan) > 2"]] <- function(d) {
  ifelse(!is.na(d$CH4_snr_allan), d$CH4_snr_allan <= 2, FALSE)
}
quality_filters[["CH4 SNR (Allan) > 3"]] <- function(d) {
  ifelse(!is.na(d$CH4_snr_allan), d$CH4_snr_allan <= 3, FALSE)
}

# Wassmann MDF thresholds (instrument-specific global precision)
quality_filters[["Wassmann 90%"]] <- function(d) {
  ifelse(!is.na(d$CH4_below_MDF_wass90), d$CH4_below_MDF_wass90, FALSE)
}
quality_filters[["Wassmann 95%"]] <- function(d) {
  ifelse(!is.na(d$CH4_below_MDF_wass95), d$CH4_below_MDF_wass95, FALSE)
}
quality_filters[["Wassmann 99%"]] <- function(d) {
  ifelse(!is.na(d$CH4_below_MDF_wass99), d$CH4_below_MDF_wass99, FALSE)
}

# Christiansen MDF (requires per-measurement Allan deviation)
quality_filters[["Christiansen 90%"]] <- function(d) {
  ifelse(!is.na(d$CH4_below_MDF_chr90), d$CH4_below_MDF_chr90, FALSE)
}
quality_filters[["Christiansen 95%"]] <- function(d) {
  ifelse(!is.na(d$CH4_below_MDF_chr95), d$CH4_below_MDF_chr95, FALSE)
}
quality_filters[["Christiansen 99%"]] <- function(d) {
  ifelse(!is.na(d$CH4_below_MDF_chr99), d$CH4_below_MDF_chr99, FALSE)
}

# CH4 R^2 thresholds
quality_filters[[paste0("CH4 R", r2_sym, " > 0.5")]] <- function(d) d$CH4_r2 <= 0.5
quality_filters[[paste0("CH4 R", r2_sym, " > 0.7")]] <- function(d) d$CH4_r2 <= 0.7
quality_filters[[paste0("CH4 R", r2_sym, " > 0.9")]] <- function(d) d$CH4_r2 <= 0.9

filter_names <- names(quality_filters)

# ============================================================
# SORT FILTERS BY STRINGENCY
# ============================================================

n_retained <- sapply(filter_names, function(fn) {
  fails <- quality_filters[[fn]](df)
  sum(!fails)
})
filter_order <- names(sort(n_retained, decreasing = TRUE))

# ============================================================
# BUILD LONG DATA
# ============================================================

# Count how many measurements have Allan deviation for annotation
n_with_allan <- sum(!is.na(df$allan_sd_CH4))

pass_labels <- sapply(filter_order, function(fn) {
  fails <- quality_filters[[fn]](df)
  n_pass <- sum(!fails)
  pct <- round(100 * n_pass / n_total, 1)
  n_neg <- sum(df$CH4_flux_nmolpm2ps[!fails] < 0)
  pct_neg <- round(100 * n_neg / n_pass, 1)
  paste0(fn, "\n(", n_pass, "/", n_total, ", ", pct, "%",
         " | ", pct_neg, "% neg)")
})

rows_list <- lapply(seq_along(filter_order), function(i) {
  fn <- filter_order[i]
  fails <- quality_filters[[fn]](df)
  data.frame(
    filter   = pass_labels[i],
    CH4_flux = df$CH4_flux_nmolpm2ps,
    passes   = !fails,
    stringsAsFactors = FALSE
  )
})

ridge_df <- do.call(rbind, rows_list)
ridge_df$filter <- factor(ridge_df$filter, levels = rev(pass_labels))

# ============================================================
# PER-FILTER SUMMARY STATS
# ============================================================

filter_stats <- ridge_df %>%
  filter(passes) %>%
  group_by(filter) %>%
  summarise(
    mean_flux   = mean(CH4_flux, na.rm = TRUE),
    median_flux = median(CH4_flux, na.rm = TRUE),
    n_neg       = sum(CH4_flux < 0),
    n_total     = n(),
    pct_neg     = round(100 * n_neg / n_total, 1),
    .groups     = "drop"
  )

cat("\n=== Filter summary (CH4 fluxes) ===\n")
print(as.data.frame(filter_stats), row.names = FALSE)

# ============================================================
# RETENTION % FOR FILL COLOR
# ============================================================

ridge_df <- ridge_df %>%
  group_by(filter) %>%
  mutate(pct_retained = 100 * sum(passes) / n()) %>%
  ungroup()

# ============================================================
# PART 9: GGRIDGES RAINFALL PLOT
# ============================================================

asinh_trans <- trans_new(
  name = "asinh",
  transform = asinh,
  inverse = sinh,
  breaks = function(x) {
    rng <- sinh(x)
    pretty(rng, n = 8)
  }
)

p_ridges <- ggplot(ridge_df, aes(x = CH4_flux, y = filter)) +
  geom_density_ridges(
    aes(fill = pct_retained, point_color = ifelse(passes, "Retained", "Excluded")),
    jittered_points = TRUE,
    point_size = 0.6, point_alpha = 0.4,
    scale = 0.9, alpha = 0.7,
    position = position_raincloud(width = 0.05, ygap = 0.05),
    bandwidth = 0.3,
    from = asinh(-5), to = asinh(30)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_point(data = filter_stats, aes(x = mean_flux, y = filter),
             shape = "|", size = 4, color = "red", inherit.aes = FALSE) +
  scale_x_continuous(trans = asinh_trans,
                     breaks = c(-5, -2, -1, 0, 0.5, 1, 2, 5, 10, 20),
                     labels = c("-5", "-2", "-1", "0", "0.5", "1", "2",
                                "5", "10", "20")) +
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                       midpoint = 70, limits = c(20, 100),
                       name = "% retained") +
  scale_discrete_manual("point_color",
                        values = c("Retained" = "grey40", "Excluded" = "red"),
                        name = "Measurement") +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = NULL,
    title = expression(CH[4]~flux~distribution~under~quality~filters),
    subtitle = paste0("n = ", n_total,
                      " measurements (", sum(df$year < 2025), " LGR 2023-24 + ",
                      sum(df$year == 2025), " LI-7810 2025)",
                      "\nAllan deviation coverage: ", n_with_allan, "/", n_total,
                      " | Labels: (n retained, % | % negative)")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7, lineheight = 0.85),
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  coord_cartesian(xlim = c(asinh(-5), asinh(25)))

# ============================================================
# SAVE
# ============================================================

out_dir <- file.path("outputs", "figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ggsave(file.path(out_dir, "ch4_flux_ridges_by_filter.png"),
       p_ridges, width = 13, height = 12, dpi = 300)
ggsave(file.path(out_dir, "ch4_flux_ridges_by_filter.pdf"),
       p_ridges, width = 13, height = 12)

message("\nSaved ridges plot to: ", out_dir)

# ============================================================
# SUMMARY TABLE
# ============================================================

cat("\n=== % Negative CH4 fluxes by filter ===\n")
neg_summary <- data.frame(
  filter = filter_order,
  stringsAsFactors = FALSE
) %>%
  rowwise() %>%
  mutate(
    n_pass = {
      fails <- quality_filters[[filter]](df)
      sum(!fails)
    },
    n_neg = {
      fails <- quality_filters[[filter]](df)
      sum(df$CH4_flux_nmolpm2ps[!fails] < 0)
    },
    pct_neg = round(100 * n_neg / n_pass, 1)
  ) %>%
  ungroup()

print(as.data.frame(neg_summary), row.names = FALSE)

# ============================================================
# ALLAN DEVIATION COMPARISON
# ============================================================

cat("\n=== Allan deviation comparison by instrument ===\n")
if (nrow(allan_df_lgr) > 0) {
  cat("LGR/UGGA (2023-24):\n")
  cat("  CH4: median =", round(median(allan_df_lgr$allan_sd_CH4, na.rm = TRUE), 4),
      "ppb, IQR =", round(quantile(allan_df_lgr$allan_sd_CH4, 0.25, na.rm = TRUE), 4),
      "-", round(quantile(allan_df_lgr$allan_sd_CH4, 0.75, na.rm = TRUE), 4), "\n")
  cat("  CO2: median =", round(median(allan_df_lgr$allan_sd_CO2, na.rm = TRUE), 3),
      "ppm, IQR =", round(quantile(allan_df_lgr$allan_sd_CO2, 0.25, na.rm = TRUE), 3),
      "-", round(quantile(allan_df_lgr$allan_sd_CO2, 0.75, na.rm = TRUE), 3), "\n")
  cat("  Manufacturer spec: CH4 =", PREC_CH4_LGR, "ppb, CO2 =", PREC_CO2_LGR, "ppm\n")
}
if (nrow(allan_df_7810) > 0) {
  cat("LI-7810 (2025):\n")
  cat("  CH4: median =", round(median(allan_df_7810$allan_sd_CH4, na.rm = TRUE), 4),
      "ppb, IQR =", round(quantile(allan_df_7810$allan_sd_CH4, 0.25, na.rm = TRUE), 4),
      "-", round(quantile(allan_df_7810$allan_sd_CH4, 0.75, na.rm = TRUE), 4), "\n")
  cat("  CO2: median =", round(median(allan_df_7810$allan_sd_CO2, na.rm = TRUE), 3),
      "ppm, IQR =", round(quantile(allan_df_7810$allan_sd_CO2, 0.25, na.rm = TRUE), 3),
      "-", round(quantile(allan_df_7810$allan_sd_CO2, 0.75, na.rm = TRUE), 3), "\n")
  cat("  Manufacturer spec: CH4 =", PREC_CH4_7810, "ppb, CO2 =", PREC_CO2_7810, "ppm\n")
}

# ============================================================
# PART 10: INSTRUMENT COMPARISON
# ============================================================

message("\n=== Instrument comparison ===")

# Label each measurement by instrument
df <- df %>%
  mutate(
    inst_label = case_when(
      year == 2025 ~ "LI-7810 (2025)",
      year < 2025  ~ "LGR/UGGA (2023-24)"
    )
  )

# --- Negative flux breakdown by instrument ---
cat("\n=== Negative CH4 fluxes by instrument ===\n")
inst_neg <- df %>%
  group_by(inst_label) %>%
  summarise(
    n_total   = n(),
    n_neg     = sum(CH4_flux_nmolpm2ps < 0),
    pct_neg   = round(100 * n_neg / n_total, 1),
    mean_flux = round(mean(CH4_flux_nmolpm2ps), 4),
    median_flux = round(median(CH4_flux_nmolpm2ps), 4),
    min_flux  = round(min(CH4_flux_nmolpm2ps), 4),
    max_flux  = round(max(CH4_flux_nmolpm2ps), 4),
    .groups   = "drop"
  )
print(as.data.frame(inst_neg), row.names = FALSE)

# --- Per-instrument filter sensitivity ---
cat("\n=== % Negative by filter, split by instrument ===\n")
inst_filter_summary <- lapply(filter_order, function(fn) {
  fails <- quality_filters[[fn]](df)
  df_pass <- df[!fails, ]
  df_pass %>%
    group_by(inst_label) %>%
    summarise(
      filter  = fn,
      n_pass  = n(),
      n_neg   = sum(CH4_flux_nmolpm2ps < 0),
      pct_neg = round(100 * n_neg / n_pass, 1),
      .groups = "drop"
    )
})
inst_filter_df <- do.call(rbind, inst_filter_summary)
# Pivot wide for readability
inst_wide <- inst_filter_df %>%
  select(filter, inst_label, n_pass, pct_neg) %>%
  pivot_wider(
    names_from  = inst_label,
    values_from = c(n_pass, pct_neg),
    names_glue  = "{inst_label}_{.value}"
  )
print(as.data.frame(inst_wide), row.names = FALSE)

# --- PLOT: Filter sensitivity by instrument (faceted ridges) ---

# Build long data with instrument facet
rows_inst <- lapply(seq_along(filter_order), function(i) {
  fn <- filter_order[i]
  fails <- quality_filters[[fn]](df)
  data.frame(
    filter     = fn,
    CH4_flux   = df$CH4_flux_nmolpm2ps,
    passes     = !fails,
    instrument = df$inst_label,
    stringsAsFactors = FALSE
  )
})
ridge_inst <- do.call(rbind, rows_inst)

# Compute per-filter per-instrument labels
ridge_inst <- ridge_inst %>%
  group_by(filter, instrument) %>%
  mutate(
    n_pass = sum(passes),
    n_inst = n(),
    pct_neg = round(100 * sum(CH4_flux[passes] < 0) / sum(passes), 1)
  ) %>%
  ungroup()

# Simplify filter labels for faceted plot
ridge_inst$filter_label <- sapply(ridge_inst$filter, function(fn) {
  paste0(fn, " (", ridge_inst$n_pass[ridge_inst$filter == fn][1], ")")
})
# Actually, make clean per-instrument per-filter labels
inst_stats <- ridge_inst %>%
  filter(passes) %>%
  group_by(filter, instrument) %>%
  summarise(
    n_pass  = n(),
    n_neg   = sum(CH4_flux < 0),
    pct_neg = round(100 * n_neg / n_pass, 1),
    mean_flux = mean(CH4_flux),
    .groups = "drop"
  )

# Order filters by total retention
filter_order_fct <- factor(ridge_inst$filter,
                            levels = rev(filter_order))
ridge_inst$filter <- filter_order_fct

p_inst_ridges <- ggplot(
  ridge_inst %>% filter(passes),
  aes(x = CH4_flux, y = filter, fill = instrument)
) +
  geom_density_ridges(
    alpha = 0.5, scale = 0.9,
    bandwidth = 0.3,
    from = asinh(-5), to = asinh(30)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red",
             linewidth = 0.5) +
  scale_x_continuous(
    trans = asinh_trans,
    breaks = c(-5, -2, -1, 0, 0.5, 1, 2, 5, 10, 20),
    labels = c("-5", "-2", "-1", "0", "0.5", "1", "2", "5", "10", "20")
  ) +
  scale_fill_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    name = "Instrument"
  ) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = NULL,
    title = expression(CH[4]~flux~distribution~by~instrument~and~quality~filter),
    subtitle = paste0("Overlaid densities: LGR/UGGA (n=", sum(df$year < 2025),
                      ") vs LI-7810 (n=", sum(df$year == 2025), ")")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  coord_cartesian(xlim = c(asinh(-5), asinh(25)))

ggsave(file.path(out_dir, "ch4_flux_ridges_by_instrument.png"),
       p_inst_ridges, width = 13, height = 12, dpi = 300)
ggsave(file.path(out_dir, "ch4_flux_ridges_by_instrument.pdf"),
       p_inst_ridges, width = 13, height = 12)
message("Saved instrument comparison ridges to: ", out_dir)

# --- PLOT: Allan deviation by instrument (violin + boxplot) ---

allan_combined <- bind_rows(
  allan_df_lgr %>%
    transmute(instrument = "LGR/UGGA (2023-24)",
              CH4_allan_sd = allan_sd_CH4,
              CO2_allan_sd = allan_sd_CO2),
  allan_df_7810 %>%
    transmute(instrument = "LI-7810 (2025)",
              CH4_allan_sd = allan_sd_CH4,
              CO2_allan_sd = allan_sd_CO2)
)

p_allan <- ggplot(allan_combined, aes(x = instrument, y = CH4_allan_sd,
                                       fill = instrument)) +
  geom_violin(alpha = 0.4, draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(width = 0.15, alpha = 0.15, size = 0.8) +
  geom_hline(yintercept = PREC_CH4_LGR,  linetype = "dashed", color = "#E69F00") +
  geom_hline(yintercept = PREC_CH4_7810, linetype = "dashed", color = "#56B4E9") +
  annotate("text", x = 2.4, y = PREC_CH4_LGR,
           label = paste0("LGR spec: ", PREC_CH4_LGR, " ppb"),
           hjust = 1, size = 3, color = "#E69F00") +
  annotate("text", x = 2.4, y = PREC_CH4_7810,
           label = paste0("LI-7810 spec: ", PREC_CH4_7810, " ppb"),
           hjust = 1, size = 3, color = "#56B4E9") +
  scale_y_log10() +
  scale_fill_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    guide = "none"
  ) +
  labs(
    x = NULL,
    y = expression(CH[4]~Allan~deviation~(ppb)),
    title = expression(Per-measurement~CH[4]~Allan~deviation~by~instrument),
    subtitle = "Dashed lines = manufacturer specs | Violin shows full distribution"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(out_dir, "ch4_allan_deviation_by_instrument.png"),
       p_allan, width = 7, height = 6, dpi = 300)
ggsave(file.path(out_dir, "ch4_allan_deviation_by_instrument.pdf"),
       p_allan, width = 7, height = 6)
message("Saved Allan deviation comparison to: ", out_dir)

# --- PLOT: Flux distributions by instrument (simple density comparison) ---

p_flux_inst <- ggplot(df, aes(x = CH4_flux_nmolpm2ps, fill = inst_label,
                                color = inst_label)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_rug(aes(color = inst_label), alpha = 0.1, sides = "b") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  scale_x_continuous(
    trans = asinh_trans,
    breaks = c(-5, -2, -1, 0, 0.5, 1, 2, 5, 10, 20),
    labels = c("-5", "-2", "-1", "0", "0.5", "1", "2", "5", "10", "20")
  ) +
  scale_fill_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    name = "Instrument"
  ) +
  scale_color_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    name = "Instrument"
  ) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = "Density",
    title = expression(CH[4]~flux~distribution~by~instrument),
    subtitle = paste0(
      "LGR/UGGA: n=", sum(df$year < 2025),
      " (", round(100 * mean(df$CH4_flux_nmolpm2ps[df$year < 2025] < 0), 1), "% neg)",
      " | LI-7810: n=", sum(df$year == 2025),
      " (", round(100 * mean(df$CH4_flux_nmolpm2ps[df$year == 2025] < 0), 1), "% neg)"
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, color = "grey40")
  ) +
  coord_cartesian(xlim = c(asinh(-5), asinh(25)))

ggsave(file.path(out_dir, "ch4_flux_density_by_instrument.png"),
       p_flux_inst, width = 9, height = 5, dpi = 300)
ggsave(file.path(out_dir, "ch4_flux_density_by_instrument.pdf"),
       p_flux_inst, width = 9, height = 5)
message("Saved flux density comparison to: ", out_dir)

# ============================================================
# PART 11: SEASONAL / MONTHLY BREAKDOWN
# ============================================================
#
# LGR covered dormant + growing season (Jun 2023–Dec 2024)
# LI-7810 covered growing season only (Apr–Oct 2025)
# Need to disentangle: which negative fluxes are real biology
# (uptake in dormant season) vs noise artifacts from noisier LGR?
# ============================================================

message("\n=== Seasonal / monthly breakdown ===")

df <- df %>%
  mutate(
    date_parsed = as.Date(date, format = "%Y-%m-%d"),
    month       = as.integer(format(date_parsed, "%m")),
    month_name  = factor(format(date_parsed, "%b"),
                          levels = c("Jan","Feb","Mar","Apr","May","Jun",
                                     "Jul","Aug","Sep","Oct","Nov","Dec")),
    year_month  = format(date_parsed, "%Y-%m"),
    # Season labels
    season = case_when(
      month %in% c(12, 1, 2) ~ "Winter (Dec-Feb)",
      month %in% c(3, 4, 5)  ~ "Spring (Mar-May)",
      month %in% c(6, 7, 8)  ~ "Summer (Jun-Aug)",
      month %in% c(9, 10, 11) ~ "Fall (Sep-Nov)"
    ),
    season = factor(season, levels = c("Spring (Mar-May)", "Summer (Jun-Aug)",
                                        "Fall (Sep-Nov)", "Winter (Dec-Feb)"))
  )

# --- TABLE 1: Monthly breakdown by instrument ---
cat("\n=== Monthly negative flux breakdown by instrument ===\n")
monthly_inst <- df %>%
  group_by(month_name, inst_label) %>%
  summarise(
    n        = n(),
    n_neg    = sum(CH4_flux_nmolpm2ps < 0),
    pct_neg  = round(100 * n_neg / n, 1),
    mean_flux  = round(mean(CH4_flux_nmolpm2ps), 4),
    median_flux = round(median(CH4_flux_nmolpm2ps), 4),
    .groups  = "drop"
  ) %>%
  arrange(month_name, inst_label)
print(as.data.frame(monthly_inst), row.names = FALSE)

# --- TABLE 2: Seasonal summary ---
cat("\n=== Seasonal summary by instrument ===\n")
seasonal_inst <- df %>%
  group_by(season, inst_label) %>%
  summarise(
    n          = n(),
    n_neg      = sum(CH4_flux_nmolpm2ps < 0),
    pct_neg    = round(100 * n_neg / n, 1),
    mean_flux  = round(mean(CH4_flux_nmolpm2ps), 4),
    median_flux = round(median(CH4_flux_nmolpm2ps), 4),
    .groups    = "drop"
  )
print(as.data.frame(seasonal_inst), row.names = FALSE)

# --- TABLE 3: Overlapping months only (Jun–Oct, when both instruments measured) ---
# This is the fairest apples-to-apples comparison
overlap_months <- c(6, 7, 8, 9, 10)  # months present in both instrument periods
cat("\n=== Overlapping months (Jun-Oct) — apples-to-apples comparison ===\n")
df_overlap <- df %>% filter(month %in% overlap_months)

overlap_summary <- df_overlap %>%
  group_by(inst_label) %>%
  summarise(
    n          = n(),
    n_neg      = sum(CH4_flux_nmolpm2ps < 0),
    pct_neg    = round(100 * n_neg / n, 1),
    mean_flux  = round(mean(CH4_flux_nmolpm2ps), 4),
    median_flux = round(median(CH4_flux_nmolpm2ps), 4),
    sd_flux    = round(sd(CH4_flux_nmolpm2ps), 4),
    min_flux   = round(min(CH4_flux_nmolpm2ps), 4),
    max_flux   = round(max(CH4_flux_nmolpm2ps), 4),
    .groups    = "drop"
  )
print(as.data.frame(overlap_summary), row.names = FALSE)

# How many LGR negatives are in dormant vs growing months?
cat("\n=== LGR negative fluxes: dormant vs growing season ===\n")
lgr_only <- df %>% filter(year < 2025)
dormant_months  <- c(11, 12, 1, 2, 3, 4)
growing_months  <- c(5, 6, 7, 8, 9, 10)
lgr_season <- lgr_only %>%
  mutate(period = if_else(month %in% dormant_months,
                           "Dormant (Nov-Apr)", "Growing (May-Oct)")) %>%
  group_by(period) %>%
  summarise(
    n          = n(),
    n_neg      = sum(CH4_flux_nmolpm2ps < 0),
    pct_neg    = round(100 * n_neg / n, 1),
    mean_flux  = round(mean(CH4_flux_nmolpm2ps), 4),
    median_flux = round(median(CH4_flux_nmolpm2ps), 4),
    .groups    = "drop"
  )
print(as.data.frame(lgr_season), row.names = FALSE)

# --- TABLE 4: Negative fluxes vs noise floor ---
# For measurements with Allan deviation: are negatives within the noise?
cat("\n=== Negative flux magnitude vs noise floor ===\n")
df_neg <- df %>%
  filter(CH4_flux_nmolpm2ps < 0, !is.na(CH4_noise_floor))
cat("Negative fluxes with Allan-derived noise floor:", nrow(df_neg), "\n")
if (nrow(df_neg) > 0) {
  df_neg <- df_neg %>%
    mutate(
      abs_flux    = abs(CH4_flux_nmolpm2ps),
      within_1sd  = abs_flux < CH4_noise_floor,
      within_2sd  = abs_flux < 2 * CH4_noise_floor,
      within_3sd  = abs_flux < 3 * CH4_noise_floor
    )
  noise_check <- df_neg %>%
    group_by(inst_label) %>%
    summarise(
      n_neg        = n(),
      within_1_noise = sum(within_1sd),
      within_2_noise = sum(within_2sd),
      within_3_noise = sum(within_3sd),
      pct_within_1   = round(100 * within_1_noise / n_neg, 1),
      pct_within_2   = round(100 * within_2_noise / n_neg, 1),
      pct_within_3   = round(100 * within_3_noise / n_neg, 1),
      median_abs_flux = round(median(abs_flux), 4),
      median_noise    = round(median(CH4_noise_floor), 4),
      .groups = "drop"
    )
  print(as.data.frame(noise_check), row.names = FALSE)
}

# --- PLOT 1: Monthly flux distributions by instrument (faceted) ---

p_monthly <- ggplot(df, aes(x = CH4_flux_nmolpm2ps, fill = inst_label)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 50, alpha = 0.6, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
  facet_wrap(~ month_name, ncol = 4, scales = "free_y") +
  scale_x_continuous(
    trans = asinh_trans,
    breaks = c(-2, -1, 0, 0.5, 1, 2, 5, 10),
    labels = c("-2", "-1", "0", "0.5", "1", "2", "5", "10")
  ) +
  scale_fill_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    name = "Instrument"
  ) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = "Density",
    title = expression(Monthly~CH[4]~flux~distributions~by~instrument),
    subtitle = "LGR covers dormant + growing season; LI-7810 covers growing season only (Apr-Oct 2025)"
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    strip.text = element_text(size = 9, face = "bold")
  ) +
  coord_cartesian(xlim = c(asinh(-3), asinh(15)))

ggsave(file.path(out_dir, "ch4_flux_monthly_by_instrument.png"),
       p_monthly, width = 12, height = 8, dpi = 300)
ggsave(file.path(out_dir, "ch4_flux_monthly_by_instrument.pdf"),
       p_monthly, width = 12, height = 8)
message("Saved monthly breakdown to: ", out_dir)

# --- PLOT 2: Negative flux fraction by month + instrument ---

monthly_neg_rate <- df %>%
  group_by(month_name, inst_label) %>%
  summarise(
    n       = n(),
    pct_neg = 100 * mean(CH4_flux_nmolpm2ps < 0),
    .groups = "drop"
  )

p_neg_month <- ggplot(monthly_neg_rate,
                       aes(x = month_name, y = pct_neg,
                           fill = inst_label, group = inst_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.8) +
  geom_text(aes(label = paste0(n)),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 2.5, color = "grey40") +
  scale_fill_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    name = "Instrument"
  ) +
  labs(
    x = "Month",
    y = "% negative fluxes",
    title = expression(Fraction~of~negative~CH[4]~fluxes~by~month~and~instrument),
    subtitle = "Numbers above bars = sample size"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  ) +
  coord_cartesian(ylim = c(0, 55))

ggsave(file.path(out_dir, "ch4_pct_negative_by_month_instrument.png"),
       p_neg_month, width = 9, height = 5, dpi = 300)
ggsave(file.path(out_dir, "ch4_pct_negative_by_month_instrument.pdf"),
       p_neg_month, width = 9, height = 5)
message("Saved negative % by month to: ", out_dir)

# --- PLOT 3: Flux vs noise floor scatter, colored by sign ---
# Shows whether negative fluxes sit within the instrument noise

if (sum(!is.na(df$CH4_noise_floor)) > 50) {
  p_noise_scatter <- ggplot(
    df %>% filter(!is.na(CH4_noise_floor)),
    aes(x = CH4_noise_floor, y = CH4_flux_nmolpm2ps, color = inst_label)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
    geom_abline(slope = -1, intercept = 0, linetype = "dotted", color = "grey50") +
    geom_point(alpha = 0.25, size = 1) +
    scale_x_log10(labels = scales::label_number()) +
    scale_y_continuous(
      trans = asinh_trans,
      breaks = c(-5, -2, -1, 0, 0.5, 1, 2, 5, 10, 20),
      labels = c("-5", "-2", "-1", "0", "0.5", "1", "2", "5", "10", "20")
    ) +
    scale_color_manual(
      values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
      name = "Instrument"
    ) +
    annotate("text", x = 0.3, y = 0.4, label = "flux = +noise floor",
             angle = 18, size = 3, color = "grey50") +
    annotate("text", x = 0.3, y = -0.4, label = "flux = \u2013noise floor",
             angle = -18, size = 3, color = "grey50") +
    labs(
      x = expression(CH[4]~noise~floor~(nmol~m^{-2}~s^{-1})),
      y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
      title = expression(CH[4]~flux~vs.~per-measurement~noise~floor),
      subtitle = "Dotted lines = \u00b1 noise floor (1\u03c3 Allan). Points within dotted wedge are indistinguishable from zero."
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 9, color = "grey40")
    )

  ggsave(file.path(out_dir, "ch4_flux_vs_noise_floor.png"),
         p_noise_scatter, width = 8, height = 6, dpi = 300)
  ggsave(file.path(out_dir, "ch4_flux_vs_noise_floor.pdf"),
         p_noise_scatter, width = 8, height = 6)
  message("Saved flux vs noise floor scatter to: ", out_dir)
}

# --- PLOT 4: Monthly time series with noise envelope ---
# Shows median flux ± noise floor by month and instrument

monthly_ts <- df %>%
  filter(!is.na(CH4_noise_floor)) %>%
  group_by(year_month, inst_label) %>%
  summarise(
    date_mid    = mean(date_parsed),
    median_flux = median(CH4_flux_nmolpm2ps),
    q25_flux    = quantile(CH4_flux_nmolpm2ps, 0.25),
    q75_flux    = quantile(CH4_flux_nmolpm2ps, 0.75),
    median_noise = median(CH4_noise_floor),
    n            = n(),
    .groups = "drop"
  )

p_ts <- ggplot(monthly_ts, aes(x = date_mid, color = inst_label, fill = inst_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 0.3) +
  # Noise envelope around zero
  geom_ribbon(aes(ymin = -median_noise, ymax = median_noise), alpha = 0.15,
              color = NA) +
  # IQR ribbon
  geom_ribbon(aes(ymin = q25_flux, ymax = q75_flux), alpha = 0.25, color = NA) +
  # Median line
  geom_line(aes(y = median_flux), linewidth = 0.8) +
  geom_point(aes(y = median_flux, size = n), alpha = 0.7) +
  scale_size_continuous(range = c(1.5, 4), name = "n measurements") +
  scale_color_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    name = "Instrument"
  ) +
  scale_fill_manual(
    values = c("LGR/UGGA (2023-24)" = "#E69F00", "LI-7810 (2025)" = "#56B4E9"),
    name = "Instrument"
  ) +
  scale_y_continuous(
    trans = asinh_trans,
    breaks = c(-2, -1, -0.5, 0, 0.5, 1, 2, 5, 10),
    labels = c("-2", "-1", "-0.5", "0", "0.5", "1", "2", "5", "10")
  ) +
  labs(
    x = "Date",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = expression(Monthly~median~CH[4]~flux~with~noise~envelope),
    subtitle = "Line = median flux | Dark ribbon = IQR | Light ribbon = \u00b1 median noise floor (1\u03c3 Allan)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave(file.path(out_dir, "ch4_flux_monthly_timeseries_noise.png"),
       p_ts, width = 11, height = 5, dpi = 300)
ggsave(file.path(out_dir, "ch4_flux_monthly_timeseries_noise.pdf"),
       p_ts, width = 11, height = 5)
message("Saved monthly time series with noise to: ", out_dir)

# ============================================================
# PART 12: BELOW-MDF TREATMENT COMPARISON
# ============================================================
#
# Three approaches for handling fluxes below detection:
#   1. REMOVE:      exclude below-MDF measurements entirely
#   2. SET TO ZERO: replace below-MDF fluxes with 0
#   3. KEEP AS-IS:  retain original values (unfiltered baseline)
#
# Stringent filtering causes a rightward shift (positive bias)
# by selectively removing near-zero fluxes. Setting to zero
# preserves sample size and avoids this bias.
# ============================================================

message("\n=== Below-MDF treatment comparison ===")

# Focus on MDF-based filters only (these have clear thresholds)
mdf_filters <- list(
  "Manufacturer MDF" = "CH4_below_MDF_manuf",
  "Wassmann 90%"     = "CH4_below_MDF_wass90",
  "Wassmann 95%"     = "CH4_below_MDF_wass95",
  "Wassmann 99%"     = "CH4_below_MDF_wass99",
  "Christiansen 90%" = "CH4_below_MDF_chr90",
  "Christiansen 95%" = "CH4_below_MDF_chr95",
  "Christiansen 99%" = "CH4_below_MDF_chr99"
)

# --- Build long data frame with three treatments per filter ---
treat_rows <- list()
k <- 0

for (filt_name in names(mdf_filters)) {
  flag_col <- mdf_filters[[filt_name]]
  below <- df[[flag_col]]
  below[is.na(below)] <- FALSE  # if no Allan data, keep measurement

  n_below <- sum(below)
  n_above <- sum(!below)

  # Treatment 1: Remove below-MDF
  k <- k + 1
  treat_rows[[k]] <- data.frame(
    filter    = filt_name,
    treatment = "Remove below MDF",
    CH4_flux  = df$CH4_flux_nmolpm2ps[!below],
    instrument = df$inst_label[!below],
    stringsAsFactors = FALSE
  )

  # Treatment 2: Set below-MDF to zero
  flux_zeroed <- df$CH4_flux_nmolpm2ps
  flux_zeroed[below] <- 0
  k <- k + 1
  treat_rows[[k]] <- data.frame(
    filter    = filt_name,
    treatment = "Set below MDF to zero",
    CH4_flux  = flux_zeroed,
    instrument = df$inst_label,
    stringsAsFactors = FALSE
  )

  # Treatment 3: Keep original (unfiltered)
  k <- k + 1
  treat_rows[[k]] <- data.frame(
    filter    = filt_name,
    treatment = "Keep all (unfiltered)",
    CH4_flux  = df$CH4_flux_nmolpm2ps,
    instrument = df$inst_label,
    stringsAsFactors = FALSE
  )
}

treat_df <- do.call(rbind, treat_rows)

# Order filters from least to most stringent
filter_order_mdf <- c("Manufacturer MDF",
                       "Wassmann 90%", "Wassmann 95%", "Wassmann 99%",
                       "Christiansen 90%", "Christiansen 95%", "Christiansen 99%")
treat_df$filter <- factor(treat_df$filter, levels = filter_order_mdf)

# Treatment ordering: unfiltered as baseline, then the two alternatives
treat_df$treatment <- factor(treat_df$treatment,
                              levels = c("Keep all (unfiltered)",
                                         "Set below MDF to zero",
                                         "Remove below MDF"))

# --- TABLE: Summary stats by filter × treatment ---
cat("\n=== Summary stats: filter × treatment ===\n")
treat_stats <- treat_df %>%
  group_by(filter, treatment) %>%
  summarise(
    n          = n(),
    mean_flux  = round(mean(CH4_flux), 4),
    median_flux = round(median(CH4_flux), 4),
    pct_neg    = round(100 * mean(CH4_flux < 0), 1),
    pct_zero   = round(100 * mean(CH4_flux == 0), 1),
    .groups    = "drop"
  )
print(as.data.frame(treat_stats), row.names = FALSE)

# --- TABLE: Same split by instrument ---
cat("\n=== Summary stats: filter × treatment × instrument ===\n")
treat_inst_stats <- treat_df %>%
  group_by(filter, treatment, instrument) %>%
  summarise(
    n           = n(),
    mean_flux   = round(mean(CH4_flux), 4),
    median_flux = round(median(CH4_flux), 4),
    pct_neg     = round(100 * mean(CH4_flux < 0), 1),
    .groups     = "drop"
  )
print(as.data.frame(treat_inst_stats), row.names = FALSE)

# --- PLOT 1: Overlaid densities, faceted by filter ---
# Shows how each treatment changes the distribution shape

p_treat_ridges <- ggplot(treat_df,
                          aes(x = CH4_flux, y = filter,
                              fill = treatment, color = treatment)) +
  geom_density_ridges(
    alpha = 0.35, scale = 0.85,
    bandwidth = 0.25,
    from = asinh(-3), to = asinh(20),
    rel_min_height = 0.005
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.4) +
  scale_x_continuous(
    trans = asinh_trans,
    breaks = c(-2, -1, 0, 0.5, 1, 2, 5, 10),
    labels = c("-2", "-1", "0", "0.5", "1", "2", "5", "10")
  ) +
  scale_fill_manual(
    values = c("Keep all (unfiltered)" = "grey70",
               "Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Below-MDF treatment"
  ) +
  scale_color_manual(
    values = c("Keep all (unfiltered)" = "grey40",
               "Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Below-MDF treatment"
  ) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = NULL,
    title = expression(Effect~of~below-MDF~treatment~on~CH[4]~flux~distribution),
    subtitle = paste0("n = ", n_total, " measurements | ",
                      "Remove: excludes below-MDF | ",
                      "Set to zero: replaces with 0 | ",
                      "Unfiltered: keeps original values")
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y  = element_text(size = 9),
    legend.position = "bottom",
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 8.5, color = "grey40")
  ) +
  coord_cartesian(xlim = c(asinh(-3), asinh(15)))

ggsave(file.path(out_dir, "ch4_mdf_treatment_ridges.png"),
       p_treat_ridges, width = 11, height = 8, dpi = 300)
ggsave(file.path(out_dir, "ch4_mdf_treatment_ridges.pdf"),
       p_treat_ridges, width = 11, height = 8)
message("Saved MDF treatment ridges to: ", out_dir)

# --- PLOT 2: Mean / median shift by treatment (dot-and-line) ---
# Directly shows the positive bias from removal

treat_summary_long <- treat_stats %>%
  select(filter, treatment, mean_flux, median_flux) %>%
  pivot_longer(cols = c(mean_flux, median_flux),
               names_to = "statistic", values_to = "value") %>%
  mutate(statistic = recode(statistic,
                             "mean_flux" = "Mean",
                             "median_flux" = "Median"))

p_bias <- ggplot(treat_summary_long,
                  aes(x = filter, y = value, color = treatment, group = treatment)) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~ statistic, ncol = 1, scales = "free_y") +
  scale_color_manual(
    values = c("Keep all (unfiltered)" = "grey50",
               "Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Below-MDF treatment"
  ) +
  scale_y_continuous(
    trans = asinh_trans,
    breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
    labels = c("0", "0.05", "0.1", "0.2", "0.5", "1", "2", "5")
  ) +
  labs(
    x = NULL,
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = expression(Positive~bias~from~removing~below-MDF~measurements),
    subtitle = "Removing below-MDF shifts mean/median upward; setting to zero preserves sample size with minimal bias"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x  = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "bottom",
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    strip.text    = element_text(size = 11, face = "bold")
  )

ggsave(file.path(out_dir, "ch4_mdf_treatment_bias.png"),
       p_bias, width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "ch4_mdf_treatment_bias.pdf"),
       p_bias, width = 10, height = 7)
message("Saved MDF treatment bias plot to: ", out_dir)

# --- PLOT 3: Same as PLOT 1 but split by instrument ---
# Shows the differential impact: removal biases LGR far more than LI-7810

p_treat_inst <- ggplot(treat_df,
                        aes(x = CH4_flux, fill = treatment, color = treatment)) +
  geom_density(alpha = 0.3, linewidth = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", linewidth = 0.3) +
  facet_grid(filter ~ instrument, scales = "free_y") +
  scale_x_continuous(
    trans = asinh_trans,
    breaks = c(-2, -1, 0, 0.5, 1, 2, 5, 10),
    labels = c("-2", "-1", "0", "0.5", "1", "2", "5", "10")
  ) +
  scale_fill_manual(
    values = c("Keep all (unfiltered)" = "grey70",
               "Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Treatment"
  ) +
  scale_color_manual(
    values = c("Keep all (unfiltered)" = "grey40",
               "Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Treatment"
  ) +
  labs(
    x = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    y = "Density",
    title = expression(Below-MDF~treatment~effect~by~instrument),
    subtitle = "Removal causes rightward shift primarily for LGR; LI-7810 is minimally affected"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "bottom",
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    strip.text.y  = element_text(size = 7, angle = 0),
    strip.text.x  = element_text(size = 9, face = "bold")
  ) +
  coord_cartesian(xlim = c(asinh(-3), asinh(15)))

ggsave(file.path(out_dir, "ch4_mdf_treatment_by_instrument.png"),
       p_treat_inst, width = 12, height = 14, dpi = 300)
ggsave(file.path(out_dir, "ch4_mdf_treatment_by_instrument.pdf"),
       p_treat_inst, width = 12, height = 14)
message("Saved MDF treatment × instrument plot to: ", out_dir)

# --- PLOT 4: Percent change in mean/median from removal vs set-to-zero ---
# relative to the unfiltered baseline

baseline <- treat_stats %>%
  filter(treatment == "Keep all (unfiltered)") %>%
  select(filter, base_mean = mean_flux, base_median = median_flux)

pct_change <- treat_stats %>%
  filter(treatment != "Keep all (unfiltered)") %>%
  left_join(baseline, by = "filter") %>%
  mutate(
    pct_change_mean   = round(100 * (mean_flux - base_mean) / abs(base_mean), 1),
    pct_change_median = round(100 * (median_flux - base_median) / abs(base_median), 1)
  )

cat("\n=== % change in mean/median relative to unfiltered baseline ===\n")
print(as.data.frame(pct_change %>%
  select(filter, treatment, mean_flux, base_mean, pct_change_mean,
         median_flux, base_median, pct_change_median)), row.names = FALSE)

pct_long <- pct_change %>%
  select(filter, treatment, pct_change_mean, pct_change_median) %>%
  pivot_longer(cols = c(pct_change_mean, pct_change_median),
               names_to = "statistic", values_to = "pct_change") %>%
  mutate(statistic = recode(statistic,
                             "pct_change_mean"   = "Mean",
                             "pct_change_median" = "Median"))

p_pct <- ggplot(pct_long,
                 aes(x = filter, y = pct_change, fill = treatment)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, alpha = 0.85) +
  geom_hline(yintercept = 0, color = "grey30") +
  facet_wrap(~ statistic, ncol = 1) +
  scale_fill_manual(
    values = c("Set below MDF to zero" = "#2166AC",
               "Remove below MDF"      = "#B2182B"),
    name = "Treatment"
  ) +
  labs(
    x = NULL,
    y = "% change relative to unfiltered",
    title = expression(Bias~introduced~by~below-MDF~treatment),
    subtitle = "Removing below-MDF inflates mean/median; setting to zero introduces much less bias"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
    legend.position = "bottom",
    plot.title    = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 9, color = "grey40"),
    strip.text    = element_text(size = 11, face = "bold")
  )

ggsave(file.path(out_dir, "ch4_mdf_treatment_pct_change.png"),
       p_pct, width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "ch4_mdf_treatment_pct_change.pdf"),
       p_pct, width = 10, height = 7)
message("Saved MDF treatment % change to: ", out_dir)

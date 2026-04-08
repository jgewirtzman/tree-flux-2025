# ============================================================
# 09_quality_flags.R
#
# Compute per-measurement quality flags using the MDF (Minimal
# Detectable Flux) framework. Three MDF approaches:
#   1. Manufacturer MDF = precision / t * flux_term
#   2. Wassmann et al. (2018) = z * global_empirical_SD / t * flux_term
#   3. Christiansen et al. (2015) = SD_allan * 3 * t_crit / t * flux_term
#
# Also computes Allan deviation (per-measurement instrument noise),
# empirical SNR, and SE-based SNR.
#
# Two instruments:
#   - LGR/UGGA GLA131 (2023-24): 0.9 ppb CH4 manufacturer spec
#   - LI-COR LI-7810 (2025):     0.6 ppb CH4 manufacturer spec
#
# Input:
#   - data/input/HF_2023-2025_tree_flux_corrected.csv
#   - data/raw/upland_wetland/processing_csvs/tree_volumes.csv
#   - data/raw/7810_Processed/Tree_Fluxes/ (LI-7810 1-Hz files)
#   - Raw LGR Data (LGR 1-Hz files, external)
#   - Field logs (measurement timestamps for LGR segmentation)
#
# Output:
#   - data/processed/flux_with_quality_flags.csv
# ============================================================

library(dplyr)
library(tidyr)

# ============================================================
# CONSTANTS
# ============================================================

SURFAREA_M2    <- pi * 0.0508^2       # m^2, collar radius = 5.08 cm
R_hPa_L       <- 83.14472            # hPa*L/(mol*K)
EXTRA_TUBE_VOL <- 0.028              # L, connecting tubing not in tree_volumes.csv

# Manufacturer precision (1-sigma, 1 sec)
PREC_CH4_LGR   <- 0.9    # ppb, GLA131 (UGGA)
PREC_CO2_LGR   <- 0.35   # ppm, GLA131
PREC_CH4_7810  <- 0.6    # ppb, LI-7810
PREC_CO2_7810  <- 3.5    # ppm, LI-7810

# Allan deviation function: isolates instrument white noise via
# first-differencing (removes smooth flux trend)
allan_sd <- function(x) {
  x <- x[!is.na(x)]
  diffs <- diff(x)
  if (length(diffs) < 2) return(NA_real_)
  sd(diffs) / sqrt(2)
}

# ============================================================
# PART 1: LOAD CORRECTED FLUX DATA
# ============================================================

message("=== Loading corrected flux data ===")

fluxes <- read.csv(file.path("data", "input",
                              "HF_2023-2025_tree_flux_corrected.csv"),
                   stringsAsFactors = FALSE)
message("Loaded: ", nrow(fluxes), " rows")

df <- fluxes %>%
  filter(
    !is.na(year),
    !is.na(CO2_flux_umolpm2ps), !is.na(CO2_r2), !is.na(CO2_SE),
    !is.na(CH4_flux_nmolpm2ps), !is.na(CH4_r2), !is.na(CH4_SE)
  ) %>%
  mutate(
    # Pre-2025 CH4_SE needs x1000 correction (flux already corrected in file)
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
tv$vol_system_tv <- tv$Volume + EXTRA_TUBE_VOL

message("tree_volumes.csv: ", nrow(tv), " trees")

df <- df %>%
  left_join(tv[, c("Tag_int", "vol_system_tv")],
            by = c("Tree" = "Tag_int")) %>%
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
    machine        = "LGR1",
    source         = "summer2024"
  )

# --- Field log 4: Sept-Dec 2024 Excel files ---
xlsx_dir <- file.path("data", "raw", "upland_wetland", "Sept2024_onwards")
fl4_std <- data.frame()
if (requireNamespace("readxl", quietly = TRUE) && dir.exists(xlsx_dir)) {
  xlsx_files <- list.files(xlsx_dir, pattern = "[.]xlsx$", full.names = TRUE)
  xlsx_files <- xlsx_files[!grepl("Template", xlsx_files)]
  xlsx_files <- xlsx_files[!grepl("2025", xlsx_files)]

  for (xf in xlsx_files) {
    tryCatch({
      xd <- as.data.frame(readxl::read_excel(xf))
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
                    field_logs$machine == ""] <- "LGR1"

# Parse dates — handle multiple formats
field_logs$date_parsed <- as.Date(field_logs$date_raw, format = "%m/%d/%Y")
bad_year <- !is.na(field_logs$date_parsed) &
            as.integer(format(field_logs$date_parsed, "%Y")) < 100
if (any(bad_year)) {
  field_logs$date_parsed[bad_year] <- as.Date(
    field_logs$date_raw[bad_year], format = "%m/%d/%y")
}
na_dates <- is.na(field_logs$date_parsed)
if (any(na_dates)) {
  field_logs$date_parsed[na_dates] <- as.Date(
    field_logs$date_raw[na_dates], format = "%m-%d-%Y")
}
na_dates <- is.na(field_logs$date_parsed)
if (any(na_dates)) {
  field_logs$date_parsed[na_dates] <- as.Date(
    field_logs$date_raw[na_dates], format = "%y-%m-%d")
}
na_dates <- is.na(field_logs$date_parsed)
if (any(na_dates)) {
  field_logs$date_parsed[na_dates] <- as.Date(
    field_logs$date_raw[na_dates], format = "%Y-%m-%d")
}

message("Field log date range: ",
        min(field_logs$date_parsed, na.rm = TRUE), " to ",
        max(field_logs$date_parsed, na.rm = TRUE))

# Build comp_start/comp_end as POSIXct
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
  ) %>%
  filter(!is.na(comp_start_posix), !is.na(comp_end_posix),
         t_sec_fl > 30, t_sec_fl < 1800)

message("After cleaning: ", nrow(field_logs), " field log entries")

# De-duplicate: prefer updated2 over summer2023, summer2024 over updated2
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

allan_df_7810$match_key <- paste(allan_df_7810$REMARK, allan_df_7810$DATE, sep = "_")

# ============================================================
# PART 5: ALLAN DEVIATION - 2023-24 LGR/UGGA
# ============================================================

message("\n=== Allan deviation: 2023-24 LGR/UGGA ===")

lgr_base <- normalizePath(
  file.path("..", "..", "..", "Matthes_Lab", "stem-CH4-flux", "Raw LGR Data"),
  mustWork = FALSE
)
if (!dir.exists(lgr_base)) {
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
    date_part <- substr(folder_name, 1, 10)
    lgr_catalog <- rbind(lgr_catalog, data.frame(
      machine = lgr_name, date_str = date_part, folder = folder,
      stringsAsFactors = FALSE
    ))
  }
}
message("LGR catalog: ", nrow(lgr_catalog), " date-folders across ",
        length(unique(lgr_catalog$machine)), " instruments")

# Parse one LGR raw file
parse_lgr_file <- function(filepath) {
  tryCatch({
    if (grepl("\\.zip$", filepath)) {
      txt_name <- sub("\\.zip$", "", basename(filepath))
      con <- unz(filepath, txt_name)
      lines <- readLines(con, warn = FALSE)
      close(con)
    } else {
      lines <- readLines(filepath, warn = FALSE)
    }
    if (length(lines) < 3) return(NULL)
    dat <- tryCatch(
      read.csv(text = paste(lines[-1], collapse = "\n"),
               stringsAsFactors = FALSE, strip.white = TRUE),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0) return(NULL)
    ts_col <- names(dat)[1]
    dat$timestamp <- as.POSIXct(dat[[ts_col]], format = "%m/%d/%Y %H:%M:%OS")
    ch4_col <- grep("CH4.*d_ppm", names(dat), value = TRUE)[1]
    co2_col <- grep("CO2.*d_ppm", names(dat), value = TRUE)[1]
    if (is.na(ch4_col) || is.na(co2_col)) return(NULL)
    data.frame(
      timestamp = dat$timestamp,
      CH4_ppm   = as.numeric(dat[[ch4_col]]),
      CO2_ppm   = as.numeric(dat[[co2_col]]),
      stringsAsFactors = FALSE
    ) %>% filter(!is.na(timestamp))
  }, error = function(e) NULL)
}

# Load and concatenate all LGR files for a given folder
load_lgr_day <- function(folder_path) {
  all_files <- list.files(folder_path, full.names = TRUE)
  data_files <- all_files[grepl("_f\\d+\\.txt(\\.zip)?$", all_files)]
  if (length(data_files) == 0) return(NULL)
  day_data <- lapply(data_files, parse_lgr_file)
  day_data <- do.call(rbind, Filter(Negate(is.null), day_data))
  if (is.null(day_data) || nrow(day_data) == 0) return(NULL)
  day_data[order(day_data$timestamp), ]
}

# Process each field log entry
message("\nProcessing LGR measurements...")

field_logs_lgr <- field_logs %>%
  filter(date_parsed < as.Date("2025-01-01"))

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

  folders <- lgr_catalog$folder[lgr_catalog$machine == m &
                                  lgr_catalog$date_str == d]
  if (length(folders) == 0) {
    folders <- lgr_catalog$folder[lgr_catalog$date_str == d]
  }
  if (length(folders) == 0) {
    n_missing <- n_missing + 1
    next
  }

  day_data <- NULL
  for (folder in folders) {
    dd <- load_lgr_day(folder)
    if (!is.null(dd)) day_data <- rbind(day_data, dd)
  }
  if (is.null(day_data) || nrow(day_data) < 10) {
    n_missing <- n_missing + 1
    next
  }

  entries <- field_logs_lgr %>% filter(date_str == d, machine == m)

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
      allan_sd_CO2 = allan_sd(segment$CO2_ppm),
      allan_sd_CH4 = allan_sd(segment$CH4_ppm) * 1000,  # ppm -> ppb
      n_pts        = nrow(segment),
      t_sec        = nrow(segment),
      instrument   = m,
      stringsAsFactors = FALSE
    )
  }

  if (i %% 10 == 0 || i == nrow(unique_day_machine)) {
    message("  Processed ", i, "/", nrow(unique_day_machine),
            " dates (", n_matched, " measurements matched)")
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

# --- Datetime fallback for unmatched 2023-24 rows ---
unmatched_idx <- which(is.na(df$allan_sd_CH4) & df$year < 2025)
message("Unmatched 2023-24 rows before datetime fallback: ", length(unmatched_idx))

if (length(unmatched_idx) > 0 && nrow(field_logs_lgr) > 0 && nrow(allan_df_lgr) > 0) {
  fl_lookup <- field_logs %>%
    filter(!is.na(date_parsed)) %>%
    mutate(
      real_start_posix = as.POSIXct(paste(date_parsed, real_start),
                                     format = "%Y-%m-%d %H:%M:%S")
    ) %>%
    filter(!is.na(real_start_posix))

  for (idx in unmatched_idx) {
    dt_str <- df$datetime_posx[idx]
    if (is.na(dt_str) || dt_str == "") next
    dt <- as.POSIXct(dt_str, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC")
    if (is.na(dt)) dt <- as.POSIXct(dt_str, format = "%Y-%m-%dT%H:%M:%S", tz = "UTC")
    if (is.na(dt)) next

    same_date <- fl_lookup$date_str == df$date[idx]
    if (!any(same_date, na.rm = TRUE)) next
    fl_sub <- fl_lookup[same_date & !is.na(same_date), ]
    if (nrow(fl_sub) == 0) next

    time_diffs <- abs(as.numeric(difftime(fl_sub$real_start_posix, dt, units = "secs")))
    best <- which.min(time_diffs)
    if (time_diffs[best] > 120) next

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
    flux_term = nmol / surfarea,
    t_sec_est = if_else(!is.na(t_sec), as.numeric(t_sec), 300),

    # Empirical noise floor in flux units
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

    # Instrument label
    inst_label = if_else(year == 2025, "LI-7810", "LGR/UGGA"),

    # Manufacturer precision (instrument-specific)
    prec_ch4 = if_else(year == 2025, PREC_CH4_7810, PREC_CH4_LGR),
    prec_co2 = if_else(year == 2025, PREC_CO2_7810, PREC_CO2_LGR),

    # 1. Manufacturer MDF
    CH4_MDF_manufacturer = ifelse(!is.na(flux_term),
                                   prec_ch4 / t_sec_est * flux_term,
                                   NA_real_),

    # 3. Christiansen MDF = allan_sd * 3 * t_crit / t * flux_term
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

# --- 2. Wassmann MDF: instrument-specific global Allan SD ---
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

# ============================================================
# CLEAN UP AND SAVE
# ============================================================

# Canonical MDF flag: Wassmann 95% is the recommended detection threshold
df <- df %>%
  mutate(CH4_below_MDF = CH4_below_MDF_wass95)

# Drop intermediate columns not needed downstream
df <- df %>%
  select(-match_key_7810, -df_meas, -t99, -t95, -t90,
         -prec_ch4, -prec_co2, -global_sd_ch4)

output_path <- file.path("data", "processed", "flux_with_quality_flags.csv")
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
write.csv(df, output_path, row.names = FALSE)

message("\n=== Output saved ===")
message("  ", output_path)
message("  Rows: ", nrow(df))
message("  Primary flag: CH4_below_MDF (Wassmann 95%)")
message("  All MDF columns: CH4_below_MDF_manuf, CH4_below_MDF_wass{90,95,99},")
message("    CH4_below_MDF_chr{90,95,99}, CH4_snr_allan, CH4_noise_floor,")
message("    allan_sd_CH4, inst_label, flux_term, t_sec")

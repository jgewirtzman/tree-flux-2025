# ============================================================
# 07_build_aligned_dataset.R
# Build an aligned hourly dataset for rolling window analysis
# 
# INPUTS (run preprocessing scripts first):
#   - data/processed/wtd_met.csv          (from 03/04 scripts - HF met/hydro)
#   - data/processed/neon_swc_hourly.csv  (from 06a - NEON soil moisture)
#   - data/processed/tower_swc_ts_hourly.csv (from 06b - tower SWC/TS)
#   - AmeriFlux tower CSVs                (for flux, met, CH4, etc.)
#   - Phenocam data                       (GCC, NDVI)
#
# OUTPUT:
#   - data/processed/aligned_hourly_dataset.csv
# ============================================================

library(tidyverse)
library(lubridate)
library(stringr)

source("scripts/helpers/find_ameriflux.R")

# ============================================================
# CONFIGURATION
# ============================================================

START_YEAR <- 2023
END_YEAR <- 2026

# File paths
PATHS <- list(
  # Preprocessed data
  wtd_met = "data/processed/wtd_met.csv",
  neon_swc = "data/processed/neon_swc_hourly.csv",
  tower_swc_ts = "data/processed/tower_swc_ts_hourly.csv",

  # AmeriFlux tower data (version-agnostic lookup)
  Ha1 = find_ameriflux("US-Ha1"),
  Ha2 = find_ameriflux("US-Ha2"),
  xHA = find_ameriflux("US-xHA"),

  # Phenocam
  phenocam = "data/raw/phenocam/harvardems2_DB_1000_ndvi_1day.txt"
)

OUTPUT_PATH <- "data/processed/aligned_hourly_dataset.csv"

# ============================================================
# HELPER FUNCTIONS
# ============================================================

to_posix_af <- function(x) {
  x_chr <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  x_chr <- str_pad(x_chr, width = 12, side = "left", pad = "0")
  ymd_hm(x_chr, tz = "UTC")
}

#' Extract and rename specified columns from a tower dataframe
extract_tower_vars <- function(df, var_map, tower_suffix = "",
                               start_year = 2022, end_year = 2026) {
  
  existing <- var_map[var_map %in% names(df)]
  missing <- var_map[!var_map %in% names(df)]
  
  if (length(missing) > 0) {
    message("    Missing: ", paste(head(missing, 5), collapse = ", "), 
            if(length(missing) > 5) paste0(" ... (", length(missing), " total)"))
  }
  
  if (length(existing) == 0) {
    warning("No specified columns found")
    return(NULL)
  }
  
  result <- df %>%
    mutate(datetime = to_posix_af(TIMESTAMP_START)) %>%
    filter(year(datetime) >= start_year, year(datetime) <= end_year) %>%
    select(datetime, all_of(unname(existing)))
  
  result <- result %>%
    mutate(across(-datetime, ~ na_if(as.numeric(.x), -9999)))
  
  # Rename columns
  new_names <- names(existing)
  old_names <- unname(existing)
  
  for (i in seq_along(existing)) {
    names(result)[names(result) == old_names[i]] <- paste0(new_names[i], tower_suffix)
  }
  
  result
}

aggregate_to_hourly <- function(df, datetime_col = "datetime") {
  value_cols <- setdiff(names(df), datetime_col)
  
  df %>%
    mutate(datetime_hour = floor_date(!!sym(datetime_col), unit = "hour")) %>%
    group_by(datetime_hour) %>%
    summarize(across(all_of(value_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    rename(datetime = datetime_hour)
}

# ============================================================
# DEFINE TOWER VARIABLES TO EXTRACT
# ============================================================

# NOTE on precipitation:
# - P_mm in wtd_met comes from HF Fisher Met station (dt10$prec)
# - P at Ha1 is tower precip
# - P_1_1_1, P_2_1_1 at xHA are two precip gauges (will average)
# - THROUGHFALL at xHA is canopy throughfall (will average 5 gauges)

# Ha1: fluxes only (met comes from Fisher)
# Variables to aggregate: FC (2 replicates -> mean), LE (2 replicates -> mean)
Ha1_vars <- c(
  # CO2 Flux (2 replicates -> will average)
  FC_1_1_1 = "FC_1_1_1", FC_1_1_2 = "FC_1_1_2",
  SC = "SC",
  # Latent heat (2 replicates -> will average)
  LE_1_1_1 = "LE_1_1_1", LE_1_1_2 = "LE_1_1_2",
  # Sensible heat
  H = "H",
  # Friction velocity
  USTAR = "USTAR"
)

# Ha2: fluxes and CO2 profile only (met comes from Fisher)
# Variables to aggregate: FC (3 replicates -> mean), LE (3 replicates -> mean),
#                         H (2 replicates -> mean), CO2 profile (8 heights -> mean)
Ha2_vars <- c(
  # CO2 Flux (3 replicates -> will average)
  FC_1_1_1 = "FC_1_1_1", FC_2_1_1 = "FC_2_1_1", FC_2_1_2 = "FC_2_1_2",
  SC_2_1_1 = "SC_2_1_1",
  # Latent heat (3 replicates -> will average)
  LE_1_1_1 = "LE_1_1_1", LE_2_1_1 = "LE_2_1_1", LE_2_1_2 = "LE_2_1_2",
  # Sensible heat (2 replicates -> will average)
  H_1_1_1 = "H_1_1_1", H_2_1_1 = "H_2_1_1",
  # CO2 profile (8 heights -> will average to single CO2_MR)
  CO2_3_1_1 = "CO2_3_1_1", CO2_3_2_1 = "CO2_3_2_1", CO2_3_3_1 = "CO2_3_3_1",
  CO2_3_4_1 = "CO2_3_4_1", CO2_3_5_1 = "CO2_3_5_1", CO2_3_6_1 = "CO2_3_6_1",
  CO2_3_7_1 = "CO2_3_7_1", CO2_3_8_1 = "CO2_3_8_1",
  # Friction velocity (2 sensors -> will average)
  USTAR_1_1_1 = "USTAR_1_1_1", USTAR_2_1_1 = "USTAR_2_1_1"
)

# xHA: CH4 and CO2 mixing ratios, ground heat, canopy temp, throughfall only
# Variables to aggregate:
#   - CH4 mixing ratio (6 heights -> mean)
#   - CO2 mixing ratio (6 heights -> mean)
#   - T_CANOPY (5 sensors -> mean)
#   - G ground heat (3 sensors -> mean)
#   - THROUGHFALL (5 gauges -> mean)
#   - REMOVED: SWC, LE, H, H2O_MR, NEE, P, TA, RH, VPD, PA, PPFD, NETRAD
xHA_vars <- c(
  # CO2 Flux (single sensor)
  FC = "FC", 
  SC = "SC",
  # CH4 mixing ratio only (6 heights -> will average)
  CH4_MR_1_1_1 = "CH4_MIXING_RATIO_1_1_1", CH4_MR_1_2_1 = "CH4_MIXING_RATIO_1_2_1",
  CH4_MR_1_3_1 = "CH4_MIXING_RATIO_1_3_1", CH4_MR_1_4_1 = "CH4_MIXING_RATIO_1_4_1",
  CH4_MR_1_5_1 = "CH4_MIXING_RATIO_1_5_1", CH4_MR_1_6_1 = "CH4_MIXING_RATIO_1_6_1",
  # CO2 mixing ratio (6 heights -> will average)
  CO2_MR_1_1_1 = "CO2_MIXING_RATIO_1_1_1", CO2_MR_1_2_2 = "CO2_MIXING_RATIO_1_2_2",
  CO2_MR_1_3_2 = "CO2_MIXING_RATIO_1_3_2", CO2_MR_1_4_2 = "CO2_MIXING_RATIO_1_4_2",
  CO2_MR_1_5_2 = "CO2_MIXING_RATIO_1_5_2", CO2_MR_1_6_2 = "CO2_MIXING_RATIO_1_6_2",
  # Canopy temperature (5 sensors -> will average)
  T_CANOPY_1_1_1 = "T_CANOPY_1_1_1", T_CANOPY_1_2_1 = "T_CANOPY_1_2_1",
  T_CANOPY_1_3_1 = "T_CANOPY_1_3_1", T_CANOPY_1_4_1 = "T_CANOPY_1_4_1",
  T_CANOPY_2_4_1 = "T_CANOPY_2_4_1",
  # Throughfall (5 gauges -> will average)
  THROUGHFALL_1_1_1 = "THROUGHFALL_1_1_1", THROUGHFALL_2_1_1 = "THROUGHFALL_2_1_1",
  THROUGHFALL_3_1_1 = "THROUGHFALL_3_1_1", THROUGHFALL_4_1_1 = "THROUGHFALL_4_1_1",
  THROUGHFALL_5_1_1 = "THROUGHFALL_5_1_1",
  # Ground heat flux (3 sensors -> will average)
  G_1_1_1 = "G_1_1_1", G_3_1_1 = "G_3_1_1", G_5_1_1 = "G_5_1_1",
  # Wind
  WS_1_1_1 = "WS_1_1_1", WD_1_1_1 = "WD_1_1_1", USTAR = "USTAR"
)

# ============================================================
# 1. LOAD PREPROCESSED DATA
# ============================================================

message("============================================================")
message("LOADING PREPROCESSED DATA")
message("============================================================")

# HF met/hydro (from scripts 03-04)
if (file.exists(PATHS$wtd_met)) {
  wtd_met <- read_csv(PATHS$wtd_met, show_col_types = FALSE) %>%
    mutate(datetime = as.POSIXct(datetime, tz = "UTC")) %>%
    # Remove ems_met variables - keep only Fisher Met and Hydro
    select(-any_of(c("vwc_top", "vwc_60cm", "LE", "t_airC_ems")))
  message("Loaded wtd_met: ", nrow(wtd_met), " rows, ", ncol(wtd_met) - 1, " variables (Fisher Met + Hydro only)")
} else {
  message("WARNING: wtd_met.csv not found - HF met/hydro data will be missing")
  message("  Run scripts 03-04 and save: write_csv(wtd_met, 'data/processed/wtd_met.csv')")
  wtd_met <- NULL
}

# NEON SWC (from script 06a)
if (file.exists(PATHS$neon_swc)) {
  neon_swc <- read_csv(PATHS$neon_swc, show_col_types = FALSE) %>%
    mutate(datetime = as.POSIXct(datetime, tz = "UTC"))
  message("Loaded NEON SWC: ", nrow(neon_swc), " rows, ", ncol(neon_swc) - 1, " variables")
} else {
  message("WARNING: neon_swc_hourly.csv not found - run 06a_process_neon_swc.R first")
  neon_swc <- NULL
}

# Tower SWC/TS (from script 06b)
if (file.exists(PATHS$tower_swc_ts)) {
  tower_swc_ts <- read_csv(PATHS$tower_swc_ts, show_col_types = FALSE) %>%
    mutate(datetime = as.POSIXct(datetime, tz = "UTC")) %>%
    # Remove xHA SWC (keep only Ha1, Ha2 SWC and all TS)
    select(-any_of("SWC_xHA"))
  message("Loaded Tower SWC/TS: ", nrow(tower_swc_ts), " rows, ", ncol(tower_swc_ts) - 1, " variables (excluding SWC_xHA)")
} else {
  message("WARNING: tower_swc_ts_hourly.csv not found - run 06b_process_tower_swc_ts.R first")
  tower_swc_ts <- NULL
}

# ============================================================
# 2. LOAD AND PROCESS AMERIFLUX TOWER DATA
# ============================================================

message("\n============================================================")
message("PROCESSING AMERIFLUX TOWER DATA")
message("============================================================")

Ha1 <- read.csv(PATHS$Ha1, header = TRUE, skip = 2)
Ha2 <- read.csv(PATHS$Ha2, header = TRUE, skip = 2)
xHA <- read.csv(PATHS$xHA, header = TRUE, skip = 2)

message("\n--- Ha1 ---")
Ha1_data <- extract_tower_vars(Ha1, Ha1_vars, tower_suffix = "_Ha1",
                               start_year = START_YEAR, end_year = END_YEAR)
Ha1_hourly <- aggregate_to_hourly(Ha1_data)

# Aggregate Ha1 variables - average FC and LE replicates
Ha1_hourly <- Ha1_hourly %>%
  rowwise() %>%
  mutate(
    # Average FC replicates
    FC_Ha1 = mean(c_across(starts_with("FC_") & ends_with("_Ha1")), na.rm = TRUE),
    # Average LE replicates
    LE_Ha1 = mean(c_across(starts_with("LE_") & ends_with("_Ha1")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(datetime, FC_Ha1, LE_Ha1, SC_Ha1, H_Ha1, USTAR_Ha1)

message("  ", nrow(Ha1_hourly), " hourly records, ", ncol(Ha1_hourly) - 1, " variables")

message("\n--- Ha2 ---")
Ha2_data <- extract_tower_vars(Ha2, Ha2_vars, tower_suffix = "_Ha2",
                               start_year = START_YEAR, end_year = END_YEAR)
Ha2_hourly <- aggregate_to_hourly(Ha2_data)

# Aggregate Ha2 variables - average replicates and profiles
Ha2_hourly <- Ha2_hourly %>%
  rowwise() %>%
  mutate(
    # Average FC replicates
    FC_Ha2 = mean(c_across(starts_with("FC_") & ends_with("_Ha2")), na.rm = TRUE),
    # Average LE replicates
    LE_Ha2 = mean(c_across(starts_with("LE_") & ends_with("_Ha2")), na.rm = TRUE),
    # Average H replicates
    H_Ha2 = mean(c_across(starts_with("H_") & ends_with("_Ha2")), na.rm = TRUE),
    # Average CO2 profile to single value
    CO2_MR_Ha2 = mean(c_across(starts_with("CO2_3_") & ends_with("_Ha2")), na.rm = TRUE),
    # Average USTAR
    USTAR_Ha2 = mean(c_across(starts_with("USTAR_") & ends_with("_Ha2")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(datetime, FC_Ha2, LE_Ha2, H_Ha2, SC_2_1_1_Ha2, CO2_MR_Ha2, USTAR_Ha2) %>%
  rename(SC_Ha2 = SC_2_1_1_Ha2)

message("  ", nrow(Ha2_hourly), " hourly records, ", ncol(Ha2_hourly) - 1, " variables")

message("\n--- xHA ---")
xHA_data <- extract_tower_vars(xHA, xHA_vars, tower_suffix = "_xHA",
                               start_year = START_YEAR, end_year = END_YEAR)
xHA_hourly <- aggregate_to_hourly(xHA_data)

# Aggregate xHA variables
xHA_hourly <- xHA_hourly %>%
  rowwise() %>%
  mutate(
    # Average CH4 mixing ratio across heights
    CH4_MR_xHA = mean(c_across(starts_with("CH4_MR_") & ends_with("_xHA")), na.rm = TRUE),
    # Average CO2 mixing ratio across heights
    CO2_MR_xHA = mean(c_across(starts_with("CO2_MR_") & ends_with("_xHA")), na.rm = TRUE),
    # Average T_CANOPY
    T_CANOPY_xHA = mean(c_across(starts_with("T_CANOPY_") & ends_with("_xHA")), na.rm = TRUE),
    # Average ground heat flux
    G_xHA = mean(c_across(starts_with("G_") & ends_with("_xHA")), na.rm = TRUE),
    # Average throughfall gauges
    THROUGHFALL_xHA = mean(c_across(starts_with("THROUGHFALL_") & ends_with("_xHA")), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  select(datetime, 
         # Fluxes (FC only, no LE/H/NEE)
         FC_xHA, SC_xHA,
         # Gas mixing ratios (aggregated, no H2O)
         CH4_MR_xHA, CO2_MR_xHA,
         # Temperature
         T_CANOPY_xHA,
         # Throughfall
         THROUGHFALL_xHA,
         # Ground heat (aggregated)
         G_xHA,
         # Wind
         WS_1_1_1_xHA, WD_1_1_1_xHA, USTAR_xHA) %>%
  rename(WS_xHA = WS_1_1_1_xHA, WD_xHA = WD_1_1_1_xHA)

message("  ", nrow(xHA_hourly), " hourly records, ", ncol(xHA_hourly) - 1, " variables")

# ============================================================
# 3. LOAD AND PROCESS PHENOCAM DATA
# ============================================================

message("\n============================================================")
message("PROCESSING PHENOCAM DATA")
message("============================================================")

phenocam_raw <- read_csv(PATHS$phenocam, comment = "#", show_col_types = FALSE)

if ("gcc_mean" %in% names(phenocam_raw)) {
  phenocam_raw$gcc <- phenocam_raw$gcc_mean
} else if ("gcc_90" %in% names(phenocam_raw)) {
  phenocam_raw$gcc <- phenocam_raw$gcc_90
} else {
  phenocam_raw$gcc <- NA_real_
}

phenocam_daily <- phenocam_raw %>%
  mutate(date = as.Date(date)) %>%
  filter(year(date) >= START_YEAR, year(date) <= END_YEAR) %>%
  select(date, gcc, any_of("ndvi_90")) %>%
  rename_with(~ ifelse(.x == "ndvi_90", "ndvi", .x))

# Expand to hourly
phenocam_hourly <- phenocam_daily %>%
  crossing(hour = 0:23) %>%
  mutate(datetime = ymd_hms(paste(date, sprintf("%02d:00:00", hour)), tz = "UTC")) %>%
  select(datetime, gcc, any_of("ndvi"))

message("Phenocam: ", nrow(phenocam_hourly), " hourly records (expanded from ", nrow(phenocam_daily), " days)")

# ============================================================
# 4. BUILD ALIGNED DATASET
# ============================================================

message("\n============================================================")
message("BUILDING ALIGNED DATASET")
message("============================================================")

# Collect all datetimes
all_datetimes <- unique(c(
  Ha1_hourly$datetime,
  Ha2_hourly$datetime,
  xHA_hourly$datetime,
  phenocam_hourly$datetime,
  if (!is.null(wtd_met)) wtd_met$datetime else NULL,
  if (!is.null(neon_swc)) neon_swc$datetime else NULL,
  if (!is.null(tower_swc_ts)) tower_swc_ts$datetime else NULL
))
all_datetimes <- sort(all_datetimes[!is.na(all_datetimes)])

# Start with datetime backbone
aligned_data <- tibble(datetime = all_datetimes)
message("Datetime backbone: ", nrow(aligned_data), " hours")

# Join preprocessed data
if (!is.null(wtd_met)) {
  aligned_data <- aligned_data %>% left_join(wtd_met, by = "datetime")
  message("  + HF met/hydro: ", ncol(wtd_met) - 1, " variables")
}

if (!is.null(neon_swc)) {
  aligned_data <- aligned_data %>% left_join(neon_swc, by = "datetime")
  message("  + NEON SWC: ", ncol(neon_swc) - 1, " variables")
}

if (!is.null(tower_swc_ts)) {
  aligned_data <- aligned_data %>% left_join(tower_swc_ts, by = "datetime")
  message("  + Tower SWC/TS: ", ncol(tower_swc_ts) - 1, " variables")
}

# Join tower flux/met data
aligned_data <- aligned_data %>%
  left_join(Ha1_hourly, by = "datetime") %>%
  left_join(Ha2_hourly, by = "datetime") %>%
  left_join(xHA_hourly, by = "datetime")
message("  + Tower flux/met: Ha1 (", ncol(Ha1_hourly) - 1, "), Ha2 (", 
        ncol(Ha2_hourly) - 1, "), xHA (", ncol(xHA_hourly) - 1, ")")

# Join phenocam
aligned_data <- aligned_data %>% left_join(phenocam_hourly, by = "datetime")
message("  + Phenocam: ", ncol(phenocam_hourly) - 1, " variables")

# Add time components
aligned_data <- aligned_data %>%
  arrange(datetime) %>%
  mutate(
    date = as.Date(datetime),
    year = year(datetime),
    month = month(datetime),
    doy = yday(datetime),
    hour = hour(datetime),
    wyear = ifelse(month >= 10, year + 1, year)
  ) %>%
  filter(year >= START_YEAR, year <= END_YEAR)

# ============================================================
# 5. SUMMARY
# ============================================================

message("\n============================================================")
message("ALIGNED DATASET SUMMARY")
message("============================================================")
message("Total rows: ", nrow(aligned_data))
message("Total columns: ", ncol(aligned_data))
message("Date range: ", min(aligned_data$datetime), " to ", max(aligned_data$datetime))

# Coverage by source
message("\nVariable counts by source:")
message("  HF met/hydro: ", sum(grepl("^(tair|RH$|VPD_kPa|PAR$|P_mm|rnet|slrr|s10t|bgs|bvs|vwc|LE$|t_air)", names(aligned_data))))
message("  NEON SWC: ", sum(grepl("^NEON_SWC", names(aligned_data))))
message("  Tower SWC: ", sum(grepl("^SWC_", names(aligned_data))))
message("  Tower TS: ", sum(grepl("^TS_", names(aligned_data))))
message("  Tower Ha1: ", sum(grepl("_Ha1$", names(aligned_data))))
message("  Tower Ha2: ", sum(grepl("_Ha2$", names(aligned_data))))
message("  Tower xHA: ", sum(grepl("_xHA$", names(aligned_data))))
message("  Phenocam: ", sum(grepl("^(gcc|ndvi)$", names(aligned_data))))

# Print final variable list
message("\nFinal variables:")
vars_to_print <- setdiff(names(aligned_data), c("datetime", "date", "year", "month", "doy", "hour", "wyear"))
for (v in vars_to_print) {
  n_valid <- sum(!is.na(aligned_data[[v]]))
  pct <- round(100 * n_valid / nrow(aligned_data), 1)
  message("  ", v, ": ", pct, "% coverage")
}

# ============================================================
# 6. SAVE
# ============================================================

message("\n--- Saving ---")
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)
write_csv(aligned_data, OUTPUT_PATH)
message("Saved to: ", OUTPUT_PATH)

# Save variable list
var_list <- tibble(
  variable = names(aligned_data),
  n_valid = sapply(aligned_data, function(x) sum(!is.na(x))),
  pct_valid = round(sapply(aligned_data, function(x) sum(!is.na(x)) / length(x) * 100), 2)
) %>%
  arrange(desc(pct_valid))

write_csv(var_list, sub(".csv", "_variables.csv", OUTPUT_PATH))
message("Variable list saved to: ", sub(".csv", "_variables.csv", OUTPUT_PATH))

message("\n============================================================")
message("DONE")
message("============================================================")
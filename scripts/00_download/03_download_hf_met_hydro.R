# ============================================================
# 03_download_hf_met_hydro.R
# Download Harvard Forest meteorological and hydrological data
# and process into wtd_met.csv
#
# Downloads from Harvard Forest LTER Data Archive:
#   - HF001: Fisher Met Station (15-min air temp, RH, precip, etc.)
#   - HF070: Stream temperature & discharge (15-min, includes BGS/BVS water table)
#
# Processes 15-min data to hourly and joins into wtd_met.csv
#
# Output: data/processed/wtd_met.csv
#
# Requires: tidyverse, lubridate, plantecophys
# ============================================================

library(tidyverse)
library(lubridate)
library(plantecophys)

# ============================================================
# CONFIGURATION
# ============================================================

OUTPUT_PATH <- "data/processed/wtd_met.csv"
dir.create(dirname(OUTPUT_PATH), recursive = TRUE, showWarnings = FALSE)

# Harvard Forest LTER Data Archive URLs
# These are stable PASTA endpoints for the datasets
FISHER_MET_URL <- "https://pasta.lternet.edu/package/data/eml/knb-lter-hfr/1/34/0b439e8fea983c9e20bb2bfaf91931e6"
HYDRO_URL      <- "https://pasta.lternet.edu/package/data/eml/knb-lter-hfr/70/36/2983b2adba6675e805d144a05087924d"

# ============================================================
# STEP 1: Download Fisher Met Station (HF001)
# ============================================================

message("=== Downloading Fisher Met Station data (HF001) ===")

fisher_file <- tempfile(fileext = ".csv")
tryCatch({
  download.file(FISHER_MET_URL, fisher_file, method = "curl", quiet = TRUE)
  if (is.na(file.size(fisher_file))) {
    download.file(FISHER_MET_URL, fisher_file, method = "auto", quiet = TRUE)
  }
  message("  Downloaded Fisher Met data")
}, error = function(e) {
  stop("Failed to download Fisher Met data: ", e$message,
       "\n  URL: ", FISHER_MET_URL)
})

dt10 <- read.csv(fisher_file, header = FALSE, skip = 1, sep = ",",
  col.names = c(
    "datetime", "jd",
    "airt", "f.airt",
    "rh", "f.rh",
    "dewp", "f.dewp",
    "prec", "f.prec",
    "slrr", "f.slrr",
    "parr", "f.parr",
    "netr", "f.netr",
    "bar", "f.bar",
    "wspd", "f.wspd",
    "wres", "f.wres",
    "wdir", "f.wdir",
    "wdev", "f.wdev",
    "gspd", "f.gspd",
    "s10t", "f.s10t"
  ), check.names = TRUE
)

unlink(fisher_file)
message("  Rows: ", nrow(dt10))

# ============================================================
# STEP 2: Download Hydro data (HF070)
# ============================================================

message("\n=== Downloading Stream/Hydro data (HF070) ===")

hydro_file <- tempfile(fileext = ".csv")
tryCatch({
  download.file(HYDRO_URL, hydro_file, method = "curl", quiet = TRUE)
  if (is.na(file.size(hydro_file))) {
    download.file(HYDRO_URL, hydro_file, method = "auto", quiet = TRUE)
  }
  message("  Downloaded Hydro data")
}, error = function(e) {
  stop("Failed to download Hydro data: ", e$message,
       "\n  URL: ", HYDRO_URL)
})

dt4 <- read.csv(hydro_file, header = FALSE, skip = 1, sep = ",",
  col.names = c(
    "datetime", "jd",
    "nb.stg", "f.nb.stg",
    "nl.stg", "f.nl.stg",
    "al.stg", "f.al.stg",
    "au.stg", "f.au.stg",
    "bgs.stg", "f.bgs.stg",
    "bvs.stg", "f.bvs.stg",
    "nb.dis", "f.nb.dis",
    "nl.dis", "f.nl.dis",
    "nt.dis", "f.nt.dis",
    "al.dis", "f.al.dis",
    "au.dis", "f.au.dis",
    "nb.wt", "f.nb.wt",
    "nl.wt", "f.nl.wt",
    "al.wt", "f.al.wt",
    "au.wt", "f.au.wt",
    "bgs.wt", "f.bgs.wt",
    "bvs.wt", "f.bvs.wt"
  ), check.names = TRUE
)

unlink(hydro_file)
message("  Rows: ", nrow(dt4))

# ============================================================
# STEP 3: Parse datetimes and filter
# ============================================================

message("\n=== Processing data ===")

# Parse hydro datetimes and filter to 2022+
hydro15 <- dt4 %>%
  mutate(
    datetime = ymd_hm(datetime),
    hour = hour(datetime),
    date = as.Date(datetime),
    year = year(datetime)
  ) %>%
  filter(year > 2022)

rm(dt4)

# Parse Fisher Met datetimes
fishermet <- dt10 %>%
  mutate(
    datetime = ymd_hm(datetime),
    hour = hour(datetime),
    date = as.Date(datetime),
    year = year(datetime),
    month = month(datetime)
  )

rm(dt10)

message("  Hydro records (post-2022): ", nrow(hydro15))
message("  Fisher Met records: ", nrow(fishermet))

# ============================================================
# STEP 4: Extract water table depth at hourly scale
# ============================================================

message("\n--- Aggregating water table to hourly ---")

wtd_hourly <- hydro15 %>%
  select(datetime, jd, bgs.stg, bvs.stg, bgs.wt, bvs.wt) %>%
  rename(
    bgs_wtd_cm = bgs.stg,
    bvs_wtd_cm = bvs.stg,
    bgs_wtemp_C = bgs.wt,
    bvs_wtemp_C = bvs.wt
  ) %>%
  mutate(
    hour = str_pad(hour(datetime), pad = "0", width = 2, side = "left"),
    date = as.Date(datetime),
    datetime = ymd_hms(paste0(date, " ", hour, ":00:00"), tz = "UTC")
  ) %>%
  group_by(datetime) %>%
  summarize(
    bgs_wtd_cm = mean(bgs_wtd_cm, na.rm = TRUE),
    bvs_wtd_cm = mean(bvs_wtd_cm, na.rm = TRUE),
    .groups = "drop"
  )

message("  Water table hourly rows: ", nrow(wtd_hourly))

# ============================================================
# STEP 5: Aggregate Fisher Met 15-min to hourly
# ============================================================

message("--- Aggregating Fisher Met to hourly ---")

met_hourly <- fishermet %>%
  group_by(date, hour) %>%
  summarize(
    tair_C = mean(airt, na.rm = TRUE),
    p_kPa = mean(bar / 10, na.rm = TRUE),
    P_mm = sum(prec, na.rm = TRUE),
    RH = mean(rh, na.rm = TRUE),
    PAR = mean(parr, na.rm = TRUE),
    rnet = mean(netr, na.rm = TRUE),
    slrr = mean(slrr, na.rm = TRUE),
    s10t = mean(s10t, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    datetime = ymd_hms(
      paste0(date, " ", str_pad(hour, pad = "0", side = "left", width = 2), ":00:00"),
      tz = "UTC"
    ),
    VPD_kPa = RHtoVPD(RH = RH, TdegC = tair_C, Pa = p_kPa * 1000)
  )

message("  Fisher Met hourly rows: ", nrow(met_hourly))

# ============================================================
# STEP 6: Join water table + met → wtd_met
# ============================================================

message("\n--- Joining water table + met ---")

wtd_met <- wtd_hourly %>%
  left_join(met_hourly, by = "datetime") %>%
  arrange(datetime)

message("  Final wtd_met rows: ", nrow(wtd_met))
message("  Columns: ", paste(names(wtd_met), collapse = ", "))
message("  Date range: ", min(wtd_met$datetime, na.rm = TRUE),
        " to ", max(wtd_met$datetime, na.rm = TRUE))

# ============================================================
# STEP 7: Save
# ============================================================

write_csv(wtd_met, OUTPUT_PATH)

message("\n============================================================")
message("DONE")
message("============================================================")
message("Output saved to: ", OUTPUT_PATH)

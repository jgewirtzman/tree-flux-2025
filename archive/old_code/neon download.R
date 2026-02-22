# ============================================================
# NEON HARV (2023-01 to 2025-12): Flux + Met + Soil clean 30-min file
#
# IMPORTANT FIX:
# - DP4.00200.001 cannot be handled by loadByProduct(); use zipsByProduct() to download
# - Do NOT assume zipsByProduct() returns a data.frame; use the filesystem as source of truth
#
# WHAT THIS SCRIPT DOES:
# A) DP4.00200.001 (eddy covariance bundle)
#    - downloads ZIPs (optional if you already downloaded)
#    - unzips
#    - reads HDF5 to extract 30-min flux series:
#        fluxCo2, fluxH2o (latent heat flux), fluxTemp (sensible heat flux)
#      and their qfFinl flags (qfqm)
#
# B) DP1 products (meteorology, radiation, precip, soils)
#    - downloads with loadByProduct()
#    - stacks to tables with stackByTable()
#    - keeps 30-min tables
#    - soils: keeps each depth DISTINCT (columns per depth); averages only across replicates at same depth+time
#
# C) Derived variables
#    - ET (mm / 30-min) from latent heat flux
#    - VPD (kPa) from RH (%) and Tair (°C) (from DP1.00098.001 if present)
#
# OUTPUTS (in out_dir):
#   HARV_30min_clean_ALL_*.parquet / .csv
#   HARV_30min_clean_QF0_*.parquet / .csv
# ============================================================

# ---------------------------
# Packages
# ---------------------------
pkgs <- c(
  "neonUtilities", "rhdf5",
  "data.table", "lubridate", "stringr", "fs",
  "arrow"
)
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ---------------------------
# User settings
# ---------------------------
site <- "HARV"
start_month <- "2023-01"
end_month   <- "2025-12"

out_dir <- "HARV_NEON_clean_30min"
dir_create(out_dir)

# If you already downloaded DP4 ZIPs, set this to your existing folder.
# Example from your console:
# dp4_zip_dir <- "/Users/jongewirtzman/Google Drive/Research/tree-flux-2025/filesToStack00200"
dp4_zip_dir <- file.path(out_dir, "dp4_zips")

dp4_unzip_dir <- file.path(out_dir, "dp4_unzipped")
dir_create(dp4_zip_dir)
dir_create(dp4_unzip_dir)

# Download provisional? Set TRUE only if you accept provisional data mixed in.
include_provisional <- TRUE

# ---------------------------
# Helpers
# ---------------------------
standardize_time <- function(dt) {
  dt <- as.data.table(dt)
  if ("startDateTime" %in% names(dt)) {
    dt[, ts := lubridate::ymd_hms(startDateTime, tz = "UTC")]
  } else if ("timeBgn" %in% names(dt)) {
    dt[, ts := lubridate::ymd_hms(timeBgn, tz = "UTC")]
  } else if ("startDate" %in% names(dt)) {
    dt[, ts := lubridate::ymd_hms(startDate, tz = "UTC")]
  } else if ("datetime" %in% names(dt)) {
    dt[, ts := lubridate::ymd_hms(datetime, tz = "UTC")]
  } else {
    stop("No recognizable timestamp column found.")
  }
  dt[]
}

pick_numeric_col <- function(dt, prefer_regex, drop_regex = "qf|flag|uncert|id|index|position|elev|lat|lon") {
  dt <- as.data.table(dt)
  num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  if (length(num_cols) == 0) return(NA_character_)
  
  drop <- names(dt)[stringr::str_detect(names(dt), stringr::regex(drop_regex, ignore_case = TRUE))]
  candidates <- setdiff(num_cols, drop)
  if (length(candidates) == 0) candidates <- num_cols
  
  hit <- candidates[stringr::str_detect(candidates, stringr::regex(prefer_regex, ignore_case = TRUE))]
  if (length(hit) > 0) return(hit[1])
  
  candidates[1]
}

locf_join <- function(x_30m, y_sparse, y_value_col, new_name) {
  x <- data.table::copy(x_30m)[, .(ts)]
  y <- data.table::copy(y_sparse)[, .(ts, val = get(y_value_col))]
  data.table::setorder(x, ts)
  data.table::setorder(y, ts)
  data.table::setkey(x, ts)
  data.table::setkey(y, ts)
  out <- y[x, roll = TRUE]
  out <- out[, .(ts, val)]
  data.table::setnames(out, "val", new_name)
  out
}

safe_loadByProduct <- function(dpID, site, startdate, enddate, package = "basic", check.size = FALSE) {
  message("\n--- Downloading ", dpID, " ---")
  tryCatch(
    neonUtilities::loadByProduct(
      dpID = dpID, site = site,
      startdate = startdate, enddate = enddate,
      package = package,
      check.size = check.size
    ),
    error = function(e) {
      message("Download failed for ", dpID, ": ", e$message)
      NULL
    }
  )
}

stack_30min_tables <- function(loadByProduct_output) {
  st <- tryCatch(neonUtilities::stackByTable(loadByProduct_output), error = function(e) NULL)
  if (is.null(st) || length(st) == 0) return(list())
  nm <- names(st)
  keep <- stringr::str_detect(nm, "30m|30min|30Min|_30")
  st[keep]
}

make_dp1_series <- function(load_obj, series_name, value_regex) {
  tabs <- stack_30min_tables(load_obj)
  if (length(tabs) == 0) return(NULL)
  
  sizes <- vapply(tabs, nrow, numeric(1))
  dt <- standardize_time(tabs[[which.max(sizes)]])
  
  val_col <- pick_numeric_col(dt, value_regex)
  if (is.na(val_col)) return(NULL)
  
  out <- dt[, .(ts, value = mean(get(val_col), na.rm = TRUE)), by = .(ts)]
  setnames(out, "value", series_name)
  out
}

soil_wide_by_depth <- function(load_obj, var_prefix, value_regex) {
  tabs <- stack_30min_tables(load_obj)
  if (length(tabs) == 0) return(NULL)
  
  sizes <- vapply(tabs, nrow, numeric(1))
  dt <- standardize_time(tabs[[which.max(sizes)]])
  
  depth_col <- names(dt)[stringr::str_detect(names(dt), stringr::regex("depth", ignore_case = TRUE))][1]
  if (is.na(depth_col)) {
    stop("No depth column found for soil variable: ", var_prefix,
         ". Inspect names(dt) for the correct depth field.")
  }
  
  val_col <- pick_numeric_col(dt, value_regex)
  if (is.na(val_col)) return(NULL)
  
  # Preserve depth labels exactly (character) to keep depths distinct
  dt[, depth_norm := as.character(get(depth_col))]
  dt[, depth_norm := stringr::str_replace_all(depth_norm, "\\s+", "")]
  dt[, depth_label := paste0("d", depth_norm)]
  
  tmp <- dt[, .(value = mean(get(val_col), na.rm = TRUE)), by = .(ts, depth_label)]
  tmp[, colname := paste0(var_prefix, "_", depth_label)]
  data.table::dcast(tmp, ts ~ colname, value.var = "value", fun.aggregate = mean)
}

# ---------------------------
# DP4 HDF5 helpers
# ---------------------------
h5_has_site_group <- function(h5_file, site_code) {
  root <- rhdf5::h5ls(h5_file, recursive = FALSE)
  any(root$name == site_code)
}

h5_list_tables <- function(h5_file, group_path) {
  lst <- tryCatch(rhdf5::h5ls(h5_file, group_path, recursive = FALSE), error = function(e) NULL)
  if (is.null(lst) || nrow(lst) == 0) return(character(0))
  lst$name
}

h5_read_table <- function(h5_file, group_path, table_name) {
  tryCatch(rhdf5::h5read(h5_file, paste0(group_path, "/", table_name)), error = function(e) NULL)
}

extract_dp4_series <- function(h5_file, site_code, var_key,
                               data_group_suffix,
                               qf_group_suffix,
                               value_prefer_regex) {
  base <- if (h5_has_site_group(h5_file, site_code)) paste0("/", site_code) else ""
  data_grp <- paste0(base, data_group_suffix)
  qf_grp   <- paste0(base, qf_group_suffix)
  
  tables <- h5_list_tables(h5_file, data_grp)
  if (length(tables) == 0) return(NULL)
  
  keep <- stringr::str_detect(tables, "30m|30min|30Min|_30")
  if (any(keep)) tables <- tables[keep]
  
  out_list <- vector("list", length(tables))
  
  for (i in seq_along(tables)) {
    tn <- tables[i]
    d <- h5_read_table(h5_file, data_grp, tn)
    if (is.null(d)) next
    dt <- standardize_time(as.data.table(d))
    
    val_col <- pick_numeric_col(dt, value_prefer_regex)
    if (is.na(val_col)) next
    
    q <- h5_read_table(h5_file, qf_grp, tn)
    qfFinl <- rep(NA_integer_, nrow(dt))
    if (!is.null(q)) {
      qdt <- as.data.table(q)
      if ("qfFinl" %in% names(qdt)) {
        qfFinl <- qdt[["qfFinl"]]
      } else {
        cand <- names(qdt)[stringr::str_detect(names(qdt), stringr::regex("qfFinl", ignore_case = TRUE))]
        if (length(cand) > 0) qfFinl <- qdt[[cand[1]]]
      }
    }
    
    res <- dt[, .(ts, value = get(val_col))]
    res[, (paste0(var_key, "_qfFinl")) := qfFinl]
    setnames(res, "value", var_key)
    
    out_list[[i]] <- res
  }
  
  out <- rbindlist(out_list, fill = TRUE)
  if (nrow(out) == 0) return(NULL)
  setorder(out, ts)
  out <- out[!duplicated(ts)]
  out
}

# ============================================================
# A) DP4.00200.001: download (optional), unzip, find HDF5, extract fluxes
# ============================================================

message("\n=== DP4.00200.001 download/unzip/extract ===")

# OPTIONAL DOWNLOAD:
# If you already downloaded DP4 ZIPs to dp4_zip_dir, you can comment this out.
# This call typically downloads into a "filesToStack..." folder if savepath is not set.
# We explicitly set savepath=dp4_zip_dir to keep it predictable.
tryCatch({
  neonUtilities::zipsByProduct(
    dpID = "DP4.00200.001",
    site = site,
    startdate = start_month,
    enddate = end_month,
    package = "basic",
    include.provisional = include_provisional,
    check.size = TRUE,
    savepath = dp4_zip_dir
  )
}, error = function(e) {
  message("zipsByProduct() message/error (continuing using whatever ZIPs are on disk): ", e$message)
})

# Source of truth: look for ZIPs on disk
zip_files <- fs::dir_ls(dp4_zip_dir, recurse = TRUE, regexp = "\\.zip$")

# If none found in dp4_zip_dir, try common neonUtilities default "filesToStack*" folder in cwd
if (length(zip_files) == 0) {
  candidates <- fs::dir_ls(getwd(), type = "directory", regexp = "filesToStack", recurse = FALSE)
  if (length(candidates) > 0) {
    zip_files <- fs::dir_ls(candidates[1], recurse = TRUE, regexp = "\\.zip$")
    if (length(zip_files) > 0) {
      message("Using ZIPs found under: ", candidates[1])
      dp4_zip_dir <- candidates[1]
    }
  }
}

if (length(zip_files) == 0) {
  stop("No DP4 ZIP files found. Set dp4_zip_dir to the folder containing your downloaded DP4 zip files.")
}

message("Found ", length(zip_files), " DP4 ZIP file(s). Unzipping to: ", dp4_unzip_dir)
for (zf in zip_files) unzip(zf, exdir = dp4_unzip_dir)

# Find HDF5 files
h5_files <- fs::dir_ls(dp4_unzip_dir, recurse = TRUE, regexp = "\\.h5$")

if (length(h5_files) == 0) {
  stop("No .h5 files found after unzip in: ", dp4_unzip_dir,
       "\nIf ZIPs contain .gz, you must gunzip first (see your unzip.neon.R utility).")
}

message("Found ", length(h5_files), " HDF5 file(s). Extracting DP4 fluxes...")

dp4_flux <- rbindlist(lapply(h5_files, function(f) {
  flux_co2 <- extract_dp4_series(
    h5_file = f, site_code = site,
    var_key = "fluxCo2",
    data_group_suffix = "/dp04/data/fluxCo2",
    qf_group_suffix   = "/dp04/qfqm/fluxCo2",
    value_prefer_regex = "fluxCo2|co2.*flux|nee|nsae"
  )
  flux_le <- extract_dp4_series(
    h5_file = f, site_code = site,
    var_key = "fluxH2o",
    data_group_suffix = "/dp04/data/fluxH2o",
    qf_group_suffix   = "/dp04/qfqm/fluxH2o",
    value_prefer_regex = "fluxH2o|latent|LE"
  )
  flux_h <- extract_dp4_series(
    h5_file = f, site_code = site,
    var_key = "fluxTemp",
    data_group_suffix = "/dp04/data/fluxTemp",
    qf_group_suffix   = "/dp04/qfqm/fluxTemp",
    value_prefer_regex = "fluxTemp|sensible|H"
  )
  
  Reduce(function(a, b) merge(a, b, by = "ts", all = TRUE),
         Filter(Negate(is.null), list(flux_co2, flux_le, flux_h)))
}), fill = TRUE)

if (is.null(dp4_flux) || nrow(dp4_flux) == 0) stop("No DP4 flux records extracted from HDF5.")

setorder(dp4_flux, ts)
dp4_flux <- dp4_flux[!duplicated(ts)]

# ============================================================
# B) DP1 products: download + stack 30-min
# ============================================================

message("\n=== DP1 download/stack ===")

dp1_list <- list(
  baroPress = "DP1.00004.001",  # barometric pressure
  irBioTemp = "DP1.00005.001",  # IR biological temperature
  rad_sw_lw = "DP1.00023.001",  # SW/LW (net radiometer product; includes in/out)
  par       = "DP1.00024.001",  # PAR
  parLine   = "DP1.00066.001",  # PAR quantum line
  relHum    = "DP1.00098.001",  # relative humidity (and usually Tair/dewpoint)
  precipDep = "DP1.00006.001",  # deprecated precipitation (older period)
  precipWgt = "DP1.00044.001",  # weighing gauge (hourly/daily)
  precipTip = "DP1.00045.001",  # tipping bucket
  precipThr = "DP1.00046.001",  # throughfall
  soilTemp  = "DP1.00041.001",  # soil temperature
  soilVWC   = "DP1.00094.001",  # soil water content
  soilHFP   = "DP1.00040.001"   # soil heat flux plate
)

dp1_downloads <- lapply(dp1_list, function(dpid) safe_loadByProduct(dpid, site, start_month, end_month))

baro   <- make_dp1_series(dp1_downloads$baroPress, "baroPress", "press|baro|station")
irbio  <- make_dp1_series(dp1_downloads$irBioTemp, "irBioTemp", "temp|ir|bio|surface")
par    <- make_dp1_series(dp1_downloads$par,       "PAR",       "par|ppfd|photosynth")
parln  <- make_dp1_series(dp1_downloads$parLine,   "PAR_line",  "par|ppfd|quantum")

prec_tip <- make_dp1_series(dp1_downloads$precipTip, "precip_tip", "precip|ppt|rain")
prec_thr <- make_dp1_series(dp1_downloads$precipThr, "precip_throughfall", "precip|ppt|rain")
prec_dep <- make_dp1_series(dp1_downloads$precipDep, "precip_deprecated", "precip|ppt|rain")

# Weighing gauge (may be hourly/daily): stack any available table, LOCF later
precip_wgt <- NULL
wgt_stacked <- tryCatch(neonUtilities::stackByTable(dp1_downloads$precipWgt), error = function(e) list())
if (length(wgt_stacked) > 0) {
  sizes <- vapply(wgt_stacked, nrow, numeric(1))
  wdt <- standardize_time(wgt_stacked[[which.max(sizes)]])
  vcol <- pick_numeric_col(wdt, "precip|ppt|rain")
  if (!is.na(vcol)) {
    precip_wgt <- wdt[, .(ts, precip_wgt_raw = get(vcol))]
    precip_wgt <- precip_wgt[, .(precip_wgt_raw = mean(precip_wgt_raw, na.rm = TRUE)), by = .(ts)]
  }
}

# Radiation: try to extract incoming SW/LW explicitly
rad <- NULL
rad_tabs <- stack_30min_tables(dp1_downloads$rad_sw_lw)
if (length(rad_tabs) > 0) {
  sizes <- vapply(rad_tabs, nrow, numeric(1))
  rdt <- standardize_time(rad_tabs[[which.max(sizes)]])
  
  sw_in_col <- names(rdt)[stringr::str_detect(names(rdt), stringr::regex("shortwave.*in|sw.*in|incoming.*shortwave|in.*sw", ignore_case = TRUE))][1]
  lw_in_col <- names(rdt)[stringr::str_detect(names(rdt), stringr::regex("longwave.*in|lw.*in|incoming.*longwave|in.*lw", ignore_case = TRUE))][1]
  
  rad <- rdt[, .(ts)]
  rad[, SW_in := if (!is.na(sw_in_col)) rdt[[sw_in_col]] else NA_real_]
  rad[, LW_in := if (!is.na(lw_in_col)) rdt[[lw_in_col]] else NA_real_]
  rad <- rad[, .(SW_in = mean(SW_in, na.rm = TRUE),
                 LW_in = mean(LW_in, na.rm = TRUE)), by = .(ts)]
}

# RH (+ temperature) to compute VPD
relhum <- NULL
rh_tabs <- stack_30min_tables(dp1_downloads$relHum)
if (length(rh_tabs) > 0) {
  sizes <- vapply(rh_tabs, nrow, numeric(1))
  hdt <- standardize_time(rh_tabs[[which.max(sizes)]])
  
  rh_col <- names(hdt)[stringr::str_detect(names(hdt), stringr::regex("^rh$|rel.*hum|relative.*humidity", ignore_case = TRUE))][1]
  t_col  <- names(hdt)[stringr::str_detect(names(hdt), stringr::regex("temp|air.*temp|temperature", ignore_case = TRUE))][1]
  
  relhum <- hdt[, .(ts)]
  relhum[, RH := if (!is.na(rh_col)) hdt[[rh_col]] else NA_real_]
  relhum[, Tair := if (!is.na(t_col)) hdt[[t_col]] else NA_real_]
  relhum <- relhum[, .(RH = mean(RH, na.rm = TRUE),
                       Tair = mean(Tair, na.rm = TRUE)), by = .(ts)]
}

# Soils (depth distinct)
soilT <- soil_wide_by_depth(dp1_downloads$soilTemp, "soilTemp", "soil.*temp|temperature")
soilV <- soil_wide_by_depth(dp1_downloads$soilVWC,  "soilVWC",  "water.*content|vwc|moist|soil.*water")
soilH <- soil_wide_by_depth(dp1_downloads$soilHFP,  "soilHFP",  "heat.*flux|soil.*heat|g.*soil")

# ============================================================
# C) Merge everything on ts
# ============================================================

message("\n=== Merge ===")
master <- dp4_flux

to_merge <- Filter(Negate(is.null), list(
  baro, irbio, rad, par, parln, relhum,
  prec_tip, prec_thr, prec_dep,
  soilT, soilV, soilH
))
for (obj in to_merge) master <- merge(master, obj, by = "ts", all = TRUE)

if (!is.null(precip_wgt) && nrow(precip_wgt) > 0) {
  wlocf <- locf_join(master, precip_wgt, "precip_wgt_raw", "precip_weighingGauge_LOCF")
  master <- merge(master, wlocf, by = "ts", all = TRUE)
}

setorder(master, ts)

# ============================================================
# D) Derived variables: ET and VPD
# ============================================================

lambda <- 2.45e6
master[, ET_mm_30min := (fluxH2o * 1800) / lambda]

es_kPa <- function(Tc) 0.6108 * exp((17.27 * Tc) / (Tc + 237.3))
if (all(c("RH", "Tair") %in% names(master))) {
  master[, VPD_kPa := pmax(es_kPa(Tair) * (1 - RH / 100), 0)]
} else {
  master[, VPD_kPa := NA_real_]
}

# ============================================================
# E) Write outputs
# ============================================================

csv_all <- file.path(out_dir, sprintf("%s_30min_clean_ALL_%s_%s.csv", site, start_month, end_month))
pq_all  <- file.path(out_dir, sprintf("%s_30min_clean_ALL_%s_%s.parquet", site, start_month, end_month))
fwrite(master, csv_all)
arrow::write_parquet(master, pq_all)

qf_cols <- intersect(c("fluxCo2_qfFinl", "fluxH2o_qfFinl", "fluxTemp_qfFinl"), names(master))
qf_ok <- copy(master)
for (qc in qf_cols) qf_ok <- qf_ok[is.na(get(qc)) | get(qc) == 0]

csv_qf <- file.path(out_dir, sprintf("%s_30min_clean_QF0_%s_%s.csv", site, start_month, end_month))
pq_qf  <- file.path(out_dir, sprintf("%s_30min_clean_QF0_%s_%s.parquet", site, start_month, end_month))
fwrite(qf_ok, csv_qf)
arrow::write_parquet(qf_ok, pq_qf)

message("\nWrote ALL:\n  ", normalizePath(csv_all), "\n  ", normalizePath(pq_all))
message("\nWrote QF0:\n  ", normalizePath(csv_qf), "\n  ", normalizePath(pq_qf))

message("\nSanity:")
message("  Rows (ALL): ", nrow(master))
message("  Rows (QF0): ", nrow(qf_ok))
message("  Soil depth cols example: ",
        paste(head(names(master)[stringr::str_detect(names(master), "^soilTemp_d|^soilVWC_d|^soilHFP_d")], 12), collapse = ", "))

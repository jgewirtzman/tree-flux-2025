library(neonUtilities)

save_dir <- "data/raw/NEON/HARV_flux_2023_2025"
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

zipsByProduct(
  dpID      = "DP4.00200.001",
  site      = "HARV",
  startdate = "2023-01",
  enddate   = "2025-12",
  package   = "basic",
  savepath  = save_dir,
  check.size = FALSE,          # avoids the y/n prompt in scripts
  include.provisional = TRUE   # optional; increases size
)

stack_path <- file.path(save_dir, "filesToStack00200")

flux_dp04 <- stackEddy(filepath = stack_path, level = "dp04")  # all fluxes
gas_dp01_30 <- stackEddy(filepath = stack_path, level = "dp01", avg = 30)# 30-min aggregated obs incl. gas conc

head(gas_dp01_30$HARV)
names(gas_dp01_30)






harv01 <- gas_dp01_30$HARV

# Keep timestamps + indices
base_cols <- c("horizontalPosition","verticalPosition","timeBgn","timeEnd")

# Keep 30-min concentration metrics (means only) for CO2, CH4, H2O
# (Dry and wet mole ratios are usually the primary “concentration” outputs.)
gas_mean_cols <- grep(
  "^(data\\.(co2|ch4|h2o)(Conc|Stor)\\.(rtioMoleDry|rtioMoleWet).+\\.mean)$",
  names(harv01),
  value = TRUE
)

# Optional: include pressure/temp/RH “environment hut” context for QA/troubleshooting
env_cols <- grep(
  "^(data\\.(co2|ch4|h2o)(Conc|Stor)\\.(pres|temp|rh).*\\.mean)$",
  names(harv01),
  value = TRUE
)

# Optional: keep final quality flags for the selected gas fields
qf_cols <- sub("^data\\.", "qfqm.", gas_mean_cols)
qf_cols <- qf_cols[qf_cols %in% names(harv01)]

harv_gas_30min <- harv01[, unique(c(base_cols, gas_mean_cols, env_cols, qf_cols))]

# Inspect
names(harv_gas_30min)
head(harv_gas_30min)

saveRDS(flux_dp04$HARV, file.path(save_dir, "HARV_dp04_flux_2023_2025.rds"))
saveRDS(harv_gas_30min, file.path(save_dir, "HARV_dp01_30min_gas_conc_2023_2025.rds"))



library(neonUtilities)

dp_list <- c(
  "DP1.00001.001", # 2D wind speed and direction
  "DP1.00004.001", # barometric pressure
  "DP1.00005.001", # IR biological temperature
  "DP1.00024.001", # PAR
  "DP1.00066.001", # PAR quantum line
  #"DP1.00044.001", # precipitation – weighing gauge
  "DP1.00098.001", # relative humidity
  "DP1.00041.001", # soil temperature
  "DP1.00094.001", # soil water content & salinity
  "DP1.00003.001", # triple aspirated air temperature
  "DP1.00002.001", # single aspirated air temperature
  #"DP1.00022.001", # shortwave radiation (primary pyranometer)
  "DP1.00014.001", # shortwave radiation (direct & diffuse)
  "DP1.00023.001"  # net radiometer (SW + LW)
)

base_dir <- "data/raw/NEON/HARV_met_soil_2023_2025_30min"
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

for (dp in dp_list) {
  dp_dir <- file.path(base_dir, dp)
  dir.create(dp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Skip if already downloaded (any filesToStack* folder exists)
  fts <- list.files(dp_dir, pattern = "^filesToStack", full.names = TRUE)
  if (length(fts) > 0) {
    message("Skipping ", dp, " (already has ", basename(fts[1]), ")")
    next
  }
  
  message("Downloading ", dp)
  
  zipsByProduct(
    dpID = dp,
    site = "HARV",
    startdate = "2023-01",
    enddate   = "2025-12",
    package   = "basic",
    release   = "current",
    timeIndex = 30,
    include.provisional = TRUE,
    check.size = FALSE,
    savepath  = dp_dir
  )
}



library(neonUtilities)

base_dir <- "data/raw/NEON/HARV_met_soil_2023_2025_30min"

dp_list <- c(
  "DP1.00001.001",
  "DP1.00004.001",
  "DP1.00005.001",
  "DP1.00024.001",
  "DP1.00066.001",
  "DP1.00098.001",
  "DP1.00041.001",
  "DP1.00094.001",
  "DP1.00003.001",
  "DP1.00002.001",
  "DP1.00014.001",
  "DP1.00023.001"
)

for (dp in dp_list) {
  dp_dir <- file.path(base_dir, dp)
  if (!dir.exists(dp_dir)) {
    message("Skipping ", dp, " (missing folder)")
    next
  }
  
  fts <- list.files(dp_dir, pattern = "^filesToStack", full.names = TRUE)
  if (length(fts) == 0) {
    message("Skipping ", dp, " (no filesToStack* found)")
    next
  }
  stack_path <- fts[order(file.info(fts)$mtime, decreasing = TRUE)][1]
  
  stacked_dir <- file.path(dp_dir, "stackedFiles")
  if (dir.exists(stacked_dir) && length(list.files(stacked_dir, recursive = TRUE)) > 0) {
    message("Skipping ", dp, " (already stacked)")
    next
  }
  
  message("Stacking ", dp, " from ", basename(stack_path), " ...")
  
  dir.create(stacked_dir, showWarnings = FALSE)
  
  tryCatch(
    {
      stackByTable(
        filepath = stack_path,
        savepath = dp_dir,
        folder   = "stackedFiles",
        saveUnzippedFiles = FALSE,
        nCores = max(1, parallel::detectCores() - 1),
        useFasttime = TRUE
      )
      message("Done stacking ", dp)
    },
    error = function(e) {
      message("ERROR stacking ", dp, ": ", conditionMessage(e))
    }
  )
}



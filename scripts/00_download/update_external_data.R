# ============================================================
# update_external_data.R  --  run interactively in R/RStudio from the project root
#
# Brings every external driver dataset up to its newest release:
#   1. AmeriFlux BASE: US-Ha1, US-Ha2 (new versions with 2025 data), US-xHA
#      -> data/raw/ameriflux/  (superseded versions moved to data/raw/ameriflux/_superseded/)
#      Needs your AmeriFlux account; running this agrees to the CC-BY-4.0 data policy.
#   2. Harvard Forest archive: Fisher met (HF001) and wetland hydrology (HF070), newest
#      revisions -> data/processed/wtd_met.csv (sources in data/processed/wtd_met_sources.txt)
#   3. NEON HARV (soil moisture, soil temperature, canopy temperature, throughfall, soil
#      heat flux, wind, eddy-covariance bundle), RELEASE-2026 + provisional
#      -> data/raw/NEON_2026/  (release of each month in release_log.csv)
#
# Usage:
#   setwd("~/My Drive/Research/tree-flux-2025")
#   source("scripts/00_download/update_external_data.R")
# Then tell Claude it is done; the import, alignment and analysis pipeline is rerun from there.
# ============================================================
stopifnot("Run from the project root" = file.exists("scripts/run_pipeline.sh"))
for (p in c("amerifluxr", "neonUtilities", "tidyverse", "lubridate", "plantecophys"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
options(timeout = 3600)

# ------------------------------------------------------------
# 1. AmeriFlux
# ------------------------------------------------------------
message("\n========== 1. AmeriFlux ==========")
AMF_DIR <- "data/raw/ameriflux"
dir.create(file.path(AMF_DIR, "_superseded"), recursive = TRUE, showWarnings = FALSE)
amf_user  <- Sys.getenv("AMERIFLUX_USER");  if (amf_user == "")  amf_user  <- readline("AmeriFlux username: ")
amf_email <- Sys.getenv("AMERIFLUX_EMAIL"); if (amf_email == "") amf_email <- readline("AmeriFlux account email: ")
version_of <- function(x) sub(".*_BASE-BADM_([0-9]+-[0-9]+).*", "\\1", x)
newer <- function(a, b) {  # is version a newer than b? ("27-5" vs "26-5")
  pa <- as.integer(strsplit(a, "-")[[1]]); pb <- as.integer(strsplit(b, "-")[[1]])
  pa[1] > pb[1] || (pa[1] == pb[1] && pa[2] > pb[2])
}
for (site in c("US-Ha1", "US-Ha2", "US-xHA")) {
  message("\n--- ", site, " ---")
  old <- list.dirs(AMF_DIR, recursive = FALSE)
  old <- old[grepl(paste0("^AMF_", site, "_BASE-BADM_"), basename(old))]
  tmp <- file.path(tempdir(), site); dir.create(tmp, showWarnings = FALSE)
  zip <- tryCatch(amerifluxr::amf_download_base(
    user_id = amf_user, user_email = amf_email, site_id = site,
    data_product = "BASE-BADM", data_policy = "CCBY4.0", agree_policy = TRUE,
    intended_use = "other_research", intended_use_text = "Tree stem CH4 flux drivers at Harvard Forest",
    out_dir = tmp, verbose = TRUE), error = function(e) { message("  download failed: ", conditionMessage(e)); NULL })
  if (is.null(zip) || !file.exists(zip[1])) next
  new_v <- version_of(basename(zip[1]))
  if (length(old) && !newer(new_v, version_of(basename(old[1])))) {
    message("  no newer version (have ", version_of(basename(old[1])), ", server ", new_v, ")"); next
  }
  dest <- file.path(AMF_DIR, sub("\\.zip$", "", basename(zip[1])))
  top <- unique(sub("/.*", "", unzip(zip[1], list = TRUE)$Name))
  if (length(top) == 1 && grepl("^AMF_", top) && !grepl("\\.", top)) unzip(zip[1], exdir = AMF_DIR) else unzip(zip[1], exdir = dest)
  for (o in old) file.rename(o, file.path(AMF_DIR, "_superseded", basename(o)))
  csv <- list.files(dest, pattern = "_BASE_.*\\.csv$", full.names = TRUE)
  last <- if (length(csv)) tail(read.csv(csv[1], skip = 2, colClasses = "character")$TIMESTAMP_END, 1) else NA
  message("  installed ", basename(dest), "; last record ", last)
}

# ------------------------------------------------------------
# 2. Harvard Forest archive (Fisher met, wetland hydrology)
# ------------------------------------------------------------
message("\n========== 2. Harvard Forest archive ==========")
source("scripts/00_download/03_download_hf_met_hydro.R", local = new.env())
w <- read.csv("data/processed/wtd_met.csv")
message("  wtd_met.csv: ", min(w$datetime), " to ", max(w$datetime))

# ------------------------------------------------------------
# 3. NEON HARV
# ------------------------------------------------------------
message("\n========== 3. NEON ==========")
source("scripts/00_download/04_download_neon_harv.R", local = new.env())

message("\nAll done. Coverage summary:")
for (f in list.files(AMF_DIR, pattern = "_BASE_.*\\.csv$", recursive = TRUE, full.names = TRUE)) {
  if (grepl("_superseded", f)) next
  message("  ", basename(f), ": last record ", tail(read.csv(f, skip = 2, colClasses = "character")$TIMESTAMP_END, 1))
}

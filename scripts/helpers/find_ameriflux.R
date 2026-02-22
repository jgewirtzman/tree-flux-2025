# ============================================================
# find_ameriflux.R
# Helper: locate AmeriFlux BASE CSV files by site ID
#
# Usage:
#   source("scripts/helpers/find_ameriflux.R")
#   Ha1 <- read.csv(find_ameriflux("US-Ha1"), header = TRUE, skip = 2)
# ============================================================

AMERIFLUX_DIR <- "data/raw/ameriflux"

#' Find the AmeriFlux BASE CSV file for a given site
#'
#' Searches data/raw/ameriflux/ for the folder and CSV matching a site ID,
#' regardless of version number. Returns the path to the CSV.
#'
#' @param site_id Character, e.g. "US-Ha1"
#' @param base_dir Directory containing AmeriFlux downloads
#' @return Character path to the CSV file
find_ameriflux <- function(site_id, base_dir = AMERIFLUX_DIR) {
  pattern <- paste0("AMF_", site_id, "_BASE_*")
  csvs <- Sys.glob(file.path(base_dir, paste0("AMF_", site_id, "_BASE-BADM_*"),
                              paste0("AMF_", site_id, "_BASE_*.csv")))
  if (length(csvs) == 0) {
    stop("No AmeriFlux CSV found for site ", site_id,
         " in ", base_dir,
         "\n  Expected pattern: AMF_", site_id, "_BASE-BADM_*/AMF_", site_id, "_BASE_*.csv",
         "\n  Run scripts/00_download/01_download_ameriflux.R first.")
  }
  if (length(csvs) > 1) {
    message("Multiple AmeriFlux versions found for ", site_id, ", using newest: ", basename(csvs[length(csvs)]))
  }
  csvs[length(csvs)]
}

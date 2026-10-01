# ============================================================
# 04_download_neon_harv.R
# NEON Harvard Forest (HARV) products for 2023-01 to 2025-12, the latest
# release (RELEASE-2026, through 2025-06) plus provisional data after it.
#
# Products:
#   DP1.00094.001  soil water content           (NEON_SWC_* drivers)
#   DP1.00041.001  soil temperature             (cross-check for tower TS)
#   DP1.00005.001  IR biological (canopy) temp  (replaces AmeriFlux T_CANOPY_xHA)
#   DP1.00046.001  precipitation - throughfall  (replaces THROUGHFALL_xHA)
#   DP1.00040.001  soil heat flux plate         (replaces G_xHA)
#   DP1.00001.001  2D wind speed and direction  (replaces WS/WD_xHA)
#   DP4.00200.001  eddy-covariance bundle       (replaces FC/SC/USTAR/CO2/CH4 _xHA)
#
# The release of every downloaded month is recorded in data/raw/NEON_2026/release_log.csv.
# Output: data/raw/NEON_2026/<dpID>/
# ============================================================
suppressPackageStartupMessages(library(neonUtilities))
# Optional NEON API token (free, from your data.neonscience.org account) raises rate limits:
#   Sys.setenv(NEON_TOKEN = "...")
TOKEN <- Sys.getenv("NEON_TOKEN", unset = NA_character_)
options(timeout = 3600)
OUT <- Sys.getenv("NEON_OUT", unset = "data/raw/NEON_2026")   # download to a local disk first if the project is on a synced drive
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
START <- "2023-01"; END <- "2025-12"

dp1 <- c("DP1.00094.001", "DP1.00041.001", "DP1.00005.001", "DP1.00046.001", "DP1.00040.001", "DP1.00001.001")
for (dp in c(dp1, "DP4.00200.001")) {
  d <- file.path(OUT, dp); dir.create(d, showWarnings = FALSE)
  if (length(list.files(d, pattern = "^filesToStack"))) { message("skip ", dp, " (already downloaded)"); next }
  message("\n=== ", dp, " ===")
  args <- list(dpID = dp, site = "HARV", startdate = START, enddate = END, package = "basic",
               release = "current", include.provisional = TRUE, check.size = FALSE, savepath = d)
  if (!is.na(TOKEN)) args$token <- TOKEN
  if (!dp %in% c("DP4.00200.001", "DP1.00046.001")) args$timeIndex <- 30   # throughfall: no timeIndex in neonUtilities 3.0
  tryCatch(do.call(zipsByProduct, args), error = function(e) message("  FAILED ", dp, ": ", conditionMessage(e)))
}

# release of each downloaded month (from the folder names NEON writes)
rl <- do.call(rbind, lapply(list.files(OUT, pattern = "^DP", full.names = TRUE), function(d) {
  f <- list.files(d, recursive = TRUE, include.dirs = TRUE, pattern = "^NEON\\.D01\\.HARV")
  f <- unique(basename(f[grepl("\\d{4}-\\d{2}\\.(basic|expanded)\\.", f)]))
  if (!length(f)) return(NULL)
  data.frame(product = basename(d), month = sub(".*\\.(\\d{4}-\\d{2})\\..*", "\\1", f),
             release = ifelse(grepl("PROVISIONAL", f), "PROVISIONAL", sub(".*\\.(RELEASE-\\d{4}).*", "\\1", f)))
}))
if (is.null(rl)) stop("No NEON files were downloaded; check the messages above.")
rl <- unique(rl[order(rl$product, rl$month), ])
write.csv(rl, file.path(OUT, "release_log.csv"), row.names = FALSE)
print(table(rl$product, rl$release))

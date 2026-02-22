# ============================================================
# 02_download_phenocam.R
# Download PhenoCam GCC/NDVI data for Harvard Forest EMS
#
# Downloads: harvardems2 camera, deciduous broadleaf ROI
# Output: data/raw/phenocam/
#
# Requires: phenocamr package
# ============================================================

library(phenocamr)

# ============================================================
# CONFIGURATION
# ============================================================

OUT_DIR <- "data/raw/phenocam"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# DOWNLOAD GCC (3-day frequency, includes NDVI)
# ============================================================

message("=== Downloading PhenoCam data for harvardems2 ===")

tryCatch({
  download_phenocam(
    site = "harvardems2$",
    veg_type = "DB",
    roi_id = "1000",
    frequency = "1",
    outlier_detection = FALSE,
    smooth = FALSE,
    out_dir = OUT_DIR
  )
  message("  Download complete")
}, error = function(e) {
  message("  ERROR: ", e$message)
  message("  Trying direct URL download as fallback...")

  # Fallback: direct download from PhenoCam server
  base_url <- "https://phenocam.nau.edu/data/archive/harvardems2/ROI"
  files <- c(
    "harvardems2_DB_1000_1day_v3.csv",
    "harvardems2_DB_1000_ndvi_1day.txt"
  )

  for (f in files) {
    url <- file.path(base_url, f)
    dest <- file.path(OUT_DIR, f)
    tryCatch({
      download.file(url, dest, mode = "wb", quiet = TRUE)
      message("  Downloaded: ", f)
    }, error = function(e2) {
      message("  Could not download: ", f)
    })
  }
})

# ============================================================
# VERIFY
# ============================================================

message("\n=== Checking downloaded files ===")
files <- list.files(OUT_DIR, pattern = "harvardems2", recursive = TRUE)
if (length(files) > 0) {
  for (f in files) message("  Found: ", f)
} else {
  message("  No files found. Download manually from:")
  message("  https://phenocam.nau.edu/webcam/sites/harvardems2/")
}

message("\nDone. Files saved to: ", OUT_DIR)

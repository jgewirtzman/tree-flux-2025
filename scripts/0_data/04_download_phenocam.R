# ============================================================
# 04_download_phenocam.R
# Download the PhenoCam 1-day NDVI summary for the Harvard Forest EMS camera
# (harvardems2, deciduous broadleaf ROI 1000), the file 1_environment/04_hourly_drivers.R
# reads (columns gcc_90 and ndvi_90):
#
#   data/external/phenocam/harvardems2_DB_1000_ndvi_1day.txt
#
# PhenoCam regenerates this file as new images arrive; its header records the creation date.
# The published results used the file created 2026-01-06 (downloaded 2026-01-07). A local
# copy is kept unless PHENOCAM_REFRESH=1, so reruns use the same file.
#
# (phenocamr::download_phenocam() only fetches the GCC summary, <site>_<veg>_<roi>_<n>day.csv,
#  not this NDVI file, so it is not used here.)
# ============================================================

OUT_DIR <- "data/external/phenocam"
FILE <- "harvardems2_DB_1000_ndvi_1day.txt"
URL <- "https://phenocam.nau.edu/data/archive/harvardems2/ROI/harvardems2_DB_1000_ndvi_1day.csv"   # served as .csv; saved under the name the drivers script reads
dest <- file.path(OUT_DIR, FILE)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (file.exists(dest) && Sys.getenv("PHENOCAM_REFRESH") != "1") {
  message("Keeping the local copy: ", dest, " (set PHENOCAM_REFRESH=1 to download again)")
} else {
  message("Downloading ", URL)
  tmp <- tempfile(fileext = ".txt")
  ok <- tryCatch(download.file(URL, tmp, mode = "wb", quiet = TRUE) == 0, error = function(e) FALSE)
  if (!ok || !file.exists(tmp) || file.size(tmp) == 0)
    stop("Could not download ", URL, "\n  Download it by hand from https://phenocam.nau.edu/webcam/sites/harvardems2/ into ", OUT_DIR)
  file.copy(tmp, dest, overwrite = TRUE)
}

# check that the file is what the drivers read, and record its version
hdr <- readLines(dest, n = 30, warn = FALSE)
cols <- strsplit(hdr[!startsWith(hdr, "#")][1], ",")[[1]]
if (!all(c("date", "gcc_90", "ndvi_90") %in% cols))
  stop(dest, " is not a PhenoCam NDVI summary (expected columns date, gcc_90, ndvi_90)")
created <- sub("^# Creation Date: *", "", grep("^# Creation Date:", hdr, value = TRUE))
message("PhenoCam file: ", dest, "; created ", if (length(created)) created else "(unknown)",
        "; downloaded ", format(file.mtime(dest), "%Y-%m-%d"))

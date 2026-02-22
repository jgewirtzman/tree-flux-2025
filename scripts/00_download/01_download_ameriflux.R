# ============================================================
# 01_download_ameriflux.R
# Download AmeriFlux BASE data for Harvard Forest towers
#
# Downloads: US-Ha1 (EMS), US-Ha2 (hemlock), US-xHA (wetland)
# Output: data/raw/ameriflux/
#
# Requires: amerifluxr package and an AmeriFlux account
#   Register at: https://ameriflux-data.lbl.gov/Pages/RequestAccount.aspx
# ============================================================

library(amerifluxr)

# ============================================================
# CONFIGURATION
# ============================================================

# AmeriFlux credentials — set these before running
# You can also set them as environment variables:
#   Sys.setenv(AMERIFLUX_USER = "your_username")
#   Sys.setenv(AMERIFLUX_EMAIL = "your_email@example.com")
USER_ID    <- Sys.getenv("AMERIFLUX_USER", unset = "")
USER_EMAIL <- Sys.getenv("AMERIFLUX_EMAIL", unset = "")

if (USER_ID == "" || USER_EMAIL == "") {
  stop(
    "AmeriFlux credentials not set.\n",
    "Either set environment variables:\n",
    "  Sys.setenv(AMERIFLUX_USER = 'your_username')\n",
    "  Sys.setenv(AMERIFLUX_EMAIL = 'your_email@example.com')\n",
    "Or edit USER_ID and USER_EMAIL in this script.\n",
    "Register at: https://ameriflux-data.lbl.gov/Pages/RequestAccount.aspx"
  )
}

# Sites to download
SITES <- c("US-Ha1", "US-Ha2", "US-xHA")

# Output directory
OUT_DIR <- "data/raw/ameriflux"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# DOWNLOAD
# ============================================================

for (site in SITES) {
  message("\n=== Downloading ", site, " ===")

  tryCatch({
    zip_path <- amf_download_base(
      user_id = USER_ID,
      user_email = USER_EMAIL,
      site_id = site,
      data_product = "BASE-BADM",
      data_policy = "CCBY4.0",
      agree_policy = TRUE,
      intended_use = "other_research",
      intended_use_text = "Tree CH4 flux analysis at Harvard Forest",
      verbose = TRUE,
      out_dir = OUT_DIR
    )

    # Unzip
    if (!is.null(zip_path) && file.exists(zip_path)) {
      unzip(zip_path, exdir = OUT_DIR)
      message("  Extracted to: ", OUT_DIR)
    }

  }, error = function(e) {
    message("  ERROR downloading ", site, ": ", e$message)
    message("  You may need to download manually from https://ameriflux.lbl.gov/")
  })
}

# ============================================================
# VERIFY
# ============================================================

message("\n=== Checking downloaded files ===")

expected_patterns <- c("AMF_US-Ha1_BASE", "AMF_US-Ha2_BASE", "AMF_US-xHA_BASE")
for (pat in expected_patterns) {
  found <- list.files(OUT_DIR, pattern = pat, recursive = TRUE)
  if (length(found) > 0) {
    message("  Found: ", found[1])
  } else {
    message("  MISSING: ", pat)
  }
}

message("\nDone. Files saved to: ", OUT_DIR)

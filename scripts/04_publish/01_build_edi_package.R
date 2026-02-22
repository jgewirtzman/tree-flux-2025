# ============================================================
# 01_build_edi_package.R
# Build EML metadata and prepare data package for EDI
#
# Creates:
#   - Zipped tomography image archives
#   - EML metadata file
#   - Ready-to-upload package in data/edi/package/
#
# Requires: EMLassemblyline, EDIutils
#   install.packages("remotes")
#   remotes::install_github("EDIorg/EMLassemblyline")
#   install.packages("EDIutils")
# ============================================================

library(EMLassemblyline)

# ============================================================
# CONFIGURATION — edit these before publishing
# ============================================================

# Package identifier: get from EDI after reserving
# For staging, use a test scope (e.g., "edi.000.1")
# For production, reserve via EDIutils::create_reservation()
PACKAGE_ID <- "edi.XXXXX.1"  # TODO: replace with reserved ID

# Paths
TEMPLATE_DIR <- "data/edi"
DATA_DIR     <- "data/input"
PACKAGE_DIR  <- "data/edi/package"
dir.create(PACKAGE_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 1: Zip tomography image directories
# ============================================================

message("=== Zipping tomography images ===")

ert_zip <- file.path(PACKAGE_DIR, "tomography_ERTs.zip")
pit_zip <- file.path(PACKAGE_DIR, "tomography_PITs.zip")

# Zip ERTs
ert_files <- list.files(file.path(DATA_DIR, "tomography/ERTs_Absolute"),
                         pattern = "\\.jpg$", full.names = TRUE)
if (length(ert_files) > 0) {
  zip(ert_zip, ert_files, flags = "-j")  # -j: junk directory paths
  message("  Created: ", ert_zip, " (", length(ert_files), " images)")
} else {
  message("  WARNING: No ERT images found")
}

# Zip PITs
pit_files <- list.files(file.path(DATA_DIR, "tomography/CH4_PITs"),
                         pattern = "\\.jpg$", full.names = TRUE)
if (length(pit_files) > 0) {
  zip(pit_zip, pit_files, flags = "-j")
  message("  Created: ", pit_zip, " (", length(pit_files), " images)")
} else {
  message("  WARNING: No PIT images found")
}

# ============================================================
# STEP 2: Copy data tables to package directory
# ============================================================

message("\n=== Copying data tables ===")

data_tables <- c(
  "HF_2023-2025_tree_flux_corrected.csv",
  "tomography_results_compiled.csv",
  "hummock_hollow.csv"
)

for (f in data_tables) {
  src <- file.path(DATA_DIR, f)
  dst <- file.path(PACKAGE_DIR, f)
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    message("  Copied: ", f)
  } else {
    message("  WARNING: ", f, " not found")
  }
}

# ============================================================
# STEP 3: Build EML metadata
# ============================================================

message("\n=== Building EML ===")

make_eml(
  path = TEMPLATE_DIR,
  data.path = PACKAGE_DIR,
  eml.path = PACKAGE_DIR,

  dataset.title = paste(
    "Tree stem methane flux, tomography, and microtopography data",
    "from Harvard Forest upland and wetland sites, 2023-2025"
  ),

  temporal.coverage = c("2023-06-29", "2025-10-31"),

  geographic.description = paste(
    "Harvard Forest, Petersham, Massachusetts, USA.",
    "Prospect Hill tract including the Environmental Measurement Station (EMS)",
    "upland eddy covariance tower footprint and Black Gum Swamp (BGS) wetland."
  ),
  geographic.coordinates = c(
    "42.5396",   # North
    "-72.1715",  # East
    "42.5310",   # South
    "-72.1800"   # West
  ),

  maintenance.description = "Completed. Data collection ended October 2025.",

  # --- Data tables (CSVs with full attribute metadata) ---
  data.table = data_tables,
  data.table.name = c(
    "Tree stem methane flux measurements",
    "Tomography results compiled",
    "Hummock-hollow microtopography classification"
  ),
  data.table.description = c(
    paste(
      "Stem CH4 and CO2 flux measurements from 60 trees (30 wetland, 30 upland)",
      "at Harvard Forest, June 2023 - October 2025. Includes 1,640 quality-controlled",
      "observations with flux rates, R-squared values, standard errors, and metadata.",
      "Fluxes calculated from closed-chamber measurements using a Los Gatos Research",
      "Ultra-Portable Greenhouse Gas Analyzer (UGGA). CH4 flux units corrected to",
      "nmol m-2 s-1. Quality control: CO2 R2 >= 0.8, positive CO2 slope, CH4 flux",
      ">= -1 nmol m-2 s-1."
    ),
    paste(
      "Non-destructive assessment of internal wood condition for 60 trees using",
      "PiCUS sonic tomography (SoT) and electrical resistance tomography (ERT).",
      "Includes mean, median, SD, CV, Gini, and entropy of resistivity; percentage",
      "sound and damaged wood from sonic tomography. Measurements taken at breast",
      "height (1.37 m)."
    ),
    paste(
      "Microtopographic classification (hummock, hollow, or intermediate) for",
      "30 trees at the Black Gum Swamp wetland site. Classification based on",
      "field assessment of tree base position relative to surrounding peat surface."
    )
  ),
  data.table.quote.character = c('"', '"', '"'),

  # --- Other entities (zipped image archives) ---
  other.entity = c("tomography_ERTs.zip", "tomography_PITs.zip"),
  other.entity.name = c(
    "Electrical resistance tomography cross-section images",
    "Sonic tomography cross-section images"
  ),
  other.entity.description = c(
    paste(
      "ZIP archive of electrical resistance tomography (ERT) cross-section images",
      "(JPG format) for study trees. Generated by PiCUS TreeTronic 3 Tomograph.",
      "Color scale represents electrical resistivity: blue = low resistance",
      "(wet/decayed), red = high resistance (dry/cavity). One image per tree."
    ),
    paste(
      "ZIP archive of sonic tomography (SoT) cross-section images",
      "(JPG format) for study trees. Generated by PiCUS Sonic Tomograph.",
      "Color scale represents sound velocity: brown = high velocity (sound wood),",
      "green = intermediate, blue/violet = low velocity (decay/cavity).",
      "One image per tree."
    )
  ),

  # --- User/package info ---
  user.id = "jgewirtzman",
  user.domain = "EDI",
  package.id = PACKAGE_ID,

  write.file = TRUE,
  return.obj = TRUE
)

message("\n=== Done ===")
message("EML written to: ", file.path(PACKAGE_DIR, paste0(PACKAGE_ID, ".xml")))
message("\nPackage contents:")
message("  ", paste(list.files(PACKAGE_DIR), collapse = "\n  "))
message("\nNext steps:")
message("  1. Review the EML file")
message("  2. Run 02_upload_edi.R to upload to staging")
message("  3. Check rendering at https://portal-s.edirepository.org/")
message("  4. When satisfied, upload to production")

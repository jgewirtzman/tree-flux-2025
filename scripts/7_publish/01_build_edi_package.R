# ============================================================
# 01_build_edi_package.R
# Assemble the data package for publication: the study's primary data (data/package/,
# the folder the pipeline reads) plus the final flux dataset, its dictionary and the
# processing log (data/final/). Tables are copied as they are; folders of raw files are
# zipped. Writes EML with EMLassemblyline.
#
# The package is published through the Harvard Forest Data Archive (knb-lter-hfr, on EDI):
# hand data/edi/package/ to the HF archivist, who assigns the package ID and DOI.
# Output: data/edi/package/
# Requires: EMLassemblyline, EDIutils. Run 00_attributes.R first.
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
PKG          <- "data/package"
FINAL        <- "data/final"
PACKAGE_DIR  <- "data/edi/package"
unlink(PACKAGE_DIR, recursive = TRUE); dir.create(PACKAGE_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 1: tables (copied as they are)
# ============================================================
tables <- c(file.path(FINAL, c("stem_ch4_flux.csv", "stem_ch4_flux_dictionary.csv", "flux_processing_log.csv")),
            file.path(PKG, c("trees.csv", "chamber_volumes.csv")),
            file.path(PKG, "tomography", c("tomography_results_compiled.csv", "ERT_application_results.csv", "Tree_ID_info.csv")),
            file.path(PKG, "qc_decisions", "manual_windows.csv"),
            file.path(PKG, "previous_processing", "HF_2023-2025_tree_flux_v1.csv"),
            file.path(PKG, "stand", "BGS_VRP_2025.csv"))
stopifnot(all(file.exists(tables)))
file.copy(tables, PACKAGE_DIR, overwrite = TRUE)
file.copy(file.path(FINAL, "flux_processing_settings.json"), PACKAGE_DIR, overwrite = TRUE)
file.copy(file.path(PKG, "stand", "Black_Gum_Swamp.kmz"), PACKAGE_DIR, overwrite = TRUE)
file.copy(file.path(PKG, "gps", c("BGS_editable.geojson", "EMS_trees_editable.geojson")), PACKAGE_DIR, overwrite = TRUE)

# ============================================================
# STEP 2: folders zipped (paths inside each zip relative to the folder)
# ============================================================
zip_dir <- function(dir, zipname) {
  old <- setwd(dir); on.exit(setwd(old))
  zip(file.path(old, PACKAGE_DIR, zipname), list.files(".", recursive = TRUE), flags = "-q")
  message("  ", zipname, ": ", length(list.files(".", recursive = TRUE)), " files")
}
zip_dir(file.path(PKG, "field_logs"), "field_logs.zip")
zip_dir(file.path(PKG, "analyzer_raw", "lgr"), "analyzer_raw_lgr.zip")
zip_dir(file.path(PKG, "analyzer_raw", "li7810"), "analyzer_raw_li7810.zip")
zip_dir(file.path(PKG, "tomography", "images"), "tomography_images.zip")
data_tables <- basename(tables)
other_entities <- c("field_logs.zip", "analyzer_raw_lgr.zip", "analyzer_raw_li7810.zip", "tomography_images.zip",
                    "flux_processing_settings.json", "Black_Gum_Swamp.kmz")

# ============================================================
# STEP 3: EML  (attribute templates: 00_attributes.R; descriptions below to finish
#               before publication -- TODO)
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
  data.table.name = c("Stem CH4 and CO2 flux, one row per measurement (final dataset)",
                      "Data dictionary of the stem flux dataset", "Flux processing log: every cleaning rule and records affected",
                      "Study trees", "Chamber volumes", "Tomography metrics (compiled)", "ERT application results",
                      "Tomography tree table", "Windows set or excluded after inspection",
                      "Earlier processed flux dataset (source of measurements without a raw record)",
                      "Prism (variable-radius) survey of Black Gum Swamp, June 2025"),
  data.table.description = data_tables,   # TODO: full descriptions
  data.table.quote.character = rep('"', length(data_tables)),
  other.entity = other_entities,
  other.entity.name = c("Field logs", "LGR/UGGA raw 1-Hz records", "LI-7810 raw 1-Hz records", "Tomography images",
                        "Flux processing settings", "Black Gum Swamp outline"),
  other.entity.description = other_entities,   # TODO: full descriptions

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

# ============================================================
# 01_build_edi_package.R
# Assemble the data package for publication: the study's primary data (data/package/,
# the folder the pipeline reads) plus the derived datasets (data/final/): the stem flux
# dataset with its dictionary, processing log and settings, the hourly environmental
# drivers with their dictionary, and the tomography decay classes. Tables are copied as
# they are; folders of raw files are zipped. Writes EML with EMLassemblyline from the
# templates written by 00_attributes.R (text in scripts/7_publish/metadata/).
#
# The package is published through the Harvard Forest Data Archive (knb-lter-hfr, on EDI):
# hand data/edi/package/ to the HF archivist, who assigns the package ID and DOI. The ID
# below is a placeholder; when the ID is issued, set it here and in
# scripts/0_data/02_download_data_package.R.
# Output: data/edi/package/
# Requires: EMLassemblyline. Run 00_attributes.R first.
# ============================================================
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))
library(EMLassemblyline)

PACKAGE_ID   <- "knb-lter-hfr.0.1"   # placeholder until the HF archive issues the ID
TEMPLATE_DIR <- "data/edi"
PKG          <- "data/package"
FINAL        <- "data/final"
META         <- "scripts/7_publish/metadata"
PACKAGE_DIR  <- "data/edi/package"
unlink(PACKAGE_DIR, recursive = TRUE); dir.create(PACKAGE_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 1: tables (copied as they are), in the order of the EML
# ============================================================
tables <- c(
  stem_ch4_flux.csv                 = file.path(FINAL, "stem_ch4_flux.csv"),
  stem_ch4_flux_dictionary.csv      = file.path(FINAL, "stem_ch4_flux_dictionary.csv"),
  flux_processing_log.csv           = file.path(FINAL, "flux_processing_log.csv"),
  environment_hourly.csv            = file.path(FINAL, "environment_hourly.csv"),
  environment_hourly_dictionary.csv = file.path(META, "dictionaries", "environment_hourly.csv"),
  tomography_classes.csv            = file.path(FINAL, "tomography_classes.csv"),
  trees.csv                         = file.path(PKG, "trees.csv"),
  chamber_volumes.csv               = file.path(PKG, "chamber_volumes.csv"),
  tomography_results_compiled.csv   = file.path(PKG, "tomography", "tomography_results_compiled.csv"),
  ERT_application_results.csv       = file.path(PKG, "tomography", "ERT_application_results.csv"),
  Tree_ID_info.csv                  = file.path(PKG, "tomography", "Tree_ID_info.csv"),
  manual_windows.csv                = file.path(PKG, "qc_decisions", "manual_windows.csv"),
  "HF_2023-2025_tree_flux_v1.csv"   = file.path(PKG, "previous_processing", "HF_2023-2025_tree_flux_v1.csv"),
  BGS_VRP_2025.csv                  = file.path(PKG, "stand", "BGS_VRP_2025.csv"))
stopifnot(all(file.exists(tables)))
stopifnot(file.copy(tables, file.path(PACKAGE_DIR, names(tables)), overwrite = TRUE))

table_info <- rbind(
  c("Stem CH4 and CO2 flux, one row per measurement",
    "Final stem flux dataset: 1,637 closed-chamber measurements on 60 trees, June 2023-October 2025, recalculated from the raw 1-Hz analyzer records with goFlux, with minimum detectable flux, detection flags and quality-control flags (fluxqc and trace checks), sampling times, chamber volume and meteorology used in the calculation. Fluxes below the detection limit are retained. Column definitions also in stem_ch4_flux_dictionary.csv."),
  c("Data dictionary of the stem flux dataset", "Column, unit and definition of every column of stem_ch4_flux.csv."),
  c("Flux processing log", "Every cleaning rule applied in building stem_ch4_flux.csv from the raw records and field logs, in order, with the number of records it affected."),
  c("Hourly environmental drivers",
    "Hourly meteorology, water level, soil temperature and moisture, ecosystem fluxes, tower gas mole fractions and canopy phenology at Harvard Forest from January 2023 (EST, UTC-5), compiled from the Harvard Forest Data Archive (HF001, HF070), AmeriFlux (US-Ha1, US-Ha2), NEON (HARV) and PhenoCam (harvardems2). Includes values derived from NEON provisional data (July-December 2025). Source of every column in environment_hourly_dictionary.csv."),
  c("Data dictionary of the hourly environmental drivers", "Column, unit, definition and source dataset of every column of environment_hourly.csv."),
  c("Tomography decay classes", "Decay class of each study tree from sonic (SoT) and electrical resistance (ERT) tomography, with the ERT principal-component score and threshold and the SoT structural loss used to assign it (classification of Thompson et al. 2026)."),
  c("Study trees", "The 60 study trees: ForestGEO tag, site, plot, species, DBH, microtopography (wetland) and GPS location."),
  c("Chamber volumes", "Collar depths (four interior points), collar, cap and tubing volumes of each tree's chamber."),
  c("Tomography metrics (compiled)", "Eight ERT resistivity metrics and SoT solid and damaged percentages for each tomogram."),
  c("ERT image analysis results", "Eight ERT resistivity metrics for each ERT image, as exported by the ERT image analysis application (Tree-Tomography repository)."),
  c("Tomography tree table", "SoT results for each tree of the tomography study: DBH, crack detected, percent solid and damaged wood."),
  c("Fitting windows set or excluded on inspection", "Decision for each of the 55 measurements inspected individually: kept as fitted, refitted over a manually chosen window (start and end given), or excluded."),
  c("Earlier processed flux dataset",
    "Earlier version of the stem flux dataset (linear fits over the logged windows). Read by the processing for curated tree IDs and meteorology, and the source of the fluxes of the 99 measurements without an archived raw record; superseded by stem_ch4_flux.csv for all other measurements."),
  c("Prism survey of Black Gum Swamp, June 2025", "Variable-radius (prism, BAF 10 ft2/acre) survey of 40 systematic points in Black Gum Swamp: species, DBH and status of each tallied stem."))
stopifnot(nrow(table_info) == length(tables))

# ============================================================
# STEP 2: other files; folders zipped (paths inside each zip relative to the folder)
# ============================================================
others <- c(file.path(FINAL, "flux_processing_settings.json"), file.path(PKG, "stand", "Black_Gum_Swamp.kmz"),
            file.path(PKG, "gps", c("BGS_editable.geojson", "EMS_trees_editable.geojson")))
stopifnot(file.copy(others, PACKAGE_DIR, overwrite = TRUE))
zip_dir <- function(dir, zipname, drop = character(0)) {
  old <- setwd(dir); on.exit(setwd(old))
  f <- list.files(".", recursive = TRUE); f <- f[!basename(f) %in% drop]
  zip(file.path(old, PACKAGE_DIR, zipname), f, flags = "-q")
  message("  ", zipname, ": ", length(f), " files")
}
zip_dir(file.path(PKG, "field_logs"), "field_logs.zip")
zip_dir(file.path(PKG, "analyzer_raw", "lgr"), "analyzer_raw_lgr.zip")
zip_dir(file.path(PKG, "analyzer_raw", "li7810"), "analyzer_raw_li7810.zip")
zip_dir(file.path(PKG, "tomography", "images"), "tomography_images.zip",
        drop = c("dummy_file.jpg", "labtest.jpg", "Icon\r", ".DS_Store"))   # not tree tomograms

other_info <- rbind(
  c("field_logs.zip", "Field logs",
    "Field datasheets and logs as recorded (LGR period, June 2023-March 2025): date, tree, analyzer, closure start and end times, notes; and the team's timing corrections for summer 2024 (Timing Updates.xlsx). CSV and Excel files."),
  c("analyzer_raw_lgr.zip", "LGR/UGGA raw 1-Hz records",
    "Raw 1-Hz records of the Los Gatos Research GLA131 analyzers (folders by analyzer unit), June 2023-March 2025, as written by the instruments."),
  c("analyzer_raw_li7810.zip", "LI-7810 raw 1-Hz records",
    "Raw 1-Hz records of the LI-COR LI-7810 analyzer, April-October 2025, with tagged remarks marking each measurement, as written by the instrument."),
  c("tomography_images.zip", "Tomography images",
    "Tomograms of each study tree at breast height: ERTs_Absolute (electrical resistance tomography, absolute resistivity scale) and CH4_PITs (sonic tomography). File names are tree_date (e.g. 412_22vii24 = tree 412, 22 July 2024). JPEG."),
  c("flux_processing_settings.json", "Flux processing settings",
    "All settings of the flux processing (deadband, analyzer volumes and precision settings, screening thresholds), as used to build stem_ch4_flux.csv. JSON."),
  c("Black_Gum_Swamp.kmz", "Black Gum Swamp outline", "Outline of Black Gum Swamp. KMZ (Google Earth)."),
  c("BGS_editable.geojson", "Wetland tree GPS points", "GPS points of the wetland study trees (April 2025). GeoJSON, WGS84."),
  c("EMS_trees_editable.geojson", "Upland tree GPS points", "GPS points of the upland study trees (January 2026). GeoJSON, WGS84."))

# ============================================================
# STEP 3: EML
# ============================================================
message("\n=== Building EML ===")
eml <- make_eml(
  path = TEMPLATE_DIR,
  data.path = PACKAGE_DIR,
  eml.path = PACKAGE_DIR,
  dataset.title = "Stem methane and carbon dioxide fluxes, environmental drivers and internal wood condition of upland and wetland trees at Harvard Forest 2023-2025",
  temporal.coverage = c("2023-06-29", "2025-10-31"),
  maintenance.description = "Completed. Data collection ended October 2025.",
  data.table = names(tables),
  data.table.name = table_info[, 1],
  data.table.description = table_info[, 2],
  data.table.quote.character = rep('"', length(tables)),
  other.entity = other_info[, 1],
  other.entity.name = other_info[, 2],
  other.entity.description = other_info[, 3],
  user.id = "jgewirtzman",
  user.domain = "EDI",
  package.id = PACKAGE_ID,
  write.file = TRUE,
  return.obj = TRUE)

message("\nEML: ", file.path(PACKAGE_DIR, paste0(PACKAGE_ID, ".xml")))
sz <- file.info(list.files(PACKAGE_DIR, full.names = TRUE))$size
message("Package: ", length(sz), " files, ", round(sum(sz) / 1e6), " MB")

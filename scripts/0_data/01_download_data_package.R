# ============================================================
# 01_download_data_package.R
# Download the two published data packages (stem flux; tomography) and unpack them into
# data/package/, the layout the pipeline reads (built for publication by
# scripts/7_publish/01_build_edi_package.R and data_package/ of the Tree-Tomography repository):
#
#   data/package/
#     trees.csv                        study trees: tag, site, species, DBH, microtopography, GPS
#     chamber_volumes.csv              measured collar + cap + tubing volume per tree
#     analyzer_raw/lgr/LGR{1,2,3}/     LGR/UGGA 1-Hz records (from analyzer_raw_lgr.zip)
#     analyzer_raw/li7810/             LI-7810 1-Hz records (from analyzer_raw_li7810.zip)
#     measurements.csv                 one row per measurement: tree, analyzer, closure times, sampling time, met
#     tomography/                      tomography metrics and tomograms (tomography package)
#     stand/                           2025 prism survey of the swamp, swamp outline
#
# The derived datasets (stem_ch4_flux.csv, environment_hourly.csv, tomography_classes.csv and
# their dictionaries and settings) are also in the package but are rebuilt by the pipeline;
# they are saved to data/package/published_final/ for comparison.
# Requires: EDIutils
# ============================================================
library(EDIutils)
PACKAGE_ID <- "knb-lter-hfr.XXX.1"   # stem flux package, Harvard Forest Data Archive (on EDI); TODO: set when issued
TOMO_ID    <- "knb-lter-hfr.YYY.1"   # tomography package; TODO: set when issued
PKG <- "data/package"
if (grepl("XXX|YYY", paste(PACKAGE_ID, TOMO_ID))) stop("Set PACKAGE_ID and TOMO_ID to the published packages first.")
dir.create(PKG, recursive = TRUE, showWarnings = FALSE)

tmp <- file.path(tempdir(), "edi_pkg"); dir.create(tmp, showWarnings = FALSE)
for (id in c(PACKAGE_ID, TOMO_ID)) {
  ent <- read_data_entity_names(id)
  message("Package ", id, ": ", nrow(ent), " entities")
  for (i in seq_len(nrow(ent))) {
    writeBin(read_data_entity(id, ent$entityId[i]), file.path(tmp, ent$entityName[i]))
    message("  ", ent$entityName[i])
  }
}
unz <- function(zip, to) { dir.create(to, recursive = TRUE, showWarnings = FALSE); unzip(zip, exdir = to) }
place <- list(
  "trees.csv" = PKG, "chamber_volumes.csv" = PKG,
  "measurements.csv" = PKG,
  "ERT_application_results.csv" = file.path(PKG, "tomography"),
  "Tree_ID_info.csv" = file.path(PKG, "tomography"),
  "BGS_editable.geojson" = file.path(PKG, "gps"), "EMS_trees_editable.geojson" = file.path(PKG, "gps"),
  "BGS_VRP_2025.csv" = file.path(PKG, "stand"), "Black_Gum_Swamp.kmz" = file.path(PKG, "stand"),
  "stem_ch4_flux.csv" = file.path(PKG, "published_final"),
  "stem_ch4_flux_dictionary.csv" = file.path(PKG, "published_final"),
  "environment_hourly.csv" = file.path(PKG, "published_final"),
  "environment_hourly_dictionary.csv" = file.path(PKG, "published_final"),
  "tomography_classes.csv" = file.path(PKG, "published_final"),
  "flux_processing_settings.json" = file.path(PKG, "published_final"))
zips <- c("analyzer_raw_lgr.zip" = file.path(PKG, "analyzer_raw", "lgr"),
          "analyzer_raw_li7810.zip" = file.path(PKG, "analyzer_raw", "li7810"),
          "tomograms_study_trees.zip" = file.path(PKG, "tomography", "images"))
for (f in names(place)) if (file.exists(file.path(tmp, f))) {
  dir.create(place[[f]], recursive = TRUE, showWarnings = FALSE); file.copy(file.path(tmp, f), place[[f]], overwrite = TRUE) }
for (z in names(zips)) if (file.exists(file.path(tmp, z))) unz(file.path(tmp, z), zips[[z]])
message("\nData package in ", PKG, ". Next: the public data (scripts/0_data/02-04), then scripts/run_pipeline.sh.")

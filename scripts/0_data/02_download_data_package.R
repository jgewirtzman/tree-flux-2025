# ============================================================
# 02_download_data_package.R
# Download the published data package and unpack it into data/package/, the layout
# the pipeline reads (built for publication by scripts/7_publish/01_build_edi_package.R):
#
#   data/package/
#     trees.csv                        study trees: tag, site, species, DBH, microtopography, GPS
#     chamber_volumes.csv              measured collar + cap + tubing volume per tree
#     field_logs/                      field logs as recorded (closure times, analyzer, notes)
#     analyzer_raw/lgr/LGR{1,2,3}/     LGR/UGGA 1-Hz records (from analyzer_raw_lgr.zip)
#     analyzer_raw/li7810/             LI-7810 1-Hz records (from analyzer_raw_li7810.zip)
#     previous_processing/             earlier processed dataset (tree IDs, met, fluxes of
#                                      measurements without a raw record)
#     qc_decisions/manual_windows.csv  windows set or excluded after inspection
#     tomography/                      tomography metrics and images (companion study)
#     stand/                           2025 prism survey of the swamp, swamp outline
#
# The final datasets (stem_ch4_flux.csv etc.) are also in the package but are rebuilt by
# the pipeline; they are saved to data/package/published_final/ for comparison.
# Requires: EDIutils
# ============================================================
library(EDIutils)
PACKAGE_ID <- "edi.XXXXX.1"   # TODO: replace with the published package ID
PKG <- "data/package"
if (grepl("XXXXX", PACKAGE_ID)) stop("Set PACKAGE_ID to the published package first.")
dir.create(PKG, recursive = TRUE, showWarnings = FALSE)

ent <- read_data_entity_names(PACKAGE_ID)
message("Package ", PACKAGE_ID, ": ", nrow(ent), " entities")
tmp <- file.path(tempdir(), "edi_pkg"); dir.create(tmp, showWarnings = FALSE)
for (i in seq_len(nrow(ent))) {
  f <- file.path(tmp, ent$entityName[i])
  writeBin(read_data_entity(PACKAGE_ID, ent$entityId[i]), f)
  message("  ", ent$entityName[i])
}
unz <- function(zip, to) { dir.create(to, recursive = TRUE, showWarnings = FALSE); unzip(zip, exdir = to) }
place <- list(
  "trees.csv" = PKG, "chamber_volumes.csv" = PKG,
  "manual_windows.csv" = file.path(PKG, "qc_decisions"),
  "HF_2023-2025_tree_flux_v1.csv" = file.path(PKG, "previous_processing"),
  "tomography_results_compiled.csv" = file.path(PKG, "tomography"),
  "ERT_application_results.csv" = file.path(PKG, "tomography"),
  "Tree_ID_info.csv" = file.path(PKG, "tomography"),
  "BGS_VRP_2025.csv" = file.path(PKG, "stand"), "Black_Gum_Swamp.kmz" = file.path(PKG, "stand"),
  "stem_ch4_flux.csv" = file.path(PKG, "published_final"),
  "stem_ch4_flux_dictionary.csv" = file.path(PKG, "published_final"),
  "flux_processing_log.csv" = file.path(PKG, "published_final"))
zips <- c("field_logs.zip" = file.path(PKG, "field_logs"),
          "analyzer_raw_lgr.zip" = file.path(PKG, "analyzer_raw", "lgr"),
          "analyzer_raw_li7810.zip" = file.path(PKG, "analyzer_raw", "li7810"),
          "tomography_images.zip" = file.path(PKG, "tomography", "images"))
for (f in names(place)) if (file.exists(file.path(tmp, f))) {
  dir.create(place[[f]], recursive = TRUE, showWarnings = FALSE); file.copy(file.path(tmp, f), place[[f]], overwrite = TRUE) }
for (z in names(zips)) if (file.exists(file.path(tmp, z))) unz(file.path(tmp, z), zips[[z]])
message("\nData package in ", PKG, ". Next: the public data (scripts/0_data/03-05), then scripts/run_pipeline.sh.")

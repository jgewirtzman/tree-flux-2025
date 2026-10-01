# ============================================================
# 01_assemble_data_package.R  (authors only; run once)
#
# Gathers the study's own primary data from the lab's working folders into
# data/package/, the folder that is published as the data package and the only
# place the pipeline reads primary data from. Anyone else gets data/package/
# from the published package instead (02_download_data_package.R).
#
# Nothing here is processed: files are copied as recorded, except trees.csv,
# which is compiled from the GPS survey, the microtopography table and the
# tree list of the earlier processed dataset.
#
# Sources (read-only):
#   Matthes_Lab/stem-CH4-flux/Raw LGR Data/LGR{1,2,3}/   LGR/UGGA 1-Hz day files
#   Matthes_Lab/stem-CH4-flux/raw-7810-data/             LI-7810 1-Hz day files (incl. Aug 2025)
#   (lab copies, moved to data/_old/ in Oct 2026) data/raw/7810_Processed/Tree_Fluxes/raw-7810-tree_flux_data/   LI-7810 day files
#   data/raw/upland_wetland/                             field logs, chamber volumes
#   data/input/                                          earlier processed dataset, tomography,
#                                                        microtopography, manual window decisions,
#                                                        GPS survey, prism survey, swamp outline
#   Tomography/Tree-Tomography/data/                     ERT application results, tree ID key
#
# Set the two lab folders with environment variables if they live elsewhere:
#   HF_LGR_DIR, HF_LI7810_DIR, HF_TOMO_DIR
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(sf) })

LGR_SRC   <- Sys.getenv("HF_LGR_DIR",    "/Users/jongewirtzman/My Drive/Matthes_Lab/stem-CH4-flux/Raw LGR Data")
LI_SRC    <- c("data/_old/raw/7810_Processed/Tree_Fluxes/raw-7810-tree_flux_data",
               Sys.getenv("HF_LI7810_DIR", "/Users/jongewirtzman/My Drive/Matthes_Lab/stem-CH4-flux/raw-7810-data"))
TOMO_SRC  <- Sys.getenv("HF_TOMO_DIR",   "/Users/jongewirtzman/My Drive/Research/Tomography/Tree-Tomography/data")
FL_SRC    <- "data/_old/raw/upland_wetland"
IN_SRC    <- "data/_old/input"
PKG       <- "data/package"
STUDY     <- as.Date(c("2023-05-01", "2025-04-30"))   # LGR period (LI-7810 from Apr 2025)

cp <- function(from, to) {
  dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
  stopifnot(file.exists(from)); invisible(file.copy(from, to, overwrite = TRUE, copy.date = TRUE))
}
cp_dir <- function(from, to) {
  f <- list.files(from, recursive = TRUE, full.names = FALSE, all.files = FALSE)
  f <- f[!grepl("(^|/)(Icon\r?|\\.DS_Store)$", f)]
  for (x in f) cp(file.path(from, x), file.path(to, x))
  length(f)
}

# ---- analyzer raw data -------------------------------------------------------
lgr <- list.files(LGR_SRC, pattern = "^micro_\\d{4}-\\d{2}-\\d{2}_f\\d+\\.txt(\\.zip)?$", recursive = TRUE)
d <- as.Date(sub("^.*micro_(\\d{4}-\\d{2}-\\d{2}).*$", "\\1", lgr))
lgr <- lgr[d >= STUDY[1] & d <= STUDY[2]]
for (x in lgr) cp(file.path(LGR_SRC, x), file.path(PKG, "analyzer_raw", "lgr", x))
li <- unlist(lapply(LI_SRC, list.files, pattern = "\\.data$", full.names = TRUE))
li <- li[!duplicated(basename(li))]
for (x in li) cp(x, file.path(PKG, "analyzer_raw", "li7810", basename(x)))
message("Analyzer files: LGR ", length(lgr), ", LI-7810 ", length(li))

# ---- field logs and chamber geometry -----------------------------------------
fl <- c(file.path("processing_csvs", "data", c("Field_Data_Monthly_Summer2023.csv", "Field_Data_Monthly_updated2.csv", "summer_2024.csv")),
        file.path("May23_Sept24", "Timing Updates.xlsx"),
        file.path("Sept2024_onwards", setdiff(list.files(file.path(FL_SRC, "Sept2024_onwards"), "\\.xlsx$"), "Template.xlsx")))
for (x in fl) cp(file.path(FL_SRC, x), file.path(PKG, "field_logs", basename(x)))
cp(file.path(FL_SRC, "processing_csvs", "tree_volumes.csv"), file.path(PKG, "chamber_volumes.csv"))
message("Field logs: ", length(fl))

# ---- earlier processed dataset and QC decisions -------------------------------
# Source of the curated tree IDs, met at each measurement and the fluxes of the 99
# measurements that have no archived 1-Hz record.
cp(file.path(IN_SRC, "HF_2023-2025_tree_flux_corrected.csv"), file.path(PKG, "previous_processing", "HF_2023-2025_tree_flux_v1.csv"))
cp(file.path(IN_SRC, "manual_windows.csv"), file.path(PKG, "qc_decisions", "manual_windows.csv"))

# ---- tomography (companion study, Thompson et al. 2026) -----------------------
cp(file.path(IN_SRC, "tomography_results_compiled.csv"), file.path(PKG, "tomography", "tomography_results_compiled.csv"))
cp(file.path(TOMO_SRC, "ERT_application_results.csv"), file.path(PKG, "tomography", "ERT_application_results.csv"))
cp(file.path(TOMO_SRC, "Tree_ID_info.csv"), file.path(PKG, "tomography", "Tree_ID_info.csv"))
n_img <- cp_dir(file.path(IN_SRC, "tomography"), file.path(PKG, "tomography", "images"))
message("Tomography images: ", n_img)

# ---- stand survey -------------------------------------------------------------
cp(file.path(IN_SRC, "spatial", "BGS_VRP_2025.csv"), file.path(PKG, "stand", "BGS_VRP_2025.csv"))
cp(file.path(IN_SRC, "spatial", "Black_Gum_Swamp.kmz"), file.path(PKG, "stand", "Black_Gum_Swamp.kmz"))

# ---- trees.csv: one row per study tree ----------------------------------------
gps <- bind_rows(
  st_read(file.path(IN_SRC, "spatial", "BGS_editable.geojson"), quiet = TRUE) %>% mutate(plot = "BGS", site = "Wetland"),
  st_read(file.path(IN_SRC, "spatial", "EMS_trees_editable.geojson"), quiet = TRUE) %>% mutate(plot = "EMS", site = "Upland")) %>%
  st_zm() %>% st_transform(4326)
xy <- st_coordinates(gps)
gps <- st_drop_geometry(gps) %>% transmute(Tree = as.integer(Name), plot, site, gps_species = tolower(species),
                                           lon = round(xy[, 1], 7), lat = round(xy[, 2], 7), gps_date = as.Date(timestamp))
v1 <- read.csv(file.path(IN_SRC, "HF_2023-2025_tree_flux_corrected.csv"), stringsAsFactors = FALSE) %>%
  filter(!is.na(PLOT)) %>% group_by(Tree) %>%
  summarise(species = first(na.omit(SPECIES)), DBH_cm = first(na.omit(DBH)), .groups = "drop")
mt <- read.csv(file.path(IN_SRC, "hummock_hollow.csv"), check.names = FALSE, fileEncoding = "UTF-8-BOM")
names(mt)[1] <- "Tree"
latin <- c(rm = "Acer rubrum", bg = "Nyssa sylvatica", hem = "Tsuga canadensis", ro = "Quercus rubra")
trees <- gps %>% left_join(v1, by = "Tree") %>%
  left_join(mt %>% transmute(Tree = as.integer(Tree), microtopography = Classification), by = "Tree") %>%
  mutate(species_name = unname(latin[species]))
mism <- trees %>% filter(species != gps_species)
if (nrow(mism)) { print(mism); stop("species in GPS survey and tree list disagree") }
trees <- trees %>% select(Tree, site, plot, species, species_name, DBH_cm, microtopography, lat, lon, gps_date)
stopifnot(nrow(trees) == 60, !anyNA(trees$species), !anyNA(trees$DBH_cm))
write.csv(trees, file.path(PKG, "trees.csv"), row.names = FALSE)
message("trees.csv: ", nrow(trees), " trees")

# ============================================================
# 00_download_edi.R
# Download study data from EDI and set up data/input/
#
# Downloads the data package for this study from the
# Environmental Data Initiative (EDI) repository and unpacks
# it into the directory structure expected by analysis scripts.
#
# EDI package contents (flat):
#   HF_2023-2025_tree_flux_corrected.csv
#   tomography_results_compiled.csv
#   hummock_hollow.csv
#   tomography_ERTs.zip
#   tomography_PITs.zip
#
# Unpacked structure (what scripts expect):
#   data/input/
#   ├── HF_2023-2025_tree_flux_corrected.csv
#   ├── tomography_results_compiled.csv
#   ├── hummock_hollow.csv
#   └── tomography/
#       ├── ERTs_Absolute/  (unzipped JPGs)
#       └── CH4_PITs/       (unzipped JPGs)
#
# Requires: EDIutils
#   install.packages("EDIutils")
# ============================================================

library(EDIutils)

# ============================================================
# CONFIGURATION
# ============================================================

# EDI package identifier — update after publication
PACKAGE_ID <- "edi.XXXXX.1"  # TODO: replace with published ID

# Output directory
OUTPUT_DIR <- "data/input"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 1: Get entity list from EDI
# ============================================================

message("=== Fetching package info from EDI ===")
message("Package: ", PACKAGE_ID)

# Parse scope, identifier, revision
parts <- strsplit(PACKAGE_ID, "\\.")[[1]]
scope <- parts[1]
identifier <- parts[2]
revision <- parts[3]

# Get list of data entities in the package
entity_names <- read_data_entity_names(
  packageId = PACKAGE_ID
)

message("Found ", nrow(entity_names), " entities:")
for (i in seq_len(nrow(entity_names))) {
  message("  ", entity_names$entityName[i])
}

# ============================================================
# STEP 2: Download each entity
# ============================================================

message("\n=== Downloading data entities ===")

for (i in seq_len(nrow(entity_names))) {
  entity_name <- entity_names$entityName[i]
  entity_id <- entity_names$entityId[i]

  # Determine output path
  out_path <- file.path(OUTPUT_DIR, entity_name)

  if (file.exists(out_path)) {
    message("  Skipping (exists): ", entity_name)
    next
  }

  message("  Downloading: ", entity_name, " ...")

  # Get the raw data bytes
  raw_data <- read_data_entity(
    packageId = PACKAGE_ID,
    entityId = entity_id
  )

  writeBin(raw_data, out_path)
  message("    Saved: ", out_path, " (", round(file.size(out_path) / 1024), " KB)")
}

# ============================================================
# STEP 3: Unzip tomography archives
# ============================================================

message("\n=== Unpacking tomography images ===")

ert_zip <- file.path(OUTPUT_DIR, "tomography_ERTs.zip")
pit_zip <- file.path(OUTPUT_DIR, "tomography_PITs.zip")

ert_dir <- file.path(OUTPUT_DIR, "tomography", "ERTs_Absolute")
pit_dir <- file.path(OUTPUT_DIR, "tomography", "CH4_PITs")

if (file.exists(ert_zip)) {
  dir.create(ert_dir, recursive = TRUE, showWarnings = FALSE)
  unzip(ert_zip, exdir = ert_dir)
  n_ert <- length(list.files(ert_dir, pattern = "\\.jpg$"))
  message("  Unpacked ERTs: ", n_ert, " images → ", ert_dir)
  # Remove zip after unpacking
  file.remove(ert_zip)
} else {
  message("  WARNING: tomography_ERTs.zip not found")
}

if (file.exists(pit_zip)) {
  dir.create(pit_dir, recursive = TRUE, showWarnings = FALSE)
  unzip(pit_zip, exdir = pit_dir)
  n_pit <- length(list.files(pit_dir, pattern = "\\.jpg$"))
  message("  Unpacked PITs: ", n_pit, " images → ", pit_dir)
  file.remove(pit_zip)
} else {
  message("  WARNING: tomography_PITs.zip not found")
}

# ============================================================
# SUMMARY
# ============================================================

message("\n=== Done ===")
message("Data ready at: ", OUTPUT_DIR)
message("\nContents:")
all_files <- list.files(OUTPUT_DIR, recursive = TRUE)
for (f in all_files) message("  ", f)
message("\nNext: run scripts in 00_download/ (01-03) for environmental data,")
message("then scripts in 01_import/ for preprocessing.")

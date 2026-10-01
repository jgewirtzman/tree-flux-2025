# ============================================================
# 02_upload_edi.R
# Upload data package to EDI repository
#
# Run 01_build_edi_package.R first to generate EML and package.
#
# Workflow:
#   1. Upload to STAGING first to test
#   2. Review at https://portal-s.edirepository.org/
#   3. When satisfied, upload to PRODUCTION
#
# Requires: EDIutils
#   install.packages("EDIutils")
#
# You need an EDI account. Request one at: support@edirepository.org
# ============================================================

library(EDIutils)

# ============================================================
# CONFIGURATION
# ============================================================

PACKAGE_DIR <- "data/edi/package"
PACKAGE_ID  <- "edi.XXXXX.1"  # TODO: must match 01_build_edi_package.R

# Environment: "staging" for testing, "production" for final
ENVIRONMENT <- "staging"

# ============================================================
# STEP 1: Authenticate
# ============================================================

message("=== Authenticating with EDI (", ENVIRONMENT, ") ===")
message("You will be prompted for your EDI username and password.")

login(environment = ENVIRONMENT)

# ============================================================
# STEP 2: Reserve package identifier (production only)
# ============================================================

# Uncomment this block when ready for production:
# message("\n=== Reserving package identifier ===")
# res <- create_reservation(scope = "edi", environment = "production")
# message("Reserved: ", res)
# # Update PACKAGE_ID in both scripts with the returned value

# ============================================================
# STEP 3: Upload package
# ============================================================

eml_file <- file.path(PACKAGE_DIR, paste0(PACKAGE_ID, ".xml"))

if (!file.exists(eml_file)) {
  stop("EML file not found: ", eml_file,
       "\nRun 01_build_edi_package.R first.")
}

message("\n=== Uploading to EDI (", ENVIRONMENT, ") ===")
message("EML: ", eml_file)

# List all data files in the package
data_files <- list.files(PACKAGE_DIR, full.names = TRUE)
data_files <- data_files[!grepl("\\.xml$", data_files)]
message("Data files:")
for (f in data_files) message("  ", basename(f))

# Create the data package
transaction <- create_data_package(
  eml = eml_file,
  env = ENVIRONMENT
)

message("\n=== Upload initiated ===")
message("Transaction ID: ", transaction)
message("\nCheck status:")
message("  check_status_create(\"", transaction, "\", env = \"", ENVIRONMENT, "\")")

if (ENVIRONMENT == "staging") {
  message("\nOnce complete, review at:")
  message("  https://portal-s.edirepository.org/nis/mapbrowse?packageid=", PACKAGE_ID)
  message("\nWhen satisfied, update ENVIRONMENT to 'production' and re-run.")
} else {
  message("\nOnce complete, view at:")
  message("  https://portal.edirepository.org/nis/mapbrowse?packageid=", PACKAGE_ID)
}

# ============================================================
# OPTIONAL: Check upload status
# ============================================================

# Uncomment to check:
# check_status_create(transaction, env = ENVIRONMENT)

# Logout when done
# logout()

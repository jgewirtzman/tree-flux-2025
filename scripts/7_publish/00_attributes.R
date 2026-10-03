# ============================================================
# 00_attributes.R
# EMLassemblyline templates for the data package, written from the tracked metadata in
# scripts/7_publish/metadata/ so that the metadata always match the files:
#   - attributes_<table>.txt for every data table, from that table's dictionary
#     (stem_ch4_flux: data/final/stem_ch4_flux_dictionary.csv; all others:
#     metadata/dictionaries/<table>.csv); each dictionary is checked against the table's columns
#   - abstract.md, methods.md, keywords.txt, personnel.txt, intellectual_rights.txt,
#     custom_units.txt, geographic_coverage.txt copied as they are
# Output: data/edi/ (templates; regenerated, not tracked)
# ============================================================
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))   # notes contain non-ASCII characters
EDI <- file.path("data", "edi"); META <- file.path("scripts", "7_publish", "metadata")
unlink(list.files(EDI, pattern = "^(attributes|catvars)_.*\\.txt$", full.names = TRUE))
dir.create(EDI, showWarnings = FALSE)
file.copy(file.path(META, c("abstract.md", "methods.md", "keywords.txt", "personnel.txt",
                            "intellectual_rights.txt", "custom_units.txt", "geographic_coverage.txt")), EDI, overwrite = TRUE)

# dictionary unit -> EML unit (standard EML units, or custom_units.txt)
units <- c("nmol m-2 s-1" = "nanomolePerMeterSquaredPerSecond", "umol m-2 s-1" = "micromolePerMeterSquaredPerSecond",
           ppb = "nanomolePerMole", ppm = "micromolePerMole", "nmol mol-1" = "nanomolePerMole",
           "umol mol-1" = "micromolePerMole", s = "second", h = "hour", d = "nominalDay", cm = "centimeter",
           mm = "millimeter", cm2 = "squareCentimeter", m2 = "squareMeter", cm3 = "cubicCentimeter", L = "liter",
           mol = "mole", degC = "celsius", kPa = "kilopascal", mbar = "millibar", "%" = "percent", deg = "degree",
           "W m-2" = "wattPerMeterSquared", "m s-1" = "meterPerSecond", "m3 m-3" = "cubicMeterPerCubicMeter",
           "ohm m" = "ohmMeter")
# date and time columns: format strings by column name (BGS_VRP_2025 date: YYYYMMDD)
DT <- c(date = "YYYY-MM-DD", DATE = "YYYY-MM-DD", gps_date = "YYYY-MM-DD", sample_time_local = "YYYY-MM-DD hh:mm:ss",
        sample_hour_est = "YYYY-MM-DD hh:mm:ss", window_start = "YYYY-MM-DD hh:mm:ss", window_end = "YYYY-MM-DD hh:mm:ss",
        start = "YYYY-MM-DD hh:mm:ss", end = "YYYY-MM-DD hh:mm:ss", datetime_posx = "YYYY-MM-DDThh:mm:ssZ",
        datetime = "YYYY-MM-DDThh:mm:ssZ")
# identifiers stored as numbers
IDS <- c("Tree", "tree", "Tag", "Old tag", "tree_id", "...1", "OLD_TAG")

write_attr <- function(csv, dict, name, date_fmt = DT) {
  d <- read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  if (!setequal(names(d), dict$column))
    stop(name, ": dictionary and table differ: ", paste(c(setdiff(names(d), dict$column), setdiff(dict$column, names(d))), collapse = ", "))
  dict <- dict[match(names(d), dict$column), ]
  cls <- vapply(names(d), function(v) {
    if (v %in% names(date_fmt)) "Date"
    else if (is.numeric(d[[v]]) && !v %in% IDS) "numeric"
    else "character"
  }, "")
  un <- ifelse(is.na(dict$unit), "", dict$unit)
  u <- ifelse(cls == "numeric", ifelse(nzchar(un), unname(units[un]), "dimensionless"), "")
  if (anyNA(u)) stop(name, ": unit not mapped: ", paste(unique(un[is.na(u)]), collapse = ", "))
  def <- dict$description
  if ("source" %in% names(dict)) def <- ifelse(nzchar(dict$source) & !is.na(dict$source), paste0(def, " [source: ", dict$source, "]"), def)
  out <- data.frame(attributeName = names(d), attributeDefinition = def, class = cls, unit = u,
                    dateTimeFormatString = ifelse(cls == "Date", unname(date_fmt[names(d)]), ""),
                    missingValueCode = "NA", missingValueCodeExplanation = "not available")
  write.table(out, file.path(EDI, paste0("attributes_", name, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  message("attributes_", name, ".txt: ", nrow(out), " columns")
}
dict <- function(name) read.csv(file.path(META, "dictionaries", paste0(name, ".csv")), stringsAsFactors = FALSE, check.names = FALSE)

TABLES <- c(
  stem_ch4_flux                = "data/final/stem_ch4_flux.csv",
  stem_ch4_flux_dictionary     = "data/final/stem_ch4_flux_dictionary.csv",
  flux_processing_log          = "data/final/flux_processing_log.csv",
  environment_hourly           = "data/final/environment_hourly.csv",
  environment_hourly_dictionary = file.path(META, "dictionaries", "environment_hourly.csv"),
  tomography_classes           = "data/final/tomography_classes.csv",
  trees                        = "data/package/trees.csv",
  chamber_volumes              = "data/package/chamber_volumes.csv",
  tomography_results_compiled  = "data/package/tomography/tomography_results_compiled.csv",
  ERT_application_results      = "data/package/tomography/ERT_application_results.csv",
  Tree_ID_info                 = "data/package/tomography/Tree_ID_info.csv",
  manual_windows               = "data/package/qc_decisions/manual_windows.csv",
  "HF_2023-2025_tree_flux_v1"  = "data/package/previous_processing/HF_2023-2025_tree_flux_v1.csv",
  BGS_VRP_2025                 = "data/package/stand/BGS_VRP_2025.csv")

meta <- dict("dictionary_meta")
for (nm in names(TABLES)) {
  dd <- switch(nm,
    stem_ch4_flux = read.csv("data/final/stem_ch4_flux_dictionary.csv", stringsAsFactors = FALSE),
    stem_ch4_flux_dictionary = meta[meta$column != "source", ],
    environment_hourly_dictionary = meta,
    dict(nm))
  write_attr(TABLES[[nm]], dd, nm,
             date_fmt = if (nm == "BGS_VRP_2025") c(date = "YYYYMMDD") else if (grepl("dictionary", nm)) character(0) else DT)
}

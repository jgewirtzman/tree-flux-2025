# ============================================================
# 00_attributes.R
# EMLassemblyline attribute templates for the tables of the data package, written
# from the data dictionaries so the metadata always match the files.
#   data/final/stem_ch4_flux.csv          <- data/final/stem_ch4_flux_dictionary.csv
#   data/final/flux_processing_log.csv
#   data/package/trees.csv
#   data/package/chamber_volumes.csv
# (tomography_results_compiled.csv keeps its hand-written template in data/edi/.)
# Categorical-variable templates: EMLassemblyline::template_categorical_variables()
# after this script, then fill in the definitions.
# Output: data/edi/attributes_<table>.txt
# ============================================================
EDI <- file.path("data", "edi"); dir.create(EDI, showWarnings = FALSE)
units <- c("nmol m-2 s-1" = "nanomolePerMeterSquaredPerSecond", "umol m-2 s-1" = "micromolePerMeterSquaredPerSecond",
           ppb = "partsPerBillion", ppm = "partsPerMillion", s = "second", cm = "centimeter", cm2 = "squareCentimeter",
           L = "liter", degC = "celsius", kPa = "kilopascal", d = "day", "%" = "percent", deg = "degree")
DT <- c(date = "YYYY-MM-DD", gps_date = "YYYY-MM-DD", sample_time_local = "YYYY-MM-DD hh:mm:ss",
        sample_hour_est = "YYYY-MM-DD hh:mm:ss", window_start = "YYYY-MM-DD hh:mm:ss", window_end = "YYYY-MM-DD hh:mm:ss",
        datetime_posx = "YYYY-MM-DDThh:mm:ssZ")

write_attr <- function(csv, dict, name) {
  d <- read.csv(csv, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(setequal(names(d), dict$column))
  dict <- dict[match(names(d), dict$column), ]
  cls <- vapply(names(d), function(v) {
    if (v %in% names(DT)) "Date"
    else if (is.numeric(d[[v]]) && !v %in% c("Tree", "year", "sampling_round")) "numeric"
    else if (is.logical(d[[v]]) || length(unique(na.omit(d[[v]]))) <= 12) "categorical"
    else "character"
  }, "")
  u <- ifelse(cls == "numeric", ifelse(nzchar(dict$unit), unname(units[dict$unit]), "dimensionless"), "")
  stopifnot(!anyNA(u))
  out <- data.frame(attributeName = names(d), attributeDefinition = dict$description, class = cls, unit = u,
                    dateTimeFormatString = ifelse(cls == "Date", unname(DT[names(d)]), ""),
                    missingValueCode = "NA", missingValueCodeExplanation = "not available")
  write.table(out, file.path(EDI, paste0("attributes_", name, ".txt")), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  message("attributes_", name, ".txt: ", nrow(out), " columns")
}

write_attr("data/final/stem_ch4_flux.csv", read.csv("data/final/stem_ch4_flux_dictionary.csv", stringsAsFactors = FALSE), "stem_ch4_flux")
write_attr("data/final/flux_processing_log.csv", data.frame(
  column = c("script", "step", "rule", "n_records", "detail"), unit = "",
  description = c("processing script (scripts/2_flux/)", "processing step", "cleaning rule applied",
                  "number of records the rule touched", "details")), "flux_processing_log")
write_attr("data/package/trees.csv", data.frame(
  column = c("Tree", "site", "plot", "species", "species_name", "DBH_cm", "microtopography", "lat", "lon", "gps_date"),
  unit = c("", "", "", "", "", "cm", "", "deg", "deg", ""),
  description = c("ForestGEO tag of the tree", "Wetland (Black Gum Swamp) or Upland (EMS)", "plot code: BGS or EMS",
                  "species code: bg, hem, rm, ro", "species name", "stem diameter at breast height",
                  "wetland trees: hummock, hollow or hollow/hummock", "latitude (WGS84; field GPS)",
                  "longitude (WGS84; field GPS)", "date of the GPS survey")), "trees")
cv <- read.csv("data/package/chamber_volumes.csv", check.names = FALSE)
write_attr("data/package/chamber_volumes.csv", data.frame(column = names(cv), unit = "",
  description = paste("chamber_volumes.csv column", names(cv), "(definition to complete)")), "chamber_volumes")

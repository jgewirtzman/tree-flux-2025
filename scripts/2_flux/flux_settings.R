# ============================================================
# flux_settings.R -- constants, paths and the processing log shared by the
# 2_flux scripts (sourced, not run on its own).
# ============================================================

# ---- inputs (data package) ----------------------------------------------------
PKG             <- "data/package"
PATH_MEASUREMENTS <- file.path(PKG, "measurements.csv")   # one row per measurement: tree, analyzer, closure times, met
PATH_VOLUMES    <- file.path(PKG, "chamber_volumes.csv")
PATH_LGR        <- file.path(PKG, "analyzer_raw", "lgr")
PATH_LI7810     <- file.path(PKG, "analyzer_raw", "li7810")
PATH_TREES      <- file.path(PKG, "trees.csv")

# ---- intermediate and final outputs -------------------------------------------
INTERIM         <- file.path("data", "interim")
FINAL           <- file.path("data", "final")
PATH_CLOSURES   <- file.path(INTERIM, "flux_closures.rds")      # 01 -> 02
PATH_FITS       <- file.path(INTERIM, "flux_fits.csv")          # 02 -> 03, 05, 06
PATH_TRACES     <- file.path(INTERIM, "flux_traces.rds")        # 02 -> 03, 04, 05
PATH_TRACE_QC   <- file.path(INTERIM, "flux_trace_qc.csv")      # 03 -> 06
PATH_FLUX       <- file.path(FINAL, "stem_ch4_flux.csv")        # 06 -> analyses
PATH_SETTINGS   <- file.path(FINAL, "flux_processing_settings.json")
PATH_LOG        <- file.path(INTERIM, "flux_processing_log.csv")   # local record of the processing steps
TAB_DIR         <- file.path("outputs", "tables", "flux_processing")
FIG_DIR         <- file.path("outputs", "figures", "flux_processing")
for (d in c(INTERIM, FINAL, TAB_DIR, FIG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- constants ------------------------------------------------------------------
TZ            <- "America/New_York"   # wall clock of both analyzers and the field logs
SURFAREA_M2   <- pi * 0.0508^2        # m2, collar inner radius 5.08 cm (all collars)
AREA_CM2      <- SURFAREA_M2 * 1e4    # goFlux wants cm2
# Analyzer internal volume added to the measured collar volume (chamber_volumes.csv already holds
# collar + cap + 28.9 cm3 of tubing): 0.028 L for both analyzers, the LI-COR LI-7810 total sample
# volume; goFlux lists 25 cm3 for the LGR microportable GLA131 cavity, within 10 %.
EXTRA_VOL_LGR_L  <- 0.028
EXTRA_VOL_7810_L <- 0.028
SHOULDER_S    <- 120                  # s of trace kept outside the window (flag = 0), for plots/QC context
DEADBAND_S    <- 20                   # s dropped after chamber closure (measurements.csv deadband_s; 04_deadband_sensitivity.R)
PREC_LGR      <- c(CO2 = 0.35, CH4 = 0.9, H2O = 100)   # datasheet 1-s precision, GLA131 (ppm, ppb, ppm)
PREC_7810     <- c(CO2 = 3.5,  CH4 = 0.6, H2O = 45)    # LI-7810 (goFlux default)

`%||%` <- function(a, b) if (is.null(a)) b else a
nz <- function(x) !is.na(x) & nchar(trimws(as.character(x))) > 0

# ---- processing log -------------------------------------------------------------
# Every cleaning rule records how many records it touched. Each script writes its own
# part to data/interim/flux_log_<script>.csv; 05_flux_dataset.R combines them into
# data/interim/flux_processing_log.csv (a local record; not published).
FLUX_LOG <- data.frame(script = character(), step = character(), rule = character(),
                       n_records = integer(), detail = character())
log_step <- function(step, rule, n, detail = "") {
  FLUX_LOG <<- rbind(FLUX_LOG, data.frame(script = SCRIPT, step = step, rule = rule,
                                          n_records = as.integer(n), detail = detail))
  message(sprintf("  [log] %s -- %s: %d%s", step, rule, as.integer(n), if (nzchar(detail)) paste0(" (", detail, ")") else ""))
}
write_log <- function() write.csv(FLUX_LOG, file.path(INTERIM, paste0("flux_log_", SCRIPT, ".csv")), row.names = FALSE)

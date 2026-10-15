# ============================================================
# 06_flux_dataset.R  (flux step 5 of 5)
#
# Compiles the final stem flux dataset used by every analysis and published with
# the data package: one row per measurement, final flux, uncertainty, detection
# limit, QC flags, tree, analyzer, window, chamber and sampling time.
#
#   detection      MDF = 1.96 * sigma / t * flux.term (1.96: benchmark multiplier, not a calibrated 95 % test);
#                  sigma = analyzer precision on that day. Retain-and-flag: below-MDF fluxes are kept and flagged.
#   sampling time  sample_time_local from measurements.csv; sample_hour_est is the hour on the EST clock
#                  used by every driver
#   QC flags       trace QC (03_trace_qc.R) and goFlux QC screens; tree metadata checked against trees.csv
#
# Inputs : data/interim/flux_fits.csv, data/interim/flux_trace_qc.csv, data/package/trees.csv
# Outputs: data/final/stem_ch4_flux.csv             the dataset
#          data/final/stem_ch4_flux_dictionary.csv  column definitions and units
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(lubridate) })
SCRIPT <- "06_flux_dataset"
source("scripts/2_flux/flux_settings.R")

g <- read.csv(PATH_FITS, stringsAsFactors = FALSE)
message("Measurements: ", nrow(g), " (fitted: ", sum(g$fitted), ")")
stopifnot(all(g$fitted))
df <- g %>% mutate(species_label = c(bg = "N. sylvatica", hem = "T. canadensis", rm = "A. rubrum", ro = "Q. rubra")[SPECIES])

# ---- detection limit ---------------------------------------------------------------
df <- df %>%
  mutate(
    closure_s = t_sec,
    CH4_sigma_ppb = CH4_sigma_mad,
    CO2_sigma_ppm = CO2_sigma_mad,
    CH4_MDF = abs(qnorm(0.975) * CH4_sigma_ppb / closure_s * flux_term),
    CO2_MDF = abs(qnorm(0.975) * CO2_sigma_ppm / closure_s * flux_term),
    CH4_below_MDF = abs(CH4_flux_nmolpm2ps) < CH4_MDF,
    CO2_below_MDF = abs(CO2_flux_umolpm2ps) < CO2_MDF,
    CH4_det_class = ifelse(CH4_flux_nmolpm2ps > CH4_MDF, "emission",
                           ifelse(CH4_flux_nmolpm2ps < -CH4_MDF, "uptake", "below detection")))
chk <- df %>% summarise(d = max(abs(CH4_MDF - g$CH4_MDF[match(closure_id, g$closure_id)]) / CH4_MDF, na.rm = TRUE))
stopifnot(chk$d < 1e-8)   # same MDF as 02_fit_fluxes.R
log_step("detection", "CH4 fluxes below the minimum detectable flux (retained, flagged)", sum(df$CH4_below_MDF, na.rm = TRUE))

# ---- sampling time --------------------------------------------------------------------
fix_t <- as.POSIXct(df$sample_time_local, format = "%Y-%m-%d %H:%M:%S", tz = "America/New_York")
df$sample_hour_est <- format(floor_date(with_tz(fix_t, "EST"), "hour"), "%Y-%m-%d %H:00:00")
# sampling campaigns: in time order, a gap of more than 7 days since the
# previous measurement starts a new campaign
df <- df %>% arrange(sample_time_local, Tree) %>%
  mutate(days_since_last = as.numeric(as.Date(date) - lag(as.Date(date))), days_since_last = coalesce(days_since_last, 0),
         new_round = days_since_last > 7, sampling_round = cumsum(new_round) + 1L)

# ---- QC flags and tree metadata --------------------------------------------------------
tq <- read.csv(PATH_TRACE_QC, stringsAsFactors = FALSE) %>%
  transmute(closure_id, trace_window_problem = window_problem, trace_data_problem = data_problem,
            trace_qc_reasons = reasons)
df <- df %>% left_join(tq, by = "closure_id")
trees <- read.csv(PATH_TREES, stringsAsFactors = FALSE)
k <- match(df$Tree, trees$Tree)
# PLOT keeps the field coding (upland trees in plot E5 are coded "E5"); site and species must agree
bad <- is.na(k) | df$SPECIES != trees$species[k] | df$location != trees$site[k]
if (any(bad)) { print(df[bad, c("closure_id", "Tree", "SPECIES", "PLOT", "location")]); stop("tree metadata disagree with trees.csv") }
log_step("QC flags", "measurements whose tree, species and site agree with trees.csv", nrow(df))

# ---- final table -------------------------------------------------------------------------
fin <- df %>% arrange(sample_time_local, Tree) %>% transmute(
  measurement_id = closure_id, Tree, location, PLOT, SPECIES, species_label, DBH,
  date, sample_time_local, sample_hour_est, year, season,
  sampling_round, new_round, days_since_last,
  instrument, inst_label, closure_start, closure_end, deadband_s, window_start, window_end, closure_s,
  end_trim_s, Area_cm2, Vtot_L, Tcham_C, Pcham_kPa,
  CH4_flux_nmolpm2ps, CH4_SE, CH4_r2, CH4_model, CH4_g_factor = round(CH4_g_factor, 4),
  CO2_flux_umolpm2ps, CO2_SE, CO2_r2, CO2_model,
  CH4_sigma_ppb, CO2_sigma_ppm, CH4_MDF, CH4_below_MDF, CH4_det_class, CO2_MDF, CO2_below_MDF,
  qc_c0, qc_co2_tracer, qc_convex, qc_noisy, qc_note,
  trace_window_problem, trace_data_problem, trace_qc_reasons)
fin <- fin %>% mutate(across(c(qc_note, trace_qc_reasons), ~ na_if(.x, "")))   # NA = no screen fired
write.csv(fin, PATH_FLUX, row.names = FALSE, na = "NA")
message("Wrote ", PATH_FLUX, ": ", nrow(fin), " measurements, ", ncol(fin), " columns")

dict_txt <- '
column,unit,description
measurement_id,,unique measurement ID: analyzer_date_tree_start (LGR/UGGA or LI-7810)
Tree,,tag number of the tree (see trees.csv)
location,,site: Wetland (Black Gum Swamp) or Upland (EMS)
PLOT,,plot code as recorded: BGS (wetland); EMS or E5 (upland trees in plot E5)
SPECIES,,species code: bg Nyssa sylvatica; hem Tsuga canadensis; rm Acer rubrum; ro Quercus rubra
species_label,,species name
DBH,cm,stem diameter at breast height
date,,measurement date (local)
sample_time_local,,sampling time (start of the measurement; America/New_York wall clock)
sample_hour_est,,sampling hour on the EST clock (UTC-5; no daylight saving) used to join the drivers
year,,calendar year
season,,winter (DJF) / spring (MAM) / summer (JJA) / fall (SON)
sampling_round,,sampling campaign (campaigns separated by more than 7 days)
new_round,,TRUE for the first measurement of a campaign
days_since_last,d,days since the previous measurement (0 within a field day)
instrument,,analyzer: LGR1 / LGR3 (LGR microportable GLA131 / UGGA) or LI-7810
inst_label,,analyzer type: LGR/UGGA or LI-7810
closure_start,,chamber closed (analyzer clock; local)
closure_end,,chamber opened (analyzer clock; local)
deadband_s,s,seconds discarded after chamber closure (0 where the start was set on inspection of the record)
window_start,,start of the fitting window (local)
window_end,,end of the fitting window (local)
closure_s,s,length of the fitting window
end_trim_s,s,seconds removed before closure_end at a detected chamber opening
Area_cm2,cm2,collar area (inner radius 5.08 cm)
Vtot_L,L,system volume: collar + cap + tubing + 0.028 L analyzer volume
Tcham_C,degC,air temperature at the time of measurement (Fisher Meteorological Station)
Pcham_kPa,kPa,air pressure at the time of measurement (Fisher Meteorological Station)
CH4_flux_nmolpm2ps,nmol m-2 s-1,stem CH4 flux (positive = emission)
CH4_SE,nmol m-2 s-1,standard error of the CH4 flux (of the selected model)
CH4_r2,,R2 of the linear CH4 fit
CH4_model,,model selected by goFlux best.flux: LM (linear) or HM (Hutchinson-Mosier)
CH4_g_factor,,ratio HM / LM flux
CO2_flux_umolpm2ps,umol m-2 s-1,stem CO2 flux
CO2_SE,umol m-2 s-1,standard error of the CO2 flux
CO2_r2,,R2 of the linear CO2 fit
CO2_model,,model selected for CO2
CH4_sigma_ppb,ppb,analyzer CH4 precision used for the MDF (that analyzer on that day; goFlux empirical.prec)
CO2_sigma_ppm,ppm,analyzer CO2 precision used for the MDF
CH4_MDF,nmol m-2 s-1,minimum detectable CH4 flux: 1.96 x analyzer precision / closure time x flux term (1.96 is a benchmark multiplier on concentration noise; not a calibrated 95 % detection test)
CH4_below_MDF,,TRUE if |CH4 flux| < CH4_MDF (retained)
CH4_det_class,,emission / uptake / below detection
CO2_MDF,umol m-2 s-1,minimum detectable CO2 flux
CO2_below_MDF,,TRUE if |CO2 flux| < CO2_MDF
qc_c0,,goFlux QC: starting concentration well above ambient
qc_co2_tracer,,goFlux QC: CO2 not clearly increasing during the closure (possible poor seal)
qc_convex,,goFlux QC: concentration curve convex (non-diffusive shape)
qc_noisy,,goFlux QC: record noisier than the analyzer precision
qc_note,,goFlux QC screens that fired
trace_window_problem,,trace QC: window problem (CO2 falling or dropping, flat start, overlap)
trace_data_problem,,trace QC: data problem (CH4 step, gap, start transient, noise)
trace_qc_reasons,,trace QC: all checks that fired
'
# column,unit,description (descriptions may contain commas: split at the first two only)
dl <- trimws(strsplit(dict_txt, "\n")[[1]]); dl <- dl[nzchar(dl)][-1]
dict <- data.frame(column = sub(",.*$", "", dl), unit = sub("^[^,]*,([^,]*),.*$", "\\1", dl),
                   description = sub("^[^,]*,[^,]*,", "", dl), stringsAsFactors = FALSE)
stopifnot(setequal(dict$column, names(fin)))
dict$unit[!nzchar(dict$unit)] <- NA   # NA = no unit
write.csv(dict[match(names(fin), dict$column), ], file.path(FINAL, "stem_ch4_flux_dictionary.csv"), row.names = FALSE, na = "NA")

write_log()
logs <- c("01_closure_table", "02_fit_fluxes", "03_trace_qc", "06_flux_dataset")   # local record (data/interim)
all_log <- bind_rows(lapply(logs, function(s) read.csv(file.path(INTERIM, paste0("flux_log_", s, ".csv")), stringsAsFactors = FALSE)))
write.csv(all_log, PATH_LOG, row.names = FALSE)
message("Wrote ", PATH_LOG, " (", nrow(all_log), " rules)")
message("Detection: ", round(100 * mean(fin$CH4_below_MDF[fin$location == "Upland"])), "% upland and ",
        round(100 * mean(fin$CH4_below_MDF[fin$location == "Wetland"])), "% wetland fluxes below the MDF")

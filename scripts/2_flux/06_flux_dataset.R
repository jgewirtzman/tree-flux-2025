# ============================================================
# 06_flux_dataset.R  (flux step 6 of 6)
#
# Compiles the final stem flux dataset used by every analysis and published with
# the data package: one row per measurement, final flux, uncertainty, detection
# limit, QC flags, tree, analyzer, window, chamber and sampling time.
#
#   9  selection      measurements of the study design (the earlier dataset's 1,640);
#                     closures excluded after inspection removed
#   10 detection      MDF = 1.96 * sigma / t * flux.term (1.96: benchmark multiplier, not a calibrated 95 % test); sigma = analyzer
#                     precision on that day; for measurements without a raw record the
#                     analyzer's median precision and the logged window length.
#                     Retain-and-flag: below-MDF fluxes are kept and flagged.
#   11 sampling time  field-log real time; closures logged without a time take the analyzer
#                     clock minus that analyzer's offset; two AM/PM errors corrected;
#                     sample_hour_est is the hour on the EST clock used by every driver
#   12 QC flags       trace QC (03_trace_qc.R) and fluxqc screens joined; tree metadata
#                     checked against data/package/trees.csv
#
# Inputs : data/interim/flux_fits.csv, data/interim/flux_trace_qc.csv, data/package/trees.csv
# Outputs: data/final/stem_ch4_flux.csv             the dataset
#          data/final/stem_ch4_flux_dictionary.csv  column definitions and units
#          data/final/flux_processing_log.csv       every cleaning rule and records affected
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(lubridate) })
SCRIPT <- "06_flux_dataset"
source("scripts/2_flux/flux_settings.R")

T_FALLBACK <- c("LGR/UGGA" = 180, "LI-7810" = 60)  # s, when no window was recorded at all

g <- read.csv(PATH_FITS, stringsAsFactors = FALSE)
message("Closures: ", nrow(g), " (fitted: ", sum(g$fitted), ")")

# ---- 9 selection -----------------------------------------------------------------
log_step("9 selection", "closures with a raw trace but outside the study design (Jan-Mar 2025 LGR days; not analysed)",
         sum(!g$in_legacy_dataset))
excl <- g$in_legacy_dataset & g$manual_exclude %in% TRUE
log_step("9 selection", "measurements excluded after inspection of the concentration record", sum(excl),
         paste(g$closure_id[excl], collapse = ", "))
df <- g %>% filter(in_legacy_dataset, !manual_exclude %in% TRUE)

# ---- 10 detection limit ---------------------------------------------------------------
pooled <- g %>% filter(fitted) %>% group_by(inst_label) %>%
  summarise(sigma_pooled_CH4 = median(CH4_sigma_mad, na.rm = TRUE),
            sigma_pooled_CO2 = median(CO2_sigma_mad, na.rm = TRUE), .groups = "drop")
df <- df %>% left_join(pooled, by = "inst_label") %>%
  mutate(
    t_sec_fl  = as.numeric(difftime(as.POSIXct(window_end), as.POSIXct(window_start), units = "secs")),
    closure_s = case_when(fitted ~ t_sec, !is.na(t_sec_fl) ~ t_sec_fl, TRUE ~ unname(T_FALLBACK[inst_label])),
    closure_s_src = case_when(fitted ~ "trace", !is.na(t_sec_fl) ~ "field-log window",
                              TRUE ~ paste0("assumed ", unname(T_FALLBACK[inst_label]), " s")),
    instrument = coalesce(instrument, ifelse(inst_label == "LI-7810", "LI-7810", "LGR1")),
    flux_term  = ifelse(fitted, flux_term, fluxqc::flux_term(Vtot_L, Pcham_kPa, Area_cm2, Tcham_C)),
    CH4_sigma_ppb = ifelse(fitted, CH4_sigma_mad, sigma_pooled_CH4),
    CO2_sigma_ppm = ifelse(fitted, CO2_sigma_mad, sigma_pooled_CO2),
    sigma_src  = ifelse(fitted, "analyzer x day (whole-day record)", "analyzer median (no raw record)"),
    CH4_MDF = abs(qnorm(0.975) * CH4_sigma_ppb / closure_s * flux_term),
    CO2_MDF = abs(qnorm(0.975) * CO2_sigma_ppm / closure_s * flux_term),
    CH4_below_MDF = abs(CH4_flux_nmolpm2ps) < CH4_MDF,
    CO2_below_MDF = abs(CO2_flux_umolpm2ps) < CO2_MDF,
    CH4_det_class = ifelse(CH4_flux_nmolpm2ps > CH4_MDF, "emission",
                           ifelse(CH4_flux_nmolpm2ps < -CH4_MDF, "uptake", "below detection")))
chk <- df %>% filter(fitted) %>% summarise(d = max(abs(CH4_MDF - g$CH4_MDF[match(closure_id, g$closure_id)]) / CH4_MDF, na.rm = TRUE))
stopifnot(chk$d < 1e-8)   # same MDF as 02_fit_fluxes.R for fitted closures
log_step("10 detection", "CH4 fluxes below the minimum detectable flux (retained, flagged)", sum(df$CH4_below_MDF, na.rm = TRUE))
log_step("10 detection", "MDF from the analyzer median precision (no raw record)", sum(!df$fitted))

# ---- 11 sampling time ------------------------------------------------------------------
# datetime_posx (= field-log real start) is local wall time stored with a "Z" label.
wall <- as.POSIXct(sub("Z$", "", sub("T", " ", df$datetime_posx)), tz = "UTC")
clo  <- as.POSIXct(df$closure_start, tz = "UTC")
hh   <- as.numeric(format(wall, "%H")) + as.numeric(format(wall, "%M")) / 60
no_time <- !is.na(wall) & format(wall, "%H:%M:%S") == "00:00:00"
ampm    <- !is.na(wall) & !no_time & hh < 6
offs <- as.numeric(difftime(wall, clo, units = "mins"))
ok_off <- !no_time & !ampm & !is.na(offs) & abs(offs) < 180
fix_t <- wall
fix_t[ampm] <- wall[ampm] + 12 * 3600
for (i in which(no_time & !is.na(clo))) {
  same <- which(ok_off & df$inst_label == df$inst_label[i] &
                  abs(as.numeric(difftime(clo, clo[i], units = "days"))) <= 60)
  fix_t[i] <- clo[i] + 60 * (if (length(same)) median(offs[same]) else 0)
}
df$sample_time_src <- ifelse(no_time, "analyzer clock - offset", ifelse(ampm, "field log (+12 h, AM/PM)", "field log real time"))
df$sample_time_local <- format(fix_t, "%Y-%m-%d %H:%M:%S")
df$sample_hour_est <- format(floor_date(with_tz(force_tz(fix_t, "America/New_York"), "EST"), "hour"), "%Y-%m-%d %H:00:00")
log_step("11 sampling time", "no time logged; analyzer clock minus that analyzer's offset to field time (median of closures within 60 d)", sum(no_time))
log_step("11 sampling time", "time logged as AM instead of PM; +12 h", sum(ampm))

# ---- 12 QC flags and tree metadata --------------------------------------------------------
tq <- read.csv(PATH_TRACE_QC, stringsAsFactors = FALSE) %>%
  transmute(closure_id, trace_window_problem = window_problem, trace_data_problem = data_problem,
            trace_qc_reasons = reasons)
df <- df %>% left_join(tq, by = "closure_id")
trees <- read.csv(PATH_TREES, stringsAsFactors = FALSE)
k <- match(df$Tree, trees$Tree)
# PLOT keeps the field coding (upland trees in the ForestGEO sub-plot E5 are coded "E5"); site and species must agree
bad <- is.na(k) | df$SPECIES != trees$species[k] | df$location != trees$site[k]
if (any(bad)) { print(df[bad, c("closure_id", "Tree", "SPECIES", "PLOT", "location")]); stop("tree metadata disagree with trees.csv") }
log_step("12 QC flags", "measurements whose tree, species and site agree with trees.csv", nrow(df))

# ---- final table -------------------------------------------------------------------------
fin <- df %>% arrange(sample_time_local, Tree) %>% transmute(
  measurement_id = closure_id, Tree, location, PLOT, SPECIES, species_label, DBH,
  date, datetime_posx, sample_time_local, sample_hour_est, sample_time_src, year, season,
  sampling_round, new_round, days_since_last,
  instrument, inst_label, window_src, window_start, window_end, closure_s, closure_s_src,
  manual_window, end_trim_s, Area_cm2, Vtot_L, Tcham_C, Pcham_kPa, met_src,
  flux_source, fitted,
  CH4_flux_nmolpm2ps, CH4_SE, CH4_r2, CH4_model, CH4_g_factor = round(CH4_g_factor, 4),
  CO2_flux_umolpm2ps, CO2_SE, CO2_r2, CO2_model,
  CH4_sigma_ppb, CO2_sigma_ppm, sigma_src, CH4_MDF, CH4_below_MDF, CH4_det_class, CO2_MDF, CO2_below_MDF,
  qc_c0, qc_co2_tracer, qc_convex, qc_noisy, qc_note,
  trace_window_problem, trace_data_problem, trace_qc_reasons,
  field_log_check = fl_check, field_log_notes = trimws(gsub("\\bNA\\b", "", fl_notes)))
write.csv(fin, PATH_FLUX, row.names = FALSE, na = "NA")
message("Wrote ", PATH_FLUX, ": ", nrow(fin), " measurements, ", ncol(fin), " columns")

dict_txt <- '
column,unit,description
measurement_id,,unique measurement ID: analyzer_date_tree_start (LGR/UGGA or LI-7810) or LEGACY_row for measurements without a recorded window
Tree,,ForestGEO tag of the tree (see trees.csv)
location,,site: Wetland (Black Gum Swamp) or Upland (EMS)
PLOT,,plot code as recorded: BGS (wetland); EMS or E5 (upland trees in ForestGEO sub-plot E5)
SPECIES,,species code: bg Nyssa sylvatica; hem Tsuga canadensis; rm Acer rubrum; ro Quercus rubra
species_label,,species name
DBH,cm,stem diameter at breast height
date,,measurement date (local)
datetime_posx,,field-log real start time as recorded (local wall clock; the trailing Z is not UTC)
sample_time_local,,corrected sampling time (America/New_York wall clock)
sample_hour_est,,sampling hour on the EST clock (UTC-5; no daylight saving) used to join the drivers
sample_time_src,,source of sample_time_local
year,,calendar year
season,,winter (DJF) / spring (MAM) / summer (JJA) / fall (SON)
sampling_round,,sampling campaign (campaigns separated by more than 7 days)
new_round,,TRUE for the first measurement of a campaign
days_since_last,d,days since the previous campaign
instrument,,analyzer: LGR1 / LGR3 (LGR microportable GLA131 / UGGA) or LI-7810
inst_label,,analyzer type: LGR/UGGA or LI-7810
window_src,,how the fitting window was defined
window_start,,start of the fitting window (local)
window_end,,end of the fitting window (local)
closure_s,s,length of the fitting window
closure_s_src,,source of closure_s
manual_window,,TRUE if the window was set by hand after inspection
end_trim_s,s,seconds removed from the logged end at a detected chamber opening
Area_cm2,cm2,collar area (inner radius 5.08 cm)
Vtot_L,L,system volume: collar + cap + tubing + 0.028 L analyzer volume
Tcham_C,degC,chamber air temperature (Fisher met station)
Pcham_kPa,kPa,air pressure (Fisher met station)
met_src,,source of Tcham_C and Pcham_kPa
flux_source,,goFlux (fitted to the raw record) or earlier processing (no raw record; linear fit)
fitted,,TRUE if fitted to the raw 1-Hz record in this processing
CH4_flux_nmolpm2ps,nmol m-2 s-1,stem CH4 flux (positive = emission)
CH4_SE,nmol m-2 s-1,standard error of the CH4 flux (of the selected model)
CH4_r2,,R2 of the linear CH4 fit
CH4_model,,model selected by goFlux best.flux: LM (linear) or HM (Hutchinson-Mosier)
CH4_g_factor,,ratio HM / LM flux
CO2_flux_umolpm2ps,umol m-2 s-1,stem CO2 flux
CO2_SE,umol m-2 s-1,standard error of the CO2 flux
CO2_r2,,R2 of the linear CO2 fit
CO2_model,,model selected for CO2
CH4_sigma_ppb,ppb,analyzer CH4 precision used for the MDF
CO2_sigma_ppm,ppm,analyzer CO2 precision used for the MDF
sigma_src,,source of the precision
CH4_MDF,nmol m-2 s-1,minimum detectable CH4 flux: 1.96 x analyzer precision / closure time x flux term (1.96 is a benchmark multiplier on concentration noise; not a calibrated 95 % detection test)
CH4_below_MDF,,TRUE if |CH4 flux| < CH4_MDF (retained)
CH4_det_class,,emission / uptake / below detection
CO2_MDF,umol m-2 s-1,minimum detectable CO2 flux
CO2_below_MDF,,TRUE if |CO2 flux| < CO2_MDF
qc_c0,,fluxqc: starting concentration well above ambient
qc_co2_tracer,,fluxqc: CO2 not increasing during the closure (possible poor seal)
qc_convex,,fluxqc: concentration curve convex (non-diffusive shape)
qc_noisy,,fluxqc: record noisier than the analyzer precision
qc_note,,fluxqc screens that fired
trace_window_problem,,trace QC: window problem (CO2 falling or dropping, flat start, overlap)
trace_data_problem,,trace QC: data problem (CH4 step, gap, start transient, noise)
trace_qc_reasons,,trace QC: all checks that fired
field_log_check,,field team check of the closure (field logs)
field_log_notes,,field team notes
'
# column,unit,description (descriptions may contain commas: split at the first two only)
dl <- trimws(strsplit(dict_txt, "\n")[[1]]); dl <- dl[nzchar(dl)][-1]
dict <- data.frame(column = sub(",.*$", "", dl), unit = sub("^[^,]*,([^,]*),.*$", "\\1", dl),
                   description = sub("^[^,]*,[^,]*,", "", dl), stringsAsFactors = FALSE)
stopifnot(setequal(dict$column, names(fin)))
write.csv(dict[match(names(fin), dict$column), ], file.path(FINAL, "stem_ch4_flux_dictionary.csv"), row.names = FALSE)

write_log()
logs <- c("01_closure_table", "02_fit_fluxes", "03_trace_qc", "06_flux_dataset")
all_log <- bind_rows(lapply(logs, function(s) read.csv(file.path(INTERIM, paste0("flux_log_", s, ".csv")), stringsAsFactors = FALSE)))
write.csv(all_log, PATH_LOG, row.names = FALSE)
message("Wrote ", PATH_LOG, " (", nrow(all_log), " rules)")
message("Detection: ", round(100 * mean(fin$CH4_below_MDF[fin$location == "Upland"])), "% upland and ",
        round(100 * mean(fin$CH4_below_MDF[fin$location == "Wetland"])), "% wetland fluxes below the MDF")

# ============================================================
# 02_fit_fluxes.R  (flux step 2 of 5)
#
# Fits CH4 and CO2 fluxes to every closure with a raw trace: goFlux linear (LM) and
# Hutchinson-Mosier (HM) models, best.flux selection with the Hueppi et al. (2018)
# criteria (g.limit = 2), and goFlux QC screens (jgewirtzman/goFlux fork, doi:10.5281/zenodo.23256675).
#   MDF = 1.96 * sigma / t * flux.term (1.96: benchmark multiplier, not a calibrated 95 % test); sigma = analyzer precision on
#   that day (goFlux empirical.prec, second differences inside the closure windows; 01_closure_table.R); t = closure
#   length in seconds. Retain-and-flag: nothing is removed here.
#
# Input : data/interim/flux_closures.rds (01_closure_table.R)
# Output: data/interim/flux_fits.csv       one row per closure
#         data/interim/flux_traces.rds     trace segments used (03-05)
#         data/final/flux_processing_settings.json
#         outputs/tables/flux_processing/  coverage, precision by day, QC review list
#         outputs/figures/flux_processing/ trace plots of every fitted closure
# ============================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(goFlux)
})
stopifnot(packageVersion("goFlux") >= "0.5.0.9002")
SCRIPT <- "02_fit_fluxes"
source("scripts/2_flux/flux_settings.R")

cl <- readRDS(PATH_CLOSURES)
closures <- cl$closures; segs <- cl$segments; day_sigma <- cl$day_sigma

# ============================================================
# PART 7: goFlux fits and QC screens
# ============================================================

message("\n=== Part 7: goFlux ===")

fit_ids <- closures$closure_id[closures$has_trace & !is.na(closures$Vtot)]
segs <- segs[fit_ids]
manID <- do.call(rbind, lapply(segs, function(s) {
  uid <- s$UniqueID[1]; c1 <- closures[closures$closure_id == uid, ]
  s$Vtot <- c1$Vtot; s$Area <- c1$Area; s$Pcham <- c1$Pcham; s$Tcham <- c1$Tcham
  s$DATE <- format(s$POSIX.time[1], "%Y-%m-%d")
  pr <- if (s$instrument[1] == "LI-7810") PREC_7810 else PREC_LGR
  s$CO2_prec <- pr[["CO2"]]; s$CH4_prec <- pr[["CH4"]]; s$H2O_prec <- pr[["H2O"]]
  s$H2O_ppm[is.na(s$H2O_ppm)] <- 0
  s$obs.length_corr <- s$obs.length; s$start.time_corr <- s$start.time; s$end.time_corr <- s$end.time
  s
}))
rownames(manID) <- NULL
manID$UniqueID <- as.character(manID$UniqueID)
message("Trace rows: ", nrow(manID), " for ", length(unique(manID$UniqueID)), " closures")
aux <- closures %>% filter(closure_id %in% fit_ids) %>% transmute(UniqueID = closure_id, instrument)

t0 <- Sys.time()
# The fit itself uses the *_prec columns set above (unchanged), so the fluxes do not depend on the
# precision estimator; det.* (goFlux's own detection columns) are kept for reference only. The MDF
# used in the dataset is computed below from the analyzer x day precision (01_closure_table.R).
res_co2 <- process.fluxes(manID, gastype = "CO2dry_ppm", auxfile = aux, by = "instrument", conf = 0.95, qc = FALSE)
message("CO2 done in ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
t0 <- Sys.time()
res_ch4 <- process.fluxes(manID, gastype = "CH4dry_ppb", auxfile = aux, by = "instrument", conf = 0.95,
                          co2.flux.result = res_co2$fluxes,
                          # qc.ambient is computed but not used: the windows begin after the closure onset
                          # (deadband), so the first seconds of a window are never at ambient by design
                          qc = list(min.obs = NULL, min.secs = 60))   # convex: quadratic test on the window (goFlux >= 0.5.0.9002)
message("CH4 done in ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
message("Campaign sigma (ppb) by analyzer (goFlux det.prec): ",
        paste(capture.output(print(res_ch4$fluxes %>% group_by(instrument) %>%
          summarise(n = n(), sigma_ppb = round(median(det.prec), 3), .groups = "drop"))), collapse = "\n"))

# closure length (goFlux closure.time: window span + one logging interval) and logging interval
cs <- manID %>% filter(flag == 1) %>% group_by(UniqueID) %>%
  summarise(t = closure.time(Etime), dt = median(diff(as.numeric(POSIX.time))), .groups = "drop")

# ============================================================
# PART 8: ASSEMBLE OUTPUT
# ============================================================

message("\n=== Part 8: output table ===")

pick_se <- function(f) ifelse(f$model == "HM" & !is.na(f$HM.SE), f$HM.SE, f$LM.SE)
gf_ch4 <- res_ch4$fluxes %>% transmute(
  closure_id = UniqueID,
  CH4_flux_goflux = best.flux, CH4_model = model, CH4_LM_flux = LM.flux, CH4_HM_flux = HM.flux,
  CH4_SE_goflux = pick_se(res_ch4$fluxes), CH4_LM_SE = LM.SE, CH4_LM_r2 = LM.r2, CH4_LM_p = LM.p.val,
  CH4_HM_r2 = HM.r2, CH4_g_factor = g.fact, CH4_HM_k = HM.k, CH4_k_max = k.max,
  CH4_C0 = C0, CH4_Ct = Ct, CH4_MAE_LM = LM.MAE, CH4_MAE_HM = HM.MAE,
  CH4_quality_check = quality.check, CH4_LM_diagnose = LM.diagnose, CH4_HM_diagnose = HM.diagnose,
  nb_obs = nb.obs, flux_term = flux.term,
  CH4_sigma_mad_record = det.prec, CH4_MDF_record = det.MDF, CH4_sigma_closure = det.prec.closure,
  qc_c0 = qc.c0, qc_c0_ratio = qc.c0.ratio,
  qc_co2_tracer = !co2.tracer,                  # goFlux co2.tracer is TRUE when CO2 rises (pass); flag = not rising
  qc_convex = qc.convex, qc_min_window = qc.min.secs, qc_noisy = qc.noisy, qc_noisy_ratio = qc.noisy.ratio) %>%
  mutate(qc_any = coalesce(qc_c0, FALSE) | coalesce(qc_co2_tracer, FALSE) | coalesce(qc_convex, FALSE) |
                  coalesce(qc_min_window, FALSE) | coalesce(qc_noisy, FALSE),
         qc_note = trimws(paste0(ifelse(coalesce(qc_c0, FALSE), "c0 ", ""), ifelse(coalesce(qc_co2_tracer, FALSE), "co2_tracer ", ""),
                                 ifelse(coalesce(qc_convex, FALSE), "convex ", ""), ifelse(coalesce(qc_min_window, FALSE), "min_window ", ""),
                                 ifelse(coalesce(qc_noisy, FALSE), "noisy", ""))))
gf_co2 <- res_co2$fluxes %>% transmute(
  closure_id = UniqueID,
  CO2_flux_goflux = best.flux, CO2_model = model, CO2_LM_flux = LM.flux, CO2_HM_flux = HM.flux,
  CO2_SE_goflux = pick_se(res_co2$fluxes), CO2_LM_r2 = LM.r2, CO2_LM_p = LM.p.val, CO2_g_factor = g.fact,
  CO2_C0 = C0, CO2_quality_check = quality.check,
  CO2_sigma_mad_record = det.prec, CO2_MDF_record = det.MDF)

write.csv(day_sigma, file.path(TAB_DIR, "goflux_precision_by_analyzer_day.csv"), row.names = FALSE)
message("Per-day precision (median CH4 ppb): ",
        paste(tapply(round(day_sigma$CH4_sigma_day, 3), day_sigma$instrument, median), collapse = " / "),
        " for ", paste(names(table(day_sigma$instrument)), collapse = " / "),
        "; days with >1 logging interval: ", sum(day_sigma$n_interval_runs > 1))
out <- closures %>%
  left_join(day_sigma %>% select(instrument, date, CH4_sigma_day, CO2_sigma_day), by = c("instrument", "date")) %>%
  left_join(cs %>% transmute(closure_id = UniqueID, t_sec = t, dt_s = dt), by = "closure_id") %>%
  left_join(gf_ch4, by = "closure_id") %>%
  left_join(gf_co2, by = "closure_id") %>%
  mutate(
    # reference detection limit: 1.96 x per-day analyzer sigma / closure seconds x flux term
    CH4_sigma_mad = coalesce(CH4_sigma_day, CH4_sigma_mad_record),
    CO2_sigma_mad = coalesce(CO2_sigma_day, CO2_sigma_mad_record),
    CH4_MDF = abs(qnorm(0.975) * CH4_sigma_mad / t_sec * flux_term),
    CO2_MDF = abs(qnorm(0.975) * CO2_sigma_mad / t_sec * flux_term),
    CH4_below_MDF = ifelse(is.na(CH4_MDF), NA, abs(CH4_flux_goflux) < CH4_MDF),
    CO2_below_MDF = ifelse(is.na(CO2_MDF), NA, abs(CO2_flux_goflux) < CO2_MDF),
    CH4_det_class = ifelse(is.na(CH4_MDF), NA_character_, ifelse(CH4_flux_goflux > CH4_MDF, "emission",
                           ifelse(CH4_flux_goflux < -CH4_MDF, "uptake", "below detection"))),
    CH4_MDF_method = ifelse(is.na(CH4_MDF), NA_character_, "1.96 x second-difference sigma (goFlux empirical.prec; analyzer x day) / t x flux term"),
    CH4_MDF_wass90 = abs(qnorm(0.95) * CH4_sigma_mad / t_sec * flux_term),
    CH4_MDF_wass95 = CH4_MDF,
    CH4_MDF_wass99 = abs(qnorm(0.995) * CH4_sigma_mad / t_sec * flux_term),
    year = as.integer(format(date, "%Y")),
    fitted = !is.na(CH4_flux_goflux),
    # flux of the analysis dataset: the goFlux fit (every measurement has a raw record)
    CH4_flux_nmolpm2ps = CH4_flux_goflux, CH4_SE = CH4_SE_goflux, CH4_r2 = CH4_LM_r2,
    CO2_flux_umolpm2ps = CO2_flux_goflux, CO2_SE = CO2_SE_goflux, CO2_r2 = CO2_LM_r2,
    inst_label = ifelse(instrument == "LI-7810", "LI-7810", "LGR/UGGA"),
    season = case_when(month(date) %in% c(12, 1, 2) ~ "winter", month(date) %in% 3:5 ~ "spring",
                       month(date) %in% 6:8 ~ "summer", TRUE ~ "fall")
  ) %>%
  rename(Tree = tree, Vtot_L = Vtot, Vcollar_L = Vcollar, Vextra_L = Vextra, Area_cm2 = Area,
         Pcham_kPa = Pcham, Tcham_C = Tcham) %>%
  arrange(window_start, Tree)

fmt_t <- function(x) format(x, "%Y-%m-%d %H:%M:%S")
out_csv <- out %>% mutate(across(where(function(v) inherits(v, "POSIXct")), fmt_t),
                          date = format(date, "%Y-%m-%d"))
write.csv(out_csv, PATH_FITS, row.names = FALSE, na = "NA")
saveRDS(list(traces = manID, closures = closures), PATH_TRACES)
writeLines(jsonlite::toJSON(list(CH4 = res_ch4$settings, CO2 = res_co2$settings,
                                 constants = list(SURFAREA_M2 = SURFAREA_M2, EXTRA_VOL_LGR_L = EXTRA_VOL_LGR_L, EXTRA_VOL_7810_L = EXTRA_VOL_7810_L,
                                                  SHOULDER_S = SHOULDER_S, DEADBAND_S = DEADBAND_S, TZ = TZ)),
                           auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null", force = TRUE), PATH_SETTINGS)
log_step("7 fits", "closures fitted (goFlux best.flux, LM or HM)", sum(out$fitted),
         sprintf("HM selected for %d", sum(out$CH4_model %in% "HM")))
log_step("7 fits", "measurements without a fit", sum(!out$fitted))
message("Wrote ", PATH_FITS, ": ", nrow(out), " rows (fitted: ", sum(out$fitted), ")")

# ============================================================
# PART 9: REPORTS
# ============================================================

message("\n=== Part 9: reports ===")

# QC review list: everything a goFlux screen or diagnostic caught
review <- out %>%
  filter(fitted) %>%
  mutate(reasons = paste0(
    ifelse(!is.na(qc_any) & qc_any, paste0("goFlux QC: ", qc_note, "; "), ""),
    ifelse(!is.na(CO2_flux_goflux) & CO2_flux_goflux <= 0, "CO2 flux <= 0; ", ""),
    "")) %>%
  filter(nz(reasons)) %>%
  mutate(n_reasons = lengths(regmatches(reasons, gregexpr("; ", reasons)))) %>%
  arrange(desc(n_reasons), date) %>%
  transmute(closure_id, date, Tree, PLOT, SPECIES, instrument, window_start, window_end,
            t_sec, CH4_flux_goflux, CH4_model, CH4_g_factor, CH4_MDF, CH4_det_class,
            CO2_flux_goflux, qc_any, n_reasons, reasons)
write.csv(review %>% mutate(across(where(function(v) inherits(v, "POSIXct")), fmt_t)),
          file.path(TAB_DIR, "goflux_qc_review.csv"), row.names = FALSE)
message("QC review list: ", nrow(review), " closures -> outputs/tables/goflux_qc_review.csv")
message("  screens fired: ", paste(
  vapply(c("qc_c0", "qc_co2_tracer", "qc_convex", "qc_min_window", "qc_noisy"),
         function(s) sprintf("%s=%d", sub("qc_", "", s), sum(out[[s]], na.rm = TRUE)), character(1)), collapse = ", "))

# 9d. goFlux trace plots for every fitted closure (one pdf per analyzer x year)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  grp <- aux %>% left_join(out %>% select(closure_id, year), by = c("UniqueID" = "closure_id")) %>%
    mutate(g = paste0(gsub("[^A-Za-z0-9]", "", instrument), "_", year))
  for (gg in unique(grp$g)) {
    ids <- grp$UniqueID[grp$g == gg]
    pl <- tryCatch(flux.plot(res_ch4$fluxes[res_ch4$fluxes$UniqueID %in% ids, ],
                             manID[manID$UniqueID %in% ids, ], "CH4dry_ppb", shoulder = SHOULDER_S),
                   error = function(e) { message("flux.plot failed for ", gg, ": ", conditionMessage(e)); NULL })
    if (!is.null(pl)) invisible(capture.output(flux2pdf(pl, outfile = file.path(FIG_DIR, paste0("flux_plots_CH4_", gg, ".pdf")))))
  }
  message("Trace plots written to ", FIG_DIR)
}

write_log()
message("\nDone: ", PATH_FITS)

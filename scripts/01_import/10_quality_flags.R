# ============================================================
# 10_quality_flags.R
#
# Build the analysis dataset (data/processed/flux_with_quality_flags.csv)
# from the goFlux-reprocessed closures written by 09_goflux_reprocess.R.
#
# Detection convention (guidelines paper / fluxqc >= 0.2.3):
#   MDF = z * sigma / t * flux.term,  z = qnorm(0.975) = 1.96 (95 %, two-sided)
#   sigma = MAD of first differences over each analyzer's whole record per
#           constant-interval run (fluxqc::flag_detection, precision = "mad")
#   t     = closure length in seconds (fluxqc::closure_seconds), NOT nb.obs
#   Retain-and-flag: nothing is deleted; CH4_below_MDF is the canonical flag.
#
# Also carried for the sensitivity analysis (07_filter_sensitivity_ridges.R):
#   CH4_MDF_manufacturer  datasheet precision / t * flux.term (goFlux form)
#   CH4_MDF_wass{90,95,99} campaign (record-wide MAD) sigma, z at 90/95/99 %
#   CH4_MDF_chr{90,95,99}  per-closure Allan sigma x 3 x t_crit / t * flux.term.
#                          This is the old "Christiansen" column, kept as a
#                          conservative comparison only; Christiansen et al.
#                          (2015) define MDF = precision / t * flux.term.
#
# Rows without a raw trace (no LGR file that day; Oct 2025 LI-7810 file not
# archived) keep their legacy flux and get an MDF from the analyzer's pooled
# sigma and the field-log window length (t_src says which).
#
# Input : data/input/HF_2023-2025_tree_flux_goflux.csv
# Output: data/processed/flux_with_quality_flags.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(fluxqc)
  library(lubridate)
})

# Closures with a raw trace but no row in the legacy (published) dataset are
# kept out of the analysis dataset by default so n stays comparable; set to
# TRUE (or GOFLUX_INCLUDE_NEW=1) to analyse every fitted closure.
INCLUDE_NEW_CLOSURES <- Sys.getenv("GOFLUX_INCLUDE_NEW", "0") == "1"

PREC_CH4 <- c("LGR/UGGA" = 0.9, "LI-7810" = 0.6)   # ppb, datasheet 1-s precision
PREC_CO2 <- c("LGR/UGGA" = 0.35, "LI-7810" = 3.5)  # ppm
T_FALLBACK <- c("LGR/UGGA" = 180, "LI-7810" = 60)  # s: usual field-log window; legacy 2025 60-s window
z90 <- qnorm(0.95); z95 <- qnorm(0.975); z99 <- qnorm(0.995)

message("=== Loading goFlux-reprocessed closures ===")
g <- read.csv(file.path("data", "input", "HF_2023-2025_tree_flux_goflux.csv"),
              stringsAsFactors = FALSE)
message("Closures: ", nrow(g), " (fitted: ", sum(g$fitted), "; in legacy dataset: ",
        sum(g$in_legacy_dataset), ")")

# Optional exclusion list for sensitivity runs (GOFLUX_EXCLUDE = path to a CSV with closure_id)
excl_ids <- character(0)
if (nzchar(Sys.getenv("GOFLUX_EXCLUDE"))) excl_ids <- read.csv(Sys.getenv("GOFLUX_EXCLUDE"), stringsAsFactors = FALSE)$closure_id
if ("manual_exclude" %in% names(g)) excl_ids <- union(excl_ids, g$closure_id[g$manual_exclude %in% TRUE])
df <- g %>% filter(in_legacy_dataset | INCLUDE_NEW_CLOSURES) %>% filter(!closure_id %in% excl_ids)
if (length(excl_ids)) message("Excluded closures (manual review / sensitivity list): ", length(excl_ids))
message("Analysis rows: ", nrow(df), if (INCLUDE_NEW_CLOSURES) " (new closures included)" else
        " (legacy rows only; new closures excluded)")

# ------------------------------------------------------------
# Pooled (campaign) sigma per analyzer, from the fitted closures
# ------------------------------------------------------------
pooled <- g %>% filter(fitted) %>% group_by(inst_label) %>%
  summarise(sigma_pooled_CH4 = median(CH4_sigma_mad, na.rm = TRUE),
            sigma_pooled_CO2 = median(CO2_sigma_mad, na.rm = TRUE), .groups = "drop")
message("Pooled MAD sigma (ppb CH4 / ppm CO2): ",
        paste(sprintf("%s = %.3f / %.3f", pooled$inst_label, pooled$sigma_pooled_CH4, pooled$sigma_pooled_CO2),
              collapse = "; "))

# ------------------------------------------------------------
# Assemble
# ------------------------------------------------------------
df <- df %>%
  left_join(pooled, by = "inst_label") %>%
  mutate(
    UniqueID   = closure_id,
    surfarea   = Area_cm2 / 1e4,
    vol_system = Vtot_L,
    # SE-based SNR (CH4_SE is the goFlux SE, or the unit-corrected legacy SE)
    CH4_SE_corr = CH4_SE,
    CO2_snr_se  = abs(CO2_flux_umolpm2ps) / CO2_SE,
    CH4_snr_se  = abs(CH4_flux_nmolpm2ps) / CH4_SE,
    # closure length: trace-based when fitted, else the field-log window, else a fallback
    t_sec_fl  = as.numeric(difftime(as.POSIXct(window_end), as.POSIXct(window_start), units = "secs")),
    t_sec_est = case_when(fitted ~ t_sec, !is.na(t_sec_fl) ~ t_sec_fl,
                          TRUE ~ unname(T_FALLBACK[inst_label])),
    t_src     = case_when(fitted ~ "trace (closure_seconds)", !is.na(t_sec_fl) ~ "field-log window",
                          TRUE ~ paste0("assumed ", unname(T_FALLBACK[inst_label]), " s")),
    dt_s      = ifelse(fitted, dt_s, NA_real_),
    n_pts     = nb_obs,
    instrument = coalesce(instrument, ifelse(inst_label == "LI-7810", "LI-7810", "LGR1")),
    # flux term for rows that were not fitted (same geometry convention)
    flux_term = ifelse(fitted, flux_term, fluxqc::flux_term(Vtot_L, Pcham_kPa, Area_cm2, Tcham_C)),
    # sigmas
    sigma_CH4   = ifelse(fitted, CH4_sigma_mad, sigma_pooled_CH4),   # campaign / pooled MAD sigma (ppb)
    sigma_src   = ifelse(fitted, "record-wide MAD per analyzer", "pooled analyzer median (no trace)"),
    allan_sd_CH4 = CH4_sigma_allan,                                  # per-closure Allan sigma (ppb)
    allan_sd_CO2 = CO2_sigma_allan,
    prec_ch4 = unname(PREC_CH4[inst_label]),
    # noise floor / SNR from the per-closure Allan sigma (diagnostic, as before)
    CH4_noise_floor = allan_sd_CH4 / t_sec_est * flux_term,
    CO2_noise_floor = allan_sd_CO2 / t_sec_est * flux_term,
    CH4_snr_allan   = ifelse(!is.na(CH4_noise_floor) & CH4_noise_floor > 0,
                             abs(CH4_flux_nmolpm2ps) / CH4_noise_floor, NA_real_),
    # MDFs (all in nmol m-2 s-1)
    CH4_MDF_manufacturer = prec_ch4 / t_sec_est * flux_term,
    CH4_MDF_wass90 = z90 * sigma_CH4 / t_sec_est * flux_term,
    CH4_MDF_wass95 = z95 * sigma_CH4 / t_sec_est * flux_term,
    CH4_MDF_wass99 = z99 * sigma_CH4 / t_sec_est * flux_term,
    df_meas = pmax(n_pts - 2, 1),
    CH4_MDF_chr90 = 3 * allan_sd_CH4 * qt(0.95,  df_meas) / t_sec_est * flux_term,
    CH4_MDF_chr95 = 3 * allan_sd_CH4 * qt(0.975, df_meas) / t_sec_est * flux_term,
    CH4_MDF_chr99 = 3 * allan_sd_CH4 * qt(0.995, df_meas) / t_sec_est * flux_term,
    across(c(CH4_MDF_manufacturer, starts_with("CH4_MDF_wass"), starts_with("CH4_MDF_chr")), abs),
    CH4_below_MDF_manuf  = abs(CH4_flux_nmolpm2ps) < CH4_MDF_manufacturer,
    CH4_below_MDF_wass90 = abs(CH4_flux_nmolpm2ps) < CH4_MDF_wass90,
    CH4_below_MDF_wass95 = abs(CH4_flux_nmolpm2ps) < CH4_MDF_wass95,
    CH4_below_MDF_wass99 = abs(CH4_flux_nmolpm2ps) < CH4_MDF_wass99,
    CH4_below_MDF_chr90  = abs(CH4_flux_nmolpm2ps) < CH4_MDF_chr90,
    CH4_below_MDF_chr95  = abs(CH4_flux_nmolpm2ps) < CH4_MDF_chr95,
    CH4_below_MDF_chr99  = abs(CH4_flux_nmolpm2ps) < CH4_MDF_chr99,
    # canonical detection flag: 95 % campaign-sigma MDF (= fluxqc MDF_emp for fitted rows)
    CH4_below_MDF = CH4_below_MDF_wass95,
    CH4_det_class = ifelse(is.na(CH4_MDF_wass95), NA_character_,
                           ifelse(CH4_flux_nmolpm2ps > CH4_MDF_wass95, "emission",
                                  ifelse(CH4_flux_nmolpm2ps < -CH4_MDF_wass95, "uptake", "below detection")))
  ) %>%
  select(-df_meas, -prec_ch4)

# sanity: for fitted rows the recomputed 95 % MDF must equal fluxqc's MDF_emp
chk <- df %>% filter(fitted) %>% summarise(max_rel_diff = max(abs(CH4_MDF_wass95 - CH4_MDF) / CH4_MDF, na.rm = TRUE))
message("Max relative difference between recomputed and fluxqc MDF (fitted rows): ",
        signif(chk$max_rel_diff, 3))

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------
report <- function(col, label) {
  v <- df[[col]]; n_eval <- sum(!is.na(v)); n_below <- sum(v, na.rm = TRUE)
  message(sprintf("  %-26s %4d/%d evaluated, %4d below (%.1f%%)", label, n_eval, nrow(df), n_below,
                  100 * n_below / max(n_eval, 1)))
}
message("\n=== MDF summary (retain-and-flag) ===")
report("CH4_below_MDF_manuf",  "Datasheet (goFlux form)")
report("CH4_below_MDF_wass90", "Campaign sigma 90%")
report("CH4_below_MDF_wass95", "Campaign sigma 95% (ref)")
report("CH4_below_MDF_wass99", "Campaign sigma 99%")
report("CH4_below_MDF_chr90",  "Per-closure sigma x3t 90%")
report("CH4_below_MDF_chr95",  "Per-closure sigma x3t 95%")
report("CH4_below_MDF_chr99",  "Per-closure sigma x3t 99%")
message("\nBy analyzer (reference MDF):")
print(as.data.frame(df %>% group_by(inst_label) %>%
  summarise(n = n(), fitted = sum(fitted), pct_below = round(100 * mean(CH4_below_MDF, na.rm = TRUE), 1),
            median_t_s = median(t_sec_est), median_sigma_ppb = round(median(sigma_CH4, na.rm = TRUE), 3),
            median_MDF = signif(median(CH4_MDF_wass95, na.rm = TRUE), 3), .groups = "drop")))
message("Flux source: ", paste(names(table(df$flux_source)), table(df$flux_source), sep = " = ", collapse = ", "))
message("t source:    ", paste(names(table(df$t_src)), table(df$t_src), sep = " = ", collapse = ", "))

# ------------------------------------------------------------
# Sampling time for joining environmental drivers
#   datetime_posx (= real_start) is the field-log real time, local wall clock (EDT in
#   summer), stored with a "Z" label. Exceptions:
#     - 146 closures (2024-03-06, Oct-Nov 2024) carry 00:00 because no real time was
#       logged; use the analyzer-clock closure start minus that analyzer's clock offset
#       (median real - analyzer difference over the nearest dated closures, +-60 d)
#     - 2 closures on 2023-06-29 were logged as 02:xx instead of 14:xx (AM/PM)
#   sample_time_local: America/New_York wall clock
#   sample_hour_est:   hour on the EST clock (UTC-5, no DST) used by every met series in
#                      aligned_hourly_dataset.csv; stored as "YYYY-MM-DD HH:00:00"
# ------------------------------------------------------------
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
loc <- force_tz(fix_t, "America/New_York")
df$sample_hour_est <- format(floor_date(with_tz(loc, "EST"), "hour"), "%Y-%m-%d %H:00:00")
message("\nSampling time source: ", paste(names(table(df$sample_time_src)), table(df$sample_time_src), sep = " = ", collapse = ", "))

output_path <- file.path("data", "processed", "flux_with_quality_flags.csv")
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
write.csv(df, output_path, row.names = FALSE, na = "NA")
message("\n=== Output saved ===\n  ", output_path, "\n  Rows: ", nrow(df),
        "\n  Primary flag: CH4_below_MDF (campaign-sigma 95 %, t in seconds)")

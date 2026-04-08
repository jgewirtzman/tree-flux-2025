# ============================================================
# 00_data_summary.R
#
# Summary statistics for the flagged flux dataset, including
# MDF detection rates, instrument comparison, and QC metrics.
#
# Run AFTER: scripts/01_import/09_quality_flags.R
#
# Outputs:
#   - Console: formatted summary
#   - outputs/tables/flux_data_summary.csv (overall stats)
#   - outputs/tables/mdf_detection_by_instrument.csv
#   - outputs/tables/flux_by_species_site.csv
# ============================================================

library(dplyr)
library(tidyr)

# ============================================================
# LOAD DATA
# ============================================================

flagged_path <- file.path("data", "processed", "flux_with_quality_flags.csv")
if (!file.exists(flagged_path)) {
  stop("flux_with_quality_flags.csv not found. Run scripts/01_import/09_quality_flags.R first.")
}

df <- read.csv(flagged_path, stringsAsFactors = FALSE)

out_dir <- file.path("outputs", "tables")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. OVERALL COUNTS
# ============================================================

cat("\n", strrep("=", 60), "\n")
cat("  FLUX DATA SUMMARY\n")
cat(strrep("=", 60), "\n\n")

n_total <- nrow(df)
n_lgr   <- sum(df$year < 2025)
n_7810  <- sum(df$year == 2025)

cat("Total measurements:     ", n_total, "\n")
cat("  LGR/UGGA (2023-24):  ", n_lgr, "\n")
cat("  LI-7810 (2025):      ", n_7810, "\n\n")

# By year
cat("By year:\n")
yr_counts <- df %>% count(year)
for (i in seq_len(nrow(yr_counts))) {
  cat(sprintf("  %d: %d\n", yr_counts$year[i], yr_counts$n[i]))
}

# By site
cat("\nBy site:\n")
if ("PLOT" %in% names(df)) {
  site_counts <- df %>%
    filter(!is.na(PLOT)) %>%
    count(PLOT)
  for (i in seq_len(nrow(site_counts))) {
    cat(sprintf("  %s: %d\n", site_counts$PLOT[i], site_counts$n[i]))
  }
}

# By species
cat("\nBy species:\n")
if ("SPECIES" %in% names(df)) {
  sp_counts <- df %>%
    filter(!is.na(SPECIES)) %>%
    count(SPECIES)
  for (i in seq_len(nrow(sp_counts))) {
    cat(sprintf("  %s: %d\n", sp_counts$SPECIES[i], sp_counts$n[i]))
  }
}

# ============================================================
# 2. CH4 FLUX DISTRIBUTION
# ============================================================

cat("\n", strrep("-", 60), "\n")
cat("  CH4 FLUX DISTRIBUTION (nmol m-2 s-1)\n")
cat(strrep("-", 60), "\n\n")

flux_stats <- function(x, label) {
  x <- x[!is.na(x)]
  cat(sprintf("%-25s n=%d  mean=%.3f  median=%.3f  SD=%.3f  [%.3f, %.3f]  IQR=%.3f\n",
              label, length(x), mean(x), median(x), sd(x),
              min(x), max(x), IQR(x)))
}

flux_stats(df$CH4_flux_nmolpm2ps, "All measurements")
flux_stats(df$CH4_flux_nmolpm2ps[df$year < 2025], "LGR/UGGA (2023-24)")
flux_stats(df$CH4_flux_nmolpm2ps[df$year == 2025], "LI-7810 (2025)")

# ============================================================
# 3. SIGN ANALYSIS
# ============================================================

cat("\n", strrep("-", 60), "\n")
cat("  SIGN ANALYSIS\n")
cat(strrep("-", 60), "\n\n")

sign_summary <- df %>%
  mutate(inst = ifelse(year < 2025, "LGR/UGGA", "LI-7810")) %>%
  group_by(inst) %>%
  summarise(
    n = n(),
    n_positive = sum(CH4_flux_nmolpm2ps > 0, na.rm = TRUE),
    n_zero     = sum(CH4_flux_nmolpm2ps == 0, na.rm = TRUE),
    n_negative = sum(CH4_flux_nmolpm2ps < 0, na.rm = TRUE),
    pct_negative = round(100 * n_negative / n, 1),
    .groups = "drop"
  )

# Overall
n_neg <- sum(df$CH4_flux_nmolpm2ps < 0, na.rm = TRUE)
cat(sprintf("Overall: %d/%d negative (%.1f%%)\n\n", n_neg, n_total,
            100 * n_neg / n_total))

cat("By instrument:\n")
for (i in seq_len(nrow(sign_summary))) {
  cat(sprintf("  %-15s %d/%d negative (%.1f%%)\n",
              sign_summary$inst[i], sign_summary$n_negative[i],
              sign_summary$n[i], sign_summary$pct_negative[i]))
}

# ============================================================
# 4. ALLAN DEVIATION COVERAGE
# ============================================================

cat("\n", strrep("-", 60), "\n")
cat("  ALLAN DEVIATION\n")
cat(strrep("-", 60), "\n\n")

n_with_allan <- sum(!is.na(df$allan_sd_CH4))
cat(sprintf("Coverage: %d/%d measurements (%.1f%%)\n",
            n_with_allan, n_total, 100 * n_with_allan / n_total))

allan_by_inst <- df %>%
  mutate(inst = ifelse(year < 2025, "LGR/UGGA", "LI-7810")) %>%
  filter(!is.na(allan_sd_CH4)) %>%
  group_by(inst) %>%
  summarise(
    n = n(),
    median_ch4 = round(median(allan_sd_CH4, na.rm = TRUE), 4),
    q25_ch4    = round(quantile(allan_sd_CH4, 0.25, na.rm = TRUE), 4),
    q75_ch4    = round(quantile(allan_sd_CH4, 0.75, na.rm = TRUE), 4),
    min_ch4    = round(min(allan_sd_CH4, na.rm = TRUE), 4),
    max_ch4    = round(max(allan_sd_CH4, na.rm = TRUE), 4),
    .groups = "drop"
  )

cat("\nCH4 Allan SD (ppb) by instrument:\n")
for (i in seq_len(nrow(allan_by_inst))) {
  cat(sprintf("  %-15s n=%d  median=%.4f  IQR=[%.4f, %.4f]  range=[%.4f, %.4f]\n",
              allan_by_inst$inst[i], allan_by_inst$n[i],
              allan_by_inst$median_ch4[i],
              allan_by_inst$q25_ch4[i], allan_by_inst$q75_ch4[i],
              allan_by_inst$min_ch4[i], allan_by_inst$max_ch4[i]))
}

cat("\nManufacturer specs for comparison:\n")
cat("  LGR/UGGA: 0.9 ppb  |  LI-7810: 0.6 ppb\n")

# ============================================================
# 5. MDF DETECTION RATES
# ============================================================

cat("\n", strrep("-", 60), "\n")
cat("  MDF DETECTION RATES (% below MDF threshold)\n")
cat(strrep("-", 60), "\n\n")

mdf_cols <- c(
  "Manufacturer MDF"  = "CH4_below_MDF_manuf",
  "Wassmann 90%"      = "CH4_below_MDF_wass90",
  "Wassmann 95%"      = "CH4_below_MDF_wass95",
  "Wassmann 99%"      = "CH4_below_MDF_wass99",
  "Christiansen 90%"  = "CH4_below_MDF_chr90",
  "Christiansen 95%"  = "CH4_below_MDF_chr95",
  "Christiansen 99%"  = "CH4_below_MDF_chr99"
)

mdf_results <- list()

for (label in names(mdf_cols)) {
  col <- mdf_cols[label]
  if (!col %in% names(df)) next

  vals <- df[[col]]
  n_eval <- sum(!is.na(vals))
  n_below <- sum(vals == TRUE, na.rm = TRUE)
  pct <- ifelse(n_eval > 0, round(100 * n_below / n_eval, 1), NA)

  # By instrument
  lgr_vals <- vals[df$year < 2025]
  li_vals  <- vals[df$year == 2025]
  lgr_eval <- sum(!is.na(lgr_vals))
  lgr_below <- sum(lgr_vals == TRUE, na.rm = TRUE)
  li_eval  <- sum(!is.na(li_vals))
  li_below <- sum(li_vals == TRUE, na.rm = TRUE)

  cat(sprintf("  %-22s Overall: %d/%d (%.1f%%)  |  LGR: %d/%d (%.1f%%)  |  LI-7810: %d/%d (%.1f%%)\n",
              label,
              n_below, n_eval, ifelse(is.na(pct), 0, pct),
              lgr_below, lgr_eval,
              ifelse(lgr_eval > 0, round(100 * lgr_below / lgr_eval, 1), 0),
              li_below, li_eval,
              ifelse(li_eval > 0, round(100 * li_below / li_eval, 1), 0)))

  mdf_results[[length(mdf_results) + 1]] <- data.frame(
    filter = label,
    n_eval_total = n_eval, n_below_total = n_below, pct_below_total = pct,
    n_eval_lgr = lgr_eval, n_below_lgr = lgr_below,
    pct_below_lgr = ifelse(lgr_eval > 0, round(100 * lgr_below / lgr_eval, 1), NA),
    n_eval_7810 = li_eval, n_below_7810 = li_below,
    pct_below_7810 = ifelse(li_eval > 0, round(100 * li_below / li_eval, 1), NA),
    stringsAsFactors = FALSE
  )
}

mdf_table <- do.call(rbind, mdf_results)

# ============================================================
# 6. SNR SUMMARY
# ============================================================

cat("\n", strrep("-", 60), "\n")
cat("  SIGNAL-TO-NOISE RATIO\n")
cat(strrep("-", 60), "\n\n")

# Allan-based SNR
if ("CH4_snr_allan" %in% names(df)) {
  snr_by_inst <- df %>%
    mutate(inst = ifelse(year < 2025, "LGR/UGGA", "LI-7810")) %>%
    filter(!is.na(CH4_snr_allan)) %>%
    group_by(inst) %>%
    summarise(
      n = n(),
      median_snr = round(median(CH4_snr_allan), 2),
      pct_gt2 = round(100 * mean(CH4_snr_allan > 2), 1),
      pct_gt3 = round(100 * mean(CH4_snr_allan > 3), 1),
      .groups = "drop"
    )

  cat("Allan-based CH4 SNR:\n")
  for (i in seq_len(nrow(snr_by_inst))) {
    cat(sprintf("  %-15s n=%d  median=%.2f  >2: %.1f%%  >3: %.1f%%\n",
                snr_by_inst$inst[i], snr_by_inst$n[i],
                snr_by_inst$median_snr[i],
                snr_by_inst$pct_gt2[i], snr_by_inst$pct_gt3[i]))
  }
}

# SE-based SNR
if ("CH4_snr_se" %in% names(df)) {
  cat("\nSE-based CH4 SNR:\n")
  snr_se_by_inst <- df %>%
    mutate(inst = ifelse(year < 2025, "LGR/UGGA", "LI-7810")) %>%
    filter(!is.na(CH4_snr_se)) %>%
    group_by(inst) %>%
    summarise(
      n = n(),
      median_snr = round(median(CH4_snr_se), 2),
      pct_gt2 = round(100 * mean(CH4_snr_se > 2), 1),
      pct_gt3 = round(100 * mean(CH4_snr_se > 3), 1),
      .groups = "drop"
    )
  for (i in seq_len(nrow(snr_se_by_inst))) {
    cat(sprintf("  %-15s n=%d  median=%.2f  >2: %.1f%%  >3: %.1f%%\n",
                snr_se_by_inst$inst[i], snr_se_by_inst$n[i],
                snr_se_by_inst$median_snr[i],
                snr_se_by_inst$pct_gt2[i], snr_se_by_inst$pct_gt3[i]))
  }
}

# ============================================================
# 7. R-SQUARED DISTRIBUTION
# ============================================================

cat("\n", strrep("-", 60), "\n")
cat("  R-SQUARED DISTRIBUTION\n")
cat(strrep("-", 60), "\n\n")

r2_summary <- df %>%
  mutate(inst = ifelse(year < 2025, "LGR/UGGA", "LI-7810")) %>%
  group_by(inst) %>%
  summarise(
    n = n(),
    median_ch4_r2 = round(median(CH4_r2, na.rm = TRUE), 3),
    pct_gt05 = round(100 * mean(CH4_r2 > 0.5, na.rm = TRUE), 1),
    pct_gt07 = round(100 * mean(CH4_r2 > 0.7, na.rm = TRUE), 1),
    pct_gt09 = round(100 * mean(CH4_r2 > 0.9, na.rm = TRUE), 1),
    median_co2_r2 = round(median(CO2_r2, na.rm = TRUE), 3),
    .groups = "drop"
  )

cat("CH4 R2 by instrument:\n")
for (i in seq_len(nrow(r2_summary))) {
  cat(sprintf("  %-15s median=%.3f  >0.5: %.1f%%  >0.7: %.1f%%  >0.9: %.1f%%\n",
              r2_summary$inst[i], r2_summary$median_ch4_r2[i],
              r2_summary$pct_gt05[i], r2_summary$pct_gt07[i],
              r2_summary$pct_gt09[i]))
}

# ============================================================
# 8. FLUX BY SPECIES AND SITE
# ============================================================

cat("\n", strrep("-", 60), "\n")
cat("  CH4 FLUX BY SPECIES AND SITE\n")
cat(strrep("-", 60), "\n\n")

flux_by_group <- df %>%
  filter(!is.na(SPECIES), !is.na(PLOT)) %>%
  mutate(location = ifelse(PLOT == "BGS", "Wetland", "Upland")) %>%
  group_by(SPECIES, location) %>%
  summarise(
    n = n(),
    mean_ch4   = round(mean(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    median_ch4 = round(median(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    sd_ch4     = round(sd(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    pct_neg    = round(100 * mean(CH4_flux_nmolpm2ps < 0, na.rm = TRUE), 1),
    pct_below_wass95 = round(100 * mean(CH4_below_MDF_wass95 == TRUE, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(location, SPECIES)

for (i in seq_len(nrow(flux_by_group))) {
  cat(sprintf("  %-5s %-8s  n=%d  mean=%.3f  median=%.3f  SD=%.3f  neg=%.1f%%  <Wass95=%.1f%%\n",
              flux_by_group$SPECIES[i], flux_by_group$location[i],
              flux_by_group$n[i], flux_by_group$mean_ch4[i],
              flux_by_group$median_ch4[i], flux_by_group$sd_ch4[i],
              flux_by_group$pct_neg[i], flux_by_group$pct_below_wass95[i]))
}

# ============================================================
# SAVE CSV TABLES
# ============================================================

# Overall summary
overall <- data.frame(
  metric = c("Total measurements", "LGR/UGGA (2023-24)", "LI-7810 (2025)",
             "CH4 mean (nmol/m2/s)", "CH4 median", "CH4 SD",
             "% negative (overall)", "% negative (LGR)", "% negative (LI-7810)",
             "Allan deviation coverage"),
  value = c(n_total, n_lgr, n_7810,
            round(mean(df$CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
            round(median(df$CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
            round(sd(df$CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
            round(100 * n_neg / n_total, 1),
            round(100 * sum(df$CH4_flux_nmolpm2ps[df$year < 2025] < 0) / n_lgr, 1),
            round(100 * sum(df$CH4_flux_nmolpm2ps[df$year == 2025] < 0) / n_7810, 1),
            round(100 * n_with_allan / n_total, 1)),
  stringsAsFactors = FALSE
)

write.csv(overall, file.path(out_dir, "flux_data_summary.csv"), row.names = FALSE)
write.csv(mdf_table, file.path(out_dir, "mdf_detection_by_instrument.csv"), row.names = FALSE)
write.csv(flux_by_group, file.path(out_dir, "flux_by_species_site.csv"), row.names = FALSE)

cat("\n", strrep("=", 60), "\n")
cat("  Tables saved to: ", out_dir, "\n")
cat("    flux_data_summary.csv\n")
cat("    mdf_detection_by_instrument.csv\n")
cat("    flux_by_species_site.csv\n")
cat(strrep("=", 60), "\n")

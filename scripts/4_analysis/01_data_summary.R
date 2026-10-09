# ============================================================
# 01_data_summary.R
#
# Summary of the final stem flux dataset: measurements by year, site, species and
# analyzer; flux distributions; detection against the MDF; analyzer precision;
# model selection and QC flags.
#
# Input : data/final/stem_ch4_flux.csv (scripts/2_flux/06_flux_dataset.R)
# Output: outputs/tables/flux_data_summary.csv     counts, detection and fit statistics
#         outputs/tables/flux_by_species_site.csv  flux distribution by site x species
#         outputs/tables/detection_by_analyzer.csv precision, MDF and detection by analyzer
# ============================================================
suppressPackageStartupMessages({ library(dplyr) })

df <- read.csv(file.path("data", "final", "stem_ch4_flux.csv"), stringsAsFactors = FALSE)
OUT <- file.path("outputs", "tables"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

say <- function(...) cat(sprintf(...), "\n", sep = "")
say("Measurements: %d (%d trees, %d dates)", nrow(df), n_distinct(df$Tree), n_distinct(df$date))
print(table(df$year, df$location))
print(table(df$inst_label, df$flux_source))

by_sp <- df %>% group_by(location, SPECIES) %>%
  summarise(n = n(), trees = n_distinct(Tree), mean = mean(CH4_flux_nmolpm2ps), median = median(CH4_flux_nmolpm2ps),
            sd = sd(CH4_flux_nmolpm2ps), q25 = quantile(CH4_flux_nmolpm2ps, 0.25), q75 = quantile(CH4_flux_nmolpm2ps, 0.75),
            pct_positive = 100 * mean(CH4_flux_nmolpm2ps > 0), pct_below_MDF = 100 * mean(CH4_below_MDF), .groups = "drop") %>%
  mutate(across(where(is.double), ~ signif(.x, 4)))
print(as.data.frame(by_sp))
write.csv(by_sp, file.path(OUT, "flux_by_species_site.csv"), row.names = FALSE)

det <- df %>% group_by(inst_label, location) %>%
  summarise(n = n(), fitted = sum(fitted), median_sigma_ppb = median(CH4_sigma_ppb), median_closure_s = median(closure_s),
            median_MDF = median(CH4_MDF), pct_below_MDF = 100 * mean(CH4_below_MDF),
            pct_uptake_detected = 100 * mean(CH4_det_class == "uptake"), .groups = "drop") %>%
  mutate(across(where(is.double), ~ signif(.x, 4)))
print(as.data.frame(det))
write.csv(det, file.path(OUT, "detection_by_analyzer.csv"), row.names = FALSE)

fit <- df %>% filter(fitted)
summ <- data.frame(
  statistic = c("measurements", "trees", "sampling dates", "fitted to the raw record", "earlier linear flux (no raw record)",
                "non-linear (HM) model selected, % of fitted", "median CH4 R2 (linear fit)",
                "below MDF, upland %", "below MDF, wetland %", "positive flux, upland %", "positive flux, wetland %",
                "goFlux QC screen fired", "trace QC window problem", "trace QC data problem", "window set by hand",
                "window ended at chamber opening"),
  value = c(nrow(df), n_distinct(df$Tree), n_distinct(df$date), sum(df$fitted), sum(!df$fitted),
            round(100 * mean(fit$CH4_model == "HM"), 1), round(median(fit$CH4_r2, na.rm = TRUE), 3),
            round(100 * mean(df$CH4_below_MDF[df$location == "Upland"]), 1), round(100 * mean(df$CH4_below_MDF[df$location == "Wetland"]), 1),
            round(100 * mean(df$CH4_flux_nmolpm2ps[df$location == "Upland"] > 0), 1), round(100 * mean(df$CH4_flux_nmolpm2ps[df$location == "Wetland"] > 0), 1),
            sum(nchar(coalesce(df$qc_note, "")) > 0), sum(df$trace_window_problem %in% TRUE), sum(df$trace_data_problem %in% TRUE),
            sum(df$manual_window %in% TRUE), sum(df$end_trim_s > 0, na.rm = TRUE)))
print(summ, row.names = FALSE)
write.csv(summ, file.path(OUT, "flux_data_summary.csv"), row.names = FALSE)

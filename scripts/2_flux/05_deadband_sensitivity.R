# ============================================================
# 05_deadband_sensitivity.R  (flux step 5 of 6; sensitivity check, not used downstream)
#
# How sensitive are the fluxes to the deadband (seconds discarded after
# chamber closure)? Refits every closure with goFlux (LM/HM, best.flux) for
# deadbands of 0, 10, 20, 30 and 45 s measured from chamber closure:
#   closure = field-log start (UGGA) or start of the tagged remark (LI-7810);
#   the adopted windows start 20 s after closure for both analyzers
# Uses the traces saved by 02_fit_fluxes.R (window +/- 120 s context).
#
# Output: outputs/tables/flux_processing/deadband_sensitivity.csv (per closure x deadband)
#         outputs/tables/flux_processing/deadband_sensitivity_summary.csv
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(goFlux) })
SCRIPT <- "05_deadband_sensitivity"
source("scripts/2_flux/flux_settings.R")
tr <- readRDS(PATH_TRACES)
traces <- tr$traces
g <- read.csv(PATH_FITS, stringsAsFactors = FALSE) %>%
  filter(fitted) %>% select(closure_id, instrument, location, SPECIES, in_legacy_dataset, date)
DB <- c(0, 10, 20, 30, 45)
traces$closure <- traces$start.time - 20   # both analyzers now use a 20-s deadband after closure
traces$closure <- as.POSIXct(traces$closure, origin = "1970-01-01", tz = attr(traces$POSIX.time, "tzone"))
res <- list()
for (db in DB) {
  d <- traces
  d$flag <- as.numeric(d$POSIX.time >= d$closure + db & d$POSIX.time <= d$end.time)
  d$Etime <- as.numeric(d$POSIX.time - (d$closure + db), units = "secs")
  keep <- tapply(d$flag, d$UniqueID, sum); keep <- names(keep)[keep >= 60]
  d <- d[d$UniqueID %in% keep, ]
  f <- suppressWarnings(suppressMessages(goFlux(d, "CH4dry_ppb", H2O_col = "H2O_ppm")))
  b <- suppressWarnings(best.flux(f))
  res[[as.character(db)]] <- data.frame(closure_id = b$UniqueID, deadband_s = db, flux = b$best.flux,
                                        LM = b$LM.flux, HM = b$HM.flux, model = b$model, nb = b$nb.obs)
  message("deadband ", db, " s: ", nrow(b), " closures")
}
res <- bind_rows(res) %>% left_join(g, by = "closure_id")
write.csv(res, file.path(TAB_DIR, "deadband_sensitivity.csv"), row.names = FALSE)
base <- res %>% filter(deadband_s == 0) %>% select(closure_id, flux0 = flux, LM0 = LM)
s <- res %>% left_join(base, by = "closure_id") %>% filter(in_legacy_dataset) %>%
  mutate(inst = ifelse(instrument == "LI-7810", "LI-7810", "LGR/UGGA")) %>%
  group_by(inst, location, deadband_s) %>%
  summarise(n = n(), mean_flux = mean(flux), median_flux = median(flux),
            median_ratio_to_db0 = median(flux / flux0, na.rm = TRUE),
            pct_change_gt25 = 100 * mean(abs(flux - flux0) > 0.25 * abs(flux0) & abs(flux0) > 0.01),
            sign_changes = sum(sign(flux) != sign(flux0)), pct_HM = 100 * mean(model == "HM"),
            .groups = "drop")
write.csv(s, file.path(TAB_DIR, "deadband_sensitivity_summary.csv"), row.names = FALSE)
options(width = 180); print(as.data.frame(s %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
bg <- res %>% filter(SPECIES == "bg", in_legacy_dataset) %>% group_by(deadband_s) %>%
  summarise(nyssa_mean = mean(flux), n = n())
print(as.data.frame(bg))

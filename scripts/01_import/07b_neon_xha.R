# ============================================================
# 07b_neon_xha.R
# Rebuild the US-xHA (NEON HARV tower) driver variables from NEON's own data
# products, so that they cover the whole flux record (the AmeriFlux US-xHA BASE
# release ends in Dec 2024). Each series is aggregated the same way as the
# AmeriFlux-based version in 08_align.R and checked against it over 2023-2024.
#
#   FC_xHA, SC_xHA      CO2 turbulent and storage flux        DP4.00200.001 dp04
#   USTAR_xHA           friction velocity                     DP4.00200.001 dp04
#   CO2_MR_xHA          CO2 dry mole fraction, profile mean   DP4.00200.001 dp01 co2Stor
#   CH4_MR_xHA          CH4 dry mole fraction, profile mean   DP4.00200.001 dp01 ch4Conc
#   T_CANOPY_xHA        IR canopy temperature, tower mean     DP1.00005.001
#   THROUGHFALL_xHA     throughfall, mean of gauges (mm per 30 min, hourly mean, as in 08_align.R)  DP1.00046.001
#   G_xHA               soil heat flux, mean of plates        DP1.00040.001
#   WS_xHA, WD_xHA      2D wind, tower top                    DP1.00001.001
#
# NEON times are UTC; output is on the EST clock (UTC-5) used by every other
# series, as hour-beginning timestamps stored with a "UTC" label.
# Outputs: data/processed/neon_xha_hourly.csv
#          outputs/tables/neon_xha_vs_ameriflux.csv (overlap agreement)
#          data/processed/neon_xha_release_by_month.csv
# ============================================================
suppressPackageStartupMessages({ library(neonUtilities); library(dplyr); library(tidyr); library(lubridate) })
NEON <- "data/raw/NEON_2026"
stack_dp <- function(dp) {
  f <- list.files(file.path(NEON, dp), pattern = "^filesToStack", full.names = TRUE)
  if (!length(f)) { message("  missing ", dp); return(NULL) }
  # stack a temporary copy: stackByTable can remove the files it unpacks
  tmp <- file.path(tempdir(), paste0("stack_", dp)); unlink(tmp, recursive = TRUE); dir.create(tmp)
  file.copy(f[1], tmp, recursive = TRUE)
  suppressMessages(stackByTable(file.path(tmp, basename(f[1])), savepath = "envt"))
}
to_est_hour <- function(t) floor_date(force_tz(with_tz(as.POSIXct(t, tz = "UTC"), "EST"), "UTC"), "hour")
rel_log <- list()
log_release <- function(dp, df) if ("release" %in% names(df) && "startDateTime" %in% names(df))
  rel_log[[dp]] <<- df %>% transmute(product = dp, month = substr(as.character(startDateTime), 1, 7), release) %>% distinct()
pick <- function(lst, pat) { nm <- names(lst)[grepl(pat, names(lst))]; if (!length(nm)) NULL else lst[[nm[1]]] }
hourly <- function(df, value, fun = mean, positions = NULL) {
  if (!is.null(positions)) df <- df %>% filter(paste(horizontalPosition, verticalPosition, sep = ".") %in% positions)
  df %>% filter(!is.na(.data[[value]])) %>% mutate(datetime = to_est_hour(startDateTime)) %>%
    group_by(datetime, horizontalPosition, verticalPosition) %>% summarise(v = fun(.data[[value]]), .groups = "drop") %>%
    group_by(datetime) %>% summarise(v = mean(v), .groups = "drop")
}
out <- list()

# ---- DP1 products ----
irbt <- stack_dp("DP1.00005.001")
if (!is.null(irbt)) { t <- pick(irbt, "IRBT_30"); log_release("DP1.00005.001", t)
  tow <- t %>% filter(horizontalPosition == "000")                    # tower levels only
  out$T_CANOPY_xHA <- hourly(tow, "bioTempMean") %>% rename(T_CANOPY_xHA = v) }
thr <- stack_dp("DP1.00046.001")
if (!is.null(thr)) { t <- pick(thr, "THRPRE_30|30min|30_min"); log_release("DP1.00046.001", t)
  vcol <- grep("^(TF)?precipBulk$", names(t), value = TRUE, ignore.case = TRUE)[1]
  out$THROUGHFALL_xHA <- hourly(t, vcol) %>% rename(THROUGHFALL_xHA = v) }
shf <- stack_dp("DP1.00040.001")
if (!is.null(shf)) { t <- pick(shf, "SHF_30"); log_release("DP1.00040.001", t)
  out$G_xHA <- hourly(t, "SHFMean") %>% rename(G_xHA = v) }
wnd <- stack_dp("DP1.00001.001")
if (!is.null(wnd)) { t <- pick(wnd, "2DWSD_30|twoDWSD_30"); log_release("DP1.00001.001", t)
  top <- max(as.integer(t$verticalPosition[t$horizontalPosition == "000"]))
  tt <- t %>% filter(horizontalPosition == "000", as.integer(verticalPosition) == top)
  out$WS_xHA <- hourly(tt, "windSpeedMean") %>% rename(WS_xHA = v)
  out$WD_xHA <- tt %>% filter(!is.na(windDirMean)) %>% mutate(datetime = to_est_hour(startDateTime)) %>%
    group_by(datetime) %>% summarise(WD_xHA = (atan2(mean(sin(windDirMean * pi / 180)), mean(cos(windDirMean * pi / 180))) * 180 / pi) %% 360, .groups = "drop") }

# ---- eddy-covariance bundle ----
ec_dir <- list.files(file.path(NEON, "DP4.00200.001"), pattern = "^filesToStack", full.names = TRUE)
if (!length(ec_dir)) {
  # fall back to the bundle already in the repo (data/raw/NEON_eddy-flux: one folder per month,
  # RELEASE-2025 through 2024-06 and PROVISIONAL after); link its H5 files into one folder
  old <- list.files("data/raw/NEON_eddy-flux", pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
  old <- old[grepl("\\.(2023|2024|2025)-\\d{2}\\.basic", old)]
  if (length(old)) {
    ec_dir <- file.path(tempdir(), "ec_h5"); dir.create(ec_dir, showWarnings = FALSE)
    for (f in old) file.symlink(normalizePath(f), file.path(ec_dir, basename(f)))
    message("  eddy-covariance bundle: using ", length(old), " monthly files from data/raw/NEON_eddy-flux")
    attr(ec_dir, "src") <- dirname(old)
  }
}
if (length(ec_dir)) {
  # stackEddy reads 36 monthly HDF5 files; cache the extracted tables so reruns are fast
  CACHE <- "data/processed/neon_ec_extract.rds"
  if (file.exists(CACHE)) { ce <- readRDS(CACHE); f4 <- ce$f4; g1 <- ce$g1 } else {
    f4 <- suppressMessages(stackEddy(ec_dir[1], level = "dp04"))$HARV
    # only the two profile concentrations used here (reading every dp01 variable takes hours)
    g1 <- suppressMessages(stackEddy(ec_dir[1], level = "dp01", avg = 30, var = c("rtioMoleDryCo2", "rtioMoleDryCh4")))$HARV
    saveRDS(list(f4 = f4, g1 = g1), CACHE)
  }
  f4$datetime <- to_est_hour(f4$timeBgn)
  # NEON final quality flags applied, as in the AmeriFlux release (unfiltered turbulent
  # CO2 flux agrees with AmeriFlux FC at r = 0.67; filtered at r = 0.999)
  out$EC <- f4 %>% mutate(fc = ifelse(qfqm.fluxCo2.turb.qfFinl == 0, data.fluxCo2.turb.flux, NA),
                          sc = ifelse(qfqm.fluxCo2.stor.qfFinl == 0, data.fluxCo2.stor.flux, NA),
                          us = ifelse(qfqm.fluxMome.turb.qfFinl == 0, data.fluxMome.turb.veloFric, NA)) %>%
    group_by(datetime) %>%
    summarise(FC_xHA = mean(fc, na.rm = TRUE), SC_xHA = mean(sc, na.rm = TRUE), USTAR_xHA = mean(us, na.rm = TRUE), .groups = "drop") %>%
    mutate(across(-datetime, ~ ifelse(is.nan(.x), NA, .x)))
  g1$datetime <- to_est_hour(g1$timeBgn)
  # tower profile heights only (the table also holds calibration-gas rows: co2Arch, co2High, ...)
  g1 <- g1 %>% filter(verticalPosition %in% sprintf("%03d", seq(10, 60, 10)))
  co2c <- grep("^data\\.co2Stor\\.rtioMoleDryCo2\\.mean$", names(g1), value = TRUE)
  ch4c <- grep("^data\\.ch4Conc\\.rtioMoleDryCh4\\.mean$", names(g1), value = TRUE)
  # profile mean over the six heights, required to have at least five heights in the hour
  # (AmeriFlux numbers profile levels from the top; NEON 010 = AmeriFlux level 6, r = 0.998)
  g1 <- g1 %>% group_by(datetime) %>%
    filter(if (length(co2c)) sum(!is.na(.data[[co2c]])) >= 5 else TRUE) %>% ungroup()
  gas <- g1 %>% group_by(datetime) %>%
    summarise(CO2_MR_xHA = if (length(co2c)) mean(.data[[co2c]], na.rm = TRUE) else NA_real_,
              CH4_MR_xHA = if (length(ch4c)) mean(.data[[ch4c]], na.rm = TRUE) * (if (length(ch4c) && median(g1[[ch4c]], na.rm = TRUE) < 10) 1000 else 1) else NA_real_,
              .groups = "drop") %>% mutate(across(-datetime, ~ ifelse(is.nan(.x), NA, .x)))
  out$GAS <- gas
  rel_log[["DP4.00200.001"]] <- data.frame(product = "DP4.00200.001",
    month = sub(".*\\.(\\d{4}-\\d{2})\\..*", "\\1", if (is.null(attr(ec_dir, "src"))) list.files(ec_dir[1], pattern = "\\.h5$") else basename(attr(ec_dir, "src"))),
    release = { src <- if (is.null(attr(ec_dir, "src"))) list.files(ec_dir[1], pattern = "\\.h5$") else basename(attr(ec_dir, "src"))
                ifelse(grepl("PROVISIONAL", src), "PROVISIONAL", sub(".*\\.(RELEASE-\\d{4}).*", "\\1", src)) }) %>% distinct()
}

xha <- Reduce(function(a, b) full_join(a, b, by = "datetime"), out) %>% arrange(datetime)
write.csv(xha %>% mutate(datetime = format(datetime, "%Y-%m-%dT%H:%M:%SZ")), "data/processed/neon_xha_hourly.csv", row.names = FALSE)
write.csv(bind_rows(rel_log), "data/processed/neon_xha_release_by_month.csv", row.names = FALSE)
message("NEON xHA variables: ", paste(setdiff(names(xha), "datetime"), collapse = ", "),
        "; ", format(min(xha$datetime)), " to ", format(max(xha$datetime)))

# ---- agreement with the AmeriFlux US-xHA series over their overlap ----
source("scripts/helpers/find_ameriflux.R")
amf <- read.csv(find_ameriflux("US-xHA"), skip = 2, na.strings = "-9999")
amf$datetime <- floor_date(ymd_hm(as.character(amf$TIMESTAMP_START), tz = "UTC"), "hour")
rm_ <- function(cols, fun = mean) { x <- amf[, intersect(cols, names(amf)), drop = FALSE]; if (!ncol(x)) return(rep(NA, nrow(amf))); apply(x, 1, function(r) if (all(is.na(r))) NA else fun(r, na.rm = TRUE)) }
ref <- data.frame(datetime = amf$datetime,
  FC_xHA = amf$FC, SC_xHA = amf$SC, USTAR_xHA = amf$USTAR,
  CO2_MR_xHA = rm_(grep("^CO2_MIXING_RATIO_1_[1-6]_[12]$", names(amf), value = TRUE)),
  CH4_MR_xHA = rm_(grep("^CH4_MIXING_RATIO_1_[1-6]_1$", names(amf), value = TRUE)),
  T_CANOPY_xHA = rm_(c("T_CANOPY_1_1_1", "T_CANOPY_1_2_1", "T_CANOPY_1_3_1", "T_CANOPY_1_4_1", "T_CANOPY_2_4_1")),
  THROUGHFALL_xHA = rm_(paste0("THROUGHFALL_", 1:5, "_1_1")),
  G_xHA = rm_(c("G_1_1_1", "G_3_1_1", "G_5_1_1")), WS_xHA = amf$WS_1_1_1) %>%
  group_by(datetime) %>% summarise(across(everything(), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
cmp <- inner_join(ref, xha, by = "datetime", suffix = c("_amf", "_neon"))
agree <- bind_rows(lapply(intersect(setdiff(names(ref), "datetime"), names(xha)), function(v) {
  a <- cmp[[paste0(v, "_amf")]]; b <- cmp[[paste0(v, "_neon")]]; ok <- is.finite(a) & is.finite(b)
  data.frame(variable = v, n_hours = sum(ok), r = if (sum(ok) > 10) cor(a[ok], b[ok]) else NA,
             mean_amf = mean(a[ok]), mean_neon = mean(b[ok]), slope = if (sum(ok) > 10) coef(lm(b[ok] ~ a[ok]))[2] else NA)
}))
write.csv(agree, "outputs/tables/neon_xha_vs_ameriflux.csv", row.names = FALSE)
print(agree, digits = 3)

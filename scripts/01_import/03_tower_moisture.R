Ha1<-read.csv(
  'data/raw/ameriflux/AMF_US-Ha1_BASE-BADM_26-5/AMF_US-Ha1_BASE_HR_26-5.csv',
  header=T, skip=2)
head(Ha1)

Ha2<-read.csv(
  'data/raw/ameriflux/AMF_US-Ha2_BASE-BADM_15-5/AMF_US-Ha2_BASE_HH_15-5.csv',
  header=T, skip=2)
head(Ha2)

xHA<-read.csv(
  'data/raw/ameriflux/AMF_US-xHA_BASE-BADM_11-5/AMF_US-xHA_BASE_HH_11-5.csv',
  header=T, skip=2)
head(xHA)


library(dplyr)
library(lubridate)
library(ggplot2)
library(stringr)

to_posix_af <- function(x) {
  x_chr <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  x_chr <- str_pad(x_chr, width = 12, side = "left", pad = "0")
  ymd_hm(x_chr, tz = "UTC")
}

ha1_pair <- Ha1 %>%
  mutate(
    datetime = to_posix_af(TIMESTAMP_START),
    FC_1_1_1 = na_if(FC_1_1_1, -9999),
    FC_1_1_2 = na_if(FC_1_1_2, -9999)
  ) %>%
  filter(
    year(datetime) >= 2022,
    year(datetime) <= 2026,
    !is.na(FC_1_1_1),
    !is.na(FC_1_1_2)
  )

ha1_monthly <- ha1_pair %>%
  mutate(month = floor_date(datetime, "month")) %>%
  group_by(month) %>%
  summarize(
    FC_1 = mean(FC_1_1_1, na.rm = TRUE),
    FC_2 = mean(FC_1_1_2, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(ha1_monthly) +
  geom_line(aes(month, FC_1, color = "FC_1_1_1")) +
  geom_line(aes(month, FC_2, color = "FC_1_1_2")) +
  labs(
    title = "Ha1 replicate EC systems – monthly mean CO₂ flux",
    x = NULL,
    y = expression("FC ("*mu*"mol m"^-2*" s"^-1*")"),
    color = "Series"
  ) +
  theme(legend.position = "bottom")


ggplot(ha1_pair, aes(FC_1_1_1, FC_1_1_2)) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() +
  labs(
    title = "Ha1 EC replicate agreement (paired half-hourly data)",
    x = expression("FC_1_1_1 ("*mu*"mol m"^-2*" s"^-1*")"),
    y = expression("FC_1_1_2 ("*mu*"mol m"^-2*" s"^-1*")")
  )




library(dplyr)
library(lubridate)
library(ggplot2)
library(stringr)

to_posix_af <- function(x) {
  x_chr <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  x_chr <- str_pad(x_chr, width = 12, side = "left", pad = "0")
  ymd_hm(x_chr, tz = "UTC")
}

ha2_pair <- Ha2 %>%
  mutate(
    datetime  = to_posix_af(TIMESTAMP_START),
    FC_2_1_1  = na_if(FC_2_1_1, -9999),
    FC_2_1_2  = na_if(FC_2_1_2, -9999)
  ) %>%
  filter(
    year(datetime) >= 2022,
    year(datetime) <= 2026,
    !is.na(FC_2_1_1),
    !is.na(FC_2_1_2)
  )


ha2_monthly <- ha2_pair %>%
  mutate(month = floor_date(datetime, "month")) %>%
  group_by(month) %>%
  summarize(
    FC_1 = mean(FC_2_1_1, na.rm = TRUE),
    FC_2 = mean(FC_2_1_2, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(ha2_monthly) +
  geom_line(aes(month, FC_1, color = "FC_2_1_1")) +
  geom_line(aes(month, FC_2, color = "FC_2_1_2")) +
  labs(
    title = "Ha2 replicate EC systems – monthly mean CO₂ flux",
    x = NULL,
    y = expression("FC ("*mu*"mol m"^-2*" s"^-1*")"),
    color = "Series"
  ) +
  theme(legend.position = "bottom")



ha2_monthly <- ha2_pair %>%
  mutate(month = floor_date(datetime, "month")) %>%
  group_by(month) %>%
  summarize(
    FC_1 = mean(FC_2_1_1, na.rm = TRUE),
    FC_2 = mean(FC_2_1_2, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(ha2_monthly) +
  geom_line(aes(month, FC_1, color = "FC_2_1_1")) +
  geom_line(aes(month, FC_2, color = "FC_2_1_2")) +
  labs(
    title = "Ha2 replicate EC systems – monthly mean CO₂ flux",
    x = NULL,
    y = expression("FC ("*mu*"mol m"^-2*" s"^-1*")"),
    color = "Series"
  ) +
  theme(legend.position = "bottom")


ggplot(ha2_pair, aes(FC_2_1_1, FC_2_1_2)) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal() +
  labs(
    title = "Ha2 EC replicate agreement (paired half-hourly data)",
    x = expression("FC_2_1_1 ("*mu*"mol m"^-2*" s"^-1*")"),
    y = expression("FC_2_1_2 ("*mu*"mol m"^-2*" s"^-1*")")
  )












library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(stringr)

to_posix_af <- function(x) {
  x_chr <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  x_chr <- str_pad(x_chr, width = 12, side = "left", pad = "0")
  ymd_hm(x_chr, tz = "UTC")  # plotting only; timezone choice won't change monthly means much
}

monthly_all_fc <- function(df, tower_name, start_year = 2022, end_year = 2026) {
  fc_cols <- names(df)[grepl("^FC($|_)", names(df))]
  if (length(fc_cols) == 0) {
    return(tibble(tower = tower_name, month = as.POSIXct(character()), sensor = character(), FC = numeric()))
  }
  
  df %>%
    mutate(datetime = to_posix_af(TIMESTAMP_START)) %>%
    filter(year(datetime) >= start_year, year(datetime) <= end_year) %>%
    select(datetime, all_of(fc_cols)) %>%
    mutate(across(all_of(fc_cols), ~ na_if(.x, -9999))) %>%
    pivot_longer(
      cols = all_of(fc_cols),
      names_to = "sensor",
      values_to = "FC"
    ) %>%
    mutate(
      month = floor_date(datetime, "month"),
      tower = tower_name
    ) %>%
    group_by(tower, sensor, month) %>%
    summarize(
      FC = mean(FC, na.rm = TRUE),
      n = sum(!is.na(FC)),
      .groups = "drop"
    ) %>%
    filter(is.finite(FC))
}

# Build monthly means for every FC sensor in each tower
all_monthly_fc <- bind_rows(
  monthly_all_fc(Ha1, "Ha1"),
  monthly_all_fc(Ha2, "Ha2"),
  monthly_all_fc(xHA, "xHA")
)

# Plot 1: Facet by tower, colored by sensor (good when there are only a few sensors)
ggplot(all_monthly_fc, aes(month, FC, color = sensor)) +
  geom_line() +
  facet_wrap(~ tower, ncol = 1, scales = "free_y") +
  labs(
    title = "Monthly mean CO2 turbulent flux (all FC sensors), 2022–2026",
    x = NULL,
    y = expression("FC ("*mu*"mol m"^-2*" s"^-1*")"),
    color = "Sensor"
  ) +
  theme(legend.position = "bottom")

# Plot 2: Facet by tower AND sensor (best when there are many sensors)
ggplot(all_monthly_fc, aes(month, FC)) +
  geom_line() +
  facet_grid(sensor ~ tower, scales = "free_y") +
  labs(
    title = "Monthly mean CO2 turbulent flux by sensor and tower, 2022–2026",
    x = NULL,
    y = expression("FC ("*mu*"mol m"^-2*" s"^-1*")")
  )

# ============================================================
# Baseline-adjusted half-hourly SWC mean for 3 towers (Ha1, Ha2, xHA)
# ============================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(stringr)
library(ggplot2)

# ---------------------------
# 0) DATA LOADING (edit as needed)
# ---------------------------
# If Ha1, Ha2, xHA already exist in your environment, you can skip this block.
# Otherwise, uncomment and set paths:
# Ha1 <- read.csv("path/to/Ha1.csv", stringsAsFactors = FALSE)
# Ha2 <- read.csv("path/to/Ha2.csv", stringsAsFactors = FALSE)
# xHA <- read.csv("path/to/xHA.csv", stringsAsFactors = FALSE)

# Put tower dfs in a named list for iteration
towers <- list(Ha1 = Ha1, Ha2 = Ha2, xHA = xHA)

# ---------------------------
# 1) TIMESTAMP PARSING
# ---------------------------
to_posix_af <- function(x) {
  x_chr <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  x_chr <- str_pad(x_chr, width = 12, side = "left", pad = "0")
  ymd_hm(x_chr, tz = "UTC")
}

# ---------------------------
# 2) LONG-FORM SWC EXTRACTION (handles multiple naming conventions)
# ---------------------------
# adjust include_pi = FALSE if you want to exclude.
swc_long <- function(df,
                     start_year = 2022, end_year = 2026,
                     offline_codes = c(-9999),
                     include_pi = FALSE) {
  
  if (!("TIMESTAMP_START" %in% names(df))) stop("Expected column TIMESTAMP_START is missing.")
  
  # Match SWC columns for Ha1/xHA: "SWC" or "SWC_*"
  swc_cols <- names(df)[str_detect(names(df), "^SWC($|_)")]
  
  # Match PI-style SWC columns (Ha2 has SWC_PI_1.. etc); include by default
  if (include_pi) {
    pi_cols <- names(df)[str_detect(names(df), "^SWC_PI($|_)")]
    swc_cols <- union(swc_cols, pi_cols)
  }
  
  if (length(swc_cols) == 0) stop("No SWC columns found (expected ^SWC($|_) and optionally ^SWC_PI).")
  
  df %>%
    mutate(datetime = to_posix_af(TIMESTAMP_START)) %>%
    filter(year(datetime) >= start_year, year(datetime) <= end_year) %>%
    select(datetime, all_of(swc_cols)) %>%
    mutate(across(all_of(swc_cols), ~ replace(as.numeric(.x), as.numeric(.x) %in% offline_codes, NA_real_))) %>%
    pivot_longer(all_of(swc_cols), names_to = "sensor", values_to = "value") %>%
    filter(!is.na(value)) %>%
    arrange(datetime, sensor)
}

# ---------------------------
# 3) BASELINE-ADJUSTED (STANDARDIZED) MEAN
# ---------------------------
# Estimand: average temporal dynamics after removing persistent sensor wet/dry offsets.
# Implementation: per-sensor baseline b_s computed over baseline window; average anomalies each time; add back grand b_bar.
half_hourly_baseline_adjusted_mean <- function(df,
                                               start_year = 2022, end_year = 2026,
                                               baseline_start = "2024-01-01",
                                               baseline_end   = "2024-12-31",
                                               baseline_fun = median,   # robust
                                               min_sensors = 8,
                                               offline_codes = c(-9999),
                                               include_pi = TRUE) {
  
  long <- swc_long(df, start_year, end_year, offline_codes, include_pi)
  
  base_tbl <- long %>%
    filter(as.Date(datetime) >= as.Date(baseline_start),
           as.Date(datetime) <= as.Date(baseline_end)) %>%
    group_by(sensor) %>%
    summarize(
      b_s = baseline_fun(value, na.rm = TRUE),
      n_base = dplyr::n(),
      .groups = "drop"
    ) %>%
    filter(!is.na(b_s))
  
  if (nrow(base_tbl) == 0) stop("No sensor baselines computed. Check baseline window or data availability.")
  
  b_bar <- mean(base_tbl$b_s, na.rm = TRUE)
  
  ts <- long %>%
    inner_join(base_tbl %>% select(sensor, b_s), by = "sensor") %>%
    mutate(anom = value - b_s) %>%
    group_by(datetime) %>%
    summarize(
      mean_anom = mean(anom, na.rm = TRUE),
      n_sensors = n_distinct(sensor),
      mean_std  = mean_anom + b_bar,
      .groups = "drop"
    ) %>%
    filter(n_sensors >= min_sensors) %>%
    arrange(datetime)
  
  list(
    half_hourly = ts,
    baselines = base_tbl,
    grand_baseline = b_bar,
    settings = list(
      start_year = start_year, end_year = end_year,
      baseline_start = baseline_start, baseline_end = baseline_end,
      baseline_fun = deparse(substitute(baseline_fun)),
      min_sensors = min_sensors,
      offline_codes = offline_codes,
      include_pi = include_pi
    )
  )
}

# ---------------------------
# 4) DIAGNOSTICS
# ---------------------------
plot_mean <- function(ts_df, tower_name) {
  ggplot(ts_df, aes(datetime, mean_std)) +
    geom_line(linewidth = 0.35) +
    labs(title = paste0(tower_name, ": half-hourly baseline-adjusted mean SWC"),
         x = NULL, y = "SWC (units)") +
    theme_bw()
}

plot_coverage <- function(ts_df, tower_name) {
  ggplot(ts_df, aes(datetime, n_sensors)) +
    geom_line(linewidth = 0.35) +
    labs(title = paste0(tower_name, ": sensor coverage contributing to mean"),
         x = NULL, y = "n sensors") +
    theme_bw()
}

plot_baseline_distribution <- function(base_tbl, tower_name) {
  ggplot(base_tbl, aes(b_s)) +
    geom_histogram(bins = 40) +
    labs(title = paste0(tower_name, ": per-sensor baseline distribution"),
         x = "baseline SWC (b_s)", y = "count") +
    theme_bw()
}

plot_baseline_n <- function(base_tbl, tower_name) {
  ggplot(base_tbl, aes(n_base)) +
    geom_histogram(bins = 40) +
    labs(title = paste0(tower_name, ": baseline window sample size per sensor"),
         x = "n observations in baseline window", y = "count") +
    theme_bw()
}

plot_mean_by_coverage_strata <- function(ts_df, tower_name) {
  
  # candidate breaks from quantiles (intended quintiles)
  brks <- as.numeric(quantile(ts_df$n_sensors, probs = seq(0, 1, 0.2), na.rm = TRUE))
  brks <- unique(brks)
  
  # If too few unique breaks, fall back to equally spaced breaks
  # across the observed range (up to 4 bins), or skip gracefully.
  if (length(brks) < 3) {
    rng <- range(ts_df$n_sensors, na.rm = TRUE)
    
    # If essentially constant coverage, this diagnostic isn't meaningful
    if (rng[1] == rng[2]) {
      message(tower_name, ": n_sensors is constant (", rng[1], "). Skipping coverage-strata diagnostic.")
      return(
        ggplot() +
          annotate("text", x = 0, y = 0,
                   label = paste0(tower_name, ": coverage constant; strata diagnostic skipped")) +
          theme_void()
      )
    }
    
    # 3–4 bins based on range
    brks <- pretty(rng, n = 4)
    brks <- unique(brks)
  }
  
  # Ensure strictly increasing breaks
  brks <- sort(brks)
  # Safety: if still not enough, skip
  if (length(brks) < 3) {
    message(tower_name, ": insufficient unique breaks for coverage strata. Skipping.")
    return(
      ggplot() +
        annotate("text", x = 0, y = 0,
                 label = paste0(tower_name, ": insufficient coverage variation; strata diagnostic skipped")) +
        theme_void()
    )
  }
  
  ts_df %>%
    mutate(
      cov_bin = cut(n_sensors, breaks = brks, include.lowest = TRUE, right = TRUE)
    ) %>%
    ggplot(aes(cov_bin, mean_std)) +
    geom_boxplot(outlier.alpha = 0.2) +
    labs(
      title = paste0(tower_name, ": mean_std by coverage strata (informativeness proxy)"),
      x = "coverage bin (n sensors)",
      y = "baseline-adjusted mean SWC"
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

# Sensitivity to baseline window choice (per tower)
baseline_sensitivity <- function(df,
                                 windows = list(
                                   c("2023-01-01", "2023-12-31"),
                                   c("2024-01-01", "2024-12-31"),
                                   c("2024-07-01", "2025-06-30")
                                 ),
                                 start_year = 2022, end_year = 2026,
                                 baseline_fun = median,
                                 min_sensors = 8,
                                 offline_codes = c(-9999),
                                 include_pi = TRUE) {
  
  runs <- lapply(seq_along(windows), function(i) {
    w <- windows[[i]]
    res <- half_hourly_baseline_adjusted_mean(
      df,
      start_year = start_year, end_year = end_year,
      baseline_start = w[1], baseline_end = w[2],
      baseline_fun = baseline_fun,
      min_sensors = min_sensors,
      offline_codes = offline_codes,
      include_pi = include_pi
    )
    res$half_hourly %>%
      select(datetime, mean_std, n_sensors) %>%
      mutate(window = paste0(w[1], " to ", w[2]))
  })
  
  bind_rows(runs)
}

plot_sensitivity <- function(sens_df, tower_name) {
  ggplot(sens_df, aes(datetime, mean_std, color = window)) +
    geom_line(linewidth = 0.35) +
    labs(title = paste0(tower_name, ": sensitivity to baseline window choice"),
         x = NULL, y = "baseline-adjusted mean SWC", color = "baseline window") +
    theme_bw() +
    theme(legend.position = "bottom")
}

# ---------------------------
# 5) RUN FOR ALL 3 TOWERS
# ---------------------------

# Global settings (edit)
start_year     <- 2022
end_year       <- 2026
baseline_start <- "2024-01-01"
baseline_end   <- "2024-12-31"
baseline_fun   <- median
min_sensors    <- 8
offline_codes  <- c(-9999)
include_pi     <- TRUE   # include SWC_PI_* where present (Ha2)

results <- lapply(names(towers), function(nm) {
  res <- half_hourly_baseline_adjusted_mean(
    towers[[nm]],
    start_year = start_year, end_year = end_year,
    baseline_start = baseline_start, baseline_end = baseline_end,
    baseline_fun = baseline_fun,
    min_sensors = min_sensors,
    offline_codes = offline_codes,
    include_pi = include_pi
  )
  res$tower <- nm
  res
})
names(results) <- names(towers)

# Combined half-hourly series for all towers (useful for faceting or comparison)
all_ts <- bind_rows(lapply(names(results), function(nm) {
  results[[nm]]$half_hourly %>%
    mutate(tower = nm)
}))

# ---------------------------
# 6) OUTPUT PLOTS (per tower) + combined
# ---------------------------

# Per tower: main series + diagnostics
for (nm in names(results)) {
  ts_df   <- results[[nm]]$half_hourly
  base_df <- results[[nm]]$baselines
  
  print(plot_mean(ts_df, nm))
  print(plot_coverage(ts_df, nm))
  print(plot_baseline_distribution(base_df, nm))
  print(plot_baseline_n(base_df, nm))
  print(plot_mean_by_coverage_strata(ts_df, nm))
  
  sens_df <- baseline_sensitivity(
    towers[[nm]],
    windows = list(
      c("2023-01-01", "2023-12-31"),
      c("2024-01-01", "2024-12-31"),
      c("2024-07-01", "2025-06-30")
    ),
    start_year = start_year, end_year = end_year,
    baseline_fun = baseline_fun,
    min_sensors = min_sensors,
    offline_codes = offline_codes,
    include_pi = include_pi
  )
  print(plot_sensitivity(sens_df, nm))
}

# Combined: compare towers on same figure
ggplot(all_ts, aes(datetime, mean_std, color = tower)) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~ tower, scales = "free_y", ncol = 1) +
  labs(title = "Baseline-adjusted SWC mean by tower (half-hourly)",
       x = NULL, y = "baseline-adjusted mean SWC", color = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

# Combined: coverage comparison
ggplot(all_ts, aes(datetime, n_sensors, color = tower)) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~ tower, scales = "free_y", ncol = 1) +
  labs(title = "Sensor coverage contributing to mean by tower",
       x = NULL, y = "n sensors", color = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

# ---------------------------
# 7) NUMERICAL SUMMARIES (for reporting)
# ---------------------------

coverage_summary <- all_ts %>%
  group_by(tower) %>%
  summarize(
    n_time = n(),
    n_sensors_min = min(n_sensors, na.rm = TRUE),
    n_sensors_p05 = quantile(n_sensors, 0.05, na.rm = TRUE),
    n_sensors_med = median(n_sensors, na.rm = TRUE),
    n_sensors_p95 = quantile(n_sensors, 0.95, na.rm = TRUE),
    n_sensors_max = max(n_sensors, na.rm = TRUE),
    .groups = "drop"
  )

baseline_summary <- bind_rows(lapply(names(results), function(nm) {
  results[[nm]]$baselines %>%
    mutate(tower = nm)
})) %>%
  group_by(tower) %>%
  summarize(
    n_sensors = n(),
    b_p05 = quantile(b_s, 0.05, na.rm = TRUE),
    b_med = median(b_s, na.rm = TRUE),
    b_p95 = quantile(b_s, 0.95, na.rm = TRUE),
    n_base_p05 = quantile(n_base, 0.05, na.rm = TRUE),
    n_base_med = median(n_base, na.rm = TRUE),
    n_base_p95 = quantile(n_base, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

coverage_summary
baseline_summary

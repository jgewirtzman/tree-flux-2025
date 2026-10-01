# ============================================================
# Baseline-adjusted half-hourly TS mean for 3 towers (Ha1, Ha2, xHA)
# + reviewer-critical diagnostics and sensitivity checks
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
# 2) LONG-FORM TS EXTRACTION (handles multiple naming conventions)
# ---------------------------
# Includes both "TS_*" and "TS_PI_*" where present (set include_pi = FALSE to exclude PI channels).
ts_long <- function(df,
                    start_year = 2022, end_year = 2026,
                    offline_codes = c(-9999),
                    include_pi = TRUE) {
  
  if (!("TIMESTAMP_START" %in% names(df))) stop("Expected column TIMESTAMP_START is missing.")
  
  # Match TS columns: "TS" or "TS_*"
  ts_cols <- names(df)[str_detect(names(df), "^TS($|_)")]
  
  # Match PI-style TS columns (e.g., TS_PI_*, TS_PI_1 etc); include by default
  if (include_pi) {
    pi_cols <- names(df)[str_detect(names(df), "^TS_PI($|_)")]
    ts_cols <- union(ts_cols, pi_cols)
  }
  
  if (length(ts_cols) == 0) stop("No TS columns found (expected ^TS($|_) and optionally ^TS_PI).")
  
  df %>%
    mutate(datetime = to_posix_af(TIMESTAMP_START)) %>%
    filter(year(datetime) >= start_year, year(datetime) <= end_year) %>%
    select(datetime, all_of(ts_cols)) %>%
    mutate(across(all_of(ts_cols), ~ replace(as.numeric(.x), as.numeric(.x) %in% offline_codes, NA_real_))) %>%
    pivot_longer(all_of(ts_cols), names_to = "sensor", values_to = "value") %>%
    filter(!is.na(value)) %>%
    arrange(datetime, sensor)
}

# ---------------------------
# 3) BASELINE-ADJUSTED (STANDARDIZED) MEAN
# ---------------------------
half_hourly_baseline_adjusted_ts <- function(df,
                                             start_year = 2022, end_year = 2026,
                                             baseline_start = "2024-01-01",
                                             baseline_end   = "2024-12-31",
                                             baseline_fun = median,
                                             min_sensors = 8,
                                             offline_codes = c(-9999),
                                             include_pi = TRUE) {
  
  long <- ts_long(df, start_year, end_year, offline_codes, include_pi)
  
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
# 4) REVIEWER-CRITICAL DIAGNOSTICS (robust coverage strata)
# ---------------------------
plot_mean <- function(ts_df, tower_name) {
  ggplot(ts_df, aes(datetime, mean_std)) +
    geom_line(linewidth = 0.35) +
    labs(title = paste0(tower_name, ": half-hourly baseline-adjusted mean TS"),
         x = NULL, y = "Soil temperature (deg C)") +
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
    labs(title = paste0(tower_name, ": per-sensor baseline distribution (TS)"),
         x = "baseline TS (b_s, deg C)", y = "count") +
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
  
  brks <- as.numeric(quantile(ts_df$n_sensors, probs = seq(0, 1, 0.2), na.rm = TRUE))
  brks <- unique(brks)
  
  if (length(brks) < 3) {
    rng <- range(ts_df$n_sensors, na.rm = TRUE)
    if (rng[1] == rng[2]) {
      message(tower_name, ": n_sensors is constant (", rng[1], "). Skipping coverage-strata diagnostic.")
      return(
        ggplot() +
          annotate("text", x = 0, y = 0,
                   label = paste0(tower_name, ": coverage constant; strata diagnostic skipped")) +
          theme_void()
      )
    }
    brks <- unique(pretty(rng, n = 4))
  }
  
  brks <- sort(brks)
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
    mutate(cov_bin = cut(n_sensors, breaks = brks, include.lowest = TRUE, right = TRUE)) %>%
    ggplot(aes(cov_bin, mean_std)) +
    geom_boxplot(outlier.alpha = 0.2) +
    labs(title = paste0(tower_name, ": mean_std by coverage strata (informativeness proxy)"),
         x = "coverage bin (n sensors)", y = "baseline-adjusted mean TS (deg C)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

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
    res <- half_hourly_baseline_adjusted_ts(
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
         x = NULL, y = "baseline-adjusted mean TS (deg C)", color = "baseline window") +
    theme_bw() +
    theme(legend.position = "bottom")
}

# ---------------------------
# 5) RUN FOR ALL 3 TOWERS
# ---------------------------

start_year     <- 2022
end_year       <- 2026
baseline_start <- "2024-01-01"
baseline_end   <- "2024-12-31"
baseline_fun   <- median
min_sensors    <- 8
offline_codes  <- c(-9999)
include_pi     <- TRUE

results <- lapply(names(towers), function(nm) {
  res <- half_hourly_baseline_adjusted_ts(
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

all_ts <- bind_rows(lapply(names(results), function(nm) {
  results[[nm]]$half_hourly %>% mutate(tower = nm)
}))

# ---------------------------
# 6) OUTPUT PLOTS (per tower) + combined
# ---------------------------

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

# Combined: compare towers (facet)
ggplot(all_ts, aes(datetime, mean_std, color = tower)) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~ tower, scales = "free_y", ncol = 1) +
  labs(title = "Baseline-adjusted soil temperature (TS) mean by tower (half-hourly)",
       x = NULL, y = "baseline-adjusted mean TS (deg C)", color = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")

# Combined: coverage comparison
ggplot(all_ts, aes(datetime, n_sensors, color = tower)) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~ tower, scales = "free_y", ncol = 1) +
  labs(title = "Sensor coverage contributing to TS mean by tower",
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
  results[[nm]]$baselines %>% mutate(tower = nm)
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

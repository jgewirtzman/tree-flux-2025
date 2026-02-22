# ============================================================
# 08_rolling_window_correlations.R
#
# Rolling-window correlation analysis between CH4 stem flux and
# environmental predictors.
#
# Approach:
#   - Anomaly filtering (DOY + hour main effects removed)
#   - Hourly resolution windows (1h to 30 days)
#   - Variable-level inference: select best window per variable,
#     then BH-FDR across variables (not across windows)
#   - Surface row plots for visualization
#
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(RcppRoll)
  library(scales)
})

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  stem_flux = "data/raw/flux_dataset.csv",
  aligned   = "data/processed/aligned_hourly_dataset.csv"
)

OUTPUT_DIR <- "figures/rolling_correlations"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "lines_by_group"), recursive = TRUE, showWarnings = FALSE)

# Date range
DATE_MIN <- as.Date("2023-01-01")
DATE_MAX <- as.Date("2026-01-01")

# Time windows (3-hour resolution, 3h to 14 days)
MAX_DAYS <- 14
WINDOW_STEP_HOURS <- 3
WINDOWS_HOURS <- seq(WINDOW_STEP_HOURS, MAX_DAYS * 24, by = WINDOW_STEP_HOURS)
WINDOWS_DAYS <- WINDOWS_HOURS / 24

message("Configuration:")
message("  Windows: ", length(WINDOWS_DAYS), " (", WINDOW_STEP_HOURS, "h to ", MAX_DAYS, "d, step = ", WINDOW_STEP_HOURS, "h)")

# Variables to drop
DROP_VARS <- c("FC_xHA", "USTAR_xHA", "TS_xHA")

# Variable groups
VAR_GROUPS <- list(
  atmospheric   = c("tair_C", "RH", "VPD_kPa", "p_kPa", "USTAR_Ha1", "USTAR_Ha2"),
  radiation     = c("PAR", "rnet", "slrr"),
  temperature   = c("s10t", "TS_Ha1", "TS_Ha2", "T_CANOPY_xHA", "G_xHA"),
  soil_moisture = c("NEON_SWC_shallow", "NEON_SWC_mid", "NEON_SWC_deep",
                    "SWC_Ha1", "SWC_Ha2"),
  water_table   = c("bgs_wtd_cm", "bvs_wtd_cm"),
  precipitation = c("P_mm", "THROUGHFALL_xHA"),
  ecosystem_fluxes = c("FC_Ha1", "FC_Ha2", "LE_Ha1", "LE_Ha2", "H_Ha1", "H_Ha2"),
  gases         = c("CH4_MR_xHA", "CO2_MR_xHA", "CO2_MR_Ha2"),
  phenology     = c("gcc", "ndvi")
)

ALL_VARS <- unique(unlist(VAR_GROUPS, use.names = FALSE))

# Display labels
VAR_LABELS <- c(
  tair_C = "Air temp", RH = "RH", VPD_kPa = "VPD", p_kPa = "Pressure",
  USTAR_Ha1 = "u* Ha1", USTAR_Ha2 = "u* Ha2",
  PAR = "PAR", rnet = "Net rad", slrr = "Solar rad",
  s10t = "Soil temp Fisher", TS_Ha1 = "Soil temp Ha1", TS_Ha2 = "Soil temp Ha2",
  T_CANOPY_xHA = "Canopy temp", G_xHA = "Ground heat flux",
  NEON_SWC_shallow = "SWC NEON shallow", NEON_SWC_mid = "SWC NEON mid", 
  NEON_SWC_deep = "SWC NEON deep", SWC_Ha1 = "SWC Ha1", SWC_Ha2 = "SWC Ha2",
  bgs_wtd_cm = "WTD BGS", bvs_wtd_cm = "WTD BVS",
  P_mm = "Precip", THROUGHFALL_xHA = "Throughfall",
  FC_Ha1 = "FCO2 Ha1", FC_Ha2 = "FCO2 Ha2",
  LE_Ha1 = "LE Ha1", LE_Ha2 = "LE Ha2",
  H_Ha1 = "Sensible heat Ha1", H_Ha2 = "Sensible heat Ha2",
  CH4_MR_xHA = "CH4 dry", CO2_MR_xHA = "CO2 dry NEON", CO2_MR_Ha2 = "CO2 dry Ha2",
  gcc = "GCC", ndvi = "NDVI"
)

# Statistical settings
FDR_ALPHA <- 0.05
MIN_OBS   <- 5

# ============================================================
# HELPER FUNCTIONS
# ============================================================

#' Rolling mean for all variables
calc_rolling_means <- function(data, vars, window_hours) {
  data <- data %>% arrange(datetime)
  result <- data %>% dplyr::select(datetime)
  
  for (var in vars) {
    if (!var %in% names(data)) next
    result[[paste0(var, "_roll")]] <- RcppRoll::roll_mean(
      data[[var]], n = window_hours, align = "right", fill = NA, na.rm = TRUE
    )
  }
  result
}

#' Safe correlation test
safe_cor_test <- function(x, y) {
  valid <- complete.cases(x, y)
  x <- x[valid]
  y <- y[valid]
  n <- length(x)
  
  if (n < MIN_OBS) {
    return(list(r = NA_real_, p = NA_real_, n = n))
  }
  
  tryCatch({
    test <- cor.test(x, y, method = "pearson")
    list(r = as.numeric(test$estimate), p = as.numeric(test$p.value), n = n)
  }, error = function(e) {
    list(r = NA_real_, p = NA_real_, n = n)
  })
}

#' Anomaly filtering: remove DOY + hour main effects
#' expected(t) = mean_DOY + mean_hour - grand_mean
#' anomaly = observed - expected
make_anomaly <- function(df, var, time_col = "datetime") {
  if (!var %in% names(df)) return(df)
  x <- df[[var]]
  if (!is.numeric(x)) return(df)
  
  df <- df %>%
    mutate(
      .doy  = yday(.data[[time_col]]),
      .hour = hour(.data[[time_col]])
    )
  
  grand <- mean(x, na.rm = TRUE)
  
  doy_means <- df %>%
    group_by(.doy) %>%
    summarize(.m_doy = mean(.data[[var]], na.rm = TRUE), .groups = "drop")
  
  hr_means <- df %>%
    group_by(.hour) %>%
    summarize(.m_hr = mean(.data[[var]], na.rm = TRUE), .groups = "drop")
  
  df %>%
    left_join(doy_means, by = ".doy") %>%
    left_join(hr_means, by = ".hour") %>%
    mutate(
      .expected = .m_doy + .m_hr - grand,
      "{var}_anom" := .data[[var]] - .expected
    ) %>%
    dplyr::select(-.doy, -.hour, -.m_doy, -.m_hr, -.expected)
}

#' Analyze all variables for one window
#' @param flux_data Data frame with flux measurements
#' @param met_roll Data frame with rolled predictor columns
#' @param vars Base variable names (without _roll suffix)
#' @param window_days Window size in days
#' @param flux_col Name of flux column
#' @param roll_suffix Suffix for rolled columns (e.g., "_roll" or "_anom_roll")
analyze_window <- function(flux_data, met_roll, vars, window_days, flux_col = "CH4_flux", roll_suffix = "_roll") {
  combined <- flux_data %>% left_join(met_roll, by = "datetime")
  results <- list()
  
  for (site_name in unique(combined$site)) {
    site_data <- combined %>% filter(site == site_name)
    
    for (var in vars) {
      var_roll <- paste0(var, roll_suffix)
      if (!var_roll %in% names(site_data)) next
      
      cor_result <- safe_cor_test(site_data[[flux_col]], site_data[[var_roll]])
      
      results[[length(results) + 1]] <- tibble(
        site = site_name,
        variable = var,
        window_days = window_days,
        window_hours = as.integer(round(window_days * 24)),
        r = cor_result$r,
        p = cor_result$p,
        n = cor_result$n
      )
    }
  }
  bind_rows(results)
}

#' Get variable group
get_var_group <- function(var) {
  case_when(
    var %in% VAR_GROUPS$atmospheric   ~ "Atmosphere",
    var %in% VAR_GROUPS$radiation     ~ "Radiation",
    var %in% VAR_GROUPS$temperature   ~ "Temperature",
    var %in% VAR_GROUPS$soil_moisture ~ "Soil Moisture",
    var %in% VAR_GROUPS$water_table   ~ "Water Table",
    var %in% VAR_GROUPS$precipitation ~ "Precipitation",
    var %in% VAR_GROUPS$ecosystem_fluxes ~ "Ecosystem Fluxes",
    var %in% VAR_GROUPS$gases         ~ "Gas Concentrations",
    var %in% VAR_GROUPS$phenology     ~ "Phenology",
    TRUE                              ~ "Other"
  )
}

#' Signed square root for color scale (enhances mid-range)
signed_sqrt <- function(x) sign(x) * sqrt(abs(x))

# ============================================================
# LOAD DATA
# ============================================================

message("\nLoading data...")

stem_flux_raw <- read_csv(PATHS$stem_flux, show_col_types = FALSE)

stem_flux <- stem_flux_raw %>%
  mutate(
    datetime = round_date(
      force_tz(as.POSIXct(real_start - 2190), tzone = "EST"), 
      "hour"
    ),
    date = as.Date(datetime),
    ID = Tree,
    site = ifelse(Plot == "BGS", "BGS", "EMS"),
    species = case_when(
      ID == 288 ~ "hem",
      ID == 153 ~ "rm",
      ID == 414 ~ "bg",
      ID == 452 ~ "bg",
      TRUE ~ species
    )
  ) %>%
  filter(date >= DATE_MIN, date <= DATE_MAX) %>%
  dplyr::select(datetime, date, ID, site, species, CH4_flux, CO2_flux)

message("  Stem flux: ", nrow(stem_flux), " observations")

aligned_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC")) %>%
  filter(as.Date(datetime) >= DATE_MIN, as.Date(datetime) <= DATE_MAX) %>%
  arrange(datetime)

message("  Environmental data: ", nrow(aligned_data), " hourly records")

# Available variables
available_vars <- ALL_VARS[ALL_VARS %in% names(aligned_data)]
available_vars <- setdiff(available_vars, DROP_VARS)
message("  Predictors: ", length(available_vars))

# ============================================================
# ANOMALY FILTERING (for predictors only)
# ============================================================

message("\nApplying anomaly filtering to predictors...")

# Anomaly-filter environmental predictors
anom_data <- aligned_data
for (v in available_vars) {
  anom_data <- make_anomaly(anom_data, v, time_col = "datetime")
}

# Prepare inputs for both analyses
# RAW: raw predictors
# ANOMALY: anomaly-filtered predictors (but still raw CH4 flux)
anom_vars <- paste0(available_vars, "_anom")
rolling_input_anom <- anom_data %>% dplyr::select(datetime, all_of(anom_vars))
rolling_input_raw <- aligned_data %>% dplyr::select(datetime, all_of(available_vars))

# ============================================================
# ROLLING WINDOW CORRELATIONS - BOTH RAW AND ANOMALY
# ============================================================

run_correlation_analysis <- function(rolling_input, flux_data, flux_col, base_vars, input_vars, roll_suffix, analysis_name) {
  # base_vars: the base variable names for output (e.g., "tair_C")
  # input_vars: the column names in rolling_input (e.g., "tair_C" or "tair_C_anom")
  # roll_suffix: suffix after rolling (e.g., "_roll" or "_anom_roll")
  
  message("\n  Running ", analysis_name, " analysis...")
  
  all_results <- vector("list", length(WINDOWS_HOURS))
  t_start <- Sys.time()
  
  for (i in seq_along(WINDOWS_HOURS)) {
    window_hours <- WINDOWS_HOURS[i]
    window_days <- WINDOWS_DAYS[i]
    
    if (i %% 40 == 0 || i == length(WINDOWS_HOURS)) {
      elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
      rate <- i / max(1, elapsed)
      remaining <- (length(WINDOWS_HOURS) - i) / max(1e-9, rate)
      message(sprintf("    %d/%d (%dh) - %.0fs elapsed, ~%.0fs remaining",
                      i, length(WINDOWS_HOURS), window_hours, elapsed, remaining))
    }
    
    # Rolling means - creates columns named {input_var}_roll
    met_roll <- calc_rolling_means(rolling_input, input_vars, window_hours)
    
    # Determine the suffix pattern for analyze_window
    # If input_vars are like "tair_C_anom", rolled columns are "tair_C_anom_roll"
    # We want to match base_vars (e.g., "tair_C") to rolled columns
    # So we need to tell analyze_window the full suffix: "_anom_roll" for anomaly, "_roll" for raw
    
    all_results[[i]] <- analyze_window(flux_data, met_roll, base_vars, window_days, flux_col, roll_suffix)
  }
  
  message("    Done in ", round(difftime(Sys.time(), t_start, units = "mins"), 1), " min")
  
  bind_rows(all_results)
}

message("\nComputing correlations across ", length(WINDOWS_HOURS), " windows...")

# Run RAW analysis (raw CH4 vs raw predictors)
# Columns in rolling_input_raw: tair_C, RH, etc.
# After rolling: tair_C_roll, RH_roll, etc.
cor_results_raw <- run_correlation_analysis(
  rolling_input = rolling_input_raw,
  flux_data = stem_flux,
  flux_col = "CH4_flux",
  base_vars = available_vars,
  input_vars = available_vars,
  roll_suffix = "_roll",
  analysis_name = "RAW"
) %>%
  mutate(
    analysis = "raw",
    var_label = VAR_LABELS[variable],
    var_label = ifelse(is.na(var_label), variable, var_label),
    var_group = get_var_group(variable)
  )

# Run ANOMALY analysis (raw CH4 vs anomaly-filtered predictors)
# Columns in rolling_input_anom: tair_C_anom, RH_anom, etc.
# After rolling: tair_C_anom_roll, RH_anom_roll, etc.
cor_results_anom <- run_correlation_analysis(
  rolling_input = rolling_input_anom,
  flux_data = stem_flux,
  flux_col = "CH4_flux",
  base_vars = available_vars,
  input_vars = anom_vars,
  roll_suffix = "_anom_roll",
  analysis_name = "ANOMALY"
) %>%
  mutate(
    analysis = "anomaly",
    var_label = VAR_LABELS[variable],
    var_label = ifelse(is.na(var_label), variable, var_label),
    var_group = get_var_group(variable)
  )

# Combine both
cor_results_all <- bind_rows(cor_results_raw, cor_results_anom)

message("\n  Total correlations: ", nrow(cor_results_all), 
        " (", nrow(cor_results_raw), " raw + ", nrow(cor_results_anom), " anomaly)")

# Save full results
write_csv(cor_results_all, file.path(OUTPUT_DIR, "correlations_all_windows.csv"))

# ============================================================
# VARIABLE-LEVEL INFERENCE - FOR BOTH ANALYSES
# ============================================================

message("\nVariable-level inference...")

process_best_windows <- function(cor_results, analysis_name) {
  best <- cor_results %>%
    filter(!is.na(r), !is.na(p)) %>%
    group_by(site, variable) %>%
    arrange(desc(abs(r)), window_hours, p) %>%
    slice(1) %>%
    ungroup() %>%
    group_by(site) %>%
    mutate(
      p_adj = p.adjust(p, method = "BH"),
      significant = p_adj < FDR_ALPHA
    ) %>%
    ungroup() %>%
    arrange(site, p_adj)
  
  sig_summary <- best %>%
    group_by(site) %>%
    summarize(n_total = n(), n_sig = sum(significant), .groups = "drop")
  
  message("  ", analysis_name, " - Significant predictors (BH-FDR < ", FDR_ALPHA, "):")
  for (i in 1:nrow(sig_summary)) {
    message("    ", sig_summary$site[i], ": ", sig_summary$n_sig[i], " / ", sig_summary$n_total[i])
  }
  
  best
}

# Process both
best_windows_raw <- process_best_windows(cor_results_raw, "RAW")
best_windows_anom <- process_best_windows(cor_results_anom, "ANOMALY")

# Add analysis labels and combine
best_windows_raw <- best_windows_raw %>% mutate(analysis = "raw")
best_windows_anom <- best_windows_anom %>% mutate(analysis = "anomaly")
best_windows_all <- bind_rows(best_windows_raw, best_windows_anom)

# Save
write_csv(best_windows_all, file.path(OUTPUT_DIR, "best_windows_by_variable.csv"))
write_csv(best_windows_all %>% filter(significant), file.path(OUTPUT_DIR, "significant_variables.csv"))

# For plotting, we'll use both - create separate objects for convenience
cor_results <- cor_results_all
best_windows <- best_windows_all

# ============================================================
# PLOTTING
# ============================================================

message("\nGenerating figures...")

# Color scale setup
r_range <- max(abs(cor_results$r), na.rm = TRUE)
r_limit <- c(-r_range, r_range)

cor_colors <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                "#F7F7F7",
                "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

# Create manual color palette for variable groups (9 groups)
group_colors <- c(
  "Atmosphere" = "#66c2a5",
  "Radiation" = "#fc8d62", 
  "Temperature" = "#8da0cb",
  "Soil Moisture" = "#e78ac3",
  "Water Table" = "#a6d854",
  "Precipitation" = "#ffd92f",
  "Ecosystem Fluxes" = "#e5c494",
  "Gas Concentrations" = "#b3b3b3",
  "Phenology" = "#1b9e77"
)

SITES <- c("EMS", "BGS")
ANALYSES <- c("raw", "anomaly")

# ------------------------------------------------------------
# FIGURE 1: Heatmaps - one per analysis type
# ------------------------------------------------------------

HEATMAP_HOURS <- c(3, 6, 12, 24, 48, 72, 120, 168, 336)

for (analysis_type in ANALYSES) {
  
  heatmap_data <- cor_results %>%
    filter(analysis == analysis_type, !is.na(r), window_hours %in% HEATMAP_HOURS) %>%
    mutate(
      window_label = case_when(
        window_hours < 24 ~ paste0(window_hours, "h"),
        TRUE ~ paste0(round(window_hours / 24), "d")
      ),
      window_label = factor(window_label, levels = unique(window_label[order(window_hours)]))
    )
  
  heatmap_stars <- best_windows %>%
    filter(analysis == analysis_type, significant) %>%
    mutate(
      window_label = case_when(
        window_hours < 24 ~ paste0(window_hours, "h"),
        TRUE ~ paste0(round(window_hours / 24), "d")
      )
    ) %>%
    filter(window_hours %in% HEATMAP_HOURS)
  
  title_suffix <- ifelse(analysis_type == "raw", "(Raw)", "(Anomaly-filtered)")
  
  p_heat <- ggplot(heatmap_data, aes(x = window_label, y = reorder(var_label, r))) +
    geom_tile(aes(fill = r), color = "white", linewidth = 0.2) +
    # White border behind stars
    geom_point(
      data = heatmap_stars,
      aes(x = window_label, y = var_label),
      shape = 8, size = 2.8, color = "white", stroke = 1.6
    ) +
    # Black stars on top
    geom_point(
      data = heatmap_stars,
      aes(x = window_label, y = var_label),
      shape = 8, size = 2, color = "black", stroke = 0.8
    ) +
    scale_fill_gradientn(
      colors = cor_colors,
      values = scales::rescale(signed_sqrt(seq(-1, 1, length.out = 11))),
      limits = r_limit,
      name = "r"
    ) +
    facet_grid(var_group ~ site, scales = "free_y", space = "free_y") +
    labs(
      title = paste("CH4 Stem Flux Correlations", title_suffix),
      subtitle = paste0("* = significant (BH-FDR < ", FDR_ALPHA, ") at selected window"),
      x = "Rolling Mean Window",
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold", size = 9),
      strip.text.y = element_text(angle = 0, hjust = 0),
      panel.grid = element_blank(),
      panel.spacing = unit(0.3, "lines"),
      legend.position = "right"
    )
  
  n_vars <- length(unique(heatmap_data$var_label))
  fig_height <- max(10, n_vars * 0.35 + 3)
  
  ggsave(file.path(OUTPUT_DIR, paste0("heatmap_", analysis_type, ".png")), p_heat, 
         width = 12, height = fig_height, dpi = 600)
  ggsave(file.path(OUTPUT_DIR, paste0("heatmap_", analysis_type, ".pdf")), p_heat, 
         width = 12, height = fig_height)
  
  message("  Saved: heatmap_", analysis_type, ".png/pdf")
}

# ------------------------------------------------------------
# FIGURE 2: Surface rows (full resolution) - per site per analysis
# ------------------------------------------------------------

for (analysis_type in ANALYSES) {
  for (site_name in SITES) {
    
    site_data <- cor_results %>%
      filter(analysis == analysis_type, site == site_name, !is.na(r))
    
    if (nrow(site_data) == 0) next
    
    # Order by extreme correlation value
    var_order <- site_data %>%
      group_by(variable, var_label, var_group) %>%
      reframe(extreme_r = r[which.max(abs(r))]) %>%
      distinct() %>%
      arrange(desc(extreme_r)) %>%
      mutate(var_num = row_number())
    
    site_data <- site_data %>%
      left_join(var_order %>% dplyr::select(variable, var_num), by = "variable")
    
    # Stars for significant variables at their best window
    site_stars <- best_windows %>%
      filter(analysis == analysis_type, site == site_name, significant) %>%
      left_join(var_order %>% dplyr::select(variable, var_num), by = "variable")
    
    # Complete grid
    site_data_full <- site_data %>%
      tidyr::complete(
        nesting(variable, var_label, var_group, var_num),
        window_hours = WINDOWS_HOURS
      )
    
    title_suffix <- ifelse(analysis_type == "raw", "(Raw)", "(Anomaly)")
    
    p_surface <- ggplot(site_data_full, aes(x = window_hours, y = 1)) +
      geom_raster(aes(fill = r), interpolate = TRUE) +
      # White border behind stars
      geom_point(
        data = site_stars,
        aes(x = window_hours, y = 1),
        shape = 8, size = 4, color = "white", stroke = 2
      ) +
      # Black stars on top
      geom_point(
        data = site_stars,
        aes(x = window_hours, y = 1),
        shape = 8, size = 3, color = "black", stroke = 1.2
      ) +
      scale_fill_gradientn(
        colors = cor_colors,
        values = scales::rescale(signed_sqrt(seq(-1, 1, length.out = 11))),
        limits = r_limit,
        name = "r",
        na.value = "grey90"
      ) +
      scale_x_continuous(
        breaks = c(3, 12, 24, 48, 72, 120, 168, 240, 336),
        labels = c("3h", "12h", "1d", "2d", "3d", "5d", "7d", "10d", "14d"),
        expand = c(0, 0)
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      facet_wrap(~ reorder(var_label, var_num), ncol = 1, strip.position = "left") +
      labs(
        title = paste("Correlation Surface -", site_name, title_suffix),
        subtitle = paste0("Ordered by extreme r; * = significant (BH-FDR < ", FDR_ALPHA, ")"),
        x = "Rolling Mean Window",
        y = NULL
      ) +
      theme_minimal(base_size = 21) +
      theme(
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 18),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        legend.position = "right"
      )
    
    n_vars_site <- n_distinct(site_data_full$variable)
    surface_height <- max(8, n_vars_site * 0.4 + 2)
    
    ggsave(file.path(OUTPUT_DIR, paste0("surface_", analysis_type, "_", site_name, ".png")),
           p_surface, width = 12, height = surface_height, dpi = 600)
    ggsave(file.path(OUTPUT_DIR, paste0("surface_", analysis_type, "_", site_name, ".pdf")),
           p_surface, width = 12, height = surface_height)
    
    message("  Saved: surface_", analysis_type, "_", site_name, ".png/pdf")
  }
}

# ------------------------------------------------------------
# FIGURE 2b: Surface rows grouped by variable group - per site per analysis
# Same aesthetics as surface plots (geom_raster with interpolation)
# ------------------------------------------------------------

message("  Creating grouped surface plots...")

for (analysis_type in ANALYSES) {
  for (site_name in SITES) {
    
    site_data <- cor_results %>%
      filter(analysis == analysis_type, site == site_name, !is.na(r))
    
    if (nrow(site_data) == 0) next
    
    # Order variables: first by group, then by extreme correlation within group
    var_order <- site_data %>%
      group_by(variable, var_label, var_group) %>%
      reframe(extreme_r = r[which.max(abs(r))]) %>%
      distinct() %>%
      # Create group order (alphabetical or custom)
      mutate(group_order = match(var_group, c("Atmosphere", "Radiation", "Temperature", 
                                              "Soil Moisture", "Water Table", "Precipitation",
                                              "Ecosystem Fluxes", "Gas Concentrations", "Phenology"))) %>%
      group_by(var_group) %>%
      arrange(desc(extreme_r), .by_group = TRUE) %>%
      mutate(var_num_in_group = row_number()) %>%
      ungroup() %>%
      arrange(group_order, var_num_in_group) %>%
      mutate(var_num = row_number())
    
    site_data <- site_data %>%
      left_join(var_order %>% dplyr::select(variable, var_num, var_group), by = c("variable", "var_group"))
    
    # Stars for significant variables - best_windows already has var_label and var_group
    # Just need to add var_num from var_order
    site_stars <- best_windows %>%
      filter(analysis == analysis_type, site == site_name, significant) %>%
      left_join(var_order %>% dplyr::select(variable, var_num), by = "variable") %>%
      mutate(facet_label = paste0(var_group, ": ", var_label))
    
    # Complete grid
    site_data_full <- site_data %>%
      tidyr::complete(
        nesting(variable, var_label, var_group, var_num),
        window_hours = WINDOWS_HOURS
      )
    
    # Create facet label with group prefix
    site_data_full <- site_data_full %>%
      mutate(facet_label = paste0(var_group, ": ", var_label))
    
    title_suffix <- ifelse(analysis_type == "raw", "(Raw)", "(Anomaly)")
    
    p_surface_grouped <- ggplot(site_data_full, aes(x = window_hours, y = 1)) +
      geom_raster(aes(fill = r), interpolate = TRUE) +
      # White border behind stars
      geom_point(
        data = site_stars,
        aes(x = window_hours, y = 1),
        shape = 8, size = 4, color = "white", stroke = 2
      ) +
      # Black stars on top
      geom_point(
        data = site_stars,
        aes(x = window_hours, y = 1),
        shape = 8, size = 3, color = "black", stroke = 1.2
      ) +
      scale_fill_gradientn(
        colors = cor_colors,
        values = scales::rescale(signed_sqrt(seq(-1, 1, length.out = 11))),
        limits = r_limit,
        name = "r",
        na.value = "grey90"
      ) +
      scale_x_continuous(
        breaks = c(3, 12, 24, 48, 72, 120, 168, 240, 336),
        labels = c("3h", "12h", "1d", "2d", "3d", "5d", "7d", "10d", "14d"),
        expand = c(0, 0)
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      facet_wrap(~ reorder(facet_label, var_num), ncol = 1, strip.position = "left") +
      labs(
        title = paste("Correlation Surface by Group -", site_name, title_suffix),
        subtitle = paste0("Ordered by group then by extreme r; * = significant (BH-FDR < ", FDR_ALPHA, ")"),
        x = "Rolling Mean Window",
        y = NULL
      ) +
      theme_minimal(base_size = 21) +
      theme(
        strip.text.y.left = element_text(angle = 0, hjust = 1, size = 18),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.1, "lines"),
        legend.position = "right"
      )
    
    n_vars <- length(unique(site_data_full$var_label))
    fig_height <- max(10, n_vars * 0.4 + 3)
    
    ggsave(file.path(OUTPUT_DIR, paste0("surface_grouped_", analysis_type, "_", site_name, ".png")),
           p_surface_grouped, width = 14, height = fig_height, dpi = 600)
    ggsave(file.path(OUTPUT_DIR, paste0("surface_grouped_", analysis_type, "_", site_name, ".pdf")),
           p_surface_grouped, width = 14, height = fig_height)
    
    message("  Saved: surface_grouped_", analysis_type, "_", site_name, ".png/pdf")
  }
}

# ------------------------------------------------------------
# FIGURE 2c: Surface rows grouped - both sites side by side, per analysis
# Same aesthetics as surface plots, grouped by var_group
# ------------------------------------------------------------

for (analysis_type in ANALYSES) {
  
  both_sites_data <- cor_results %>%
    filter(analysis == analysis_type, site %in% SITES, !is.na(r))
  
  # Order variables: first by group, then by extreme correlation within group
  var_order_both <- both_sites_data %>%
    group_by(variable, var_label, var_group) %>%
    reframe(extreme_r = r[which.max(abs(r))]) %>%
    distinct() %>%
    mutate(group_order = match(var_group, c("Atmosphere", "Radiation", "Temperature", 
                                            "Soil Moisture", "Water Table", "Precipitation",
                                            "Ecosystem Fluxes", "Gas Concentrations", "Phenology"))) %>%
    group_by(var_group) %>%
    arrange(desc(extreme_r), .by_group = TRUE) %>%
    mutate(var_num_in_group = row_number()) %>%
    ungroup() %>%
    arrange(group_order, var_num_in_group) %>%
    mutate(var_num = row_number())
  
  both_sites_data <- both_sites_data %>%
    left_join(var_order_both %>% dplyr::select(variable, var_num, var_num_in_group), by = "variable")
  
  both_stars <- best_windows %>%
    filter(analysis == analysis_type, site %in% SITES, significant) %>%
    left_join(var_order_both %>% dplyr::select(variable, var_num, var_num_in_group), by = "variable")
  
  # Complete grid
  both_sites_full <- both_sites_data %>%
    tidyr::complete(
      nesting(variable, var_label, var_group, var_num, var_num_in_group),
      site = SITES,
      window_hours = WINDOWS_HOURS
    )
  
  title_suffix <- ifelse(analysis_type == "raw", "(Raw)", "(Anomaly-filtered)")
  
  # Use facet_grid with var_group on rows and site on columns
  # y-axis shows individual variables within each group
  p_surface_both <- ggplot(both_sites_full, 
                           aes(x = window_hours, y = reorder(var_label, -var_num_in_group))) +
    geom_tile(aes(fill = r), color = NA) +
    # White border behind stars
    geom_point(
      data = both_stars,
      aes(x = window_hours, y = var_label),
      shape = 8, size = 2.8, color = "white", stroke = 1.6
    ) +
    # Black stars on top
    geom_point(
      data = both_stars,
      aes(x = window_hours, y = var_label),
      shape = 8, size = 2, color = "black", stroke = 0.8
    ) +
    scale_fill_gradientn(
      colors = cor_colors,
      values = scales::rescale(signed_sqrt(seq(-1, 1, length.out = 11))),
      limits = r_limit,
      name = "r",
      na.value = "grey90"
    ) +
    scale_x_continuous(
      breaks = c(3, 12, 24, 48, 96, 144, 192, 240, 288, 336),
      labels = c("3h", "12h", "1d", "2d", "4d", "6d", "8d", "10d", "12d", "14d"),
      expand = c(0, 0),
      guide = guide_axis(n.dodge = 2)
    ) +
    facet_grid(var_group ~ site, scales = "free_y", space = "free_y") +
    labs(
      x = "Rolling Mean Window",
      y = NULL
    ) +
    theme_minimal(base_size = 22) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),
      axis.text.y = element_text(size = 16),
      axis.ticks.x = element_line(color = "black"),
      axis.ticks.length.x = unit(0.15, "cm"),
      strip.text = element_text(face = "bold", size = 20),
      strip.text.y = element_text(angle = 0, hjust = 0),
      panel.grid = element_blank(),
      panel.spacing.x = unit(2, "lines"),
      panel.spacing.y = unit(0.3, "lines"),
      legend.position = "right"
    )
  
  n_vars <- length(unique(both_sites_full$var_label))
  fig_height <- max(10, n_vars * 0.35 + 3)
  
  ggsave(file.path(OUTPUT_DIR, paste0("surface_grouped_both_", analysis_type, ".png")),
         p_surface_both, width = 20, height = fig_height, dpi = 600)
  ggsave(file.path(OUTPUT_DIR, paste0("surface_grouped_both_", analysis_type, ".pdf")),
         p_surface_both, width = 20, height = fig_height)
  
  message("  Saved: surface_grouped_both_", analysis_type, ".png/pdf")
}

# ------------------------------------------------------------
# FIGURE 3: Bar plot of best correlations - per analysis
# ------------------------------------------------------------

for (analysis_type in ANALYSES) {
  
  best_for_plot <- best_windows %>%
    filter(analysis == analysis_type)
  
  title_suffix <- ifelse(analysis_type == "raw", "(Raw)", "(Anomaly)")
  
  p_bars <- ggplot(best_for_plot,
                   aes(x = reorder(var_label, r), y = r, fill = var_group)) +
    geom_col(width = 0.7) +
    geom_point(aes(shape = significant), size = 2, color = "black") +
    scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 1), guide = "none") +
    geom_text(
      aes(label = ifelse(
        window_hours < 24, 
        paste0(window_hours, "h"),
        paste0(round(window_days, 1), "d")
      )),
      hjust = ifelse(best_for_plot$r > 0, -0.1, 1.1),
      size = 2.5, color = "gray30"
    ) +
    scale_fill_manual(values = group_colors, name = "Group") +
    facet_wrap(~ site, scales = "free_y") +
    coord_flip(ylim = r_limit) +
    labs(
      title = paste("Best Window Correlations with CH4 Stem Flux", title_suffix),
      subtitle = paste0("Labels = selected window; * = BH-FDR < ", FDR_ALPHA),
      x = NULL, y = "Pearson r"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    ) +
    guides(fill = guide_legend(nrow = 2))
  
  ggsave(file.path(OUTPUT_DIR, paste0("barplot_", analysis_type, ".png")),
         p_bars, width = 12, height = max(8, nrow(best_for_plot) * 0.12 + 3), dpi = 600)
  
  message("  Saved: barplot_", analysis_type, ".png")
}

# ------------------------------------------------------------
# FIGURE 4: Line plots by variable group - per analysis
# ------------------------------------------------------------

for (analysis_type in ANALYSES) {
  for (grp in names(VAR_GROUPS)) {
    grp_vars <- intersect(VAR_GROUPS[[grp]], available_vars)
    if (length(grp_vars) == 0) next
    
    grp_data <- cor_results %>%
      filter(analysis == analysis_type, variable %in% grp_vars, !is.na(r))
    
    if (nrow(grp_data) == 0) next
    
    # Mark significant variables - include var_label for color mapping
    grp_sig <- best_windows %>%
      filter(analysis == analysis_type, variable %in% grp_vars, significant) %>%
      dplyr::select(site, variable, var_label, window_hours, r)
    
    title_suffix <- ifelse(analysis_type == "raw", "(Raw)", "(Anomaly)")
    
    p_lines <- ggplot(grp_data, aes(x = window_days, y = r, color = var_label)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_line(linewidth = 0.7) +
      # White border behind stars
      geom_point(
        data = grp_sig,
        aes(x = window_hours / 24, y = r),
        shape = 8, size = 4, color = "white", stroke = 1.8,
        inherit.aes = FALSE
      ) +
      # Colored stars on top
      geom_point(
        data = grp_sig,
        aes(x = window_hours / 24, y = r, color = var_label),
        shape = 8, size = 3, stroke = 1,
        inherit.aes = FALSE
      ) +
      scale_x_log10(
        breaks = c(3/24, 6/24, 1, 3, 7, 14),
        labels = c("3h", "6h", "1d", "3d", "7d", "14d")
      ) +
      scale_color_brewer(palette = "Set1", name = NULL) +
      facet_wrap(~ site, ncol = 2) +
      coord_cartesian(ylim = r_limit) +
      labs(
        title = paste("CH4 Correlations:", tools::toTitleCase(gsub("_", " ", grp)), title_suffix),
        subtitle = paste0("* = significant at selected window"),
        x = "Rolling Mean Window",
        y = "Pearson r"
      ) +
      theme_minimal(base_size = 21) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold")
      )
    
    ggsave(file.path(OUTPUT_DIR, "lines_by_group", paste0("lines_", grp, "_", analysis_type, ".png")),
           p_lines, width = 10, height = 5, dpi = 600)
  }
}

message("  Saved: line plots by group")

# ============================================================
# SUMMARY
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("RESULTS SUMMARY")
message(paste(rep("=", 60), collapse = ""))

for (analysis_type in ANALYSES) {
  message("\n--- ", toupper(analysis_type), " ANALYSIS ---")
  message("\nSignificant variables by site:\n")
  
  sig_table <- best_windows %>%
    filter(analysis == analysis_type, significant) %>%
    arrange(site, p_adj) %>%
    dplyr::select(site, variable, window_hours, r, p, p_adj) %>%
    mutate(
      r = round(r, 3),
      p = formatC(p, format = "e", digits = 2),
      p_adj = formatC(p_adj, format = "e", digits = 2)
    )
  
  print(sig_table, n = 50)
}

# Comparison table
message("\n--- COMPARISON: RAW vs ANOMALY ---\n")

comparison <- best_windows_raw %>%
  dplyr::select(site, variable, var_label, 
                window_hours_raw = window_hours, r_raw = r, sig_raw = significant) %>%
  left_join(
    best_windows_anom %>%
      dplyr::select(site, variable, 
                    window_hours_anom = window_hours, r_anom = r, sig_anom = significant),
    by = c("site", "variable")
  ) %>%
  mutate(
    r_diff = r_anom - r_raw,
    same_sign = sign(r_raw) == sign(r_anom)
  ) %>%
  arrange(site, desc(abs(r_diff)))

message("Variables with largest difference between raw and anomaly:")
print(comparison %>% head(20), n = 20)

write_csv(comparison, file.path(OUTPUT_DIR, "comparison_raw_vs_anomaly.csv"))

# ============================================================
# FIGURE 6: Time series of all variables grouped by variable group
# ============================================================

message("\nGenerating time series plots...")

# Prepare data - pivot to long format
ts_data <- aligned_data %>%
  dplyr::select(datetime, all_of(available_vars)) %>%
  pivot_longer(
    cols = -datetime,
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    var_label = VAR_LABELS[variable],
    var_label = ifelse(is.na(var_label), variable, var_label),
    var_group = get_var_group(variable)
  ) %>%
  filter(!is.na(value))

# Order variables within groups
var_order_ts <- ts_data %>%
  group_by(variable, var_label, var_group) %>%
  summarize(mean_val = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(group_order = match(var_group, c("Atmosphere", "Radiation", "Temperature", 
                                          "Soil Moisture", "Water Table", "Precipitation",
                                          "Ecosystem Fluxes", "Gas Concentrations", "Phenology"))) %>%
  group_by(var_group) %>%
  arrange(var_label, .by_group = TRUE) %>%
  mutate(var_num_in_group = row_number()) %>%
  ungroup() %>%
  arrange(group_order, var_num_in_group) %>%
  mutate(var_num = row_number())

ts_data <- ts_data %>%
  left_join(var_order_ts %>% dplyr::select(variable, var_num), by = "variable")

# Create facet label
ts_data <- ts_data %>%
  mutate(facet_label = paste0(var_group, ": ", var_label))

# Normalize values within each variable for comparable visualization
ts_data <- ts_data %>%
  group_by(variable) %>%
  mutate(
    value_scaled = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE)
  ) %>%
  ungroup()

# Plot - one panel per variable, grouped by var_group
p_timeseries <- ggplot(ts_data, aes(x = datetime, y = value_scaled)) +
  geom_line(linewidth = 0.3, alpha = 0.7, color = "steelblue") +
  facet_wrap(~ reorder(facet_label, var_num), ncol = 1, scales = "free_y", strip.position = "left") +
  labs(
    title = "Environmental Variable Time Series",
    subtitle = "Z-scored within each variable for visual comparison",
    x = "Date",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.spacing = unit(0.1, "lines"),
    legend.position = "none"
  )

n_vars_ts <- n_distinct(ts_data$variable)
ts_height <- max(12, n_vars_ts * 0.5 + 2)

ggsave(file.path(OUTPUT_DIR, "timeseries_all_variables.png"),
       p_timeseries, width = 14, height = ts_height, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "timeseries_all_variables.pdf"),
       p_timeseries, width = 14, height = ts_height)

message("  Saved: timeseries_all_variables.png/pdf")

# Also create a version with raw (unscaled) values, one plot per group
dir.create(file.path(OUTPUT_DIR, "timeseries_by_group"), recursive = TRUE, showWarnings = FALSE)

for (grp in unique(ts_data$var_group)) {
  grp_data <- ts_data %>%
    filter(var_group == grp)
  
  if (nrow(grp_data) == 0) next
  
  p_grp <- ggplot(grp_data, aes(x = datetime, y = value, color = var_label)) +
    geom_line(linewidth = 0.4, alpha = 0.8) +
    scale_color_brewer(palette = "Set2", name = NULL) +
    labs(
      title = paste("Time Series:", grp),
      x = "Date",
      y = "Value (original units)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(nrow = 2))
  
  grp_filename <- tolower(gsub(" ", "_", grp))
  ggsave(file.path(OUTPUT_DIR, "timeseries_by_group", paste0("timeseries_", grp_filename, ".png")),
         p_grp, width = 12, height = 6, dpi = 300)
}

message("  Saved: timeseries by group plots")

message("\nOutputs saved to: ", OUTPUT_DIR)
message("Done!")
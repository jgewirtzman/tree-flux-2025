# ============================================================
# rolling_correlations.R
#
# Rolling-window correlation analysis between CH4 stem flux and
# environmental predictors.
#
# Run AFTER: timeseries.R (which creates the corrected flux file)
#
# Approach:
#   - Uses asinh-transformed CH4 flux (handles skewed distribution)
#   - Both raw and anomaly-filtered correlations
#   - Anomaly filtering removes DOY + hour main effects
#   - 3-hour resolution windows (3h to 14 days)
#   - BH-FDR correction across variables (not windows)
#
# Outputs:
#   - CSV: correlations_all_windows.csv, best_windows_by_variable.csv
#   - Figures: heatmaps, surface plots, line plots by group, bar plots
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
  stem_flux = "data/HF_2023-2025_tree_flux_corrected.csv",
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

# Variables to drop (insufficient data or redundant)
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

# Use asinh transform on already-corrected CH4 flux
stem_flux <- stem_flux_raw %>%
  mutate(
    datetime = round_date(
      force_tz(as.POSIXct(datetime_posx), tzone = "EST"), 
      "hour"
    ),
    date = as.Date(datetime),
    ID = Tree,
    site = factor(ifelse(PLOT == "BGS", "Wetland", "Upland"), levels = c("Wetland", "Upland")),
    CH4_flux = asinh(CH4_flux_nmolpm2ps)
  ) %>%
  filter(date >= DATE_MIN, date <= DATE_MAX) %>%
  dplyr::select(datetime, date, ID, site, species = SPECIES, CH4_flux)

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

aligned_anom <- aligned_data
for (var in available_vars) {
  aligned_anom <- make_anomaly(aligned_anom, var)
}

anom_vars <- paste0(available_vars, "_anom")
message("  Created ", length(anom_vars), " anomaly variables")

# ============================================================
# RUN ROLLING WINDOW ANALYSIS
# ============================================================

ANALYSES <- c("raw", "anomaly")

message("\nRunning rolling window analysis...")
t_start <- Sys.time()

all_results <- list()

for (i in seq_along(WINDOWS_DAYS)) {
  window_days <- WINDOWS_DAYS[i]
  window_hours <- as.integer(round(window_days * 24))
  
  if (i %% 20 == 0 || i == length(WINDOWS_DAYS)) {
    elapsed <- difftime(Sys.time(), t_start, units = "secs")
    message(sprintf("  Window %d/%d (%.1f days) - %.0fs elapsed",
                    i, length(WINDOWS_DAYS), window_days, elapsed))
  }
  
  # Raw analysis
  met_roll_raw <- calc_rolling_means(aligned_data, available_vars, window_hours)
  results_raw <- analyze_window(stem_flux, met_roll_raw, available_vars, window_days,
                                flux_col = "CH4_flux", roll_suffix = "_roll") %>%
    mutate(analysis = "raw")
  
  # Anomaly analysis
  met_roll_anom <- calc_rolling_means(aligned_anom, anom_vars, window_hours)
  names(met_roll_anom) <- gsub("_anom_roll", "_roll", names(met_roll_anom))
  results_anom <- analyze_window(stem_flux, met_roll_anom, available_vars, window_days,
                                 flux_col = "CH4_flux", roll_suffix = "_roll") %>%
    mutate(analysis = "anomaly")
  
  all_results[[i]] <- bind_rows(results_raw, results_anom)
}

message("  Completed in ", round(difftime(Sys.time(), t_start, units = "mins"), 1), " minutes")

# Combine results
cor_results_all <- bind_rows(all_results) %>%
  mutate(
    var_label = VAR_LABELS[variable],
    var_label = ifelse(is.na(var_label), variable, var_label),
    var_group = get_var_group(variable)
  )

message("\nTotal correlations: ", nrow(cor_results_all))

write_csv(cor_results_all, file.path(OUTPUT_DIR, "correlations_all_windows.csv"))

# ============================================================
# SELECT BEST WINDOW PER VARIABLE + FDR CORRECTION
# ============================================================

message("\nSelecting optimal windows and applying FDR correction...")

best_windows_all <- cor_results_all %>%
  filter(!is.na(r)) %>%
  group_by(analysis, site, variable) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(analysis, site) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(
    significant = p_adj < FDR_ALPHA,
    var_label = VAR_LABELS[variable],
    var_label = ifelse(is.na(var_label), variable, var_label),
    var_group = get_var_group(variable)
  ) %>%
  arrange(analysis, site, p_adj)

write_csv(best_windows_all, file.path(OUTPUT_DIR, "best_windows_by_variable.csv"))
write_csv(best_windows_all %>% filter(significant), file.path(OUTPUT_DIR, "significant_variables.csv"))

# Split by analysis type
best_windows_raw <- best_windows_all %>% filter(analysis == "raw")
best_windows_anom <- best_windows_all %>% filter(analysis == "anomaly")

message("  Significant variables (raw):     ", sum(best_windows_raw$significant))
message("  Significant variables (anomaly): ", sum(best_windows_anom$significant))

# ============================================================
# VISUALIZATION
# ============================================================

message("\nGenerating figures...")

# Dynamic color scale limits (symmetric, based on max |r|, rounded to nearest 0.05)
r_range <- ceiling(max(abs(cor_results_all$r), na.rm = TRUE) * 20) / 20
r_limit <- c(-r_range, r_range)
message("  Color scale r_limit: ", r_limit[1], " to ", r_limit[2])

# ------------------------------------------------------------
# FIGURE 1: Heatmaps
# ------------------------------------------------------------

for (analysis_type in ANALYSES) {
  
  cor_data <- cor_results_all %>%
    filter(analysis == analysis_type, !is.na(r))
  
  best_for_heatmap <- best_windows_all %>%
    filter(analysis == analysis_type)
  
  # Order variables by mean |r|
  var_order <- cor_data %>%
    group_by(variable, var_label, var_group) %>%
    summarize(mean_abs_r = mean(abs(r), na.rm = TRUE), .groups = "drop") %>%
    arrange(var_group, desc(mean_abs_r)) %>%
    mutate(var_num = row_number())
  
  cor_data <- cor_data %>%
    left_join(var_order %>% dplyr::select(variable, var_num), by = "variable")
  
  title_suffix <- ifelse(analysis_type == "raw", "(Raw)", "(Anomaly-filtered)")
  
  p_heat <- ggplot(cor_data, aes(x = window_hours, y = reorder(var_label, -var_num))) +
    geom_tile(aes(fill = r), color = NA) +
    geom_point(
      data = best_for_heatmap %>% filter(significant),
      aes(x = window_hours, y = var_label),
      shape = 8, size = 2, color = "black"
    ) +
    scale_fill_gradientn(
      colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                 "#F7F7F7",
                 "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
      values = scales::rescale(signed_sqrt(seq(-1, 1, length.out = 11))),
      limits = r_limit,
      name = "r"
    ) +
    scale_x_continuous(
      breaks = c(3, 12, 24, 72, 168, 336),
      labels = c("3h", "12h", "1d", "3d", "7d", "14d")
    ) +
    facet_grid(var_group ~ site, scales = "free_y", space = "free_y") +
    labs(
      title = paste("CH4 Flux Correlations", title_suffix),
      subtitle = "✳ = significant after FDR correction",
      x = "Rolling Mean Window",
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.y = element_text(size = 8),
      strip.text = element_text(face = "bold", size = 9),
      strip.text.y = element_text(angle = 0, hjust = 0),
      panel.grid = element_blank(),
      panel.spacing = unit(0.3, "lines"),
      legend.position = "right"
    )
  
  n_vars <- length(unique(cor_data$var_label))
  fig_height <- max(10, n_vars * 0.35 + 3)
  
  print(p_heat)
  ggsave(file.path(OUTPUT_DIR, paste0("heatmap_", analysis_type, ".png")), p_heat, 
         width = 12, height = fig_height, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, paste0("heatmap_", analysis_type, ".pdf")), p_heat, 
         width = 12, height = fig_height)
  
  message("  Saved: heatmap_", analysis_type, ".png/pdf")
}

# ------------------------------------------------------------
# FIGURE 2: Bar plots of optimal window correlations
# ------------------------------------------------------------

for (analysis_type in ANALYSES) {
  
  best_for_plot <- best_windows_all %>%
    filter(analysis == analysis_type)
  
  title_suffix <- ifelse(analysis_type == "raw", "(Raw)", "(Anomaly)")
  
  p_bars <- ggplot(best_for_plot, 
                   aes(x = reorder(var_label, r), y = r, fill = var_group)) +
    geom_col(width = 0.7) +
    geom_point(
      data = best_for_plot %>% filter(significant),
      aes(x = var_label, y = r),
      shape = 8, size = 2, color = "black"
    ) +
    geom_text(
      aes(label = paste0(window_hours, "h")),
      hjust = ifelse(best_for_plot$r > 0, -0.1, 1.1),
      size = 2.5, color = "gray30"
    ) +
    scale_fill_brewer(palette = "Set2", name = "Group") +
    facet_wrap(~ site, scales = "free_y") +
    coord_flip(ylim = r_limit) +
    labs(
      title = paste("Optimal Window Correlations with asinh(CH4)", title_suffix),
      subtitle = "Labels show optimal window; ✳ = significant (FDR < 0.05)",
      x = NULL, y = "Pearson r"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    ) +
    guides(fill = guide_legend(nrow = 2))
  
  print(p_bars)
  ggsave(file.path(OUTPUT_DIR, paste0("barplot_", analysis_type, ".png")),
         p_bars, width = 12, height = max(8, nrow(best_for_plot) * 0.12 + 3), dpi = 300)
  
  message("  Saved: barplot_", analysis_type, ".png")
}

# ------------------------------------------------------------
# FIGURE 3: Line plots by variable group
# ------------------------------------------------------------

for (analysis_type in ANALYSES) {
  for (grp in names(VAR_GROUPS)) {
    grp_vars <- intersect(VAR_GROUPS[[grp]], available_vars)
    if (length(grp_vars) == 0) next
    
    grp_data <- cor_results_all %>%
      filter(analysis == analysis_type, variable %in% grp_vars, !is.na(r))
    
    if (nrow(grp_data) == 0) next
    
    grp_sig <- best_windows_all %>%
      filter(analysis == analysis_type, variable %in% grp_vars, significant) %>%
      dplyr::select(site, variable, var_label, window_hours, r)
    
    title_suffix <- ifelse(analysis_type == "raw", "(Raw)", "(Anomaly)")
    
    p_lines <- ggplot(grp_data, aes(x = window_days, y = r, color = var_label)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      geom_line(linewidth = 0.7) +
      geom_point(
        data = grp_sig,
        aes(x = window_hours / 24, y = r),
        shape = 8, size = 4, color = "white", stroke = 1.8,
        inherit.aes = FALSE
      ) +
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
        subtitle = "✳ = significant at optimal window",
        x = "Rolling Mean Window",
        y = "Pearson r"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold")
      )
    
    print(p_lines)
    ggsave(file.path(OUTPUT_DIR, "lines_by_group", paste0("lines_", grp, "_", analysis_type, ".png")),
           p_lines, width = 10, height = 5, dpi = 300)
  }
}

message("  Saved: line plots by group")

# ------------------------------------------------------------
# FIGURE 4: Surface plots grouped by variable group (per site)
# ------------------------------------------------------------

message("  Creating grouped surface plots...")

# Color palette
cor_colors <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                "#F7F7F7",
                "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")

SITES <- c("Wetland", "Upland")

for (analysis_type in ANALYSES) {
  for (site_name in SITES) {
    
    site_data <- cor_results_all %>%
      filter(analysis == analysis_type, site == site_name, !is.na(r))
    
    if (nrow(site_data) == 0) next
    
    # Order variables: first by group, then by extreme correlation within group
    var_order <- site_data %>%
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
    
    site_data <- site_data %>%
      left_join(var_order %>% dplyr::select(variable, var_num, var_group), by = c("variable", "var_group"))
    
    # Stars for significant variables
    site_stars <- best_windows_all %>%
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
    
    # Use lowercase for filenames
    site_code <- tolower(site_name)
    
    print(p_surface_grouped)
    ggsave(file.path(OUTPUT_DIR, paste0("surface_grouped_", analysis_type, "_", site_code, ".png")),
           p_surface_grouped, width = 14, height = fig_height, dpi = 600)
    ggsave(file.path(OUTPUT_DIR, paste0("surface_grouped_", analysis_type, "_", site_code, ".pdf")),
           p_surface_grouped, width = 14, height = fig_height)
    
    message("  Saved: surface_grouped_", analysis_type, "_", site_code, ".png/pdf")
  }
}

# ------------------------------------------------------------
# FIGURE 5: Surface plots - both sites side by side
# ------------------------------------------------------------

for (analysis_type in ANALYSES) {
  
  both_sites_data <- cor_results_all %>%
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
  
  both_stars <- best_windows_all %>%
    filter(analysis == analysis_type, site %in% SITES, significant) %>%
    left_join(var_order_both %>% dplyr::select(variable, var_num, var_num_in_group), by = "variable")
  
  both_stars <- both_stars %>%
    mutate(site = factor(site, levels = c("Wetland", "Upland")))
  
  # Complete grid
  both_sites_full <- both_sites_data %>%
    tidyr::complete(
      nesting(variable, var_label, var_group, var_num, var_num_in_group),
      site = SITES,
      window_hours = WINDOWS_HOURS
    )
  
  both_sites_full <- both_sites_full %>%
    mutate(site = factor(site, levels = c("Wetland", "Upland")))
  
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
  
  print(p_surface_both)
  ggsave(file.path(OUTPUT_DIR, paste0("surface_grouped_both_", analysis_type, ".png")),
         p_surface_both, width = 20, height = fig_height, dpi = 600)
  ggsave(file.path(OUTPUT_DIR, paste0("surface_grouped_both_", analysis_type, ".pdf")),
         p_surface_both, width = 20, height = fig_height)
  
  message("  Saved: surface_grouped_both_", analysis_type, ".png/pdf")
}

# ============================================================
# SUMMARY
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("RESULTS SUMMARY")
message(paste(rep("=", 60), collapse = ""))

for (analysis_type in ANALYSES) {
  message("\n--- ", toupper(analysis_type), " ANALYSIS ---")
  message("\nSignificant variables (FDR < 0.05):\n")
  
  sig_table <- best_windows_all %>%
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

message("\nOutputs saved to: ", OUTPUT_DIR)
message("Done!")
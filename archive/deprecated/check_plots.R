# ============================================================
# 07b_plot_timeseries.R
# Plot timeseries of all variables in the aligned dataset
# ============================================================

library(tidyverse)
library(lubridate)
library(patchwork)

# ============================================================
# CONFIGURATION
# ============================================================

DATA_PATH <- "data/processed/aligned_hourly_dataset.csv"
OUTPUT_DIR <- "figures/timeseries"

# Date range for all plots
DATE_MIN <- as.Date("2023-01-01")
DATE_MAX <- as.Date("2026-01-01")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD DATA
# ============================================================

message("Loading aligned dataset...")
aligned_data <- read_csv(DATA_PATH, show_col_types = FALSE)
aligned_data <- aligned_data %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC"))

message("  ", nrow(aligned_data), " rows, ", ncol(aligned_data), " columns")

# ============================================================
# DEFINE VARIABLE GROUPS FOR PLOTTING
# ============================================================

# Get all numeric columns (excluding time columns)
time_cols <- c("datetime", "date", "year", "month", "doy", "hour", "wyear")
all_vars <- names(aligned_data)[!names(aligned_data) %in% time_cols]
all_vars <- all_vars[sapply(aligned_data[all_vars], is.numeric)]

message("  ", length(all_vars), " variables to plot")

# Define variable groups based on simplified structure
var_groups <- list(
  # Fluxes (averaged within tower)
  "CO2_Flux_FC" = all_vars[grepl("^FC_", all_vars)],
  "CO2_Storage_SC" = all_vars[grepl("^SC_", all_vars)],
  "Latent_Heat_LE" = all_vars[grepl("^LE_", all_vars)],
  "Sensible_Heat_H" = all_vars[grepl("^H_", all_vars)],
  "Ground_Heat_G" = all_vars[grepl("^G_", all_vars)],
  
  # Gas mixing ratios (from xHA and Ha2)
  "CH4_MixingRatio" = all_vars[grepl("^CH4_MR", all_vars)],
  "CO2_MixingRatio" = all_vars[grepl("^CO2_MR", all_vars)],
  
  # Radiation - Fisher only (PAR, slrr, rnet)
  "Radiation_Fisher" = all_vars[grepl("^PAR$|^slrr$|^rnet$", all_vars)],
  
  # Temperature - Fisher air temp, xHA canopy temp, tower soil temp
  "Air_Temp_Fisher" = all_vars[grepl("^tair_C$", all_vars)],
  "Canopy_Temp_xHA" = all_vars[grepl("^T_CANOPY", all_vars)],
  "Soil_Temp_Fisher" = all_vars[grepl("^s10t$", all_vars)],
  "Soil_Temp_Tower" = all_vars[grepl("^TS_", all_vars)],
  
  # Atmospheric moisture - Fisher only
  "RH_Fisher" = all_vars[grepl("^RH$", all_vars)],
  "VPD_Fisher" = all_vars[grepl("^VPD_kPa$", all_vars)],
  
  # Soil moisture (NEON + Ha1/Ha2 towers)
  "SWC_NEON" = all_vars[grepl("^NEON_SWC", all_vars)],
  "SWC_Tower" = all_vars[grepl("^SWC_", all_vars)],
  "Water_Table_HF" = all_vars[grepl("wtd", all_vars, ignore.case = TRUE)],
  
  # Precipitation - Fisher only
  "Precip_Fisher" = all_vars[grepl("^P_mm$", all_vars)],
  # Throughfall - xHA
  "Throughfall_xHA" = all_vars[grepl("^THROUGHFALL", all_vars)],
  
  # Pressure - Fisher only
  "Pressure_Fisher" = all_vars[grepl("^p_kPa$", all_vars)],
  
  # Wind - xHA
  "Wind_xHA" = all_vars[grepl("^WS_|^WD_|^USTAR", all_vars)],
  
  # Phenology
  "Phenocam" = all_vars[grepl("^gcc|^ndvi", all_vars, ignore.case = TRUE)]
)

# Remove empty groups
var_groups <- var_groups[sapply(var_groups, length) > 0]

message("\nVariable groups:")
for (grp in names(var_groups)) {
  message("  ", grp, ": ", length(var_groups[[grp]]), " variables")
}

# ============================================================
# PLOTTING FUNCTION
# ============================================================

plot_timeseries_group <- function(data, vars, group_name, 
                                  resample = "daily",
                                  show_legend = TRUE) {
  
  if (length(vars) == 0) return(NULL)
  
  # Check which vars exist and have data
  vars <- vars[vars %in% names(data)]
  vars <- vars[sapply(data[vars], function(x) sum(!is.na(x)) > 0)]
  
  if (length(vars) == 0) {
    message("  No valid data for ", group_name)
    return(NULL)
  }
  
  # Resample if requested
  if (resample == "daily") {
    plot_data <- data %>%
      mutate(date = as.Date(datetime)) %>%
      filter(date >= DATE_MIN, date <= DATE_MAX) %>%
      group_by(date) %>%
      summarize(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      mutate(datetime = as.POSIXct(date))
  } else if (resample == "weekly") {
    plot_data <- data %>%
      filter(as.Date(datetime) >= DATE_MIN, as.Date(datetime) <= DATE_MAX) %>%
      mutate(week = floor_date(datetime, "week")) %>%
      group_by(week) %>%
      summarize(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
      rename(datetime = week)
  } else {
    plot_data <- data %>%
      filter(as.Date(datetime) >= DATE_MIN, as.Date(datetime) <= DATE_MAX)
  }
  
  # Pivot to long format
  plot_long <- plot_data %>%
    select(datetime, all_of(vars)) %>%
    pivot_longer(-datetime, names_to = "variable", values_to = "value") %>%
    filter(!is.na(value))
  
  # Extract tower suffix for coloring
  plot_long <- plot_long %>%
    mutate(
      tower = case_when(
        grepl("_Ha1$", variable) ~ "Ha1",
        grepl("_Ha2$", variable) ~ "Ha2",
        grepl("_xHA$", variable) ~ "xHA",
        grepl("^gcc|^ndvi", variable) ~ "Phenocam",
        TRUE ~ "HF"
      ),
      var_short = gsub("_(Ha1|Ha2|xHA)$", "", variable)
    )
  
  # Create plot with fixed x-axis limits
  p <- ggplot(plot_long, aes(x = datetime, y = value, color = variable)) +
    geom_line(alpha = 0.7, linewidth = 0.4) +
    scale_x_datetime(limits = c(as.POSIXct(DATE_MIN), as.POSIXct(DATE_MAX)),
                     date_breaks = "6 months", date_labels = "%b %Y") +
    labs(
      title = gsub("_", " ", group_name),
      x = NULL,
      y = NULL,
      color = NULL
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 12, face = "bold"),
      legend.position = if (show_legend && length(vars) <= 12) "right" else "bottom",
      legend.text = element_text(size = 7),
      axis.text = element_text(size = 8),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Adjust legend for many variables
  if (length(vars) > 12) {
    p <- p + guides(color = guide_legend(ncol = 4, override.aes = list(linewidth = 1.5)))
  }
  
  p
}

# ============================================================
# PLOT ALL GROUPS - DAILY RESOLUTION
# ============================================================

message("\n--- Generating daily timeseries plots ---")

plots_daily <- list()

for (grp in names(var_groups)) {
  message("  Plotting: ", grp)
  p <- plot_timeseries_group(aligned_data, var_groups[[grp]], grp, resample = "daily")
  if (!is.null(p)) {
    plots_daily[[grp]] <- p
    
    # Save individual plot
    ggsave(
      file.path(OUTPUT_DIR, paste0("daily_", grp, ".png")),
      p, width = 12, height = 5, dpi = 150
    )
  }
}

# ============================================================
# COMBINED MULTI-PANEL FIGURES
# ============================================================

message("\n--- Creating combined figures ---")

# Figure 1: CO2 Fluxes (FC and SC)
flux_groups <- c("CO2_Flux_FC", "CO2_Storage_SC")
flux_plots <- plots_daily[flux_groups[flux_groups %in% names(plots_daily)]]

if (length(flux_plots) > 0) {
  p_flux <- wrap_plots(flux_plots, ncol = 1) +
    plot_annotation(title = "CO2 Fluxes - Daily Means",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  ggsave(file.path(OUTPUT_DIR, "combined_CO2_fluxes.png"), p_flux, 
         width = 14, height = 3 * length(flux_plots), dpi = 150)
  message("  Saved: combined_CO2_fluxes.png")
}

# Figure 2: Energy Fluxes (LE, H, G)
energy_groups <- c("Latent_Heat_LE", "Sensible_Heat_H", "Ground_Heat_G")
energy_plots <- plots_daily[energy_groups[energy_groups %in% names(plots_daily)]]

if (length(energy_plots) > 0) {
  p_energy <- wrap_plots(energy_plots, ncol = 1) +
    plot_annotation(title = "Energy Fluxes - Daily Means",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  ggsave(file.path(OUTPUT_DIR, "combined_energy_fluxes.png"), p_energy,
         width = 14, height = 3 * length(energy_plots), dpi = 150)
  message("  Saved: combined_energy_fluxes.png")
}

# Figure 3: Gas Mixing Ratios (CH4, CO2 only - no H2O)
gas_groups <- c("CH4_MixingRatio", "CO2_MixingRatio")
gas_plots <- plots_daily[gas_groups[gas_groups %in% names(plots_daily)]]

if (length(gas_plots) > 0) {
  p_gas <- wrap_plots(gas_plots, ncol = 1) +
    plot_annotation(title = "Gas Mixing Ratios - Daily Means",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  ggsave(file.path(OUTPUT_DIR, "combined_gas_mixing_ratios.png"), p_gas,
         width = 14, height = 3 * length(gas_plots), dpi = 150)
  message("  Saved: combined_gas_mixing_ratios.png")
}

# Figure 4: Radiation (Fisher only - PAR, slrr, rnet)
rad_groups <- c("Radiation_Fisher")
rad_plots <- plots_daily[rad_groups[rad_groups %in% names(plots_daily)]]

if (length(rad_plots) > 0) {
  p_rad <- wrap_plots(rad_plots, ncol = 1) +
    plot_annotation(title = "Radiation (Fisher Met) - Daily Means",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  ggsave(file.path(OUTPUT_DIR, "combined_radiation.png"), p_rad,
         width = 14, height = 4, dpi = 150)
  message("  Saved: combined_radiation.png")
}

# Figure 5: Temperature (Fisher air, xHA canopy, Fisher/Tower soil)
temp_groups <- c("Air_Temp_Fisher", "Canopy_Temp_xHA", "Soil_Temp_Fisher", "Soil_Temp_Tower")
temp_plots <- plots_daily[temp_groups[temp_groups %in% names(plots_daily)]]

if (length(temp_plots) > 0) {
  p_temp <- wrap_plots(temp_plots, ncol = 1) +
    plot_annotation(title = "Temperature - Daily Means",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  ggsave(file.path(OUTPUT_DIR, "combined_temperature.png"), p_temp,
         width = 14, height = 3 * length(temp_plots), dpi = 150)
  message("  Saved: combined_temperature.png")
}

# Figure 6: Atmospheric Moisture (Fisher RH, VPD)
atm_moist_groups <- c("RH_Fisher", "VPD_Fisher")
atm_moist_plots <- plots_daily[atm_moist_groups[atm_moist_groups %in% names(plots_daily)]]

if (length(atm_moist_plots) > 0) {
  p_atm <- wrap_plots(atm_moist_plots, ncol = 1) +
    plot_annotation(title = "Atmospheric Moisture (Fisher Met) - Daily Means",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  ggsave(file.path(OUTPUT_DIR, "combined_atm_moisture.png"), p_atm,
         width = 14, height = 3 * length(atm_moist_plots), dpi = 150)
  message("  Saved: combined_atm_moisture.png")
}

# Figure 7: Soil Moisture (NEON, Ha1/Ha2 towers, Water Table)
soil_moist_groups <- c("SWC_NEON", "SWC_Tower", "Water_Table_HF")
soil_moist_plots <- plots_daily[soil_moist_groups[soil_moist_groups %in% names(plots_daily)]]

if (length(soil_moist_plots) > 0) {
  p_soil <- wrap_plots(soil_moist_plots, ncol = 1) +
    plot_annotation(title = "Soil Moisture & Water Table - Daily Means",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  ggsave(file.path(OUTPUT_DIR, "combined_soil_moisture.png"), p_soil,
         width = 14, height = 3 * length(soil_moist_plots), dpi = 150)
  message("  Saved: combined_soil_moisture.png")
}

# Figure 8: Precipitation & Throughfall (Fisher P_mm, xHA throughfall)
precip_groups <- c("Precip_Fisher", "Throughfall_xHA")
precip_plots <- plots_daily[precip_groups[precip_groups %in% names(plots_daily)]]

if (length(precip_plots) > 0) {
  p_precip <- wrap_plots(precip_plots, ncol = 1) +
    plot_annotation(title = "Precipitation & Throughfall - Daily Means",
                    theme = theme(plot.title = element_text(size = 14, face = "bold")))
  ggsave(file.path(OUTPUT_DIR, "combined_precipitation.png"), p_precip,
         width = 14, height = 3 * length(precip_plots), dpi = 150)
  message("  Saved: combined_precipitation.png")
}

# Figure 9: Phenology
if ("Phenocam" %in% names(plots_daily)) {
  ggsave(file.path(OUTPUT_DIR, "phenocam.png"), plots_daily[["Phenocam"]],
         width = 12, height = 4, dpi = 150)
  message("  Saved: phenocam.png")
}

# ============================================================
# FACETED PLOTS BY TOWER
# ============================================================

message("\n--- Creating tower comparison plots ---")

# Function for faceted tower comparison
plot_tower_comparison <- function(data, vars, title, y_label = NULL) {
  
  vars <- vars[vars %in% names(data)]
  vars <- vars[sapply(data[vars], function(x) sum(!is.na(x)) > 0)]
  
  if (length(vars) == 0) return(NULL)
  
  # Daily means with date filter
  plot_data <- data %>%
    mutate(date = as.Date(datetime)) %>%
    filter(date >= DATE_MIN, date <= DATE_MAX) %>%
    group_by(date) %>%
    summarize(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")
  
  # Long format with tower extraction
  plot_long <- plot_data %>%
    pivot_longer(-date, names_to = "variable", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(
      tower = case_when(
        grepl("_Ha1$", variable) ~ "Ha1",
        grepl("_Ha2$", variable) ~ "Ha2",
        grepl("_xHA$", variable) ~ "xHA",
        TRUE ~ "Other"
      ),
      var_base = gsub("_(Ha1|Ha2|xHA)$", "", variable)
    )
  
  ggplot(plot_long, aes(x = date, y = value, color = var_base)) +
    geom_line(alpha = 0.7, linewidth = 0.5) +
    scale_x_date(limits = c(DATE_MIN, DATE_MAX),
                 date_breaks = "6 months", date_labels = "%b %Y") +
    facet_wrap(~ tower, ncol = 1, scales = "free_y") +
    labs(title = title, x = NULL, y = y_label, color = NULL) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    guides(color = guide_legend(ncol = 4))
}

# FC comparison across towers (with replicates)
p_fc_compare <- plot_tower_comparison(
  aligned_data, 
  var_groups$CO2_Flux_FC,
  "CO2 Flux (FC) by Tower - Including Replicates",
  expression("FC ("*mu*"mol m"^-2*" s"^-1*")")
)
if (!is.null(p_fc_compare)) {
  ggsave(file.path(OUTPUT_DIR, "tower_compare_FC.png"), p_fc_compare,
         width = 12, height = 8, dpi = 150)
  message("  Saved: tower_compare_FC.png")
}

# LE comparison across towers (with replicates)
p_le_compare <- plot_tower_comparison(
  aligned_data,
  var_groups$Latent_Heat_LE,
  "Latent Heat Flux (LE) by Tower - Including Replicates",
  expression("LE (W m"^-2*")")
)
if (!is.null(p_le_compare)) {
  ggsave(file.path(OUTPUT_DIR, "tower_compare_LE.png"), p_le_compare,
         width = 12, height = 8, dpi = 150)
  message("  Saved: tower_compare_LE.png")
}

# Soil moisture comparison (NEON vs Tower)
swc_vars <- c(var_groups$SWC_NEON, var_groups$SWC_Tower)
if (length(swc_vars) > 0) {
  swc_data <- aligned_data %>%
    mutate(date = as.Date(datetime)) %>%
    filter(date >= DATE_MIN, date <= DATE_MAX) %>%
    group_by(date) %>%
    summarize(across(any_of(swc_vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    pivot_longer(-date, names_to = "variable", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(
      source = case_when(
        grepl("^NEON", variable) ~ "NEON",
        grepl("_Ha1$", variable) ~ "Ha1",
        grepl("_Ha2$", variable) ~ "Ha2",
        grepl("_xHA$", variable) ~ "xHA",
        TRUE ~ "Other"
      )
    )
  
  p_swc <- ggplot(swc_data, aes(x = date, y = value, color = variable)) +
    geom_line(alpha = 0.7, linewidth = 0.5) +
    scale_x_date(limits = c(DATE_MIN, DATE_MAX),
                 date_breaks = "6 months", date_labels = "%b %Y") +
    facet_wrap(~ source, ncol = 1, scales = "free_y") +
    labs(
      title = "Soil Water Content Comparison",
      subtitle = "NEON (by depth) vs AmeriFlux Towers (baseline-adjusted)",
      x = NULL,
      y = "SWC",
      color = NULL
    ) +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(OUTPUT_DIR, "soil_moisture_comparison.png"), p_swc,
         width = 12, height = 8, dpi = 150)
  message("  Saved: soil_moisture_comparison.png")
}

# ============================================================
# SUMMARY
# ============================================================

message("\n============================================================")
message("PLOTTING COMPLETE")
message("============================================================")
message("Output directory: ", OUTPUT_DIR)
message("Individual plots: ", length(list.files(OUTPUT_DIR, pattern = "daily_.*\\.png")))
message("Combined plots: ", length(list.files(OUTPUT_DIR, pattern = "combined_.*\\.png")))
message("Comparison plots: ", length(list.files(OUTPUT_DIR, pattern = "tower_.*\\.png|soil_moisture.*\\.png")))
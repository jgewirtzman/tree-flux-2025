# ============================================================
# 08_rolling_window_correlations.R
# 
# Exploratory analysis: correlations between CH4 stem flux and 
# environmental variables across multiple time integration windows.
#
# Purpose: Variable dplyr::selection for mixed models
# ============================================================

library(tidyverse)
library(lubridate)
library(RcppRoll)
library(patchwork)
library(viridis)

# ============================================================
# CONFIGURATION
# ============================================================

# File paths
PATHS <- list(
  stem_flux = "data/raw/flux_dataset.csv",
  aligned = "data/processed/aligned_hourly_dataset.csv"
)

OUTPUT_DIR <- "figures/rolling_correlations"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Date range
DATE_MIN <- as.Date("2023-01-01")
DATE_MAX <- as.Date("2026-01-01")

# Time windows to analyze (in days)
# Finer resolution for better surface visualization
WINDOWS_DAYS <- c(
  # Sub-daily (hours)
  1/24, 2/24, 3/24, 4/24, 6/24, 8/24, 12/24, 18/24,
  # Daily to weekly (finer spacing)
  1, 1.5, 2, 2.5, 3, 4, 5, 6, 7,
  # Weekly to monthly
  8, 10, 12, 14, 17, 21, 30
)

# Variables to correlate with CH4 flux
VAR_GROUPS <- list(
  atmospheric = c("tair_C", "RH", "VPD_kPa", "p_kPa"),
  radiation = c("PAR", "rnet", "slrr"),
  soil_temp = c("s10t", "TS_Ha1", "TS_Ha2", "TS_xHA"),
  soil_moisture = c("NEON_SWC_shallow", "NEON_SWC_mid", "NEON_SWC_deep",
                    "SWC_Ha1", "SWC_Ha2"),
  water_table = c("bgs_wtd_cm", "bvs_wtd_cm"),
  precipitation = c("P_mm", "THROUGHFALL_xHA"), 
  tower_fluxes = c("FC_Ha1", "FC_Ha2", "FC_xHA", 
                   "LE_Ha1", "LE_Ha2",
                   "H_Ha1", "H_Ha2"),
  gases = c("CH4_MR_xHA", "CO2_MR_xHA", "CO2_MR_Ha2"),
  canopy = c("T_CANOPY_xHA", "gcc", "ndvi"),
  other = c("G_xHA", "USTAR_Ha1", "USTAR_Ha2", "USTAR_xHA")
)

ALL_VARS <- unlist(VAR_GROUPS, use.names = FALSE)

# Variable display names
VAR_LABELS <- c(
  tair_C = "Air temp", RH = "RH", VPD_kPa = "VPD", p_kPa = "Pressure",
  PAR = "PAR", rnet = "Net rad", slrr = "Solar rad",
  s10t = "Soil temp 10cm", TS_Ha1 = "Soil temp Ha1", TS_Ha2 = "Soil temp Ha2", TS_xHA = "Soil temp xHA",
  NEON_SWC_shallow = "NEON SWC shallow", NEON_SWC_mid = "NEON SWC mid", NEON_SWC_deep = "NEON SWC deep",
  SWC_Ha1 = "SWC Ha1", SWC_Ha2 = "SWC Ha2",
  bgs_wtd_cm = "WTD BGS", bvs_wtd_cm = "WTD BVS",
  P_mm = "Precip", THROUGHFALL_xHA = "Throughfall",
  FC_Ha1 = "FC Ha1", FC_Ha2 = "FC Ha2", FC_xHA = "FC xHA",
  LE_Ha1 = "LE Ha1", LE_Ha2 = "LE Ha2",
  H_Ha1 = "H Ha1", H_Ha2 = "H Ha2",
  CH4_MR_xHA = "CH4 MR", CO2_MR_xHA = "CO2 MR xHA", CO2_MR_Ha2 = "CO2 MR Ha2",
  T_CANOPY_xHA = "Canopy temp", gcc = "GCC", ndvi = "NDVI",
  G_xHA = "Ground heat", USTAR_Ha1 = "USTAR Ha1", USTAR_Ha2 = "USTAR Ha2", USTAR_xHA = "USTAR xHA"
)

# Variables to sum (precipitation) vs mean (everything else)
SUM_VARS <- c("P_mm", "THROUGHFALL_xHA")

# ============================================================
# HELPER FUNCTIONS
# ============================================================

calc_rolling_stats <- function(data, vars, window_hours) {
  data <- data %>% arrange(datetime)
  result <- data %>% dplyr::select(datetime)
  
  for (var in vars) {
    if (!var %in% names(data)) next
    x <- data[[var]]
    
    if (var %in% SUM_VARS) {
      result[[paste0(var, "_roll")]] <- roll_sum(x, n = window_hours, 
                                                 align = "right", fill = NA, na.rm = TRUE)
    } else {
      result[[paste0(var, "_roll")]] <- roll_mean(x, n = window_hours, 
                                                  align = "right", fill = NA, na.rm = TRUE)
    }
  }
  result
}

safe_cor_test <- function(x, y) {
  valid <- complete.cases(x, y)
  x <- x[valid]; y <- y[valid]
  n <- length(x)
  
  if (n < 5) return(list(r = NA_real_, p = NA_real_, n = n))
  
  tryCatch({
    test <- cor.test(x, y, method = "pearson")
    list(r = as.numeric(test$estimate), p = as.numeric(test$p.value), n = n)
  }, error = function(e) list(r = NA_real_, p = NA_real_, n = n))
}

analyze_window <- function(flux_data, met_data, vars, window_days) {
  combined <- flux_data %>% left_join(met_data, by = "datetime")
  results <- list()
  
  for (site_name in unique(combined$site)) {
    site_data <- combined %>% filter(site == site_name)
    
    for (var in vars) {
      var_roll <- paste0(var, "_roll")
      if (!var_roll %in% names(site_data)) next
      
      cor_result <- safe_cor_test(site_data$CH4_flux, site_data[[var_roll]])
      
      results[[length(results) + 1]] <- tibble(
        site = site_name, variable = var, window_days = window_days,
        r = cor_result$r, p = cor_result$p, n = cor_result$n
      )
    }
  }
  bind_rows(results)
}

# ============================================================
# LOAD DATA
# ============================================================

message("Loading data...")

stem_flux_raw <- read_csv(PATHS$stem_flux, show_col_types = FALSE)

stem_flux <- stem_flux_raw %>%
  mutate(
    datetime = round_date(force_tz(as.POSIXct(real_start - 2190), tzone = "EST"), "hour"),
    date = as.Date(datetime),
    ID = Tree,
    site = ifelse(Plot == "BGS", "BGS", "EMS"),
    species = case_when(
      ID == 288 ~ "hem", ID == 153 ~ "rm", ID == 414 ~ "bg", ID == 452 ~ "bg",
      TRUE ~ species
    )
  ) %>%
  filter(date >= DATE_MIN, date <= DATE_MAX) %>%
  dplyr::select(datetime, date, ID, site, species, CH4_flux, CO2_flux)

message("  Stem flux: ", nrow(stem_flux), " measurements")

aligned_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC")) %>%
  filter(as.Date(datetime) >= DATE_MIN, as.Date(datetime) <= DATE_MAX)

message("  Aligned data: ", nrow(aligned_data), " hourly records")

available_vars <- ALL_VARS[ALL_VARS %in% names(aligned_data)]
message("  Available variables: ", length(available_vars))

# ============================================================
# RUN ANALYSIS
# ============================================================

message("\nRunning rolling window analysis...")

all_results <- list()

for (i in seq_along(WINDOWS_DAYS)) {
  window_days <- WINDOWS_DAYS[i]
  window_hours <- round(window_days * 24)
  
  message(sprintf("  Window %2d/%d: %.2f days...", i, length(WINDOWS_DAYS), window_days))
  
  met_roll <- calc_rolling_stats(aligned_data, available_vars, window_hours)
  results <- analyze_window(stem_flux, met_roll, available_vars, window_days)
  all_results[[i]] <- results
}

cor_results <- bind_rows(all_results) %>%
  mutate(
    significant = p < 0.05,
    var_label = VAR_LABELS[variable],
    var_label = ifelse(is.na(var_label), variable, var_label),
    var_group = case_when(
      variable %in% VAR_GROUPS$atmospheric ~ "Atmospheric",
      variable %in% VAR_GROUPS$radiation ~ "Radiation",
      variable %in% VAR_GROUPS$soil_temp ~ "Soil Temperature",
      variable %in% VAR_GROUPS$soil_moisture ~ "Soil Moisture",
      variable %in% VAR_GROUPS$water_table ~ "Water Table",
      variable %in% VAR_GROUPS$precipitation ~ "Precipitation",
      variable %in% VAR_GROUPS$tower_fluxes ~ "Tower Fluxes",
      variable %in% VAR_GROUPS$gases ~ "Gas Concentrations",
      variable %in% VAR_GROUPS$canopy ~ "Canopy/Phenology",
      TRUE ~ "Other"
    )
  )

message("\nTotal correlations: ", nrow(cor_results))

write_csv(cor_results, file.path(OUTPUT_DIR, "rolling_correlations.csv"))

# ============================================================
# COMPREHENSIVE HEATMAP
# ============================================================

message("\nGenerating figures...")

# Window labels
heatmap_data <- cor_results %>%
  filter(!is.na(r)) %>%
  mutate(
    window_label = case_when(
      window_days < 1 ~ paste0(round(window_days * 24), "h"),
      TRUE ~ paste0(round(window_days, 1), "d")
    ),
    window_label = factor(window_label, levels = unique(window_label[order(window_days)]))
  )

# Custom color scale function - signed square root to enhance mid-range differences
# This compresses extreme values and expands differences near zero
signed_sqrt <- function(x) sign(x) * sqrt(abs(x))
signed_sqrt_inv <- function(x) sign(x) * x^2

# Combined heatmap - both sites
p_heat <- ggplot(heatmap_data, 
                 aes(x = window_label, y = reorder(var_label, r))) +
  geom_tile(aes(fill = r), color = "white", size = 0.2) +
  geom_point(data = heatmap_data %>% filter(significant),
             shape = 16, size = 0.8, color = "black") +
  scale_fill_gradientn(
    colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
               "#F7F7F7",
               "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
    values = scales::rescale(signed_sqrt(seq(-1, 1, length.out = 11))),
    limits = c(-1, 1),
    name = "r"
  ) +
  facet_grid(var_group ~ site, scales = "free_y", space = "free_y") +
  labs(
    title = "CH4 Stem Flux Correlations with Environmental Variables",
    subtitle = "Dots indicate p < 0.05. Color scale enhanced for mid-range visibility.",
    x = "Integration Window",
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
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12)
  )

# Calculate figure height based on number of variables
n_vars <- length(unique(heatmap_data$var_label))
fig_height <- max(10, n_vars * 0.35 + 3)
fig_width <- 14  # Wider to accommodate more windows

ggsave(file.path(OUTPUT_DIR, "heatmap_comprehensive.png"),
       p_heat, width = fig_width, height = fig_height, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "heatmap_comprehensive.pdf"),
       p_heat, width = fig_width, height = fig_height)

message("  Saved: heatmap_comprehensive.png/pdf")

# ============================================================
# INTERPOLATED SURFACE PLOT
# ============================================================

message("  Creating interpolated surface...")

library(akima)  # For interpolation

for (site_name in unique(cor_results$site)) {
  
  site_data <- cor_results %>%
    filter(site == site_name, !is.na(r))
  
  # Create numeric y-axis (variable order by mean correlation)
  var_order <- site_data %>%
    group_by(variable, var_label) %>%
    summarize(mean_r = mean(r, na.rm = TRUE), .groups = "drop") %>%
    arrange(mean_r) %>%
    mutate(var_num = row_number())
  
  site_data <- site_data %>%
    left_join(var_order %>% dplyr::select(variable, var_num), by = "variable")
  
  # Interpolate to finer grid
  # Use log scale for window to give more weight to short windows
  interp_result <- tryCatch({
    akima::interp(
      x = log10(site_data$window_days),
      y = site_data$var_num,
      z = site_data$r,
      xo = seq(log10(min(WINDOWS_DAYS)), log10(max(WINDOWS_DAYS)), length.out = 100),
      yo = seq(1, max(site_data$var_num), length.out = 100),
      linear = FALSE,  # Use spline interpolation
      extrap = FALSE
    )
  }, error = function(e) NULL)
  
  if (!is.null(interp_result)) {
    # Convert to data frame for ggplot
    interp_df <- expand.grid(
      log_window = interp_result$x,
      var_num = interp_result$y
    ) %>%
      mutate(
        r = as.vector(interp_result$z),
        window_days = 10^log_window
      ) %>%
      filter(!is.na(r))
    
    # Create y-axis labels
    y_breaks <- var_order$var_num
    y_labels <- var_order$var_label
    
    p_surface <- ggplot(interp_df, aes(x = window_days, y = var_num)) +
      geom_raster(aes(fill = r), interpolate = TRUE) +
      geom_contour(aes(z = r), color = "white", alpha = 0.4, size = 0.3,
                   breaks = seq(-0.8, 0.8, by = 0.2)) +
      # Add original data points
      geom_point(data = site_data, aes(x = window_days, y = var_num),
                 shape = 16, size = 0.5, alpha = 0.5) +
      scale_fill_gradientn(
        colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                   "#F7F7F7",
                   "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
        values = scales::rescale(signed_sqrt(seq(-1, 1, length.out = 11))),
        limits = c(-1, 1),
        name = "r"
      ) +
      scale_x_log10(
        breaks = c(1/24, 6/24, 1, 7, 30),
        labels = c("1h", "6h", "1d", "7d", "30d")
      ) +
      scale_y_continuous(breaks = y_breaks, labels = y_labels) +
      labs(
        title = paste("Correlation Surface -", site_name),
        subtitle = "Spline-interpolated. Contours at r = ±0.2, ±0.4, ±0.6, ±0.8",
        x = "Integration Window (log scale)",
        y = NULL
      ) +
      theme_minimal(base_size = 10) +
      theme(
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold")
      )
    
    ggsave(file.path(OUTPUT_DIR, paste0("surface_", site_name, ".png")),
           p_surface, width = 10, height = fig_height, dpi = 300)
    ggsave(file.path(OUTPUT_DIR, paste0("surface_", site_name, ".pdf")),
           p_surface, width = 10, height = fig_height)
    
    message("  Saved: surface_", site_name, ".png/pdf")
  }
}

# ============================================================
# LINE PLOTS BY VARIABLE GROUP
# ============================================================

for (grp in names(VAR_GROUPS)) {
  grp_data <- cor_results %>%
    filter(variable %in% VAR_GROUPS[[grp]], !is.na(r))
  
  if (nrow(grp_data) == 0) next
  
  p_lines <- ggplot(grp_data, aes(x = window_days, y = r, color = var_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(size = 0.8) +
    geom_point(aes(shape = significant), size = 2) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), 
                       name = "p < 0.05", guide = "none") +
    scale_x_log10(breaks = c(0.1, 0.5, 1, 5, 10, 30),
                  labels = c("2h", "12h", "1d", "5d", "10d", "30d")) +
    scale_color_brewer(palette = "Set1", name = NULL) +
    facet_wrap(~ site, ncol = 2) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(
      title = paste("CH4 Correlations:", tools::toTitleCase(gsub("_", " ", grp))),
      x = "Integration Window",
      y = "Pearson r"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(file.path(OUTPUT_DIR, paste0("lines_", grp, ".png")),
         p_lines, width = 10, height = 5, dpi = 300)
}

message("  Saved: line plots by variable group")

# ============================================================
# OPTIMAL WINDOW SUMMARY
# ============================================================

optimal_windows <- cor_results %>%
  filter(!is.na(r)) %>%
  group_by(site, variable, var_label, var_group) %>%
  slice_max(abs(r), n = 1) %>%
  ungroup() %>%
  arrange(site, desc(abs(r)))

write_csv(optimal_windows, file.path(OUTPUT_DIR, "optimal_windows.csv"))

# Bar plot
p_optimal <- ggplot(optimal_windows, 
                    aes(x = reorder(var_label, r), y = r, fill = var_group)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(round(window_days, 1), "d")),
            hjust = ifelse(optimal_windows$r > 0, -0.1, 1.1),
            size = 2.5, color = "gray30") +
  scale_fill_brewer(palette = "Set2", name = "Group") +
  facet_wrap(~ site, scales = "free_y") +
  coord_flip(ylim = c(-1, 1)) +
  labs(
    title = "Optimal Window Correlations with CH4 Stem Flux",
    subtitle = "Labels show optimal integration window",
    x = NULL, y = "Pearson r"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8)
  ) +
  guides(fill = guide_legend(nrow = 2))

ggsave(file.path(OUTPUT_DIR, "optimal_windows.png"),
       p_optimal, width = 12, height = fig_height * 0.8, dpi = 300)

message("  Saved: optimal_windows.png")

# ============================================================
# SUMMARY
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("TOP 15 CORRELATIONS BY SITE")
message(paste(rep("=", 60), collapse = ""))

top_cors <- cor_results %>%
  filter(!is.na(r)) %>%
  group_by(site) %>%
  slice_max(abs(r), n = 15) %>%
  dplyr::select(site, variable, window_days, r, p, n) %>%
  arrange(site, desc(abs(r)))

print(top_cors, n = 30)

message("\nOutputs saved to: ", OUTPUT_DIR)
message("Done!")
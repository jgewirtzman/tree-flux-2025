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
# Hourly resolution: 1 hour to 14 days (336 hours)
WINDOWS_DAYS <- (1:336) / 24  # 1h, 2h, 3h, ..., 336h (14 days)

message("Total windows to analyze: ", length(WINDOWS_DAYS), " (hourly resolution, 1h to 14d)")

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
message("  This may take a few minutes with ", length(WINDOWS_DAYS), " windows...")

all_results <- list()

# Progress tracking
t_start <- Sys.time()
report_interval <- 50  # Report every 50 windows

for (i in seq_along(WINDOWS_DAYS)) {
  window_days <- WINDOWS_DAYS[i]
  window_hours <- max(1, round(window_days * 24))  # Minimum 1 hour
  
  # Progress report
  if (i %% report_interval == 0 || i == length(WINDOWS_DAYS)) {
    elapsed <- difftime(Sys.time(), t_start, units = "secs")
    rate <- i / as.numeric(elapsed)
    remaining <- (length(WINDOWS_DAYS) - i) / rate
    message(sprintf("  Window %d/%d (%.1f days) - elapsed: %.0fs, remaining: ~%.0fs",
                    i, length(WINDOWS_DAYS), window_days, elapsed, remaining))
  }
  
  met_roll <- calc_rolling_stats(aligned_data, available_vars, window_hours)
  results <- analyze_window(stem_flux, met_roll, available_vars, window_days)
  all_results[[i]] <- results
}

message("  Completed in ", round(difftime(Sys.time(), t_start, units = "mins"), 1), " minutes")

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

# Discrete time windows for readable heatmap
HEATMAP_WINDOWS <- c(
  1/24, 2/24, 4/24, 6/24, 12/24,      # Sub-daily
  1, 2, 3, 5, 7,                       # Daily to weekly
  10, 14                               # Up to 2 weeks
)

# Window labels
heatmap_data <- cor_results %>%
  filter(!is.na(r), window_days %in% HEATMAP_WINDOWS) %>%
  mutate(
    window_label = case_when(
      window_days < 1 ~ paste0(round(window_days * 24), "h"),
      TRUE ~ paste0(round(window_days, 1), "d")
    ),
    window_label = factor(window_label, levels = unique(window_label[order(window_days)]))
  )

# Custom color scale function - signed square root to enhance mid-range differences
signed_sqrt <- function(x) sign(x) * sqrt(abs(x))

# Get global min/max correlation values for consistent color scale across sites
r_min <- min(cor_results$r, na.rm = TRUE)
r_max <- max(cor_results$r, na.rm = TRUE)
r_abs_max <- max(abs(r_min), abs(r_max))
r_limit <- c(-r_abs_max, r_abs_max)

message("  Correlation range: ", round(r_min, 3), " to ", round(r_max, 3))
message("  Color scale limits: ", round(r_limit[1], 3), " to ", round(r_limit[2], 3))

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
    limits = r_limit,
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

message("  Creating correlation surface...")

# With hourly resolution, we have enough data that we can just plot directly
# No interpolation needed for the 1h-14d range

for (site_name in unique(cor_results$site)) {
  
  site_data <- cor_results %>%
    filter(site == site_name, !is.na(r))
  
  # Create variable order (by mean correlation within group)
  var_order <- site_data %>%
    group_by(variable, var_label, var_group) %>%
    summarize(mean_r = mean(r, na.rm = TRUE), .groups = "drop") %>%
    arrange(var_group, mean_r) %>%
    mutate(var_num = row_number())
  
  site_data <- site_data %>%
    left_join(var_order %>% dplyr::select(variable, var_num), by = "variable")
  
  # Y-axis labels
  y_breaks <- var_order$var_num
  y_labels <- var_order$var_label
  
  # Surface plot with faceted rows (one per variable, horizontal interpolation only)
  # Linear time scale with geom_raster
  p_surface_rows <- ggplot(site_data, aes(x = window_days * 24, y = 1)) +
    geom_raster(aes(fill = r), interpolate = TRUE) +
    scale_fill_gradientn(
      colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                 "#F7F7F7",
                 "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
      values = scales::rescale(signed_sqrt(seq(-1, 1, length.out = 11))),
      limits = r_limit,
      name = "r"
    ) +
    scale_x_continuous(
      breaks = c(1, 6, 12, 24, 48, 72, 120, 168, 240, 336),
      labels = c("1h", "6h", "12h", "1d", "2d", "3d", "5d", "7d", "10d", "14d"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~ reorder(var_label, var_num), ncol = 1, strip.position = "left") +
    labs(
      title = paste("CH4 Stem Flux Correlation Surface -", site_name),
      subtitle = "Hourly resolution, 1h to 14 days",
      x = "Integration Window",
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "right"
    )
  
  n_vars <- length(unique(site_data$variable))
  surface_height <- max(8, n_vars * 0.4 + 2)
  
  ggsave(file.path(OUTPUT_DIR, paste0("surface_rows_", site_name, ".png")),
         p_surface_rows, width = 12, height = surface_height, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, paste0("surface_rows_", site_name, ".pdf")),
         p_surface_rows, width = 12, height = surface_height)
  
  message("  Saved: surface_rows_", site_name, ".png/pdf")
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
    scale_x_log10(breaks = c(1/24, 6/24, 1, 3, 7, 14),
                  labels = c("1h", "6h", "1d", "3d", "7d", "14d")) +
    scale_color_brewer(palette = "Set1", name = NULL) +
    facet_wrap(~ site, ncol = 2) +
    coord_cartesian(ylim = r_limit) +
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
  coord_flip(ylim = r_limit) +
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










# ============================================================
# INTERPOLATED SURFACE PLOT - ORDERED BY CORRELATION VALUE
# ============================================================

message("  Creating correlation surface...")

for (site_name in unique(cor_results$site)) {
  
  site_data <- cor_results %>%
    filter(site == site_name, !is.na(r))
  
  var_order <- site_data %>%
    group_by(variable, var_label, var_group) %>%
    summarize(
      # Get the r value with the largest absolute value (keeping sign)
      extreme_r = r[which.max(abs(r))],
      .groups = "drop"
    ) %>%
    arrange(desc(extreme_r)) %>%  # Most positive first, most negative last
    mutate(var_num = row_number())
  
  site_data <- site_data %>%
    left_join(var_order %>% dplyr::select(variable, var_num, extreme_r), by = "variable")
  
  # Surface plot ordered by correlation value
  p_surface_rows <- ggplot(site_data, aes(x = window_days * 24, y = 1)) +
    geom_raster(aes(fill = r), interpolate = TRUE) +
    scale_fill_gradientn(
      colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                 "#F7F7F7",
                 "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
      values = scales::rescale(signed_sqrt(seq(-1, 1, length.out = 11))),
      limits = r_limit,
      name = "r"
    ) +
    scale_x_continuous(
      breaks = c(1, 6, 12, 24, 48, 72, 120, 168, 240, 336),
      labels = c("1h", "6h", "12h", "1d", "2d", "3d", "5d", "7d", "10d", "14d"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~ reorder(var_label, var_num), ncol = 1, strip.position = "left") +
    labs(
      title = paste("CH4 Stem Flux Correlation Surface -", site_name),
      subtitle = "Ordered by max r (positive to negative), hourly resolution 1h to 14 days",
      x = "Integration Window",
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "right"
    )
  
  n_vars <- length(unique(site_data$variable))
  surface_height <- max(8, n_vars * 0.4 + 2)
  
  ggsave(file.path(OUTPUT_DIR, paste0("surface_rows_ranked_", site_name, ".png")),
         p_surface_rows, width = 12, height = surface_height, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, paste0("surface_rows_ranked_", site_name, ".pdf")),
         p_surface_rows, width = 12, height = surface_height)
  
  message("  Saved: surface_rows_ranked_", site_name, ".png/pdf")
}










# ============================================================
# INTERPOLATED SURFACE PLOT - ORDERED BY CORRELATION VALUE
# ============================================================

message("  Creating correlation surface...")

for (site_name in unique(cor_results$site)) {
  
  site_data <- cor_results %>%
    filter(site == site_name, !is.na(r))
  
  # Create variable order by most extreme correlation (keeping sign)
  var_order <- site_data %>%
    group_by(variable, var_label, var_group) %>%
    summarize(
      extreme_r = r[which.max(abs(r))],
      .groups = "drop"
    ) %>%
    arrange(desc(extreme_r)) %>%
    mutate(var_num = row_number())
  
  site_data <- site_data %>%
    left_join(var_order %>% dplyr::select(variable, var_num, extreme_r), by = "variable")
  
  # Find strongest significant correlation for each variable
  best_significant <- site_data %>%
    filter(significant) %>%
    group_by(variable, var_label, var_num) %>%
    slice_max(abs(r), n = 1) %>%
    ungroup() %>%
    mutate(x_pos = window_days * 24)
  
  # Surface plot ordered by correlation value
  p_surface_rows <- ggplot(site_data, aes(x = window_days * 24, y = 1)) +
    geom_raster(aes(fill = r), interpolate = TRUE) +
    # Add star at strongest significant correlation
    geom_point(data = best_significant, 
               aes(x = x_pos, y = 1),
               shape = 8, size = 3, color = "black", stroke = 1.2) +
    scale_fill_gradientn(
      colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
                 "#F7F7F7",
                 "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
      values = scales::rescale(signed_sqrt(seq(-1, 1, length.out = 11))),
      limits = r_limit,
      name = "r"
    ) +
    scale_x_continuous(
      breaks = c(1, 6, 12, 24, 48, 72, 120, 168, 240, 336),
      labels = c("1h", "6h", "12h", "1d", "2d", "3d", "5d", "7d", "10d", "14d"),
      expand = c(0, 0)
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_wrap(~ reorder(var_label, var_num), ncol = 1, strip.position = "left") +
    labs(
      title = paste("CH4 Stem Flux Correlation Surface -", site_name),
      subtitle = "Ordered by extreme r; ✳ = strongest significant correlation",
      x = "Integration Window",
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      panel.spacing = unit(0.1, "lines"),
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "right"
    )
  
  n_vars <- length(unique(site_data$variable))
  surface_height <- max(8, n_vars * 0.4 + 2)
  
  ggsave(file.path(OUTPUT_DIR, paste0("surface_rows_ranked_", site_name, ".png")),
         p_surface_rows, width = 12, height = surface_height, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, paste0("surface_rows_ranked_", site_name, ".pdf")),
         p_surface_rows, width = 12, height = surface_height)
  
  message("  Saved: surface_rows_ranked_", site_name, ".png/pdf")
}
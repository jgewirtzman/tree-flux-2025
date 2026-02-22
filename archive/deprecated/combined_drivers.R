# ============================================================
# 10_combined_flux_drivers_plot.R
# 
# Combined timeseries of CH4 stem flux with key drivers:
# - Water table (BVS)
# - Soil temperature (Fisher s10t)
#
# Uses the same flux plotting style as fig_temporal
# ============================================================

library(tidyverse)
library(lubridate)
library(patchwork)

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  flux = "data/HF_2023-2025_tree_flux.csv",  # UPDATE PATH AS NEEDED
  aligned = "data/processed/aligned_hourly_dataset.csv"
)

OUTPUT_DIR <- "figures"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD AND PROCESS FLUX DATA (same as your analysis script)
# ============================================================

message("Loading flux data...")

fluxes <- read.csv(PATHS$flux)

fluxes <- fluxes %>%
  # Fix units for pre-2025 data
  mutate(CH4_flux_nmolpm2ps = ifelse(year < 2025, CH4_flux_nmolpm2ps * 1000, CH4_flux_nmolpm2ps)) %>%
  # Fill missing species from same tree
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  # Filter bad data
  filter(CH4_flux_nmolpm2ps >= -1,
         !is.na(PLOT)) %>%
  # Add derived columns
  mutate(
    date = as.Date(datetime_posx),
    location = ifelse(PLOT == "BGS", "wetland", "upland")
  )

# Add sampling rounds (>7 day gap = new round)
fluxes <- fluxes %>%
  arrange(date) %>%
  mutate(
    days_since_last = as.numeric(date - lag(date), units = "days"),
    days_since_last = replace_na(days_since_last, 0),
    new_round = days_since_last > 4,
    sampling_round = cumsum(new_round) + 1
  )

message("  Flux data: ", nrow(fluxes), " observations")
message("  Date range: ", min(fluxes$date), " to ", max(fluxes$date))

# Create temporal_round summary (matching fig_temporal style)
temporal_round <- fluxes %>%
  group_by(sampling_round, location) %>%
  summarise(
    date = mean(date),
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# ============================================================
# LOAD ENVIRONMENTAL DATA
# ============================================================

message("Loading aligned environmental data...")

aligned_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC"))

# Daily environmental data
env_daily <- aligned_data %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarize(
    bvs_wtd_cm = mean(bvs_wtd_cm, na.rm = TRUE),
    s10t = mean(s10t, na.rm = TRUE),
    NEON_SWC_shallow = mean(NEON_SWC_shallow, na.rm = TRUE),
    NEON_SWC_mid = mean(NEON_SWC_mid, na.rm = TRUE),
    NEON_SWC_deep = mean(NEON_SWC_deep, na.rm = TRUE),
    .groups = "drop"
  )

# Pivot SWC for plotting
env_swc_long <- env_daily %>%
  dplyr::select(date, NEON_SWC_shallow, NEON_SWC_mid, NEON_SWC_deep) %>%
  pivot_longer(cols = starts_with("NEON_SWC"),
               names_to = "depth",
               values_to = "SWC") %>%
  mutate(depth = factor(depth, 
                        levels = c("NEON_SWC_shallow", "NEON_SWC_mid", "NEON_SWC_deep"),
                        labels = c("Shallow", "Mid", "Deep")))

message("  Env data range: ", min(env_daily$date), " to ", max(env_daily$date))

# ============================================================
# CREATE COMBINED PLOT
# ============================================================

message("Creating combined plot...")

# Get date range from flux data
DATE_MIN <- min(temporal_round$date)
DATE_MAX <- max(temporal_round$date)

message("  Flux date range: ", DATE_MIN, " to ", DATE_MAX)

# Panel 1: CH4 flux (fig_temporal style)
fig_temporal <- ggplot(temporal_round, aes(x = date, y = mean_CH4, color = location)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 5) +
  scale_color_manual(values = c("upland" = "#C2703D", "wetland" = "#3D7C9C"),
                     labels = c("Upland", "Wetland")) +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       color = "Location") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "top"
  )

# Panel 2: Water table depth (BVS)
p_wtd <- ggplot(env_daily, aes(x = date, y = bvs_wtd_cm)) +
  geom_line(color = "#4682B4", linewidth = 0.6) +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(y = "Water table\nBVS (cm)") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

# Panel 3: Soil temperature (Fisher s10t)
p_ts <- ggplot(env_daily, aes(x = date, y = s10t)) +
  geom_line(color = "#B2182B", linewidth = 0.6) +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(y = "Soil temp\n(°C)") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

# Panel 4: NEON soil moisture (all 3 depths)
p_swc <- ggplot(env_swc_long, aes(x = date, y = SWC, color = depth)) +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("Shallow" = "#90EE90", "Mid" = "#228B22", "Deep" = "#006400"),
                     name = "Depth") +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(x = "Date",
       y = "Soil moisture\nNEON (%)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Combine panels
p_combined <- fig_temporal / p_wtd / p_ts / p_swc +
  plot_layout(heights = c(2, 1, 1, 1))

# Save
ggsave(file.path(OUTPUT_DIR, "combined_flux_wtd_ts.png"),
       p_combined, width = 12, height = 8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "combined_flux_wtd_ts.pdf"),
       p_combined, width = 12, height = 8)

message("Saved: combined_flux_wtd_ts.png/pdf")
message("\nDone!")


# ============================================================
# 10_combined_flux_drivers_plot.R
# 
# Combined timeseries of CH4 stem flux with key drivers:
# - Water table (BVS)
# - Soil temperature (Fisher s10t)
#
# Uses the same flux plotting style as fig_temporal
# ============================================================

library(tidyverse)
library(lubridate)
library(patchwork)

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  flux = "data/HF_2023-2025_tree_flux.csv",  # UPDATE PATH AS NEEDED
  aligned = "data/processed/aligned_hourly_dataset.csv"
)

OUTPUT_DIR <- "figures"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD AND PROCESS FLUX DATA (same as your analysis script)
# ============================================================

message("Loading flux data...")

fluxes <- read.csv(PATHS$flux)

fluxes <- fluxes %>%
  # Fix units for pre-2025 data
  mutate(CH4_flux_nmolpm2ps = ifelse(year < 2025, CH4_flux_nmolpm2ps * 1000, CH4_flux_nmolpm2ps)) %>%
  # Fill missing species from same tree
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  # Filter bad data
  filter(CH4_flux_nmolpm2ps >= -1,
         !is.na(PLOT)) %>%
  # Add derived columns
  mutate(
    date = as.Date(datetime_posx),
    location = ifelse(PLOT == "BGS", "wetland", "upland")
  )

# Add sampling rounds (>7 day gap = new round)
fluxes <- fluxes %>%
  arrange(date) %>%
  mutate(
    days_since_last = as.numeric(date - lag(date), units = "days"),
    days_since_last = replace_na(days_since_last, 0),
    new_round = days_since_last > 7,
    sampling_round = cumsum(new_round) + 1
  )

message("  Flux data: ", nrow(fluxes), " observations")
message("  Date range: ", min(fluxes$date), " to ", max(fluxes$date))

# Create temporal_round summary (matching fig_temporal style)
temporal_round <- fluxes %>%
  group_by(sampling_round, location) %>%
  summarise(
    date = mean(date),
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# ============================================================
# LOAD ENVIRONMENTAL DATA
# ============================================================

message("Loading aligned environmental data...")

aligned_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC"))

# Daily environmental data
env_daily <- aligned_data %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarize(
    bvs_wtd_cm = mean(bvs_wtd_cm, na.rm = TRUE),
    s10t = mean(s10t, na.rm = TRUE),
    NEON_SWC_shallow = mean(NEON_SWC_shallow, na.rm = TRUE),
    NEON_SWC_mid = mean(NEON_SWC_mid, na.rm = TRUE),
    NEON_SWC_deep = mean(NEON_SWC_deep, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(date) %>%
  # Convert to z-scores
  mutate(
    wtd_z = (bvs_wtd_cm - mean(bvs_wtd_cm, na.rm = TRUE)) / sd(bvs_wtd_cm, na.rm = TRUE),
    ts_z = (s10t - mean(s10t, na.rm = TRUE)) / sd(s10t, na.rm = TRUE),
    swc_shallow_z = (NEON_SWC_shallow - mean(NEON_SWC_shallow, na.rm = TRUE)) / sd(NEON_SWC_shallow, na.rm = TRUE),
    swc_mid_z = (NEON_SWC_mid - mean(NEON_SWC_mid, na.rm = TRUE)) / sd(NEON_SWC_mid, na.rm = TRUE),
    swc_deep_z = (NEON_SWC_deep - mean(NEON_SWC_deep, na.rm = TRUE)) / sd(NEON_SWC_deep, na.rm = TRUE)
  ) %>%
  # Rolling 5-day mean of z-scores
  mutate(
    wtd_z_5d = zoo::rollmean(wtd_z, k = 5, fill = NA, align = "right"),
    ts_z_5d = zoo::rollmean(ts_z, k = 5, fill = NA, align = "right")
  ) %>%
  # Classify conditions based on rolling 5-day mean
  mutate(
    condition = case_when(
      wtd_z_5d > 0.5 & ts_z_5d > 0.5 ~ "high_both",
      wtd_z_5d < -0.5 & ts_z_5d < -0.5 ~ "low_both",
      TRUE ~ "mixed"
    )
  )

# Create shading regions (find contiguous periods)
env_daily <- env_daily %>%
  mutate(
    condition_change = condition != lag(condition, default = first(condition)),
    condition_group = cumsum(condition_change)
  )

# Summarize shading regions
shading_regions <- env_daily %>%
  filter(condition != "mixed") %>%
  group_by(condition_group, condition) %>%
  summarize(
    xmin = min(date),
    xmax = max(date),
    duration = as.numeric(xmax - xmin) + 1,
    .groups = "drop"
  )

# Pivot SWC z-scores for plotting
env_swc_long <- env_daily %>%
  dplyr::select(date, swc_shallow_z, swc_mid_z, swc_deep_z) %>%
  pivot_longer(cols = starts_with("swc_"),
               names_to = "depth",
               values_to = "SWC_z") %>%
  mutate(depth = factor(depth, 
                        levels = c("swc_shallow_z", "swc_mid_z", "swc_deep_z"),
                        labels = c("Shallow", "Mid", "Deep")))

message("  Env data range: ", min(env_daily$date), " to ", max(env_daily$date))

# ============================================================
# CREATE COMBINED PLOT
# ============================================================

message("Creating combined plot...")

# Get date range from flux data
DATE_MIN <- min(temporal_round$date)
DATE_MAX <- max(temporal_round$date)

message("  Flux date range: ", DATE_MIN, " to ", DATE_MAX)

# Filter shading to flux date range
shading_regions <- shading_regions %>%
  filter(xmax >= DATE_MIN, xmin <= DATE_MAX) %>%
  mutate(
    xmin = pmax(xmin, DATE_MIN),
    xmax = pmin(xmax, DATE_MAX)
  )

# Panel 1: CH4 flux (fig_temporal style)
fig_temporal <- ggplot(temporal_round, aes(x = date, y = mean_CH4, color = location)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = condition),
            inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_manual(values = c("high_both" = "#FF6B6B", "low_both" = "#4ECDC4"),
                    labels = c("high_both" = "Warm & Wet", "low_both" = "Cold & Dry"),
                    name = "Conditions") +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 5) +
  scale_color_manual(values = c("upland" = "#D2691E", "wetland" = "#4682B4"),
                     labels = c("Upland", "Wetland")) +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       color = "Location") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "top"
  ) +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))

# Panel 2: Water table depth (z-score)
p_wtd <- ggplot(env_daily, aes(x = date, y = wtd_z)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = condition),
            inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_manual(values = c("high_both" = "#FF6B6B", "low_both" = "#4ECDC4"), guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(color = "#4682B4", linewidth = 0.6) +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(y = "Water table\n(z-score)") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

# Panel 3: Soil temperature (z-score)
p_ts <- ggplot(env_daily, aes(x = date, y = ts_z)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = condition),
            inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_manual(values = c("high_both" = "#FF6B6B", "low_both" = "#4ECDC4"), guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(color = "#B2182B", linewidth = 0.6) +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(y = "Soil temp\n(z-score)") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

# Panel 4: NEON soil moisture (all 3 depths, z-scores)
p_swc <- ggplot(env_swc_long, aes(x = date, y = SWC_z, color = depth)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = condition),
            inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_manual(values = c("high_both" = "#FF6B6B", "low_both" = "#4ECDC4"), guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("Shallow" = "#90EE90", "Mid" = "#228B22", "Deep" = "#006400"),
                     name = "Depth") +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(x = "Date",
       y = "Soil moisture\n(z-score)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Combine panels
p_combined <- fig_temporal / p_wtd / p_ts / p_swc +
  plot_layout(heights = c(2, 1, 1, 1))

# Save
ggsave(file.path(OUTPUT_DIR, "combined_flux_wtd_ts.png"),
       p_combined, width = 12, height = 8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "combined_flux_wtd_ts.pdf"),
       p_combined, width = 12, height = 8)

message("Saved: combined_flux_wtd_ts.png/pdf")
message("\nDone!")



# ============================================================
# 10_combined_flux_drivers_plot.R
# 
# Combined timeseries of CH4 stem flux with key drivers:
# - Water table (BVS)
# - Soil temperature (Fisher s10t)
#
# Uses the same flux plotting style as fig_temporal
# ============================================================

library(tidyverse)
library(lubridate)
library(patchwork)

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  flux = "data/HF_2023-2025_tree_flux.csv",  # UPDATE PATH AS NEEDED
  aligned = "data/processed/aligned_hourly_dataset.csv"
)

OUTPUT_DIR <- "figures"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD AND PROCESS FLUX DATA (same as your analysis script)
# ============================================================

message("Loading flux data...")

fluxes <- read.csv(PATHS$flux)

fluxes <- fluxes %>%
  # Fix units for pre-2025 data
  mutate(CH4_flux_nmolpm2ps = ifelse(year < 2025, CH4_flux_nmolpm2ps * 1000, CH4_flux_nmolpm2ps)) %>%
  # Fill missing species from same tree
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  # Filter bad data
  filter(CH4_flux_nmolpm2ps >= -1,
         !is.na(PLOT)) %>%
  # Add derived columns
  mutate(
    date = as.Date(datetime_posx),
    location = ifelse(PLOT == "BGS", "wetland", "upland")
  )

# Add sampling rounds (>7 day gap = new round)
fluxes <- fluxes %>%
  arrange(date) %>%
  mutate(
    days_since_last = as.numeric(date - lag(date), units = "days"),
    days_since_last = replace_na(days_since_last, 0),
    new_round = days_since_last > 7,
    sampling_round = cumsum(new_round) + 1
  )

message("  Flux data: ", nrow(fluxes), " observations")
message("  Date range: ", min(fluxes$date), " to ", max(fluxes$date))

# Create temporal_round summary (matching fig_temporal style)
temporal_round <- fluxes %>%
  group_by(sampling_round, location) %>%
  summarise(
    date = mean(date),
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# ============================================================
# LOAD ENVIRONMENTAL DATA
# ============================================================

message("Loading aligned environmental data...")

aligned_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC"))

# Daily environmental data
env_daily <- aligned_data %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarize(
    bvs_wtd_cm = mean(bvs_wtd_cm, na.rm = TRUE),
    s10t = mean(s10t, na.rm = TRUE),
    NEON_SWC_shallow = mean(NEON_SWC_shallow, na.rm = TRUE),
    NEON_SWC_mid = mean(NEON_SWC_mid, na.rm = TRUE),
    NEON_SWC_deep = mean(NEON_SWC_deep, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(date) %>%
  # Convert to z-scores
  mutate(
    wtd_z = (bvs_wtd_cm - mean(bvs_wtd_cm, na.rm = TRUE)) / sd(bvs_wtd_cm, na.rm = TRUE),
    ts_z = (s10t - mean(s10t, na.rm = TRUE)) / sd(s10t, na.rm = TRUE),
    swc_shallow_z = (NEON_SWC_shallow - mean(NEON_SWC_shallow, na.rm = TRUE)) / sd(NEON_SWC_shallow, na.rm = TRUE),
    swc_mid_z = (NEON_SWC_mid - mean(NEON_SWC_mid, na.rm = TRUE)) / sd(NEON_SWC_mid, na.rm = TRUE),
    swc_deep_z = (NEON_SWC_deep - mean(NEON_SWC_deep, na.rm = TRUE)) / sd(NEON_SWC_deep, na.rm = TRUE)
  ) %>%
  # Classify conditions (both > 0 = "high", both < 0 = "low")
  mutate(
    condition = case_when(
      wtd_z > 0 & ts_z > 0 ~ "high_both",
      wtd_z < 0 & ts_z < 0 ~ "low_both",
      TRUE ~ "mixed"
    )
  )

# Create shading regions (find contiguous periods)
env_daily <- env_daily %>%
  mutate(
    condition_change = condition != lag(condition, default = first(condition)),
    condition_group = cumsum(condition_change)
  )

# Summarize shading regions
shading_regions <- env_daily %>%
  filter(condition != "mixed") %>%
  group_by(condition_group, condition) %>%
  summarize(
    xmin = min(date),
    xmax = max(date),
    .groups = "drop"
  )

# Pivot SWC z-scores for plotting
env_swc_long <- env_daily %>%
  dplyr::select(date, swc_shallow_z, swc_mid_z, swc_deep_z) %>%
  pivot_longer(cols = starts_with("swc_"),
               names_to = "depth",
               values_to = "SWC_z") %>%
  mutate(depth = factor(depth, 
                        levels = c("swc_shallow_z", "swc_mid_z", "swc_deep_z"),
                        labels = c("Shallow", "Mid", "Deep")))

message("  Env data range: ", min(env_daily$date), " to ", max(env_daily$date))

# ============================================================
# CREATE COMBINED PLOT
# ============================================================

message("Creating combined plot...")

# Get date range from flux data
DATE_MIN <- min(temporal_round$date)
DATE_MAX <- max(temporal_round$date)

message("  Flux date range: ", DATE_MIN, " to ", DATE_MAX)

# Filter shading to flux date range
shading_regions <- shading_regions %>%
  filter(xmax >= DATE_MIN, xmin <= DATE_MAX) %>%
  mutate(
    xmin = pmax(xmin, DATE_MIN),
    xmax = pmin(xmax, DATE_MAX)
  )

# Panel 1: CH4 flux (fig_temporal style)
fig_temporal <- ggplot(temporal_round, aes(x = date, y = mean_CH4, color = location)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = condition),
            inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_manual(values = c("high_both" = "#FF6B6B", "low_both" = "#4ECDC4"),
                    labels = c("high_both" = "Warm & Wet", "low_both" = "Cold & Dry"),
                    name = "Conditions") +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 5) +
  scale_color_manual(values = c("upland" = "#D2691E", "wetland" = "#4682B4"),
                     labels = c("Upland", "Wetland")) +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       color = "Location") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "top"
  ) +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))

# Panel 2: Water table depth (z-score)
p_wtd <- ggplot(env_daily, aes(x = date, y = wtd_z)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = condition),
            inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_manual(values = c("high_both" = "#FF6B6B", "low_both" = "#4ECDC4"), guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(color = "#4682B4", linewidth = 0.6) +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(y = "Water table\n(z-score)") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

# Panel 3: Soil temperature (z-score)
p_ts <- ggplot(env_daily, aes(x = date, y = ts_z)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = condition),
            inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_manual(values = c("high_both" = "#FF6B6B", "low_both" = "#4ECDC4"), guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(color = "#B2182B", linewidth = 0.6) +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(y = "Soil temp\n(z-score)") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

# Panel 4: NEON soil moisture (all 3 depths, z-scores)
p_swc <- ggplot(env_swc_long, aes(x = date, y = SWC_z, color = depth)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = condition),
            inherit.aes = FALSE, alpha = 0.3) +
  scale_fill_manual(values = c("high_both" = "#FF6B6B", "low_both" = "#4ECDC4"), guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(linewidth = 0.6) +
  scale_color_manual(values = c("Shallow" = "#90EE90", "Mid" = "#228B22", "Deep" = "#006400"),
                     name = "Depth") +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX)) +
  labs(x = "Date",
       y = "Soil moisture\n(z-score)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Combine panels
p_combined <- fig_temporal / p_wtd / p_ts / p_swc +
  plot_layout(heights = c(2, 1, 1, 1))

# Save
ggsave(file.path(OUTPUT_DIR, "combined_flux_wtd_ts.png"),
       p_combined, width = 12, height = 8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "combined_flux_wtd_ts.pdf"),
       p_combined, width = 12, height = 8)

message("Saved: combined_flux_wtd_ts.png/pdf")
message("\nDone!")





# ============================================================
# 10_combined_flux_drivers_plot.R
# 
# Combined timeseries of CH4 stem flux with key drivers:
# - Water table (BVS)
# - Soil temperature (Fisher s10t)
#
# Uses the same flux plotting style as fig_temporal
# ============================================================

library(tidyverse)
library(lubridate)
library(patchwork)

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  flux = "data/HF_2023-2025_tree_flux.csv",  # UPDATE PATH AS NEEDED
  aligned = "data/processed/aligned_hourly_dataset.csv"
)

OUTPUT_DIR <- "figures"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD AND PROCESS FLUX DATA (same as your analysis script)
# ============================================================

message("Loading flux data...")

fluxes <- read.csv(PATHS$flux)

fluxes <- fluxes %>%
  # Fix units for pre-2025 data
  mutate(CH4_flux_nmolpm2ps = ifelse(year < 2025, CH4_flux_nmolpm2ps * 1000, CH4_flux_nmolpm2ps)) %>%
  # Fill missing species from same tree
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  # Filter bad data
  filter(CH4_flux_nmolpm2ps >= -1,
         !is.na(PLOT)) %>%
  # Add derived columns
  mutate(
    date = as.Date(datetime_posx),
    location = ifelse(PLOT == "BGS", "wetland", "upland")
  )

# Add sampling rounds (>7 day gap = new round)
fluxes <- fluxes %>%
  arrange(date) %>%
  mutate(
    days_since_last = as.numeric(date - lag(date), units = "days"),
    days_since_last = replace_na(days_since_last, 0),
    new_round = days_since_last > 7,
    sampling_round = cumsum(new_round) + 1
  )

message("  Flux data: ", nrow(fluxes), " observations")
message("  Date range: ", min(fluxes$date), " to ", max(fluxes$date))

# Create temporal_round summary (matching fig_temporal style)
temporal_round <- fluxes %>%
  group_by(sampling_round, location) %>%
  summarise(
    date = mean(date),
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

# ============================================================
# LOAD ENVIRONMENTAL DATA
# ============================================================

message("Loading aligned environmental data...")

aligned_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC"))

# Daily environmental data
env_daily <- aligned_data %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarize(
    bvs_wtd_cm = mean(bvs_wtd_cm, na.rm = TRUE),
    s10t = mean(s10t, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(date) %>%
  # Convert to z-scores
  mutate(
    wtd_z = (bvs_wtd_cm - mean(bvs_wtd_cm, na.rm = TRUE)) / sd(bvs_wtd_cm, na.rm = TRUE),
    ts_z = (s10t - mean(s10t, na.rm = TRUE)) / sd(s10t, na.rm = TRUE)
  ) %>%
  # 7-day antecedent rolling mean of z-scores
  mutate(
    wtd_z_7d = zoo::rollmean(wtd_z, k = 7, fill = NA, align = "right"),
    ts_z_7d = zoo::rollmean(ts_z, k = 7, fill = NA, align = "right")
  ) %>%
  # Classify conditions based on 7-day antecedent (only highlight when both are high)
  mutate(
    condition = case_when(
      wtd_z_7d > 0.5 & ts_z_7d > 0.5 ~ "high_both",
      TRUE ~ "mixed"
    )
  )

# Create shading regions (find contiguous periods)
env_daily <- env_daily %>%
  mutate(
    condition_change = condition != lag(condition, default = first(condition)),
    condition_group = cumsum(condition_change)
  )

# Summarize shading regions
shading_regions <- env_daily %>%
  filter(condition != "mixed") %>%
  group_by(condition_group, condition) %>%
  summarize(
    xmin = min(date),
    xmax = max(date),
    .groups = "drop"
  )

message("  Env data range: ", min(env_daily$date), " to ", max(env_daily$date))

# ============================================================
# CREATE COMBINED PLOT
# ============================================================

message("Creating combined plot...")

# Get date range from flux data
DATE_MIN <- min(temporal_round$date)
DATE_MAX <- max(temporal_round$date)

message("  Flux date range: ", DATE_MIN, " to ", DATE_MAX)

# Filter shading to flux date range
shading_regions <- shading_regions %>%
  filter(xmax >= DATE_MIN, xmin <= DATE_MAX) %>%
  mutate(
    xmin = pmax(xmin, DATE_MIN),
    xmax = pmin(xmax, DATE_MAX)
  )

# Panel 1: CH4 flux (fig_temporal style)
fig_temporal <- ggplot(temporal_round, aes(x = date, y = mean_CH4, color = location)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "gray80", alpha = 0.5) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 5) +
  scale_color_manual(values = c("upland" = "#C2703D", "wetland" = "#3D7C9C"),
                     labels = c("Upland", "Wetland")) +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX),
               date_breaks = "3 months", date_labels = "%b %Y") +
  labs(y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       color = "Location") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "top"
  )

# Panel 2: Water table depth (actual values, colored by z-score)
p_wtd <- ggplot(env_daily, aes(x = date, y = bvs_wtd_cm)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "gray80", alpha = 0.5) +
  geom_segment(aes(xend = lead(date), yend = lead(bvs_wtd_cm), color = wtd_z), 
               linewidth = 0.8) +
  scale_color_gradientn(colors = c("#2166AC", "gray80", "#B2182B"),
                        values = scales::rescale(c(min(env_daily$wtd_z, na.rm = TRUE), 
                                                   0, 
                                                   max(env_daily$wtd_z, na.rm = TRUE))),
                        guide = "none") +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX),
               date_breaks = "3 months", date_labels = "%b %Y") +
  labs(y = "Water table\nBVS (cm)") +
  theme_classic(base_size = 12) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank()
  )

# Panel 3: Soil temperature (actual values, colored by z-score)
p_ts <- ggplot(env_daily, aes(x = date, y = s10t)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "gray80", alpha = 0.5) +
  geom_segment(aes(xend = lead(date), yend = lead(s10t), color = ts_z), 
               linewidth = 0.8) +
  scale_color_gradientn(colors = c("#2166AC", "gray80", "#B2182B"),
                        values = scales::rescale(c(min(env_daily$ts_z, na.rm = TRUE), 
                                                   0, 
                                                   max(env_daily$ts_z, na.rm = TRUE))),
                        guide = "none") +
  scale_x_date(limits = c(DATE_MIN, DATE_MAX),
               date_breaks = "3 months", date_labels = "%b %Y") +
  labs(x = "Date",
       y = "Soil temp\n(°C)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Combine panels (equal heights)
p_combined <- fig_temporal / p_wtd / p_ts +
  plot_layout(heights = c(2, 1, 1))

# Save
ggsave(file.path(OUTPUT_DIR, "combined_flux_wtd_ts.png"),
       p_combined, width = 12, height = 8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "combined_flux_wtd_ts.pdf"),
       p_combined, width = 12, height = 8)

message("Saved: combined_flux_wtd_ts.png/pdf")
message("\nDone!")

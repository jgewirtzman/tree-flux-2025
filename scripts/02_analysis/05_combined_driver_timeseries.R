# ============================================================
# combined_drivers.R
# 
# Combined timeseries of CH4 stem flux with key environmental drivers:
# - Water table depth (BVS)
# - Soil temperature (Fisher s10t)
#
# Shaded regions indicate periods when both drivers are elevated
# (7-day rolling z-score > 0.5 for both).
#
# Run AFTER: timeseries.R (which creates the corrected flux file)
#
# Outputs:
#   - combined_flux_drivers.png/pdf
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(patchwork)
  library(zoo)
})

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  flux = "data/processed/flux_with_quality_flags.csv",
  aligned = "data/processed/aligned_hourly_dataset.csv"
)

OUTPUT_DIR <- "outputs/figures"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD AND PROCESS FLUX DATA
# ============================================================

message("Loading flux data...")

fluxes <- read_csv(PATHS$flux, show_col_types = FALSE) %>%
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  filter(!is.na(PLOT)) %>%
  mutate(
    date = as.Date(datetime_posx),
    location = factor(
      ifelse(PLOT == "BGS", "Wetland", "Upland"),
      levels = c("Wetland", "Upland")
    )
  )

# Add sampling rounds (>4 day gap = new round)
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

# Create temporal_round summary
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

# Daily environmental data with z-scores and rolling means
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
    wtd_z_7d = rollmean(wtd_z, k = 7, fill = NA, align = "right"),
    ts_z_7d = rollmean(ts_z, k = 7, fill = NA, align = "right")
  ) %>%
  # Classify conditions (highlight when both are elevated)
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

message("  Plot date range: ", DATE_MIN, " to ", DATE_MAX)

# Filter shading to flux date range
shading_regions <- shading_regions %>%
  filter(xmax >= DATE_MIN, xmin <= DATE_MAX) %>%
  mutate(
    xmin = pmax(xmin, DATE_MIN),
    xmax = pmin(xmax, DATE_MAX)
  )

# Panel 1: CH4 flux
p_flux <- ggplot(temporal_round, aes(x = date, y = mean_CH4, color = location)) +
  geom_rect(data = shading_regions,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "gray80", alpha = 0.5) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 5) +
  scale_color_manual(values = c("Wetland" = "#2A7F7A", "Upland" = "#6E8B3D")) +
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

# Panel 2: Water table depth (colored by z-score)
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

# Panel 3: Soil temperature (colored by z-score)
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

# Combine panels
p_combined <- p_flux / p_wtd / p_ts +
  plot_layout(heights = c(2, 1, 1))

print(p_combined)

# ============================================================
# SAVE
# ============================================================

ggsave(file.path(OUTPUT_DIR, "combined_flux_drivers.png"),
       p_combined, width = 12, height = 8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "combined_flux_drivers.pdf"),
       p_combined, width = 12, height = 8)

message("\nSaved: combined_flux_drivers.png/pdf")
message("Done!")
# ============================================================
# Temporal flux plot by location with species colors
# ============================================================

library(tidyverse)

# ============================================================
# CONFIG
# ============================================================

PATHS <- list(
  flux = "data/processed/flux_with_quality_flags.csv"
)

OUTPUT_DIR <- "outputs/figures"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD AND PROCESS DATA
# ============================================================

fluxes <- read_csv(PATHS$flux, show_col_types = FALSE) %>%
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  filter(!is.na(PLOT)) %>%
  mutate(
    date = as.Date(datetime_posx),
    location = ifelse(PLOT == "BGS", "Wetland", "Upland")
  )

# Define sampling rounds (>7 day gap = new round)
fluxes <- fluxes %>%
  arrange(date) %>%
  mutate(
    days_since_last = as.numeric(date - lag(date), units = "days"),
    days_since_last = replace_na(days_since_last, 0),
    new_round = days_since_last > 7,
    sampling_round = cumsum(new_round) + 1
  )

# Species labels
species_labels <- c(
  "bg" = "N. sylvatica",
  "rm" = "A. rubrum", 
  "ro" = "Q. rubra",
  "hem" = "T. canadensis"
)

fluxes <- fluxes %>%
  mutate(species_label = species_labels[SPECIES])

# Calculate round means per species and location
round_means <- fluxes %>%
  group_by(sampling_round, location, species_label) %>%
  summarize(
    date = mean(date),
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    .groups = "drop"
  )

# ============================================================
# PLOT
# ============================================================

# Color palette (teal + olive, ecological)
# Nyssa = dark teal (dominant wetland), Quercus = dark olive (upland oak)
# Acer = light teal, Tsuga = light olive
species_colors <- c(
  "N. sylvatica" = "#2A7F7A",
  "Q. rubra" = "#6E8B3D",
  "A. rubrum" = "#A7DAD1",
  "T. canadensis" = "#C9D6A4"
)

# Custom asinh transform for better outlier handling
asinh_trans <- function() {
  scales::trans_new(
    name = "asinh",
    transform = asinh,
    inverse = sinh
  )
}

p <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "gray40") +
  geom_point(
    data = fluxes,
    aes(x = date, y = CH4_flux_nmolpm2ps, fill = species_label,
        shape = ifelse(CH4_below_MDF == TRUE, "Below MDF", "Above MDF")),
    position = position_jitter(width = 2, height = 0),
    size = 1.5, alpha = 0.4, color = "gray30", stroke = 0.2
  ) +
  scale_shape_manual(
    values = c("Above MDF" = 21, "Below MDF" = 1),
    name = NULL,
    guide = guide_legend(override.aes = list(alpha = 0.8, size = 2))
  ) +
  geom_line(
    data = round_means,
    aes(x = date, y = mean_CH4, color = species_label),
    linewidth = 0.8, alpha = 0.8
  ) +
  geom_point(
    data = round_means,
    aes(x = date, y = mean_CH4, fill = species_label),
    size = 2.5, shape = 21, color = "gray30", stroke = 0.3
  ) +
  scale_y_continuous(
    trans = asinh_trans(),
    breaks = c(-1, 0, 1, 2, 5, 10, 20, 50, 100)
  ) +
  scale_fill_manual(values = species_colors) +
  scale_color_manual(values = species_colors) +
  scale_x_date(date_breaks = "3 months", date_labels = "%b '%y") +
  facet_wrap(~ location, ncol = 1, scales = "free_y") +
  geom_blank(data = data.frame(location = "Upland", y = -1), aes(y = y)) +
  labs(
    x = NULL,
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    fill = NULL,
    color = NULL
  ) +
  guides(color = "none") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "top",
    legend.text = element_text(face = "italic", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines"),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# Print and save
print(p)

ggsave(file.path(OUTPUT_DIR, "flux_temporal_species.png"), p, 
       width = 8, height = 8, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "flux_temporal_species.pdf"), p, 
       width = 8, height = 8)

message("Saved: flux_temporal_species.png/pdf")






# ============================================================
# TEMPORAL SUMMARY STATISTICS FOR RESULTS
# ============================================================

library(lubridate)

message("\n", paste(rep("=", 60), collapse = ""))
message("TEMPORAL PATTERNS FOR RESULTS")
message(paste(rep("=", 60), collapse = ""))

# Add temporal variables
fluxes <- fluxes %>%
  mutate(
    year = year(date),
    month = month(date),
    month_name = month(date, label = TRUE)
  )

# Sampling overview
message("\n--- Sampling Overview ---")
message("Date range: ", min(fluxes$date), " to ", max(fluxes$date))
message("Total observations: ", nrow(fluxes))
message("Sampling rounds: ", max(fluxes$sampling_round))

# Observations by site and species
message("\n--- Observations by Site × Species ---")
fluxes %>%
  group_by(location, species_label) %>%
  summarise(
    n_obs = n(),
    n_trees = n_distinct(Tree),
    .groups = "drop"
  ) %>%
  print()

# Flux summary by site (mean ± SE, median, IQR)
message("\n--- CH4 Flux by Site (nmol m-2 s-1) ---")
fluxes %>%
  group_by(location) %>%
  summarise(
    n = n(),
    mean = round(mean(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    se = round(sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()), 3),
    median = round(median(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    q25 = round(quantile(CH4_flux_nmolpm2ps, 0.25, na.rm = TRUE), 3),
    q75 = round(quantile(CH4_flux_nmolpm2ps, 0.75, na.rm = TRUE), 3),
    min = round(min(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    max = round(max(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

# Flux summary by species × site
message("\n--- CH4 Flux by Site × Species (nmol m-2 s-1) ---")
fluxes %>%
  group_by(location, species_label) %>%
  summarise(
    n = n(),
    mean = round(mean(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    se = round(sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()), 3),
    median = round(median(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    q25 = round(quantile(CH4_flux_nmolpm2ps, 0.25, na.rm = TRUE), 3),
    q75 = round(quantile(CH4_flux_nmolpm2ps, 0.75, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

# Negative fluxes (uptake)
message("\n--- Negative Flux (CH4 Uptake) ---")
fluxes %>%
  group_by(location) %>%
  summarise(
    n_negative = sum(CH4_flux_nmolpm2ps < 0),
    pct_negative = round(100 * mean(CH4_flux_nmolpm2ps < 0), 1),
    .groups = "drop"
  ) %>%
  print()

# Peak and minimum months by site
message("\n--- Peak and Minimum Months by Site ---")
monthly_site <- fluxes %>%
  group_by(month_name, location) %>%
  summarise(
    n = n(),
    mean = round(mean(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    se = round(sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()), 3),
    .groups = "drop"
  )

monthly_site %>%
  group_by(location) %>%
  summarise(
    peak_month = month_name[which.max(mean)],
    peak_mean = mean[which.max(mean)],
    peak_se = se[which.max(mean)],
    low_month = month_name[which.min(mean)],
    low_mean = mean[which.min(mean)],
    low_se = se[which.min(mean)],
    .groups = "drop"
  ) %>%
  print()

# Summer means by year × site (June-August)
message("\n--- Summer (Jun-Aug) Mean Flux by Year × Site ---")
fluxes %>%
  filter(month %in% c(6, 7, 8)) %>%
  group_by(year, location) %>%
  summarise(
    n = n(),
    mean = round(mean(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    se = round(sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()), 3),
    median = round(median(CH4_flux_nmolpm2ps, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(location, year) %>%
  print()

# Summer means by year × site × species
message("\n--- Summer (Jun-Aug) Mean Flux by Year × Site × Species ---")
fluxes %>%
  filter(month %in% c(6, 7, 8)) %>%
  group_by(year, location, species_label) %>%
  summarise(
    mean = round(mean(CH4_flux_nmolpm2ps, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = year, values_from = mean) %>%
  arrange(location, species_label) %>%
  print()
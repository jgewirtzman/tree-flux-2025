# ============================================================
# 12_temp_moisture_interaction.R
# 
# Response of CH4 flux to soil temperature at different
# soil moisture percentiles - showing temp × moisture interaction
# ============================================================

library(tidyverse)
library(lubridate)

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  flux = "data/HF_2023-2025_tree_flux.csv",
  aligned = "data/processed/aligned_hourly_dataset.csv"
)

OUTPUT_DIR <- "figures"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD AND PROCESS DATA
# ============================================================

message("Loading flux data...")

fluxes <- read.csv(PATHS$flux)

fluxes <- fluxes %>%
  mutate(CH4_flux_nmolpm2ps = ifelse(year < 2025, CH4_flux_nmolpm2ps * 1000, CH4_flux_nmolpm2ps)) %>%
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  filter(CH4_flux_nmolpm2ps >= -1,
         !is.na(PLOT)) %>%
  mutate(
    datetime = as.POSIXct(datetime_posx),
    date = as.Date(datetime_posx),
    location = ifelse(PLOT == "BGS", "wetland", "upland"),
    species_full = case_when(
      SPECIES == "bg"  ~ "Black gum",
      SPECIES == "hem" ~ "Eastern hemlock",
      SPECIES == "rm"  ~ "Red maple",
      SPECIES == "ro"  ~ "Red oak",
      TRUE ~ SPECIES
    )
  )

message("  Flux data: ", nrow(fluxes), " observations")

# Load aligned environmental data
message("Loading aligned environmental data...")
aligned_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC"))

# Join flux data with environmental data
fluxes <- fluxes %>%
  mutate(datetime_hour = floor_date(datetime, "hour"))

env_hourly <- aligned_data %>%
  dplyr::select(datetime, s10t, bvs_wtd_cm)

flux_env <- fluxes %>%
  left_join(env_hourly, by = c("datetime_hour" = "datetime")) %>%
  filter(!is.na(s10t), !is.na(bvs_wtd_cm), !is.na(CH4_flux_nmolpm2ps))

message("  Joined data: ", nrow(flux_env), " observations with complete env data")

# ============================================================
# FIT INTERACTION MODEL AND PREDICT
# ============================================================

# Define water table percentiles for prediction (10 bins)
wtd_percentiles <- seq(10, 90, by = 10)

# Function to fit model and generate predictions for a dataset
fit_interaction_model <- function(data, group_name) {
  
  if (nrow(data) < 30) {
    message("  Skipping ", group_name, " - too few observations")
    return(NULL)
  }
  
  # Fit model with interaction
  # Using log-linear model for exponential response
  # Add small constant to handle zeros
  data <- data %>%
    mutate(log_CH4 = log(CH4_flux_nmolpm2ps + 0.01))
  
  model <- lm(log_CH4 ~ s10t * bvs_wtd_cm, data = data)
  
  # Get water table percentile values from data
  wtd_values <- quantile(data$bvs_wtd_cm, probs = wtd_percentiles/100, na.rm = TRUE)
  
  # Generate predictions at each water table percentile
  # Only predict within temperature range observed for similar water table values
  pred_list <- lapply(seq_along(wtd_percentiles), function(i) {
    
    # Get temperature range for observations near this water table level
    # Define "near" as within +/- 0.5 percentile bands
    wtd_lo <- ifelse(i == 1, -Inf, (wtd_values[i] + wtd_values[i-1]) / 2)
    wtd_hi <- ifelse(i == length(wtd_percentiles), Inf, (wtd_values[i] + wtd_values[i+1]) / 2)
    
    nearby_data <- data %>%
      filter(bvs_wtd_cm >= wtd_lo, bvs_wtd_cm <= wtd_hi)
    
    # If not enough nearby data, use full temperature range but be conservative
    if (nrow(nearby_data) < 10) {
      temp_min <- quantile(data$s10t, 0.05, na.rm = TRUE)
      temp_max <- quantile(data$s10t, 0.95, na.rm = TRUE)
    } else {
      temp_min <- quantile(nearby_data$s10t, 0.05, na.rm = TRUE)
      temp_max <- quantile(nearby_data$s10t, 0.95, na.rm = TRUE)
    }
    
    temp_range <- seq(temp_min, temp_max, length.out = 100)
    
    pred_data <- data.frame(
      s10t = temp_range,
      bvs_wtd_cm = wtd_values[i]
    )
    
    pred <- predict(model, newdata = pred_data, se.fit = TRUE)
    
    pred_data %>%
      mutate(
        log_fit = pred$fit,
        log_se = pred$se.fit,
        fit = exp(log_fit) - 0.01,  # back-transform
        lwr = exp(log_fit - 1.96 * log_se) - 0.01,
        upr = exp(log_fit + 1.96 * log_se) - 0.01,
        wtd_percentile = wtd_percentiles[i],
        wtd_value = wtd_values[i],
        group = group_name
      )
  })
  
  predictions <- bind_rows(pred_list)
  
  return(list(
    model = model,
    predictions = predictions,
    data = data,
    wtd_quantiles = wtd_values
  ))
}

# ============================================================
# FIT BY SPECIES
# ============================================================

message("\nFitting interaction models by species...")

species_list <- c("Eastern hemlock", "Red maple", "Black gum", "Red oak")

results_species <- list()

for (sp in species_list) {
  message("  ", sp, "...")
  sp_data <- flux_env %>% filter(species_full == sp)
  results_species[[sp]] <- fit_interaction_model(sp_data, sp)
}

# Remove NULLs
results_species <- results_species[!sapply(results_species, is.null)]

# ============================================================
# FIT BY LOCATION
# ============================================================

message("\nFitting interaction models by location...")

results_location <- list()

for (loc in c("wetland", "upland")) {
  message("  ", loc, "...")
  loc_data <- flux_env %>% filter(location == loc)
  results_location[[loc]] <- fit_interaction_model(loc_data, loc)
}

results_location <- results_location[!sapply(results_location, is.null)]

# ============================================================
# CREATE SPECIES PLOT
# ============================================================

message("\nCreating species plot...")

# Combine predictions
pred_species <- bind_rows(lapply(results_species, function(x) x$predictions))

# Color palette - gradient for water table (low = brown/dry, high = blue/wet)
wtd_colors <- colorRampPalette(c("#8B4513", "#CD853F", "#DEB887", "#B0C4DE", "#6495ED", "#4169E1", "#0000CD"))(length(wtd_percentiles))
names(wtd_colors) <- as.character(wtd_percentiles)

# Ensure factor ordering
pred_species <- pred_species %>%
  mutate(
    wtd_percentile = factor(wtd_percentile, levels = wtd_percentiles),
    group = factor(group, levels = species_list)
  )

p_species <- ggplot(pred_species, aes(x = s10t, y = fit, color = wtd_percentile)) +
  geom_line(linewidth = 1) +
  facet_wrap(~ group, scales = "free_y", ncol = 4) +
  scale_color_manual(
    values = wtd_colors,
    labels = paste0(wtd_percentiles, "%"),
    name = "Water table\npercentile"
  ) +
  labs(
    x = "Soil temperature (°C)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}))
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    legend.position = "right"
  )

print(p_species)

ggsave(file.path(OUTPUT_DIR, "temp_wtd_interaction_species.png"),
       p_species, width = 14, height = 4, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "temp_wtd_interaction_species.pdf"),
       p_species, width = 14, height = 4)

message("Saved: temp_wtd_interaction_species.png/pdf")

# ============================================================
# CREATE LOCATION PLOT
# ============================================================

message("\nCreating location plot...")

pred_location <- bind_rows(lapply(results_location, function(x) x$predictions))

# Same color scheme for both locations
wtd_colors_loc <- colorRampPalette(c("#8B4513", "#CD853F", "#DEB887", "#B0C4DE", "#6495ED", "#4169E1", "#0000CD"))(length(wtd_percentiles))

pred_location <- pred_location %>%
  mutate(
    wtd_percentile = factor(wtd_percentile, levels = wtd_percentiles),
    group = factor(group, levels = c("wetland", "upland"))
  )

# Create separate color mapping
pred_location <- pred_location %>%
  mutate(
    color_group = paste(group, wtd_percentile, sep = "_")
  )

# Build color vector
all_colors <- c(
  setNames(wtd_colors_loc, paste("wetland", wtd_percentiles, sep = "_")),
  setNames(wtd_colors_loc, paste("upland", wtd_percentiles, sep = "_"))
)

p_location <- ggplot(pred_location, aes(x = s10t, y = fit, color = color_group, fill = color_group)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~ group, scales = "fixed", 
             labeller = labeller(group = c("wetland" = "Wetland", "upland" = "Upland"))) +
  scale_color_manual(values = all_colors, guide = "none") +
  scale_fill_manual(values = all_colors, guide = "none") +
  labs(
    x = "Soil temperature (°C)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}))
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  )

# Add a manual legend
# Create a dummy plot for legend
legend_data <- data.frame(
  x = rep(1:length(wtd_percentiles), 1),
  y = rep(1, length(wtd_percentiles)),
  percentile = factor(wtd_percentiles, levels = wtd_percentiles)
)

p_legend <- ggplot(legend_data, aes(x = x, y = y, fill = percentile)) +
  geom_tile() +
  scale_fill_manual(
    values = wtd_colors_loc,
    labels = paste0(wtd_percentiles, "%"),
    name = "Water table\npercentile"
  ) +
  theme_void() +
  theme(legend.position = "right")

# Extract legend
library(cowplot)
legend <- get_legend(p_legend)

# Combine
p_location_final <- plot_grid(p_location, legend, rel_widths = c(1, 0.2))

print(p_location_final)

ggsave(file.path(OUTPUT_DIR, "temp_wtd_interaction_location.png"),
       p_location_final, width = 10, height = 4, dpi = 300)

message("Saved: temp_wtd_interaction_location.png")

# ============================================================
# ALTERNATIVE: SINGLE PANEL WITH BETTER LEGEND
# ============================================================

# Wetland only (where the interaction is likely strongest)
message("\nCreating wetland-only plot...")

wetland_colors <- colorRampPalette(c("#8B4513", "#CD853F", "#DEB887", "#B0C4DE", "#6495ED", "#4169E1", "#0000CD"))(length(wtd_percentiles))

pred_wetland <- results_location[["wetland"]]$predictions %>%
  mutate(wtd_percentile = factor(wtd_percentile, levels = wtd_percentiles))

p_wetland <- ggplot(pred_wetland, aes(x = s10t, y = fit, color = wtd_percentile)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr, fill = wtd_percentile), 
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(
    values = wetland_colors,
    labels = paste0(wtd_percentiles, "%"),
    name = "Water table\npercentile"
  ) +
  scale_fill_manual(
    values = wetland_colors,
    guide = "none"
  ) +
  labs(
    x = "Soil temperature (°C)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Wetland: CH4 flux response to temperature × water table interaction"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right"
  )

print(p_wetland)

ggsave(file.path(OUTPUT_DIR, "temp_wtd_interaction_wetland.png"),
       p_wetland, width = 8, height = 5, dpi = 300)

message("Saved: temp_wtd_interaction_wetland.png")

message("\nDone!")
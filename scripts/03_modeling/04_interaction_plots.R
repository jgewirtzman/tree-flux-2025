# ============================================================
# 107_interaction_plot_by_site.R
# 
# Temperature × Water Table interaction plots by site
# Using CORE models (TS × WTD × species only) from:
#   - BGS (Wetland): 104_bgs_model_complete.R
#   - EMS (Upland): 105b_ems_model_bgs_drivers.R
#
# Both use: asinh(CH4*1000) ~ TS_Ha2 * bvs_wtd_cm * species + (1|Tree)
# ============================================================

library(tidyverse)
library(lme4)
library(cowplot)

# ============================================================
# LOAD MODELS AND DATA
# ============================================================

message("Loading models...")

m_bgs <- readRDS("outputs/models/bgs_final/m_core_asinh.rds")
m_ems <- readRDS("outputs/models/ems_bgs_drivers/m_core_asinh.rds")

# Extract model frames for scaling parameters
frame_bgs <- m_bgs@frame
frame_ems <- m_ems@frame

# Check model formulas
cat("BGS model formula:\n")
print(formula(m_bgs))
cat("\nEMS model formula:\n")
print(formula(m_ems))

# Get all predictor columns from each model frame
bgs_vars <- names(frame_bgs)
ems_vars <- names(frame_ems)

cat("\nBGS model variables:", paste(bgs_vars, collapse = ", "), "\n")
cat("EMS model variables:", paste(ems_vars, collapse = ", "), "\n")

# Identify the core predictor columns for each model
# BGS may have different predictors than EMS
ts_pred_bgs <- bgs_vars[grepl("^TS_Ha", bgs_vars) & grepl("_raw_", bgs_vars)][1]
wtd_pred_bgs <- bgs_vars[grepl("^bvs_wtd_cm_raw_", bgs_vars)][1]

ts_pred_ems <- ems_vars[grepl("^TS_Ha", ems_vars) & grepl("_raw_", ems_vars)][1]
wtd_pred_ems <- ems_vars[grepl("^bvs_wtd_cm_raw_", ems_vars)][1]

cat("\nBGS - Temperature predictor:", ts_pred_bgs, "\n")
cat("BGS - Water table predictor:", wtd_pred_bgs, "\n")
cat("EMS - Temperature predictor:", ts_pred_ems, "\n")
cat("EMS - Water table predictor:", wtd_pred_ems, "\n")

# Get all additional predictors in each model (for setting to mean values)
response_var <- "CH4_flux_asinh"
exclude_vars <- c(response_var, "Tree", "species")

additional_preds_bgs <- setdiff(bgs_vars, c(exclude_vars, ts_pred_bgs, wtd_pred_bgs))
additional_preds_ems <- setdiff(ems_vars, c(exclude_vars, ts_pred_ems, wtd_pred_ems))

cat("\nBGS additional predictors:", paste(additional_preds_bgs, collapse = ", "), "\n")
cat("EMS additional predictors:", paste(additional_preds_ems, collapse = ", "), "\n")

# ============================================================
# LOAD UNSCALED DATA FOR AXIS LABELS
# ============================================================

# We need the original unscaled data to convert predictions back
# Load flux and environmental data

PATHS <- list(
  flux = "data/input/HF_2023-2025_tree_flux_corrected.csv",
  aligned = "data/processed/aligned_hourly_dataset.csv"
)

message("Loading raw data for scaling parameters...")

aligned_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC"))

flux_data <- read_csv(PATHS$flux, show_col_types = FALSE)

# Build the same features as in the model scripts
library(RcppRoll)

roll_mean <- function(x, n) {
  RcppRoll::roll_mean(x, n = n, align = "right", fill = NA, na.rm = TRUE)
}

met <- aligned_data %>% arrange(datetime)

# Extract window lengths from predictor names
extract_window <- function(pred_name) {
  as.numeric(str_extract(pred_name, "\\d+(?=h$)"))
}

# Get windows for each model's predictors
ts_window_bgs <- extract_window(ts_pred_bgs)
wtd_window_bgs <- extract_window(wtd_pred_bgs)
ts_window_ems <- extract_window(ts_pred_ems)
wtd_window_ems <- extract_window(wtd_pred_ems)

cat("\nBuilding features with windows:\n")
cat("  BGS - TS window:", ts_window_bgs, "h, WTD window:", wtd_window_bgs, "h\n")
cat("  EMS - TS window:", ts_window_ems, "h, WTD window:", wtd_window_ems, "h\n")

# Extract base variable names
ts_var_bgs <- str_extract(ts_pred_bgs, "^[^_]+_[^_]+")
wtd_var_bgs <- str_extract(wtd_pred_bgs, "^[^_]+_[^_]+_[^_]+")

ts_var_ems <- str_extract(ts_pred_ems, "^[^_]+_[^_]+")
wtd_var_ems <- str_extract(wtd_pred_ems, "^[^_]+_[^_]+_[^_]+")

# Build features for BGS predictors
if (!is.na(ts_window_bgs) && ts_var_bgs %in% names(met)) {
  met[[ts_pred_bgs]] <- roll_mean(met[[ts_var_bgs]], ts_window_bgs)
}
if (!is.na(wtd_window_bgs) && wtd_var_bgs %in% names(met)) {
  met[[wtd_pred_bgs]] <- roll_mean(met[[wtd_var_bgs]], wtd_window_bgs)
}

# Build features for EMS predictors (if different)
if (!ts_pred_ems %in% names(met) && !is.na(ts_window_ems) && ts_var_ems %in% names(met)) {
  met[[ts_pred_ems]] <- roll_mean(met[[ts_var_ems]], ts_window_ems)
}
if (!wtd_pred_ems %in% names(met) && !is.na(wtd_window_ems) && wtd_var_ems %in% names(met)) {
  met[[wtd_pred_ems]] <- roll_mean(met[[wtd_var_ems]], wtd_window_ems)
}

# Process flux data
flux_processed <- flux_data %>%
  mutate(
    datetime = round_date(as.POSIXct(datetime_posx, tz = "EST"), "hour"),
    site = location,
    Tree = as.factor(Tree),
    species = as.factor(SPECIES),
    CH4_flux = CH4_flux_nmolpm2ps / 1000
  )

# Select features to join
features_to_join <- unique(c(ts_pred_bgs, wtd_pred_bgs, ts_pred_ems, wtd_pred_ems))
features_to_join <- features_to_join[features_to_join %in% names(met)]

# Join with environmental features
model_data_full <- flux_processed %>%
  left_join(met %>% dplyr::select(datetime, all_of(features_to_join)), 
            by = "datetime")

# Split by site and filter for complete data
data_bgs <- model_data_full %>% 
  filter(site == "Wetland") %>%
  filter(!is.na(.data[[ts_pred_bgs]]), !is.na(.data[[wtd_pred_bgs]]))

data_ems <- model_data_full %>% 
  filter(site == "Upland") %>%
  filter(!is.na(.data[[ts_pred_ems]]), !is.na(.data[[wtd_pred_ems]]))

cat("\nData after filtering:\n")
cat("  BGS:", nrow(data_bgs), "observations\n")
cat("  EMS:", nrow(data_ems), "observations\n")

# Get scaling parameters for each site using their respective predictors
scaling_bgs <- list(
  ts_mean = mean(data_bgs[[ts_pred_bgs]], na.rm = TRUE),
  ts_sd = sd(data_bgs[[ts_pred_bgs]], na.rm = TRUE),
  wtd_mean = mean(data_bgs[[wtd_pred_bgs]], na.rm = TRUE),
  wtd_sd = sd(data_bgs[[wtd_pred_bgs]], na.rm = TRUE)
)

scaling_ems <- list(
  ts_mean = mean(data_ems[[ts_pred_ems]], na.rm = TRUE),
  ts_sd = sd(data_ems[[ts_pred_ems]], na.rm = TRUE),
  wtd_mean = mean(data_ems[[wtd_pred_ems]], na.rm = TRUE),
  wtd_sd = sd(data_ems[[wtd_pred_ems]], na.rm = TRUE)
)

cat("\nBGS scaling - TS: mean =", round(scaling_bgs$ts_mean, 2), 
    ", sd =", round(scaling_bgs$ts_sd, 2), "\n")
cat("BGS scaling - WTD: mean =", round(scaling_bgs$wtd_mean, 2), 
    ", sd =", round(scaling_bgs$wtd_sd, 2), "\n")
cat("EMS scaling - TS: mean =", round(scaling_ems$ts_mean, 2), 
    ", sd =", round(scaling_ems$ts_sd, 2), "\n")
cat("EMS scaling - WTD: mean =", round(scaling_ems$wtd_mean, 2), 
    ", sd =", round(scaling_ems$wtd_sd, 2), "\n")

# ============================================================
# GENERATE PREDICTIONS
# ============================================================

message("\nGenerating predictions...")

# Water table percentiles for prediction (every 10th)
wtd_percentiles <- seq(10, 90, by = 10)

# Function to generate predictions for a model
generate_predictions <- function(model, scaling, site_data, site_name, 
                                 ts_pred, wtd_pred, additional_preds,
                                 limit_range = TRUE) {
  
  # Get species levels from the model
  species_levels <- levels(model@frame$species)
  
  # Get water table percentile values (unscaled)
  wtd_quantiles_unscaled <- quantile(site_data[[wtd_pred]], 
                                     probs = wtd_percentiles/100, na.rm = TRUE)
  
  # Scale to match model input
  wtd_quantiles_scaled <- (wtd_quantiles_unscaled - scaling$wtd_mean) / scaling$wtd_sd
  
  # Full temperature range (used if limit_range = FALSE)
  full_temp_min <- quantile(site_data[[ts_pred]], 0.05, na.rm = TRUE)
  full_temp_max <- quantile(site_data[[ts_pred]], 0.95, na.rm = TRUE)
  
  # Generate predictions for each species and wtd percentile
  pred_list <- list()
  
  for (sp in species_levels) {
    for (i in seq_along(wtd_percentiles)) {
      
      if (limit_range) {
        # Get temperature range for observations near this water table level
        # Define "near" as within +/- 0.5 percentile bands
        wtd_lo <- ifelse(i == 1, -Inf, (wtd_quantiles_unscaled[i] + wtd_quantiles_unscaled[i-1]) / 2)
        wtd_hi <- ifelse(i == length(wtd_percentiles), Inf, (wtd_quantiles_unscaled[i] + wtd_quantiles_unscaled[i+1]) / 2)
        
        nearby_data <- site_data %>%
          filter(.data[[wtd_pred]] >= wtd_lo, .data[[wtd_pred]] <= wtd_hi)
        
        # If not enough nearby data, use full temperature range but be conservative
        if (nrow(nearby_data) < 10) {
          temp_min <- full_temp_min
          temp_max <- full_temp_max
        } else {
          temp_min <- quantile(nearby_data[[ts_pred]], 0.05, na.rm = TRUE)
          temp_max <- quantile(nearby_data[[ts_pred]], 0.95, na.rm = TRUE)
        }
      } else {
        # Use full range
        temp_min <- full_temp_min
        temp_max <- full_temp_max
      }
      
      temp_range_unscaled <- seq(temp_min, temp_max, length.out = 100)
      temp_range_scaled <- (temp_range_unscaled - scaling$ts_mean) / scaling$ts_sd
      
      n_pts <- length(temp_range_scaled)
      
      # Build prediction dataframe
      pred_df <- data.frame(
        species = factor(rep(sp, n_pts), levels = species_levels),
        Tree = factor(rep(model@frame$Tree[1], n_pts), levels = levels(model@frame$Tree))
      )
      
      # Add core predictors
      pred_df[[ts_pred]] <- temp_range_scaled
      pred_df[[wtd_pred]] <- rep(wtd_quantiles_scaled[i], n_pts)
      
      # Add any additional predictors at their mean value (0 for scaled data)
      for (ap in additional_preds) {
        if (!ap %in% names(pred_df)) {
          pred_df[[ap]] <- rep(0, n_pts)
        }
      }
      
      # Add dummy response variable (needed for model.matrix)
      pred_df$CH4_flux_asinh <- 0
      
      # Predict with SE (fixed effects only)
      # Build model matrix for predictions
      mm <- model.matrix(formula(model, fixed.only = TRUE)[-2], pred_df)
      
      # Get variance-covariance matrix of fixed effects
      vcov_mat <- as.matrix(vcov(model))
      
      # Predicted values on asinh scale
      pred_df$fit_asinh <- as.vector(mm %*% fixef(model))
      
      # SE on asinh scale
      pred_df$se_asinh <- sqrt(diag(mm %*% vcov_mat %*% t(mm)))
      
      # Back-transform: sinh(asinh_value) gives nmol directly
      pred_df$fit_nmol <- sinh(pred_df$fit_asinh)
      
      # Back-transform CI bounds (delta method approximation)
      pred_df$lwr_nmol <- sinh(pred_df$fit_asinh - 1.96 * pred_df$se_asinh)
      pred_df$upr_nmol <- sinh(pred_df$fit_asinh + 1.96 * pred_df$se_asinh)
      
      pred_df$temp_unscaled <- temp_range_unscaled
      pred_df$wtd_unscaled <- wtd_quantiles_unscaled[i]
      pred_df$wtd_percentile <- wtd_percentiles[i]
      pred_df$species_name <- sp
      pred_df$site <- site_name
      
      pred_list[[length(pred_list) + 1]] <- pred_df
    }
  }
  
  bind_rows(pred_list)
}

# Generate predictions for both sites - LIMITED RANGE
pred_bgs_limited <- generate_predictions(m_bgs, scaling_bgs, data_bgs, "Wetland",
                                         ts_pred_bgs, wtd_pred_bgs, additional_preds_bgs,
                                         limit_range = TRUE)
pred_ems_limited <- generate_predictions(m_ems, scaling_ems, data_ems, "Upland",
                                         ts_pred_ems, wtd_pred_ems, additional_preds_ems,
                                         limit_range = TRUE)

# Generate predictions for both sites - FULL RANGE
pred_bgs_full <- generate_predictions(m_bgs, scaling_bgs, data_bgs, "Wetland",
                                      ts_pred_bgs, wtd_pred_bgs, additional_preds_bgs,
                                      limit_range = FALSE)
pred_ems_full <- generate_predictions(m_ems, scaling_ems, data_ems, "Upland",
                                      ts_pred_ems, wtd_pred_ems, additional_preds_ems,
                                      limit_range = FALSE)

# ============================================================
# AVERAGE ACROSS SPECIES (for site-level plot)
# ============================================================

message("Averaging across species...")

# Helper function to average predictions
average_predictions <- function(pred_data) {
  pred_data %>%
    group_by(site, wtd_percentile, temp_unscaled) %>%
    summarize(
      fit_nmol = mean(fit_nmol, na.rm = TRUE),
      lwr_nmol = mean(lwr_nmol, na.rm = TRUE),
      upr_nmol = mean(upr_nmol, na.rm = TRUE),
      wtd_unscaled = first(wtd_unscaled),
      .groups = "drop"
    )
}

# Average for limited range
pred_combined_limited <- bind_rows(
  average_predictions(pred_bgs_limited),
  average_predictions(pred_ems_limited)
) %>%
  mutate(
    wtd_percentile = factor(wtd_percentile, levels = wtd_percentiles),
    site = factor(site, levels = c("Wetland", "Upland"))
  )

# Average for full range
pred_combined_full <- bind_rows(
  average_predictions(pred_bgs_full),
  average_predictions(pred_ems_full)
) %>%
  mutate(
    wtd_percentile = factor(wtd_percentile, levels = wtd_percentiles),
    site = factor(site, levels = c("Wetland", "Upland"))
  )

# ============================================================
# CREATE PLOTS
# ============================================================

message("Creating plots...")

# Color palette - gradient for water table (low = brown/dry, high = blue/wet)
wtd_colors <- colorRampPalette(c("#8B4513", "#CD853F", "#DEB887", "#B0C4DE", "#6495ED", "#4169E1", "#0000CD"))(length(wtd_percentiles))
names(wtd_colors) <- as.character(wtd_percentiles)

# Helper function to create plot
make_plot <- function(pred_data, y_scales = "fixed") {
  ggplot(pred_data, 
         aes(x = temp_unscaled, y = fit_nmol, color = wtd_percentile, 
             fill = wtd_percentile)) +
    geom_ribbon(aes(ymin = lwr_nmol, ymax = upr_nmol), alpha = 0.15, color = NA) +
    geom_line(linewidth = 1) +
    facet_wrap(~ site, scales = y_scales) +
    scale_color_manual(
      values = wtd_colors,
      labels = paste0(wtd_percentiles, "%"),
      name = "Water table\npercentile"
    ) +
    scale_fill_manual(
      values = wtd_colors,
      guide = "none"
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
}

# ---- LIMITED RANGE plots ----
p_limited_fixed <- make_plot(pred_combined_limited, "fixed")
p_limited_free <- make_plot(pred_combined_limited, "free_y")

ggsave("outputs/figures/interaction_by_site_limited_fixed.png", p_limited_fixed,
       width = 10, height = 5, dpi = 300)
ggsave("outputs/figures/interaction_by_site_limited_free.png", p_limited_free,
       width = 10, height = 5, dpi = 300)

message("Saved: interaction_by_site_limited_fixed.png")
message("Saved: interaction_by_site_limited_free.png")

# ---- FULL RANGE plots ----
p_full_fixed <- make_plot(pred_combined_full, "fixed")
p_full_free <- make_plot(pred_combined_full, "free_y")

ggsave("outputs/figures/interaction_by_site_full_fixed.png", p_full_fixed,
       width = 10, height = 5, dpi = 300)
ggsave("outputs/figures/interaction_by_site_full_free.png", p_full_free,
       width = 10, height = 5, dpi = 300)

message("Saved: interaction_by_site_full_fixed.png")
message("Saved: interaction_by_site_full_free.png")

# Print one for display
print(p_full_free)

# ---- INSET VERSION: Fixed y-axis with upland inset on own scale ----
message("Creating inset version...")

# Wetland panel (full scale)
p_wetland <- ggplot(pred_combined_full %>% filter(site == "Wetland"), 
                    aes(x = temp_unscaled, y = fit_nmol, color = wtd_percentile, 
                        fill = wtd_percentile)) +
  geom_ribbon(aes(ymin = lwr_nmol, ymax = upr_nmol), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = wtd_colors,
    labels = paste0(wtd_percentiles, "%"),
    name = "Water table\npercentile"
  ) +
  scale_fill_manual(
    values = wtd_colors,
    guide = "none"
  ) +
  labs(
    x = "Soil temperature (°C)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Wetland"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    legend.position = "none"
  )

# Upland panel (same y-axis as wetland for main plot)
wetland_ylim <- range(pred_combined_full %>% filter(site == "Wetland") %>% 
                        pull(fit_nmol), na.rm = TRUE)

p_upland_main <- ggplot(pred_combined_full %>% filter(site == "Upland"), 
                        aes(x = temp_unscaled, y = fit_nmol, color = wtd_percentile, 
                            fill = wtd_percentile)) +
  geom_ribbon(aes(ymin = lwr_nmol, ymax = upr_nmol), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(
    values = wtd_colors,
    labels = paste0(wtd_percentiles, "%"),
    name = "Water table\npercentile"
  ) +
  scale_fill_manual(
    values = wtd_colors,
    guide = "none"
  ) +
  coord_cartesian(ylim = wetland_ylim) +
  labs(
    x = "Soil temperature (°C)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Upland"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
    legend.position = "none"
  )

# Upland inset (own y-axis scale)
p_upland_inset <- ggplot(pred_combined_full %>% filter(site == "Upland"), 
                         aes(x = temp_unscaled, y = fit_nmol, color = wtd_percentile, 
                             fill = wtd_percentile)) +
  geom_ribbon(aes(ymin = lwr_nmol, ymax = upr_nmol), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(values = wtd_colors, guide = "none") +
  scale_fill_manual(values = wtd_colors, guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    plot.background = element_rect(fill = "white", color = "grey50", linewidth = 0.5),
    axis.text = element_text(size = 6),
    plot.margin = margin(2, 2, 2, 2)
  )

# Combine upland main + inset
p_upland_with_inset <- p_upland_main + 
  annotation_custom(
    grob = ggplotGrob(p_upland_inset),
    xmin = max(pred_combined_full$temp_unscaled, na.rm = TRUE) - 
      0.45 * diff(range(pred_combined_full$temp_unscaled, na.rm = TRUE)),
    xmax = max(pred_combined_full$temp_unscaled, na.rm = TRUE) - 0.02,
    ymin = wetland_ylim[2] * 0.55,
    ymax = wetland_ylim[2] * 0.98
  )

# Extract legend from one of the main plots
legend_plot <- ggplot(pred_combined_full %>% filter(site == "Wetland"), 
                      aes(x = temp_unscaled, y = fit_nmol, color = wtd_percentile)) +
  geom_line() +
  scale_color_manual(
    values = wtd_colors,
    labels = paste0(wtd_percentiles, "%"),
    name = "Water table\npercentile"
  ) +
  theme(legend.position = "right")

legend <- cowplot::get_legend(legend_plot)

# Combine all panels
p_inset <- cowplot::plot_grid(
  p_wetland, 
  p_upland_with_inset, 
  legend,
  nrow = 1,
  rel_widths = c(1, 1, 0.25)
)

print(p_inset)

ggsave("outputs/figures/interaction_by_site_full_inset.png", p_inset,
       width = 12, height = 5, dpi = 300)

message("Saved: interaction_by_site_full_inset.png")

# ============================================================
# SUMMARY STATS
# ============================================================

cat("\n══════════════════════════════════════════════════════════════\n")
cat("                    SUMMARY\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# Model R² values
r2_bgs <- var(predict(m_bgs, re.form = NA)) / var(m_bgs@frame[[1]])
r2_ems <- var(predict(m_ems, re.form = NA)) / var(m_ems@frame[[1]])

cat("Model fit:\n")
cat(sprintf("  Wetland (BGS): R² = %.1f%%\n", r2_bgs * 100))
cat(sprintf("  Upland (EMS):  R² = %.1f%%\n", r2_ems * 100))

cat("\nPrediction ranges (nmol m⁻² s⁻¹):\n")
cat(sprintf("  Wetland: %.2f to %.2f\n", 
            min(pred_bgs_avg$fit_nmol), max(pred_bgs_avg$fit_nmol)))
cat(sprintf("  Upland:  %.2f to %.2f\n", 
            min(pred_ems_avg$fit_nmol), max(pred_ems_avg$fit_nmol)))

cat("\n══════════════════════════════════════════════════════════════\n")

message("\nDone!")




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
  flux = "data/input/HF_2023-2025_tree_flux_corrected.csv",
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
  filter(CH4_flux_nmolpm2ps >= -1,
         !is.na(PLOT)) %>%
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







# ============================================================
# INTERACTION PLOT RESULTS SUMMARY
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("INTERACTION PLOT - SUMMARY STATS")
message(paste(rep("=", 60), collapse = ""))

# Wetland flux range across conditions
message("\n--- WETLAND PREDICTIONS ---")
wetland_preds <- pred_combined_full %>% filter(site == "Wetland")

message("\nFlux range by water table percentile:")
wetland_summary <- wetland_preds %>%
  group_by(wtd_percentile) %>%
  summarise(
    min_flux = round(min(fit_nmol), 2),
    max_flux = round(max(fit_nmol), 2),
    range = round(max(fit_nmol) - min(fit_nmol), 2),
    .groups = "drop"
  )
print(wetland_summary)

message("\nDry (10th percentile) vs Wet (90th percentile):")
dry_wet <- wetland_preds %>%
  filter(wtd_percentile %in% c(10, 90)) %>%
  group_by(wtd_percentile) %>%
  summarise(
    temp_min = round(min(temp_unscaled), 1),
    temp_max = round(max(temp_unscaled), 1),
    flux_min = round(min(fit_nmol), 2),
    flux_max = round(max(fit_nmol), 2),
    .groups = "drop"
  )
print(dry_wet)

# Upland flux range
message("\n--- UPLAND PREDICTIONS ---")
upland_preds <- pred_combined_full %>% filter(site == "Upland")

message("\nFlux range by water table percentile:")
upland_summary <- upland_preds %>%
  group_by(wtd_percentile) %>%
  summarise(
    min_flux = round(min(fit_nmol), 3),
    max_flux = round(max(fit_nmol), 3),
    range = round(max(fit_nmol) - min(fit_nmol), 3),
    .groups = "drop"
  )
print(upland_summary)

message("\nOverall upland range:")
message("  Min flux: ", round(min(upland_preds$fit_nmol), 3), " nmol m⁻² s⁻¹")
message("  Max flux: ", round(max(upland_preds$fit_nmol), 3), " nmol m⁻² s⁻¹")

# Wetland:Upland ratio
message("\n--- WETLAND:UPLAND COMPARISON ---")
message("Maximum flux ratio: ", 
        round(max(wetland_preds$fit_nmol) / max(upland_preds$fit_nmol), 1), 
        "× higher at wetland")



# ============================================================
# COMBINED FIGURE: INTERACTION + TIMESERIES
# ============================================================

library(cowplot)

# Create combined figure with interaction on top, timeseries below
p_combined_fig <- plot_grid(
  p_inset,           # Interaction plot with inset (from interaction_updated.R)
  p_combined,        # Timeseries plot (from combined_driver_timeseries.R)
  ncol = 1,
  rel_heights = c(1, 1.5),
  labels = c("A", "B"),
  label_size = 14
)

print(p_combined_fig)

# Save
ggsave("outputs/figures/interaction_timeseries_combined.png", p_combined_fig,
       width = 12, height = 11, dpi = 300)
ggsave("outputs/figures/interaction_timeseries_combined.pdf", p_combined_fig,
       width = 12, height = 11)

message("Saved: interaction_timeseries_combined.png/pdf")

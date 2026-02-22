# Check the distribution of CH4 flux
summary(model_data_scaled$CH4_flux)
hist(model_data_scaled$CH4_flux, breaks = 50, main = "Raw CH4 flux")

# Check for negative values
sum(model_data_scaled$CH4_flux < 0)

# Create asinh-transformed response
model_data_scaled <- model_data_scaled %>%
  mutate(CH4_flux_asinh = asinh(CH4_flux * 1000))  # scale up first since values are small

# Check transformed distribution
hist(model_data_scaled$CH4_flux_asinh, breaks = 50, main = "asinh(CH4 flux)")

# Refit the final model with transformed response
m_final_asinh <- lmer(CH4_flux_asinh ~ TS_Ha2_raw_78h * bvs_wtd_cm_raw_189h * species + 
                        SWC_Ha2_raw_189h + (1|Tree), 
                      data = model_data_scaled, REML = FALSE)

summary(m_final_asinh)

# Compare R² 
var_total_asinh <- var(model_data_scaled$CH4_flux_asinh)
r2_asinh <- var(predict(m_final_asinh, re.form = NA)) / var_total_asinh

cat("\nR² with asinh transform:", round(r2_asinh, 4), "\n")
cat("R² with raw (previous): ", round(r2_simple, 4), "\n")



# Species-specific slopes from asinh model
cat("=== SPECIES-SPECIFIC SLOPES (asinh model) ===\n\n")

# Temperature main effect
cat("Temperature slope:\n")
cat("  Black gum:", round(0.58495, 3), "\n")
cat("  Hemlock:  ", round(0.58495 - 0.49351, 3), "\n")
cat("  Red maple:", round(0.58495 - 0.47256, 3), "\n")

# Water table main effect  
cat("\nWater table slope:\n")
cat("  Black gum:", round(0.57434, 3), "\n")
cat("  Hemlock:  ", round(0.57434 - 0.45464, 3), "\n")
cat("  Red maple:", round(0.57434 - 0.30047, 3), "\n")

# Temperature × Water table INTERACTION
cat("\nTemp × WTD interaction:\n")
cat("  Black gum:", round(0.40042, 3), "\n")
cat("  Hemlock:  ", round(0.40042 - 0.41145, 3), "\n")
cat("  Red maple:", round(0.40042 - 0.26332, 3), "\n")

# Check if we should keep SWC (marginal, t = -2.31)
# Try without it
m_no_swc <- lmer(CH4_flux_asinh ~ TS_Ha2_raw_78h * bvs_wtd_cm_raw_189h * species + (1|Tree), 
                 data = model_data_scaled, REML = FALSE)

cat("\n\n=== MODEL COMPARISON ===\n")
cat("With SWC - AIC:", AIC(m_final_asinh), "BIC:", BIC(m_final_asinh), "\n")
cat("No SWC   - AIC:", AIC(m_no_swc), "BIC:", BIC(m_no_swc), "\n")

anova(m_no_swc, m_final_asinh)

# VIF check
cat("\n\nVIF:\n")
vif(m_final_asinh)







# Update visualization with asinh model
# Back-transform predictions: sinh(predicted) / 1000 to get original units

# Quick test of back-transformation
pred_test <- predict(m_final_asinh, re.form = NA)
ch4_backtransformed <- sinh(pred_test) / 1000

cat("Back-transformed predictions range:\n")
cat("Min:", min(ch4_backtransformed), "\n")
cat("Max:", max(ch4_backtransformed), "\n")
cat("Mean:", mean(ch4_backtransformed), "\n")

# Compare to original
cat("\nOriginal CH4 flux range:\n")
cat("Min:", min(model_data_scaled$CH4_flux), "\n")
cat("Max:", max(model_data_scaled$CH4_flux), "\n")
cat("Mean:", mean(model_data_scaled$CH4_flux), "\n")











# ============================================================
# 17_final_model_visualization_asinh.R
# 
# Final model with asinh-transformed CH4 flux:
# asinh(CH4*1000) ~ TS_Ha2 * bvs_wtd_cm * species + SWC_Ha2 + (1|Tree)
#
# Creates:
# 1. Effect size plot (coefficient plot)
# 2. Species-specific interaction plots with back-transformed predictions
# ============================================================

library(tidyverse)
library(lme4)
library(RcppRoll)
library(cowplot)
library(car)

set.seed(42)

# ============================================================
# CONFIGURATION
# ============================================================

OUTPUT_DIR <- "figures/final_model"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

roll_mean <- function(x, n) RcppRoll::roll_mean(x, n = n, align = "right", fill = NA, na.rm = TRUE)

# ============================================================
# LOAD AND PREPARE DATA
# ============================================================

message("Loading data...")

aligned_data <- read_csv("data/processed/aligned_hourly_dataset.csv", show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC"))

stem_flux_raw <- read_csv("data/raw/flux_dataset.csv", show_col_types = FALSE)

stem_flux <- stem_flux_raw %>%
  mutate(
    datetime = round_date(force_tz(as.POSIXct(real_start - 2190), tzone = "EST"), "hour"),
    site = ifelse(Plot == "BGS", "BGS", "EMS"),
    Tree = as.factor(Tree),
    species = as.factor(species)
  ) %>%
  filter(site == "BGS") %>%
  dplyr::select(datetime, Tree, species, CH4_flux)

# Build features
met <- aligned_data %>% arrange(datetime)
features <- met %>% dplyr::select(datetime)

features$TS_Ha2_raw_78h <- roll_mean(met$TS_Ha2, 78)
features$bvs_wtd_cm_raw_189h <- roll_mean(met$bvs_wtd_cm, 189)
features$SWC_Ha2_raw_189h <- roll_mean(met$SWC_Ha2, 189)

# Prepare model data
model_data <- stem_flux %>%
  left_join(features, by = "datetime") %>%
  drop_na()

pred_names <- c("TS_Ha2_raw_78h", "bvs_wtd_cm_raw_189h", "SWC_Ha2_raw_189h")

# Store unscaled for plotting
model_data_unscaled <- model_data

# Get scaling parameters BEFORE scaling
scaling_params <- list(
  TS_mean = mean(model_data$TS_Ha2_raw_78h),
  TS_sd = sd(model_data$TS_Ha2_raw_78h),
  WTD_mean = mean(model_data$bvs_wtd_cm_raw_189h),
  WTD_sd = sd(model_data$bvs_wtd_cm_raw_189h),
  SWC_mean = mean(model_data$SWC_Ha2_raw_189h),
  SWC_sd = sd(model_data$SWC_Ha2_raw_189h)
)

# Scale predictors and create asinh-transformed response
model_data_scaled <- model_data %>%
  mutate(
    across(all_of(pred_names), ~ scale(.)[,1]),
    CH4_flux_asinh = asinh(CH4_flux * 1000)
  )

var_total_asinh <- var(model_data_scaled$CH4_flux_asinh)

cat("Model data:", nrow(model_data_scaled), "observations\n")
cat("N trees:", n_distinct(model_data_scaled$Tree), "\n")
cat("N species:", n_distinct(model_data_scaled$species), "\n")

# ============================================================
# FIT FINAL MODEL
# ============================================================

message("\nFitting final model (asinh-transformed)...")

m_final <- lmer(CH4_flux_asinh ~ TS_Ha2_raw_78h * bvs_wtd_cm_raw_189h * species + 
                  SWC_Ha2_raw_189h + (1|Tree), 
                data = model_data_scaled, REML = FALSE)

cat("\nModel summary:\n")
print(summary(m_final))

# R²
r2_final <- var(predict(m_final, re.form = NA)) / var_total_asinh
cat("\nR² (fixed effects):", round(r2_final, 4), "\n")

# Save model
saveRDS(m_final, file.path(OUTPUT_DIR, "m_final_asinh.rds"))

# ============================================================
# 1. EFFECT SIZE PLOT
# ============================================================

message("\n========== CREATING EFFECT SIZE PLOT ==========")

# Extract coefficients and SEs
coef_df <- as.data.frame(summary(m_final)$coefficients) %>%
  rownames_to_column("term") %>%
  rename(estimate = Estimate, se = `Std. Error`, t_value = `t value`) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    significant = abs(t_value) > 2
  )

# Create readable labels
coef_df <- coef_df %>%
  mutate(
    term_clean = case_when(
      term == "TS_Ha2_raw_78h" ~ "Temperature (black gum)",
      term == "bvs_wtd_cm_raw_189h" ~ "Water table (black gum)",
      term == "SWC_Ha2_raw_189h" ~ "Soil moisture (all species)",
      term == "specieshem" ~ "Hemlock (intercept diff)",
      term == "speciesrm" ~ "Red maple (intercept diff)",
      term == "TS_Ha2_raw_78h:bvs_wtd_cm_raw_189h" ~ "Temp × WTD (black gum)",
      term == "TS_Ha2_raw_78h:specieshem" ~ "Temp × Hemlock",
      term == "TS_Ha2_raw_78h:speciesrm" ~ "Temp × Red maple",
      term == "bvs_wtd_cm_raw_189h:specieshem" ~ "WTD × Hemlock",
      term == "bvs_wtd_cm_raw_189h:speciesrm" ~ "WTD × Red maple",
      term == "TS_Ha2_raw_78h:bvs_wtd_cm_raw_189h:specieshem" ~ "Temp × WTD × Hemlock",
      term == "TS_Ha2_raw_78h:bvs_wtd_cm_raw_189h:speciesrm" ~ "Temp × WTD × Red maple",
      TRUE ~ term
    ),
    effect_type = case_when(
      grepl("× .* ×", term_clean) ~ "3-way interaction",
      grepl("×", term_clean) ~ "2-way interaction",
      grepl("intercept", term_clean, ignore.case = TRUE) ~ "Species intercept",
      TRUE ~ "Main effect"
    )
  )

# Order terms logically
term_order <- c(
  "Temperature (black gum)",
  "Water table (black gum)",
  "Soil moisture (all species)",
  "Temp × WTD (black gum)",
  "Hemlock (intercept diff)",
  "Red maple (intercept diff)",
  "Temp × Hemlock",
  "Temp × Red maple",
  "WTD × Hemlock",
  "WTD × Red maple",
  "Temp × WTD × Hemlock",
  "Temp × WTD × Red maple"
)

coef_df <- coef_df %>%
  mutate(term_clean = factor(term_clean, levels = rev(term_order)))

# Color by effect type
effect_colors <- c(
  "Main effect" = "#2166AC",
  "Species intercept" = "#4DAF4A", 
  "2-way interaction" = "#FF7F00",
  "3-way interaction" = "#E31A1C"
)

p_effects <- ggplot(coef_df, aes(x = estimate, y = term_clean, color = effect_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 0.8) +
  geom_point(aes(shape = significant), size = 3) +
  scale_color_manual(values = effect_colors, name = "Effect type") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), 
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                     name = "Significance") +
  labs(
    x = "Coefficient (± 95% CI) on asinh scale",
    y = NULL,
    title = "Effect sizes for CH₄ flux model",
    subtitle = paste0("asinh(CH₄) ~ Temp × WTD × Species + SWC + (1|Tree)  |  R² = ", 
                      round(r2_final * 100, 1), "%")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_effects)

ggsave(file.path(OUTPUT_DIR, "effect_sizes_asinh.png"), p_effects, 
       width = 10, height = 7, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "effect_sizes_asinh.pdf"), p_effects, 
       width = 10, height = 7)

message("Saved: effect_sizes_asinh.png/pdf")

# ============================================================
# 2. SPECIES-SPECIFIC INTERACTION PLOTS (back-transformed)
# ============================================================

message("\n========== CREATING INTERACTION PLOTS ==========")

# Species-specific slopes for reference
species_effects <- tibble(
  species = c("Black gum", "Hemlock", "Red maple"),
  species_code = c("bg", "hem", "rm"),
  temp_slope = c(0.585, 0.585 - 0.494, 0.585 - 0.473),
  wtd_slope = c(0.574, 0.574 - 0.455, 0.574 - 0.300),
  interaction = c(0.400, 0.400 - 0.411, 0.400 - 0.263)
)

cat("\nSpecies-specific slopes:\n")
print(species_effects)

# Create prediction grid
# Temperature range (scaled): -2 to 2 SD
# Water table: 10th, 25th, 50th, 75th, 90th percentiles

temp_range_scaled <- seq(-2, 2, length.out = 100)
wtd_percentiles <- c(10, 25, 50, 75, 90)

# Get WTD values at percentiles (unscaled then scale)
wtd_quantiles_unscaled <- quantile(model_data_unscaled$bvs_wtd_cm_raw_189h, 
                                   probs = wtd_percentiles/100)
wtd_quantiles_scaled <- (wtd_quantiles_unscaled - scaling_params$WTD_mean) / scaling_params$WTD_sd

# Convert temp range to unscaled for plotting
temp_unscaled <- temp_range_scaled * scaling_params$TS_sd + scaling_params$TS_mean

# Generate predictions for each species and WTD level
pred_list <- list()

for (sp in c("bg", "hem", "rm")) {
  for (i in seq_along(wtd_percentiles)) {
    pred_df <- data.frame(
      TS_Ha2_raw_78h = temp_range_scaled,
      bvs_wtd_cm_raw_189h = wtd_quantiles_scaled[i],
      SWC_Ha2_raw_189h = 0,  # median
      species = factor(sp, levels = c("bg", "hem", "rm")),
      Tree = factor(model_data_scaled$Tree[1])  # dummy
    )
    
    # Predict on asinh scale (fixed effects only)
    pred_df$fit_asinh <- predict(m_final, newdata = pred_df, re.form = NA)
    
    # Back-transform to original units: sinh(pred) / 1000
    pred_df$fit <- sinh(pred_df$fit_asinh) / 1000
    
    pred_df$temp_unscaled <- temp_unscaled
    pred_df$wtd_percentile <- wtd_percentiles[i]
    pred_df$wtd_unscaled <- wtd_quantiles_unscaled[i]
    pred_df$species_name <- case_when(
      sp == "bg" ~ "Black gum",
      sp == "hem" ~ "Hemlock",
      sp == "rm" ~ "Red maple"
    )
    
    pred_list[[length(pred_list) + 1]] <- pred_df
  }
}

pred_all <- bind_rows(pred_list) %>%
  mutate(
    species_name = factor(species_name, levels = c("Black gum", "Hemlock", "Red maple")),
    wtd_percentile = factor(wtd_percentile)
  )

# Color palette for water table (brown = low/dry, blue = high/wet)
wtd_colors <- c("10" = "#8B4513", "25" = "#CD853F", "50" = "#808080", 
                "75" = "#6495ED", "90" = "#0000CD")

# Panel plot by species
p_interaction <- ggplot(pred_all, aes(x = temp_unscaled, y = fit * 1000, 
                                      color = wtd_percentile)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ species_name, scales = "free_y") +
  scale_color_manual(
    values = wtd_colors,
    labels = paste0(wtd_percentiles, "th"),
    name = "Water table\npercentile"
  ) +
  labs(
    x = "Soil temperature (°C, 78h mean)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Species-specific temperature × water table interactions",
    subtitle = "Back-transformed predictions (SWC held at median)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p_interaction)

ggsave(file.path(OUTPUT_DIR, "interaction_by_species_asinh.png"), p_interaction,
       width = 12, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "interaction_by_species_asinh.pdf"), p_interaction,
       width = 12, height = 5)

message("Saved: interaction_by_species_asinh.png/pdf")

# ============================================================
# 3. SINGLE PANEL COMPARISON (wet vs dry)
# ============================================================

message("\nCreating comparison plot...")

# Just show 10th and 90th percentile for clarity
pred_extremes <- pred_all %>%
  filter(wtd_percentile %in% c("10", "90")) %>%
  mutate(
    condition = paste(species_name, ifelse(wtd_percentile == "90", "(wet)", "(dry)")),
    linetype = ifelse(wtd_percentile == "90", "solid", "dashed")
  )

species_colors <- c("Black gum" = "#1B9E77", "Hemlock" = "#D95F02", "Red maple" = "#7570B3")

p_comparison <- ggplot(pred_extremes, aes(x = temp_unscaled, y = fit * 1000, 
                                          color = species_name, linetype = wtd_percentile)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = species_colors, name = "Species") +
  scale_linetype_manual(values = c("10" = "dashed", "90" = "solid"),
                        labels = c("10" = "Dry (10th %ile)", "90" = "Wet (90th %ile)"),
                        name = "Water table") +
  labs(
    x = "Soil temperature (°C, 78h mean)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Temperature response differs by species and water table",
    subtitle = "Black gum shows strong synergy between warm and wet conditions"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p_comparison)

ggsave(file.path(OUTPUT_DIR, "interaction_comparison_asinh.png"), p_comparison,
       width = 9, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "interaction_comparison_asinh.pdf"), p_comparison,
       width = 9, height = 6)

message("Saved: interaction_comparison_asinh.png/pdf")

# ============================================================
# 4. SLOPE COMPARISON BAR PLOT
# ============================================================

message("\nCreating slope comparison plot...")

slopes_long <- species_effects %>%
  pivot_longer(cols = c(temp_slope, wtd_slope, interaction),
               names_to = "effect", values_to = "slope") %>%
  mutate(
    effect_label = case_when(
      effect == "temp_slope" ~ "Temperature",
      effect == "wtd_slope" ~ "Water table", 
      effect == "interaction" ~ "Temp × WTD"
    ),
    effect_label = factor(effect_label, levels = c("Temperature", "Water table", "Temp × WTD"))
  )

p_slopes <- ggplot(slopes_long, aes(x = species, y = slope, fill = species)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ effect_label, scales = "free_y") +
  scale_fill_manual(values = species_colors, guide = "none") +
  labs(
    x = NULL,
    y = "Slope (asinh scale, per SD)",
    title = "Species-specific environmental sensitivities",
    subtitle = "Black gum responds most strongly to all drivers"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(p_slopes)

ggsave(file.path(OUTPUT_DIR, "slopes_by_species_asinh.png"), p_slopes,
       width = 9, height = 5, dpi = 300)

message("Saved: slopes_by_species_asinh.png")

# ============================================================
# 5. COMBINED FIGURE FOR PAPER
# ============================================================

message("\n========== CREATING COMBINED FIGURE ==========")

p_combined <- plot_grid(
  p_effects + theme(legend.position = "bottom", legend.box = "horizontal"),
  p_interaction + theme(legend.position = "bottom"),
  ncol = 1,
  rel_heights = c(1, 0.8),
  labels = c("A", "B"),
  label_size = 14
)

ggsave(file.path(OUTPUT_DIR, "figure_combined_asinh.png"), p_combined,
       width = 12, height = 12, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "figure_combined_asinh.pdf"), p_combined,
       width = 12, height = 12)

message("Saved: figure_combined_asinh.png/pdf")

# ============================================================
# SAVE SUMMARY STATISTICS
# ============================================================

message("\n========== SAVING SUMMARY STATS ==========")

summary_stats <- tibble(
  metric = c("R² (fixed effects)", "AIC", "BIC", "N observations", "N trees", "N species",
             "Residual SD", "Tree SD (random intercept)"),
  value = c(
    round(r2_final, 4),
    round(AIC(m_final), 1),
    round(BIC(m_final), 1),
    nrow(model_data_scaled),
    n_distinct(model_data_scaled$Tree),
    n_distinct(model_data_scaled$species),
    round(sigma(m_final), 4),
    round(as.data.frame(VarCorr(m_final))$sdcor[1], 4)
  )
)

write_csv(summary_stats, file.path(OUTPUT_DIR, "model_summary_stats_asinh.csv"))
write_csv(coef_df, file.path(OUTPUT_DIR, "coefficients_asinh.csv"))
write_csv(species_effects, file.path(OUTPUT_DIR, "species_specific_slopes_asinh.csv"))

# VIF table
vif_df <- as.data.frame(vif(m_final)) %>%
  rownames_to_column("term")
write_csv(vif_df, file.path(OUTPUT_DIR, "vif_asinh.csv"))

message("\n========== DONE ==========")
message("Output directory: ", OUTPUT_DIR)
message("Final model R² = ", round(r2_final * 100, 1), "%")

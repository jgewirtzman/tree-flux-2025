# Full set of interactions to test
# 2-way: all pairs among env predictors + species interactions
# 3-way: key combinations

env_preds <- c("TS_Ha2_raw_78h", "bvs_wtd_cm_raw_189h", "SWC_Ha2_raw_189h", "gcc_raw_189h", "LE_Ha1_raw_123h")

# Generate all 2-way interactions among env predictors
two_way_env <- combn(env_preds, 2, function(x) paste(x, collapse = ":"))

# 2-way interactions with species
two_way_species <- paste0(env_preds, ":species")

# 3-way interactions (env × env × species, and key env × env × env)
three_way_species <- combn(env_preds, 2, function(x) paste(c(x, "species"), collapse = ":"))

# Key 3-way env interactions (focus on temperature interactions)
three_way_env <- c(
  "TS_Ha2_raw_78h:bvs_wtd_cm_raw_189h:SWC_Ha2_raw_189h",
  "TS_Ha2_raw_78h:bvs_wtd_cm_raw_189h:gcc_raw_189h",
  "TS_Ha2_raw_78h:bvs_wtd_cm_raw_189h:LE_Ha1_raw_123h",
  "TS_Ha2_raw_78h:SWC_Ha2_raw_189h:gcc_raw_189h",
  "TS_Ha2_raw_78h:SWC_Ha2_raw_189h:LE_Ha1_raw_123h"
)

all_interactions <- c(two_way_env, two_way_species, three_way_species, three_way_env)

cat("Testing", length(all_interactions), "interactions...\n\n")

base_formula <- "CH4_flux ~ TS_Ha2_raw_78h + bvs_wtd_cm_raw_189h + SWC_Ha2_raw_189h + gcc_raw_189h + LE_Ha1_raw_123h + species + (1|Tree)"
m_base <- lmer(as.formula(base_formula), data = model_data_scaled, REML = FALSE)

results <- list()

for (int in all_interactions) {
  formula_int <- paste(base_formula, "+", int)
  
  tryCatch({
    m_int <- lmer(as.formula(formula_int), data = model_data_scaled, REML = FALSE)
    
    # LRT
    lrt <- anova(m_base, m_int)
    
    results[[int]] <- tibble(
      interaction = int,
      aic = AIC(m_int),
      delta_aic = AIC(m_int) - AIC(m_base),
      df = lrt$Df[2],
      lrt_chisq = lrt$Chisq[2],
      lrt_p = lrt$`Pr(>Chisq)`[2]
    )
  }, error = function(e) {
    message("Error with ", int, ": ", e$message)
  }, warning = function(w) NULL)
}

int_results <- bind_rows(results) %>%
  arrange(lrt_p)

cat("All interactions tested (sorted by p-value):\n\n")
print(int_results %>% 
        mutate(
          delta_aic = round(delta_aic, 2),
          lrt_chisq = round(lrt_chisq, 2),
          lrt_p = format.pval(lrt_p, digits = 3),
          sig = ifelse(as.numeric(lrt_p) < 0.05, "*", "")
        ), n = 40)

# Show significant ones
cat("\n\nSignificant interactions (p < 0.05):\n")
sig_int <- int_results %>% filter(lrt_p < 0.05)
print(sig_int %>% 
        mutate(
          delta_aic = round(delta_aic, 2),
          lrt_chisq = round(lrt_chisq, 2),
          lrt_p = format.pval(lrt_p, digits = 3)
        ))











# Model with top species interactions
m_species_int <- lmer(CH4_flux ~ TS_Ha2_raw_78h * species + bvs_wtd_cm_raw_189h * species + 
                        SWC_Ha2_raw_189h + gcc_raw_189h + LE_Ha1_raw_123h + (1|Tree), 
                      data = model_data_scaled, REML = FALSE)

summary(m_species_int)

# Compare to base
cat("\n\nModel comparison:\n")
cat("Base AIC:", AIC(m_base), "\n")
cat("Species interactions AIC:", AIC(m_species_int), "\n")
cat("ΔAIC:", AIC(m_species_int) - AIC(m_base), "\n")

# R²
r2_species_int <- var(predict(m_species_int, re.form = NA)) / var_total
cat("\nR² base:", round(var(predict(m_base, re.form = NA)) / var_total, 4), "\n")
cat("R² species int:", round(r2_species_int, 4), "\n")

# What are the species-specific slopes?
cat("\n\nSpecies-specific effects:\n")
coef_table <- summary(m_species_int)$coefficients
print(round(coef_table, 5))








# Calculate species-specific slopes
cat("Species-specific slopes:\n\n")

cat("Temperature effect (TS_Ha2):\n")
cat("  Black gum:", round(0.00686, 5), "\n")
cat("  Hemlock:  ", round(0.00686 - 0.00465, 5), "\n")
cat("  Red maple:", round(0.00686 - 0.00447, 5), "\n")

cat("\nWater table effect (bvs_wtd_cm):\n")
cat("  Black gum:", round(0.00828, 5), "\n")
cat("  Hemlock:  ", round(0.00828 - 0.00627, 5), "\n")
cat("  Red maple:", round(0.00828 - 0.00572, 5), "\n")

# Should we add the 3-way interaction (TS × WTD × species)?
m_3way <- lmer(CH4_flux ~ TS_Ha2_raw_78h * bvs_wtd_cm_raw_189h * species + 
                 SWC_Ha2_raw_189h + gcc_raw_189h + LE_Ha1_raw_123h + (1|Tree), 
               data = model_data_scaled, REML = FALSE)

cat("\n\n3-way interaction test:\n")
anova(m_species_int, m_3way)

cat("\nAIC comparison:\n")
cat("2-way species int:", AIC(m_species_int), "\n")
cat("3-way species int:", AIC(m_3way), "\n")

# Also check VIF
cat("\n\nVIF for species interaction model:\n")
vif(m_species_int)









# Examine the 3-way interaction model
summary(m_3way)

# R²
r2_3way <- var(predict(m_3way, re.form = NA)) / var_total
cat("\n\nR² comparison:\n")
cat("Base (no interactions):", round(var(predict(m_base, re.form = NA)) / var_total, 4), "\n")
cat("2-way species int:     ", round(r2_species_int, 4), "\n")
cat("3-way species int:     ", round(r2_3way, 4), "\n")

# VIF for 3-way model
cat("\n\nVIF for 3-way model:\n")
vif(m_3way)

# Extract species-specific slopes for TS, WTD, and TS×WTD interaction
cat("\n\nSpecies-specific coefficients from 3-way model:\n")
coef_3way <- summary(m_3way)$coefficients
print(round(coef_3way, 5))







# Calculate complete species-specific effects
cat("=== SPECIES-SPECIFIC SLOPES ===\n\n")

# Temperature main effect
cat("Temperature slope:\n")
cat("  Black gum:", round(0.00354, 5), "\n")
cat("  Hemlock:  ", round(0.00354 - 0.00153, 5), "\n")
cat("  Red maple:", round(0.00354 - 0.00161, 5), "\n")

# Water table main effect  
cat("\nWater table slope:\n")
cat("  Black gum:", round(0.00919, 5), "\n")
cat("  Hemlock:  ", round(0.00919 - 0.00799, 5), "\n")
cat("  Red maple:", round(0.00919 - 0.00734, 5), "\n")

# Temperature × Water table INTERACTION
cat("\nTemp × WTD interaction (synergy):\n")
cat("  Black gum:", round(0.00606, 5), "(strong synergy)\n")
cat("  Hemlock:  ", round(0.00606 - 0.00711, 5), "(NEGATIVE - no synergy!)\n")
cat("  Red maple:", round(0.00606 - 0.00657, 5), "(near zero - no synergy)\n")

# Now check if we can simplify by dropping non-significant terms
cat("\n\n=== SIGNIFICANCE CHECK ===\n")
cat("gcc is marginal (t = -1.31)\n")
cat("LE_Ha1 is marginal (t = -1.94)\n")
cat("TS×species interactions marginal (t ~ -1.9)\n")

# Try dropping gcc and LE
m_3way_simple <- lmer(CH4_flux ~ TS_Ha2_raw_78h * bvs_wtd_cm_raw_189h * species + 
                        SWC_Ha2_raw_189h + (1|Tree), 
                      data = model_data_scaled, REML = FALSE)

cat("\n\nSimplified model (drop gcc, LE):\n")
cat("AIC full:    ", AIC(m_3way), "\n")
cat("AIC simple:  ", AIC(m_3way_simple), "\n")
cat("BIC full:    ", BIC(m_3way), "\n")
cat("BIC simple:  ", BIC(m_3way_simple), "\n")

summary(m_3way_simple)$coefficients






# Final comparison
r2_simple <- var(predict(m_3way_simple, re.form = NA)) / var_total

cat("=== FINAL MODEL COMPARISON ===\n\n")
cat("R² simple model:", round(r2_simple, 4), "\n")
cat("R² full model:  ", round(r2_3way, 4), "\n")

# VIF check
cat("\n\nVIF for simplified model:\n")
vif(m_3way_simple)

# Summary table for paper
cat("\n\n=== FOR PAPER ===\n")
cat("Final model explains", round(r2_simple * 100, 1), "% of variance\n")
cat("Species alone:", round(r2_species * 100, 1), "%\n")
cat("Environment adds:", round((r2_simple - r2_species) * 100, 1), "%\n")












# ============================================================
# 17_final_model_visualization.R
# 
# Visualizations for final model:
# CH4_flux ~ TS_Ha2 * bvs_wtd_cm * species + SWC_Ha2 + (1|Tree)
#
# 1. Effect size plot (coefficient plot)
# 2. Species-specific interaction plots (Temp × WTD)
# ============================================================

library(tidyverse)
library(lme4)
library(RcppRoll)
library(cowplot)

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

# Scale predictors
model_data_scaled <- model_data %>%
  mutate(across(all_of(pred_names), ~ scale(.)[,1]))

# Get scaling parameters
scaling_params <- model_data_unscaled %>%
  summarize(
    TS_mean = mean(TS_Ha2_raw_78h),
    TS_sd = sd(TS_Ha2_raw_78h),
    WTD_mean = mean(bvs_wtd_cm_raw_189h),
    WTD_sd = sd(bvs_wtd_cm_raw_189h),
    SWC_mean = mean(SWC_Ha2_raw_189h),
    SWC_sd = sd(SWC_Ha2_raw_189h)
  )

var_total <- var(model_data_scaled$CH4_flux)

cat("Model data:", nrow(model_data_scaled), "observations\n")

# ============================================================
# FIT FINAL MODEL
# ============================================================

message("Fitting final model...")

m_final <- lmer(CH4_flux ~ TS_Ha2_raw_78h * bvs_wtd_cm_raw_189h * species + 
                  SWC_Ha2_raw_189h + (1|Tree), 
                data = model_data_scaled, REML = FALSE)

cat("\nModel summary:\n")
print(summary(m_final))

# Save model
saveRDS(m_final, file.path(OUTPUT_DIR, "m_final_3way.rds"))

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
      term == "specieshem" ~ "Hemlock intercept",
      term == "speciesrm" ~ "Red maple intercept",
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
      grepl("×.*×", term_clean) ~ "3-way interaction",
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
  "Hemlock intercept",
  "Red maple intercept",
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
    x = "Standardized coefficient (± 95% CI)",
    y = NULL,
    title = "Effect sizes for CH₄ flux model",
    subtitle = "Coefficients from: CH4 ~ Temp × WTD × Species + SWC + (1|Tree)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(size = 10)
  )

print(p_effects)

ggsave(file.path(OUTPUT_DIR, "effect_sizes.png"), p_effects, 
       width = 10, height = 7, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "effect_sizes.pdf"), p_effects, 
       width = 10, height = 7)

message("Saved: effect_sizes.png/pdf")

# ============================================================
# 2. SPECIES-SPECIFIC INTERACTION PLOTS
# ============================================================

message("\n========== CREATING INTERACTION PLOTS ==========")

# Calculate species-specific slopes
species_effects <- tibble(
  species = c("Black gum", "Hemlock", "Red maple"),
  species_code = c("bg", "hem", "rm"),
  temp_slope = c(0.002217, 0.002217 - 0.001581, 0.002217 - 0.001623),
  wtd_slope = c(0.009241, 0.009241 - 0.007974, 0.009241 - 0.007393),
  interaction = c(0.006447, 0.006447 - 0.007142, 0.006447 - 0.006632)
)

cat("\nSpecies-specific slopes:\n")
print(species_effects)

# Create prediction grid for each species
# Vary temperature and water table, hold SWC at median (0 in scaled)

temp_range <- seq(-2, 2, length.out = 50)  # scaled
wtd_percentiles <- c(10, 25, 50, 75, 90)

# Get WTD values at percentiles (in scaled units)
wtd_quantiles_unscaled <- quantile(model_data_unscaled$bvs_wtd_cm_raw_189h, 
                                   probs = wtd_percentiles/100)
wtd_quantiles_scaled <- (wtd_quantiles_unscaled - scaling_params$WTD_mean) / scaling_params$WTD_sd

# Convert temp range back to unscaled for plotting
temp_unscaled <- temp_range * scaling_params$TS_sd + scaling_params$TS_mean

# Generate predictions for each species and WTD level
pred_list <- list()

for (sp in c("bg", "hem", "rm")) {
  for (i in seq_along(wtd_percentiles)) {
    pred_df <- data.frame(
      TS_Ha2_raw_78h = temp_range,
      bvs_wtd_cm_raw_189h = wtd_quantiles_scaled[i],
      SWC_Ha2_raw_189h = 0,  # median
      species = factor(sp, levels = c("bg", "hem", "rm")),
      Tree = factor(model_data_scaled$Tree[1])  # dummy for prediction
    )
    
    # Predict (fixed effects only)
    pred_df$fit <- predict(m_final, newdata = pred_df, re.form = NA)
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
    labels = paste0(wtd_percentiles, "th %ile"),
    name = "Water table\npercentile"
  ) +
  labs(
    x = "Soil temperature (°C, 78h mean)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Species-specific temperature × water table interactions",
    subtitle = "Lines show predicted CH₄ flux at different water table levels (SWC at median)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p_interaction)

ggsave(file.path(OUTPUT_DIR, "interaction_by_species.png"), p_interaction,
       width = 12, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "interaction_by_species.pdf"), p_interaction,
       width = 12, height = 5)

message("Saved: interaction_by_species.png/pdf")

# ============================================================
# 3. SINGLE PANEL WITH ALL SPECIES (for comparison)
# ============================================================

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

ggsave(file.path(OUTPUT_DIR, "interaction_comparison.png"), p_comparison,
       width = 9, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "interaction_comparison.pdf"), p_comparison,
       width = 9, height = 6)

message("Saved: interaction_comparison.png/pdf")

# ============================================================
# 4. SLOPE COMPARISON BAR PLOT
# ============================================================

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

p_slopes <- ggplot(slopes_long, aes(x = species, y = slope * 1000, fill = species)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ effect_label, scales = "free_y") +
  scale_fill_manual(values = species_colors, guide = "none") +
  labs(
    x = NULL,
    y = "Slope (nmol m⁻² s⁻¹ per SD)",
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

ggsave(file.path(OUTPUT_DIR, "slopes_by_species.png"), p_slopes,
       width = 9, height = 5, dpi = 300)

message("Saved: slopes_by_species.png")

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

ggsave(file.path(OUTPUT_DIR, "figure_combined.png"), p_combined,
       width = 12, height = 12, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "figure_combined.pdf"), p_combined,
       width = 12, height = 12)

message("Saved: figure_combined.png/pdf")

# ============================================================
# SAVE SUMMARY STATISTICS
# ============================================================

message("\n========== SAVING SUMMARY STATS ==========")

# Model comparison table
r2_final <- var(predict(m_final, re.form = NA)) / var_total

summary_stats <- tibble(
  metric = c("R² (total)", "R² (fixed effects)", "AIC", "BIC", "N observations", "N trees"),
  value = c(
    round(r2_final, 3),
    round(r2_final, 3),
    round(AIC(m_final), 1),
    round(BIC(m_final), 1),
    nrow(model_data_scaled),
    n_distinct(model_data_scaled$Tree)
  )
)

write_csv(summary_stats, file.path(OUTPUT_DIR, "model_summary_stats.csv"))
write_csv(coef_df, file.path(OUTPUT_DIR, "coefficients_with_ci.csv"))
write_csv(species_effects, file.path(OUTPUT_DIR, "species_specific_slopes.csv"))

message("\nDone!")
message("Output directory: ", OUTPUT_DIR)


# ============================================
# 16l_theory_driven_randomized.R
#
# Theory-driven model building:
# 1. Core drivers (forced): soil temp, water table, soil moisture (all raw)
# 2. Candidate pool: 
#    - ALL other significant raw predictors
#    - PLUS specific LE/FC anomalies only
# 3. 100 randomized forward selection runs
#
# Model structure:
#   CH4_flux ~ predictors + species + (1|Tree)
# ============================================

library(tidyverse)
library(RcppRoll)
library(lme4)

set.seed(42)

# ============================================
# CONFIG
# ============================================

OUTPUT_DIR <- "figures/predictor_selection"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

SUM_VARS <- c("P_mm", "THROUGHFALL_xHA")
N_RANDOM_RUNS <- 100
MAX_PREDICTORS <- 6

# Core drivers (forced in, all raw) - defined by variable name
CORE_VARIABLES <- c("TS_Ha2", "bvs_wtd_cm", "NEON_SWC_shallow")

# Specific anomaly variables to include (LE and FC only)
SPECIFIC_ANOM_VARIABLES <- c("LE_Ha1", "LE_Ha2", "FC_Ha1", "FC_Ha2")

# ============================================
# HELPER FUNCTIONS
# ============================================

roll_mean <- function(x, n) RcppRoll::roll_mean(x, n = n, align = "right", fill = NA, na.rm = TRUE)
roll_sum  <- function(x, n) RcppRoll::roll_sum(x, n = n, align = "right", fill = NA, na.rm = TRUE)

make_anomaly <- function(df, var, time_col = "datetime") {
  if (!var %in% names(df)) return(df)
  x <- df[[var]]
  if (!is.numeric(x)) return(df)
  
  df <- df %>%
    mutate(
      .doy  = yday(.data[[time_col]]),
      .hour = hour(.data[[time_col]])
    )
  
  grand <- mean(x, na.rm = TRUE)
  
  doy_means <- df %>%
    group_by(.doy) %>%
    summarize(.m_doy = mean(.data[[var]], na.rm = TRUE), .groups = "drop")
  
  hr_means <- df %>%
    group_by(.hour) %>%
    summarize(.m_hr = mean(.data[[var]], na.rm = TRUE), .groups = "drop")
  
  df %>%
    left_join(doy_means, by = ".doy") %>%
    left_join(hr_means, by = ".hour") %>%
    mutate(
      .expected = .m_doy + .m_hr - grand,
      "{var}_anom" := .data[[var]] - .expected
    ) %>%
    dplyr::select(-.doy, -.hour, -.m_doy, -.m_hr, -.expected)
}

# ============================================
# LOAD DATA
# ============================================

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

bw <- read_csv("figures/rolling_correlations/best_windows_by_variable.csv", show_col_types = FALSE)

# ============================================
# BUILD PREDICTOR SETS FROM bw
# ============================================

message("\nBuilding predictor sets...")

# All significant RAW predictors
bgs_raw <- bw %>%
  filter(analysis == "raw", significant, site == "BGS") %>%
  dplyr::select(variable, var_group, window_hours, r) %>%
  mutate(
    predictor = paste0(variable, "_raw_", window_hours, "h"),
    mode = "raw"
  )

cat("Significant raw predictors:", nrow(bgs_raw), "\n")

# Specific LE/FC ANOMALY predictors only
bgs_anom_subset <- bw %>%
  filter(analysis == "anomaly", significant, site == "BGS",
         variable %in% SPECIFIC_ANOM_VARIABLES) %>%
  dplyr::select(variable, var_group, window_hours, r) %>%
  mutate(
    predictor = paste0(variable, "_anom_", window_hours, "h"),
    mode = "anom"
  )

cat("Specific anomaly predictors (LE/FC):", nrow(bgs_anom_subset), "\n")

# Core drivers (subset of raw)
core_info <- bgs_raw %>%
  filter(variable %in% CORE_VARIABLES)

cat("\nCore drivers:\n")
print(core_info %>% dplyr::select(predictor, variable, r))

# Candidates = all raw (except core) + specific anomalies
candidate_info <- bind_rows(
  bgs_raw %>% filter(!variable %in% CORE_VARIABLES),
  bgs_anom_subset
)

cat("\nCandidate predictors:", nrow(candidate_info), "\n")
print(candidate_info %>% dplyr::select(predictor, variable, mode, r) %>% arrange(mode, desc(abs(r))))

# ============================================
# BUILD FEATURES
# ============================================

message("\nBuilding features...")

met <- aligned_data %>% arrange(datetime)

# Get all variables we need
all_vars <- unique(c(bgs_raw$variable, bgs_anom_subset$variable))

# Create anomaly versions for the specific variables
for (v in SPECIFIC_ANOM_VARIABLES) {
  if (v %in% names(met)) {
    met <- make_anomaly(met, v, time_col = "datetime")
  }
}

# Build features
features <- met %>% dplyr::select(datetime)

# Raw features
for (i in seq_len(nrow(bgs_raw))) {
  v <- bgs_raw$variable[i]
  w <- bgs_raw$window_hours[i]
  if (!v %in% names(met)) next
  feat_name <- paste0(v, "_raw_", w, "h")
  if (v %in% SUM_VARS) {
    features[[feat_name]] <- roll_sum(met[[v]], w)
  } else {
    features[[feat_name]] <- roll_mean(met[[v]], w)
  }
}

# Specific anomaly features
for (i in seq_len(nrow(bgs_anom_subset))) {
  v <- bgs_anom_subset$variable[i]
  w <- bgs_anom_subset$window_hours[i]
  src <- paste0(v, "_anom")
  if (!src %in% names(met)) next
  feat_name <- paste0(v, "_anom_", w, "h")
  features[[feat_name]] <- roll_mean(met[[src]], w)
}

cat("Total features built:", ncol(features) - 1, "\n")

# ============================================
# PREPARE MODEL DATA
# ============================================

model_data <- stem_flux %>%
  left_join(features, by = "datetime") %>%
  drop_na()

core_names <- core_info$predictor
candidate_names <- candidate_info$predictor

# Verify all predictors exist
core_names <- core_names[core_names %in% names(model_data)]
candidate_names <- candidate_names[candidate_names %in% names(model_data)]

all_pred_names <- c(core_names, candidate_names)

# Scale predictors
model_data_scaled <- model_data %>%
  mutate(across(all_of(all_pred_names), ~ scale(.)[,1]))

cat("\nModel data:", nrow(model_data_scaled), "observations\n")
cat("Core predictors:", length(core_names), "\n")
cat("Candidate predictors:", length(candidate_names), "\n\n")

# Predictor info for later joining
pred_info <- bind_rows(core_info, candidate_info) %>%
  dplyr::select(predictor, variable, var_group, mode, r)

# ============================================
# 1. RANDOMIZED FORWARD SELECTION
# ============================================

message("========== 1. RANDOMIZED FORWARD SELECTION ==========\n")
cat("Running", N_RANDOM_RUNS, "iterations...\n")
cat("Core drivers (always first):", paste(core_names, collapse = ", "), "\n")
cat("Candidates:", length(candidate_names), "predictors\n\n")

var_total <- var(model_data_scaled$CH4_flux)

step_records <- list()

for (run in 1:N_RANDOM_RUNS) {
  if (run %% 20 == 0) cat("  Run", run, "/", N_RANDOM_RUNS, "\n")
  
  # Start with core drivers (always included)
  selected <- core_names
  
  # Shuffle candidates
  remaining <- sample(candidate_names)
  
  for (step in 1:MAX_PREDICTORS) {
    best_r2 <- -Inf
    best_pred <- NULL
    
    for (pred in remaining) {
      current_preds <- c(selected, pred)
      formula_str <- paste("CH4_flux ~", paste(current_preds, collapse = " + "), 
                           "+ species + (1|Tree)")
      
      tryCatch({
        fit <- lmer(as.formula(formula_str), data = model_data_scaled, REML = FALSE)
        var_fixed <- var(predict(fit, re.form = NA))
        r2 <- var_fixed / var_total
        
        if (r2 > best_r2) {
          best_r2 <- r2
          best_pred <- pred
        }
      }, error = function(e) NULL, warning = function(w) NULL)
    }
    
    if (is.null(best_pred)) break
    
    selected <- c(selected, best_pred)
    remaining <- setdiff(remaining, best_pred)
    
    step_records[[length(step_records) + 1]] <- tibble(
      run = run, 
      step = step, 
      predictor = best_pred, 
      r2 = best_r2
    )
  }
}

all_steps <- bind_rows(step_records)

# ============================================
# 2. SUMMARIZE SELECTION
# ============================================

message("\n========== 2. SELECTION FREQUENCY ==========\n")

selection_freq <- all_steps %>%
  group_by(predictor) %>%
  summarize(
    n_selected = n(),
    avg_step = mean(step),
    min_step = min(step),
    max_step = max(step),
    .groups = "drop"
  ) %>%
  mutate(pct_selected = n_selected / N_RANDOM_RUNS * 100) %>%
  left_join(pred_info, by = "predictor") %>%
  arrange(desc(n_selected))

cat("Selection frequency across", N_RANDOM_RUNS, "runs:\n\n")
print(selection_freq %>% 
        dplyr::select(predictor, variable, var_group, mode, r, pct_selected, avg_step) %>%
        mutate(r = round(r, 3), pct_selected = round(pct_selected, 1), avg_step = round(avg_step, 2)),
      n = 20)

# First step selection
first_step <- all_steps %>%
  filter(step == 1) %>%
  count(predictor) %>%
  mutate(pct_first = n / N_RANDOM_RUNS * 100) %>%
  left_join(pred_info, by = "predictor") %>%
  arrange(desc(n))

cat("\n\nFirst candidate selected (most important beyond core):\n")
print(first_step %>% 
        dplyr::select(predictor, variable, mode, n, pct_first) %>%
        head(10))

# ============================================
# 3. CONSENSUS PREDICTORS
# ============================================

message("\n========== 3. CONSENSUS ==========\n")

strong_consensus <- selection_freq %>% filter(pct_selected >= 80)
cat("Selected in ≥80% of runs:", nrow(strong_consensus), "predictors\n")
if (nrow(strong_consensus) > 0) {
  print(strong_consensus %>% 
          dplyr::select(predictor, variable, var_group, mode, r, pct_selected))
}

moderate_consensus <- selection_freq %>% filter(pct_selected >= 50, pct_selected < 80)
cat("\nSelected in 50-80% of runs:", nrow(moderate_consensus), "predictors\n")
if (nrow(moderate_consensus) > 0) {
  print(moderate_consensus %>% 
          dplyr::select(predictor, variable, var_group, mode, r, pct_selected))
}

# ============================================
# 4. FIT FINAL MODELS
# ============================================

message("\n========== 4. FINAL MODELS ==========\n")

# Species only
species_model <- lmer(CH4_flux ~ species + (1|Tree), data = model_data_scaled, REML = FALSE)
r2_species <- var(predict(species_model, re.form = NA)) / var_total

# Core only
core_formula <- paste("CH4_flux ~", paste(core_names, collapse = " + "), "+ species + (1|Tree)")
core_model <- lmer(as.formula(core_formula), data = model_data_scaled, REML = FALSE)
r2_core <- var(predict(core_model, re.form = NA)) / var_total

# Core + strong consensus
if (nrow(strong_consensus) > 0) {
  strong_preds <- c(core_names, strong_consensus$predictor)
  strong_formula <- paste("CH4_flux ~", paste(strong_preds, collapse = " + "), "+ species + (1|Tree)")
  strong_model <- lmer(as.formula(strong_formula), data = model_data_scaled, REML = FALSE)
  r2_strong <- var(predict(strong_model, re.form = NA)) / var_total
} else {
  strong_preds <- core_names
  strong_model <- core_model
  r2_strong <- r2_core
}

# Comparison
comparison <- tibble(
  model = c("Species only", "Core drivers", "Core + consensus (≥80%)"),
  n_env_preds = c(0, length(core_names), length(strong_preds)),
  aic = c(AIC(species_model), AIC(core_model), AIC(strong_model)),
  bic = c(BIC(species_model), BIC(core_model), BIC(strong_model)),
  r2_total = c(r2_species, r2_core, r2_strong),
  r2_env = c(0, r2_core - r2_species, r2_strong - r2_species)
)

cat("Model comparison:\n\n")
print(comparison %>% mutate(
  aic = round(aic, 2),
  bic = round(bic, 2),
  r2_total = round(r2_total, 4),
  r2_env = round(r2_env, 4)
))

# ============================================
# 5. BEST MODEL DETAILS
# ============================================

message("\n========== 5. BEST MODEL DETAILS ==========\n")

cat("Formula:", strong_formula, "\n\n")

cat("Fixed effects:\n")
print(summary(strong_model)$coefficients %>% round(5))

cat("\nRandom effects:\n")
print(VarCorr(strong_model))

cat("\nModel fit:\n")
cat("  AIC:", round(AIC(strong_model), 2), "\n")
cat("  BIC:", round(BIC(strong_model), 2), "\n")
cat("  R² (total):", round(r2_strong, 4), "\n")
cat("  R² (species):", round(r2_species, 4), "\n")
cat("  R² (environment):", round(r2_strong - r2_species, 4), "\n")

# ============================================
# 6. COLLINEARITY CHECK
# ============================================

message("\n========== 6. COLLINEARITY CHECK ==========\n")

cor_matrix <- model_data_scaled %>%
  dplyr::select(all_of(strong_preds)) %>%
  cor(use = "pairwise.complete.obs")

cat("Correlation matrix of final predictors:\n")
print(round(cor_matrix, 3))

# ============================================
# SAVE OUTPUTS
# ============================================

write_csv(selection_freq, file.path(OUTPUT_DIR, "theory_driven_selection_freq.csv"))
write_csv(first_step, file.path(OUTPUT_DIR, "theory_driven_first_step.csv"))
write_csv(comparison, file.path(OUTPUT_DIR, "theory_driven_comparison.csv"))

# Plot
png(file.path(OUTPUT_DIR, "theory_driven_selection_stability.png"), 
    width = 10, height = 8, units = "in", res = 300)

selection_freq %>%
  head(20) %>%
  mutate(
    predictor = fct_reorder(predictor, pct_selected),
    mode_label = ifelse(mode == "raw", "raw", "anomaly")
  ) %>%
  ggplot(aes(x = predictor, y = pct_selected, fill = var_group, alpha = mode)) +
  geom_col() +
  geom_hline(yintercept = 80, linetype = "dashed", color = "red") +
  geom_hline(yintercept = 50, linetype = "dotted", color = "orange") +
  scale_alpha_manual(values = c("raw" = 1, "anom" = 0.6)) +
  coord_flip() +
  labs(
    title = "Selection Stability: Theory-Driven Model",
    subtitle = paste0(N_RANDOM_RUNS, " randomized runs | Core: ", 
                      paste(CORE_VARIABLES, collapse = ", "),
                      " | Candidates: all sig. raw + LE/FC anomaly"),
    x = NULL,
    y = "% of runs where predictor was selected",
    fill = "Process",
    alpha = "Mode"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

dev.off()
message("Saved: theory_driven_selection_stability.png")

message("\n\nDone!")

# Install if needed: install.packages("car")
library(car)

# VIF for the full model
vif(strong_model)


# Compare models with each LE variant (dropping tair in all)
base_preds <- setdiff(strong_preds, c("tair_C_raw_285h", "LE_Ha1_raw_123h", "LE_Ha2_raw_120h", "LE_Ha1_anom_33h"))

# Model with LE_Ha1_raw only
m1_preds <- c(base_preds, "LE_Ha1_raw_123h")
m1 <- lmer(as.formula(paste("CH4_flux ~", paste(m1_preds, collapse = " + "), "+ species + (1|Tree)")), 
           data = model_data_scaled, REML = FALSE)

# Model with LE_Ha2_raw only
m2_preds <- c(base_preds, "LE_Ha2_raw_120h")
m2 <- lmer(as.formula(paste("CH4_flux ~", paste(m2_preds, collapse = " + "), "+ species + (1|Tree)")), 
           data = model_data_scaled, REML = FALSE)

# Model with LE_Ha1_anom only
m3_preds <- c(base_preds, "LE_Ha1_anom_33h")
m3 <- lmer(as.formula(paste("CH4_flux ~", paste(m3_preds, collapse = " + "), "+ species + (1|Tree)")), 
           data = model_data_scaled, REML = FALSE)

# Model with no LE
m0 <- lmer(as.formula(paste("CH4_flux ~", paste(base_preds, collapse = " + "), "+ species + (1|Tree)")), 
           data = model_data_scaled, REML = FALSE)

# Compare
data.frame(
  model = c("No LE", "LE_Ha1_raw", "LE_Ha2_raw", "LE_Ha1_anom"),
  aic = c(AIC(m0), AIC(m1), AIC(m2), AIC(m3)),
  bic = c(BIC(m0), BIC(m1), BIC(m2), BIC(m3))
) %>% 
  mutate(delta_aic = aic - min(aic)) %>%
  arrange(aic)




# Final reduced model
final_preds <- c(base_preds, "LE_Ha1_raw_123h")
final_formula <- paste("CH4_flux ~", paste(final_preds, collapse = " + "), "+ species + (1|Tree)")
final_model <- lmer(as.formula(final_formula), data = model_data_scaled, REML = FALSE)

cat("Final predictors:", paste(final_preds, collapse = ", "), "\n\n")

# VIF check
cat("VIF:\n")
vif(final_model)

# Summary
cat("\n\nCoefficients:\n")
print(summary(final_model)$coefficients %>% round(5))

# Fit
r2_final <- var(predict(final_model, re.form = NA)) / var_total
cat("\n\nModel fit:\n")
cat("AIC:", round(AIC(final_model), 2), "\n")
cat("BIC:", round(BIC(final_model), 2), "\n")
cat("R² (total):", round(r2_final, 4), "\n")
cat("R² (env):", round(r2_final - r2_species, 4), "\n")



cor_matrix <- model_data_scaled %>%
  dplyr::select(all_of(final_preds)) %>%
  cor(use = "pairwise.complete.obs")

round(cor_matrix, 2)







# Model 1: Minimal (just TS_Ha2 + water table)
m_minimal <- lmer(CH4_flux ~ TS_Ha2_raw_78h + bvs_wtd_cm_raw_189h + species + (1|Tree), 
                  data = model_data_scaled, REML = FALSE)

# Model 2: Full with SWC_Ha2 (drop NEON_SWC_shallow)
preds_swc_ha2 <- c("TS_Ha2_raw_78h", "bvs_wtd_cm_raw_189h", "SWC_Ha2_raw_189h", "gcc_raw_189h", "LE_Ha1_raw_123h")
m_swc_ha2 <- lmer(as.formula(paste("CH4_flux ~", paste(preds_swc_ha2, collapse = " + "), "+ species + (1|Tree)")), 
                  data = model_data_scaled, REML = FALSE)

# Model 3: Full with NEON_SWC_shallow (drop SWC_Ha2)
preds_neon <- c("TS_Ha2_raw_78h", "bvs_wtd_cm_raw_189h", "NEON_SWC_shallow_raw_63h", "gcc_raw_189h", "LE_Ha1_raw_123h")
m_neon <- lmer(as.formula(paste("CH4_flux ~", paste(preds_neon, collapse = " + "), "+ species + (1|Tree)")), 
               data = model_data_scaled, REML = FALSE)

# Compare
cat("=== MODEL COMPARISON ===\n\n")

models <- list(minimal = m_minimal, swc_ha2 = m_swc_ha2, neon_swc = m_neon)

comparison <- data.frame(
  model = names(models),
  n_preds = c(2, 5, 5),
  aic = sapply(models, AIC),
  bic = sapply(models, BIC),
  r2_total = sapply(models, function(m) var(predict(m, re.form = NA)) / var_total),
  r2_env = sapply(models, function(m) var(predict(m, re.form = NA)) / var_total - r2_species)
)

print(comparison %>% mutate(across(c(aic, bic), round, 2), 
                            across(c(r2_total, r2_env), round, 4)))

# VIF for each
cat("\n\n=== VIF ===\n")
cat("\nMinimal model:\n")
print(vif(m_minimal))

cat("\nSWC_Ha2 model:\n")
print(vif(m_swc_ha2))

cat("\nNEON_SWC model:\n")
print(vif(m_neon))

# Coefficients
cat("\n\n=== COEFFICIENTS ===\n")
cat("\nMinimal model:\n")
print(round(summary(m_minimal)$coefficients, 5))

cat("\nSWC_Ha2 model:\n")
print(round(summary(m_swc_ha2)$coefficients, 5))

cat("\nNEON_SWC model:\n")
print(round(summary(m_neon)$coefficients, 5))
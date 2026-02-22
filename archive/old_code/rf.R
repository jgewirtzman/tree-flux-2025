# ============================================================
# 09_random_forest_CH4.R
# 
# Random forest model of CH4 stem flux with variable importance
# and SHAP value analysis
#
# Uses optimal integration windows from correlation analysis

library(tidyverse)
library(lubridate)
library(RcppRoll)
library(ranger)       # Faster RF implementation
library(fastshap)     # SHAP values via Monte Carlo
library(patchwork)

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  stem_flux = "data/raw/flux_dataset.csv",
  aligned = "data/processed/aligned_hourly_dataset.csv",
  optimal_windows = "figures/rolling_correlations/optimal_windows.csv"
)

OUTPUT_DIR <- "figures/random_forest"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

DATE_MIN <- as.Date("2023-01-01")
DATE_MAX <- as.Date("2026-01-01")

# Variables to sum (precipitation) vs mean (everything else)
SUM_VARS <- c("P_mm", "THROUGHFALL_xHA")

# ============================================================
# LOAD DATA
# ============================================================

message("Loading data...")

# Load stem flux
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
  select(datetime, date, ID, site, species, CH4_flux, CO2_flux)

message("  Stem flux: ", nrow(stem_flux), " observations")

# Load aligned environmental data
aligned_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC")) %>%
  filter(as.Date(datetime) >= DATE_MIN, as.Date(datetime) <= DATE_MAX) %>%
  arrange(datetime)

message("  Aligned data: ", nrow(aligned_data), " hourly records")

# Load optimal windows
optimal_windows <- read_csv(PATHS$optimal_windows, show_col_types = FALSE)

message("  Optimal windows: ", nrow(optimal_windows), " variable-site combinations")

# ============================================================
# CALCULATE ROLLING STATISTICS AT OPTIMAL WINDOWS
# ============================================================

message("\n============================================================")
message("CALCULATING ROLLING STATISTICS AT OPTIMAL WINDOWS")
message("============================================================")

# Function to calculate rolling stat for one variable at one window
calc_rolling_var <- function(data, var, window_hours, is_sum = FALSE) {
  x <- data[[var]]
  if (is_sum) {
    roll_sum(x, n = window_hours, align = "right", fill = NA, na.rm = TRUE)
  } else {
    roll_mean(x, n = window_hours, align = "right", fill = NA, na.rm = TRUE)
  }
}

# Process each site separately (different optimal windows)
process_site_windows <- function(site_name) {
  
  message("\n--- ", site_name, " ---")
  
  # Get optimal windows for this site
  site_windows <- optimal_windows %>%
    filter(site == site_name) %>%
    select(variable, window_days, r) %>%
    filter(variable %in% names(aligned_data))
  
  message("  Variables with optimal windows: ", nrow(site_windows))
  
  # Start with datetime
  result <- aligned_data %>% select(datetime)
  
  # Calculate rolling stat for each variable at its optimal window
  for (i in seq_len(nrow(site_windows))) {
    var <- site_windows$variable[i]
    window_days <- site_windows$window_days[i]
    window_hours <- max(1, round(window_days * 24))
    
    is_sum <- var %in% SUM_VARS
    
    # Create column name with window info
    col_name <- paste0(var, "_", window_hours, "h")
    
    result[[col_name]] <- calc_rolling_var(aligned_data, var, window_hours, is_sum)
    
    message(sprintf("    %s: %d hour window (r = %.3f)", 
                    var, window_hours, site_windows$r[i]))
  }
  
  result
}

# Process each site
bgs_env <- process_site_windows("BGS")
ems_env <- process_site_windows("EMS")

# ============================================================
# JOIN FLUX DATA WITH OPTIMAL-WINDOW ENVIRONMENT DATA
# ============================================================

message("\n============================================================")
message("PREPARING MODEL DATA")
message("============================================================")

# Join flux with site-specific environmental data
bgs_flux <- stem_flux %>% filter(site == "BGS")
ems_flux <- stem_flux %>% filter(site == "EMS")

bgs_model_data <- bgs_flux %>%
  left_join(bgs_env, by = "datetime") %>%
  filter(!is.na(CH4_flux))

ems_model_data <- ems_flux %>%
  left_join(ems_env, by = "datetime") %>%
  filter(!is.na(CH4_flux))

message("  BGS observations: ", nrow(bgs_model_data))
message("  EMS observations: ", nrow(ems_model_data))

# ============================================================
# PREPARE DATA FOR MODELING
# ============================================================

prepare_model_data <- function(data, site_name) {
  
  # Exclude non-predictor columns
  exclude_cols <- c("datetime", "date", "year", "month", "doy", "hour", "wyear",
                    "ID", "site", "species", "CH4_flux", "CO2_flux")
  
  predictor_cols <- names(data)[!names(data) %in% exclude_cols]
  predictor_cols <- predictor_cols[sapply(data[predictor_cols], is.numeric)]
  
  model_df <- data %>%
    select(CH4_flux, all_of(predictor_cols))
  
  # Remove columns with all NA
  na_cols <- sapply(model_df, function(x) all(is.na(x)))
  model_df <- model_df[, !na_cols]
  
  # Complete cases only
  model_df <- model_df[complete.cases(model_df), ]
  
  predictors <- names(model_df)[names(model_df) != "CH4_flux"]
  
  message("  ", site_name, ": ", nrow(model_df), " complete cases, ", 
          length(predictors), " predictors")
  
  list(
    data = model_df,
    predictors = predictors
  )
}

bgs_prep <- prepare_model_data(bgs_model_data, "BGS")
ems_prep <- prepare_model_data(ems_model_data, "EMS")

# ============================================================
# FIT RANDOM FOREST MODELS
# ============================================================

message("\n============================================================")
message("FITTING RANDOM FOREST MODELS")
message("============================================================")

set.seed(42)

fit_rf_model <- function(prep_data, site_name) {
  
  message("\n--- ", site_name, " ---")
  
  X <- as.data.frame(prep_data$data[, prep_data$predictors])
  y <- prep_data$data$CH4_flux
  
  # Fit ranger model
  rf_model <- ranger(
    x = X,
    y = y,
    num.trees = 500,
    importance = "permutation",
    mtry = floor(sqrt(ncol(X))),
    min.node.size = 5,
    seed = 42
  )
  
  message("  R² (OOB): ", round(rf_model$r.squared, 3))
  message("  MSE (OOB): ", round(rf_model$prediction.error, 4))
  
  # Variable importance
  var_imp <- tibble(
    variable = names(rf_model$variable.importance),
    importance = rf_model$variable.importance
  ) %>%
    arrange(desc(importance)) %>%
    # Extract base variable name and window
    mutate(
      base_var = str_remove(variable, "_\\d+h$"),
      window_hours = as.numeric(str_extract(variable, "\\d+(?=h$)"))
    )
  
  message("\n  Top 10 variables:")
  print(head(var_imp %>% select(variable, importance), 10))
  
  list(
    model = rf_model,
    X = X,
    y = y,
    var_imp = var_imp,
    r_squared = rf_model$r.squared
  )
}

bgs_rf <- fit_rf_model(bgs_prep, "BGS")
ems_rf <- fit_rf_model(ems_prep, "EMS")

# ============================================================
# VARIABLE IMPORTANCE PLOTS
# ============================================================

message("\n============================================================")
message("GENERATING VARIABLE IMPORTANCE PLOTS")
message("============================================================")

# Combined importance plot
var_imp_combined <- bind_rows(
  bgs_rf$var_imp %>% mutate(site = "BGS"),
  ems_rf$var_imp %>% mutate(site = "EMS")
)

# Top 20 per site
top_vars <- var_imp_combined %>%
  group_by(site) %>%
  slice_max(importance, n = 20) %>%
  pull(variable) %>%
  unique()

p_imp <- ggplot(var_imp_combined %>% filter(variable %in% top_vars),
                aes(x = importance, y = reorder(variable, importance), fill = site)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("BGS" = "#1f77b4", "EMS" = "#ff7f0e")) +
  labs(
    title = "Random Forest Variable Importance",
    subtitle = paste0("BGS R² = ", round(bgs_rf$r_squared, 3), 
                      ", EMS R² = ", round(ems_rf$r_squared, 3),
                      "\nVariables at optimal integration windows"),
    x = "Permutation Importance",
    y = NULL,
    fill = "Site"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(OUTPUT_DIR, "variable_importance.png"), p_imp,
       width = 10, height = 10, dpi = 300)
message("  Saved: variable_importance.png")

# Separate plots per site
for (site_name in c("BGS", "EMS")) {
  rf_result <- if (site_name == "BGS") bgs_rf else ems_rf
  
  p_site <- ggplot(rf_result$var_imp %>% head(20),
                   aes(x = importance, y = reorder(variable, importance))) +
    geom_col(fill = ifelse(site_name == "BGS", "#1f77b4", "#ff7f0e")) +
    labs(
      title = paste("Variable Importance -", site_name),
      subtitle = paste0("R² = ", round(rf_result$r_squared, 3)),
      x = "Permutation Importance",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(OUTPUT_DIR, paste0("importance_", site_name, ".png")), p_site,
         width = 8, height = 8, dpi = 300)
}

# ============================================================
# SHAP VALUES (using fastshap)
# ============================================================

message("\n============================================================")
message("CALCULATING SHAP VALUES")
message("============================================================")

calc_shap <- function(rf_result, site_name, nsim = 50, max_obs = 300) {
  
  message("\n--- ", site_name, " ---")
  
  X <- rf_result$X
  
  # Sample if too many observations
  n_obs <- nrow(X)
  if (n_obs > max_obs) {
    set.seed(42)
    idx <- sample(1:n_obs, max_obs)
    X_sample <- X[idx, ]
    message("  Sampled ", max_obs, " of ", n_obs, " observations for SHAP")
  } else {
    X_sample <- X
    idx <- 1:n_obs
  }
  
  # Prediction wrapper for ranger
  pfun <- function(object, newdata) {
    predict(object, data = newdata)$predictions
  }
  
  # Calculate SHAP values using fastshap
  message("  Calculating SHAP values (nsim = ", nsim, ")...")
  shap_values <- fastshap::explain(
    rf_result$model,
    X = X,
    pred_wrapper = pfun,
    newdata = X_sample,
    nsim = nsim,
    .progress = "text"
  )
  
  shap_df <- as.data.frame(shap_values)
  
  message("  SHAP values calculated")
  
  list(
    shap = shap_df,
    X = X_sample,
    y = rf_result$y[idx]
  )
}

bgs_shap <- calc_shap(bgs_rf, "BGS")
ems_shap <- calc_shap(ems_rf, "EMS")

# ============================================================
# SHAP PLOTS
# ============================================================

message("\n--- Generating SHAP plots ---")

create_shap_plots <- function(shap_result, site_name, output_dir) {
  
  shap_df <- shap_result$shap
  X <- shap_result$X
  
  # 1. Mean |SHAP| bar plot
  mean_shap <- colMeans(abs(shap_df))
  shap_imp <- tibble(
    variable = names(mean_shap),
    mean_abs_shap = mean_shap
  ) %>%
    arrange(desc(mean_abs_shap))
  
  p_bar <- ggplot(shap_imp %>% head(20),
                  aes(x = mean_abs_shap, y = reorder(variable, mean_abs_shap))) +
    geom_col(fill = ifelse(site_name == "BGS", "#1f77b4", "#ff7f0e")) +
    labs(
      title = paste("Mean |SHAP| -", site_name),
      x = "Mean |SHAP value|",
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  
  ggsave(file.path(output_dir, paste0("shap_bar_", site_name, ".png")),
         p_bar, width = 8, height = 8, dpi = 300)
  
  # 2. Beeswarm-style summary plot
  top_vars <- shap_imp$variable[1:min(20, nrow(shap_imp))]
  
  # Prepare long format data for beeswarm
  shap_long <- shap_df %>%
    select(all_of(top_vars)) %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(-row_id, names_to = "variable", values_to = "shap_value")
  
  # Add scaled feature values
  X_scaled <- as.data.frame(scale(X[, top_vars]))
  X_long <- X_scaled %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(-row_id, names_to = "variable", values_to = "feature_value")
  
  shap_long <- shap_long %>%
    left_join(X_long, by = c("row_id", "variable")) %>%
    mutate(variable = factor(variable, levels = rev(top_vars)))
  
  p_summary <- ggplot(shap_long, aes(x = shap_value, y = variable, color = feature_value)) +
    geom_jitter(height = 0.2, width = 0, alpha = 0.6, size = 1) +
    scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                          midpoint = 0, name = "Feature\nvalue\n(scaled)") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      title = paste("SHAP Summary -", site_name),
      x = "SHAP value (impact on CH4 flux prediction)",
      y = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave(file.path(output_dir, paste0("shap_summary_", site_name, ".png")),
         p_summary, width = 10, height = 10, dpi = 300)
  
  # 3. Dependence plots for top 9 variables
  top9 <- shap_imp$variable[1:min(9, nrow(shap_imp))]
  
  dep_plots <- list()
  for (var in top9) {
    plot_df <- tibble(
      feature = X[[var]],
      shap = shap_df[[var]]
    )
    
    p <- ggplot(plot_df, aes(x = feature, y = shap)) +
      geom_point(alpha = 0.5, size = 1, color = ifelse(site_name == "BGS", "#1f77b4", "#ff7f0e")) +
      geom_smooth(method = "loess", se = FALSE, color = "black", size = 0.8) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
      labs(title = var, x = "Feature value", y = "SHAP value") +
      theme_minimal(base_size = 9)
    
    dep_plots[[var]] <- p
  }
  
  p_dep_combined <- wrap_plots(dep_plots, ncol = 3) +
    plot_annotation(
      title = paste("SHAP Dependence Plots -", site_name),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
  
  ggsave(file.path(output_dir, paste0("shap_dependence_", site_name, ".png")),
         p_dep_combined, width = 12, height = 10, dpi = 300)
  
  message("  Saved SHAP plots for ", site_name)
  
  return(shap_imp)
}

bgs_shap_imp <- create_shap_plots(bgs_shap, "BGS", OUTPUT_DIR)
ems_shap_imp <- create_shap_plots(ems_shap, "EMS", OUTPUT_DIR)

# ============================================================
# COMBINED SHAP IMPORTANCE COMPARISON
# ============================================================

shap_imp_combined <- bind_rows(
  bgs_shap_imp %>% mutate(site = "BGS"),
  ems_shap_imp %>% mutate(site = "EMS")
)

top_shap_vars <- shap_imp_combined %>%
  group_by(variable) %>%
  summarize(max_shap = max(mean_abs_shap)) %>%
  slice_max(max_shap, n = 20) %>%
  pull(variable)

p_shap_compare <- ggplot(shap_imp_combined %>% filter(variable %in% top_shap_vars),
                         aes(x = mean_abs_shap, y = reorder(variable, mean_abs_shap), fill = site)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("BGS" = "#1f77b4", "EMS" = "#ff7f0e")) +
  labs(
    title = "SHAP Variable Importance Comparison",
    subtitle = "Variables at optimal integration windows",
    x = "Mean |SHAP value|",
    y = NULL,
    fill = "Site"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(OUTPUT_DIR, "shap_importance_comparison.png"), p_shap_compare,
       width = 10, height = 10, dpi = 300)

# ============================================================
# SAVE RESULTS
# ============================================================

write_csv(var_imp_combined, file.path(OUTPUT_DIR, "variable_importance.csv"))
write_csv(shap_imp_combined, file.path(OUTPUT_DIR, "shap_importance.csv"))

# Save window info used
windows_used <- bind_rows(
  optimal_windows %>% filter(site == "BGS") %>% mutate(site = "BGS"),
  optimal_windows %>% filter(site == "EMS") %>% mutate(site = "EMS")
) %>%
  select(site, variable, window_days, r)

write_csv(windows_used, file.path(OUTPUT_DIR, "windows_used.csv"))

message("\n============================================================")
message("DONE")
message("============================================================")
message("Outputs saved to: ", OUTPUT_DIR)

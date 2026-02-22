# ============================================================
# 03_base_model_correlations.R
#
# APPROACH B: Correlations on residuals from base model
#
# 1. Find optimal windows for core drivers (T, WTD/SWC)
# 2. Fit base model: T × WTD × species (Wetland) or T × SWC × species (Upland)
# 3. Correlate RESIDUALS with all other predictors at rolling windows
# 4. Check for remaining temporal structure in residuals
#
# This approach asks: "After accounting for temperature, moisture,
# and species, what else explains CH4 variance?"
#
# Inputs:
#   - data/processed/flux_cleaned.csv
#   - data/processed/aligned_hourly_dataset.csv
#
# Outputs:
#   - results/approach_b/base_model_{site}.rds
#   - results/approach_b/best_windows.csv
#   - figures/03_base_model/residual_correlations_{site}.png
#   - figures/03_base_model/temporal_structure_{site}.png
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(RcppRoll)
  library(lme4)
  library(patchwork)
})

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  flux = "data/processed/flux_cleaned.csv",
  env = "data/processed/aligned_hourly_dataset.csv"
)

OUTPUT_DIR <- "figures/03_base_model"
RESULTS_DIR <- "results/approach_b"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Core drivers for base model (site-specific)
CORE_DRIVERS <- list(
  Wetland = c("TS_Ha2", "bvs_wtd_cm"),
  Upland = c("TS_Ha2", "SWC_Ha2")
)

# Candidate variables to test on residuals (everything else)
CANDIDATE_VARS <- c(
  "LE_Ha1", "LE_Ha2", "FC_Ha1", "FC_Ha2", "H_Ha1", "H_Ha2",
  "VPD_kPa", "PAR", "rnet", "gcc", "ndvi",
  "tair_C", "RH", "USTAR_Ha1", "USTAR_Ha2",
  "P_mm", "THROUGHFALL_xHA"
)

# Windows for core driver optimization
CORE_WINDOWS <- c(3, 6, 12, 24, 48, 72, 120, 168, 240, 336)

# Windows for candidate variable correlations
CANDIDATE_WINDOWS <- c(6, 12, 24, 48, 72, 120, 168, 240, 336)

# FDR threshold
FDR_ALPHA <- 0.05

# Variables to sum
SUM_VARS <- c("P_mm", "THROUGHFALL_xHA")

# ============================================================
# HELPER FUNCTIONS
# ============================================================

roll_mean <- function(x, n) {
  if (n == 1) return(x)
  RcppRoll::roll_mean(x, n = n, align = "right", fill = NA, na.rm = TRUE)
}

roll_sum <- function(x, n) {
  if (n == 1) return(x)
  RcppRoll::roll_sum(x, n = n, align = "right", fill = NA, na.rm = TRUE)
}

# ============================================================
# LOAD DATA
# ============================================================

message("Loading data...")

flux <- read_csv(PATHS$flux, show_col_types = FALSE) %>%
  mutate(
    datetime = as.POSIXct(datetime, tz = "UTC"),
    site = factor(site, levels = c("Wetland", "Upland")),
    Tree = as.factor(Tree),
    species = as.factor(species)
  )

env <- read_csv(PATHS$env, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC")) %>%
  arrange(datetime)

message("  Flux observations: ", nrow(flux))
message("  Environmental records: ", nrow(env))

# ============================================================
# PROCESS EACH SITE
# ============================================================

results_all <- list()
base_models <- list()

for (site_name in c("Wetland", "Upland")) {
  
  message("\n", paste(rep("=", 60), collapse = ""))
  message("  ", site_name)
  message(paste(rep("=", 60), collapse = ""))
  
  site_flux <- flux %>% filter(site == site_name)
  core_vars <- CORE_DRIVERS[[site_name]]
  
  message("\n  Observations: ", nrow(site_flux))
  message("  Core drivers: ", paste(core_vars, collapse = ", "))
  
  # ============================================================
  # STEP 1: Find optimal windows for core drivers
  # ============================================================
  
  message("\n  --- Step 1: Optimizing core driver windows ---")
  
  # Test all window combinations for core drivers
  # Use AIC of the full interaction model as criterion
  
  window_combos <- expand.grid(
    w1 = CORE_WINDOWS,
    w2 = CORE_WINDOWS
  )
  
  combo_results <- list()
  
  for (i in 1:nrow(window_combos)) {
    w1 <- window_combos$w1[i]
    w2 <- window_combos$w2[i]
    
    # Build features at these windows
    env_roll <- env %>% dplyr::select(datetime)
    env_roll[[paste0(core_vars[1], "_", w1, "h")]] <- roll_mean(env[[core_vars[1]]], w1)
    env_roll[[paste0(core_vars[2], "_", w2, "h")]] <- roll_mean(env[[core_vars[2]]], w2)
    
    # Join and prepare model data
    model_data <- site_flux %>%
      left_join(env_roll, by = "datetime") %>%
      drop_na()
    
    if (nrow(model_data) < 100) next
    
    # Scale predictors
    v1_name <- paste0(core_vars[1], "_", w1, "h")
    v2_name <- paste0(core_vars[2], "_", w2, "h")
    
    model_data <- model_data %>%
      mutate(
        v1_scaled = scale(.data[[v1_name]])[,1],
        v2_scaled = scale(.data[[v2_name]])[,1]
      )
    
    # Fit model with 3-way interaction
    tryCatch({
      m <- lmer(CH4_flux_asinh ~ v1_scaled * v2_scaled * species + (1|Tree),
                data = model_data, REML = FALSE)
      
      var_total <- var(model_data$CH4_flux_asinh)
      r2 <- var(predict(m, re.form = NA)) / var_total
      
      combo_results[[length(combo_results) + 1]] <- tibble(
        w1 = w1, w2 = w2,
        aic = AIC(m),
        bic = BIC(m),
        r2 = r2,
        n = nrow(model_data)
      )
    }, error = function(e) NULL)
  }
  
  combo_df <- bind_rows(combo_results) %>%
    arrange(aic)
  
  best_combo <- combo_df[1, ]
  
  cat("\n  Best window combination (by AIC):\n")
  cat("    ", core_vars[1], ": ", best_combo$w1, "h\n")
  cat("    ", core_vars[2], ": ", best_combo$w2, "h\n")
  cat("    AIC: ", round(best_combo$aic, 2), "\n")
  cat("    R²: ", round(best_combo$r2, 4), "\n")
  
  # ============================================================
  # STEP 2: Fit base model at optimal windows
  # ============================================================
  
  message("\n  --- Step 2: Fitting base model ---")
  
  w1 <- best_combo$w1
  w2 <- best_combo$w2
  v1_name <- paste0(core_vars[1], "_", w1, "h")
  v2_name <- paste0(core_vars[2], "_", w2, "h")
  
  # Build features
  env_roll <- env %>% dplyr::select(datetime)
  env_roll[[v1_name]] <- roll_mean(env[[core_vars[1]]], w1)
  env_roll[[v2_name]] <- roll_mean(env[[core_vars[2]]], w2)
  
  model_data <- site_flux %>%
    left_join(env_roll, by = "datetime") %>%
    drop_na()
  
  # Scale and store scaling params
  v1_mean <- mean(model_data[[v1_name]], na.rm = TRUE)
  v1_sd <- sd(model_data[[v1_name]], na.rm = TRUE)
  v2_mean <- mean(model_data[[v2_name]], na.rm = TRUE)
  v2_sd <- sd(model_data[[v2_name]], na.rm = TRUE)
  
  model_data <- model_data %>%
    mutate(
      v1_scaled = (.data[[v1_name]] - v1_mean) / v1_sd,
      v2_scaled = (.data[[v2_name]] - v2_mean) / v2_sd
    )
  
  # Fit base model
  m_base <- lmer(CH4_flux_asinh ~ v1_scaled * v2_scaled * species + (1|Tree),
                 data = model_data, REML = FALSE)
  
  var_total <- var(model_data$CH4_flux_asinh)
  r2_base <- var(predict(m_base, re.form = NA)) / var_total
  
  cat("\n  Base model summary:\n")
  cat("    Formula: CH4 ~ ", core_vars[1], " * ", core_vars[2], " * species + (1|Tree)\n")
  cat("    R² (fixed): ", round(r2_base, 4), "\n")
  cat("    AIC: ", round(AIC(m_base), 2), "\n")
  cat("    Residual variance: ", round(var(residuals(m_base)), 4), "\n")
  
  # Extract residuals
  model_data$residuals <- residuals(m_base)
  
  # Save base model
  base_models[[site_name]] <- list(
    model = m_base,
    core_windows = c(w1, w2),
    core_vars = core_vars,
    scaling = list(v1_mean = v1_mean, v1_sd = v1_sd, v2_mean = v2_mean, v2_sd = v2_sd),
    r2 = r2_base
  )
  
  saveRDS(base_models[[site_name]], file.path(RESULTS_DIR, paste0("base_model_", tolower(site_name), ".rds")))
  
  # ============================================================
  # STEP 3: Check temporal structure in residuals
  # ============================================================
  
  message("\n  --- Step 3: Checking residual temporal structure ---")
  
  model_data <- model_data %>%
    mutate(
      doy = yday(datetime),
      hour = hour(datetime),
      month = month(datetime)
    )
  
  # Month effect
  month_aov <- aov(residuals ~ factor(month), data = model_data)
  month_ss <- summary(month_aov)[[1]]["factor(month)", "Sum Sq"]
  total_ss <- sum(summary(month_aov)[[1]][, "Sum Sq"])
  month_r2 <- month_ss / total_ss
  month_p <- summary(month_aov)[[1]]["factor(month)", "Pr(>F)"]
  
  # Hour effect
  hour_aov <- aov(residuals ~ factor(hour), data = model_data)
  hour_ss <- summary(hour_aov)[[1]]["factor(hour)", "Sum Sq"]
  total_ss_h <- sum(summary(hour_aov)[[1]][, "Sum Sq"])
  hour_r2 <- hour_ss / total_ss_h
  hour_p <- summary(hour_aov)[[1]]["factor(hour)", "Pr(>F)"]
  
  cat("\n  Temporal structure in residuals:\n")
  cat("    Month R²: ", round(month_r2, 4), " (p = ", format.pval(month_p, digits = 2), ")\n")
  cat("    Hour R²:  ", round(hour_r2, 4), " (p = ", format.pval(hour_p, digits = 2), ")\n")
  
  if (month_r2 > 0.05 || hour_r2 > 0.05) {
    cat("    ⚠ WARNING: Substantial temporal structure remains (>5%)\n")
  } else {
    cat("    ✓ Minimal temporal structure\n")
  }
  
  # Plot temporal structure
  monthly_resid <- model_data %>%
    group_by(month) %>%
    summarize(mean_r = mean(residuals), se_r = sd(residuals)/sqrt(n()), .groups = "drop")
  
  hourly_resid <- model_data %>%
    group_by(hour) %>%
    summarize(mean_r = mean(residuals), se_r = sd(residuals)/sqrt(n()), .groups = "drop")
  
  p_month <- ggplot(monthly_resid, aes(x = month, y = mean_r)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(ymin = mean_r - 1.96*se_r, ymax = mean_r + 1.96*se_r)) +
    scale_x_continuous(breaks = 1:12, labels = month.abb) +
    labs(title = "Mean Residual by Month", 
         subtitle = paste0("R² = ", round(month_r2, 3)),
         x = "Month", y = "Mean residual") +
    theme_minimal()
  
  p_hour <- ggplot(hourly_resid, aes(x = hour, y = mean_r)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_pointrange(aes(ymin = mean_r - 1.96*se_r, ymax = mean_r + 1.96*se_r)) +
    scale_x_continuous(breaks = seq(0, 23, 3)) +
    labs(title = "Mean Residual by Hour",
         subtitle = paste0("R² = ", round(hour_r2, 3)),
         x = "Hour", y = "Mean residual") +
    theme_minimal()
  
  p_temporal <- p_month + p_hour +
    plot_annotation(title = paste0(site_name, " - Residual Temporal Structure"))
  
  ggsave(file.path(OUTPUT_DIR, paste0("temporal_structure_", tolower(site_name), ".png")),
         p_temporal, width = 10, height = 4, dpi = 300)
  
  # ============================================================
  # STEP 4: Correlate residuals with candidate variables
  # ============================================================
  
  message("\n  --- Step 4: Correlating residuals with candidates ---")
  
  # Build candidate features at all windows
  available_candidates <- intersect(CANDIDATE_VARS, names(env))
  
  site_results <- list()
  
  for (v in available_candidates) {
    for (w in CANDIDATE_WINDOWS) {
      
      # Calculate rolling feature
      if (v %in% SUM_VARS) {
        env_feat <- roll_sum(env[[v]], w)
      } else {
        env_feat <- roll_mean(env[[v]], w)
      }
      
      env_roll <- tibble(datetime = env$datetime, feat = env_feat)
      
      # Join with residuals
      resid_with_feat <- model_data %>%
        dplyr::select(datetime, residuals) %>%
        left_join(env_roll, by = "datetime")
      
      valid <- complete.cases(resid_with_feat$residuals, resid_with_feat$feat)
      n_valid <- sum(valid)
      
      if (n_valid < 30) next
      
      ct <- cor.test(resid_with_feat$residuals[valid], resid_with_feat$feat[valid])
      
      site_results[[length(site_results) + 1]] <- tibble(
        site = site_name,
        variable = v,
        window_hours = w,
        window_days = w / 24,
        n = n_valid,
        r = ct$estimate,
        p = ct$p.value
      )
    }
  }
  
  results_all[[site_name]] <- bind_rows(site_results)
  
  # Best window per variable
  site_best <- results_all[[site_name]] %>%
    group_by(variable) %>%
    slice_max(abs(r), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p, method = "BH"),
           significant = p_adj < FDR_ALPHA) %>%
    arrange(p_adj)
  
  n_sig <- sum(site_best$significant)
  
  cat("\n  Significant residual correlations: ", n_sig, "/", nrow(site_best), "\n")
  
  if (n_sig > 0) {
    cat("  Top significant:\n")
    site_best %>%
      filter(significant) %>%
      head(5) %>%
      mutate(r = round(r, 3), p_adj = signif(p_adj, 2)) %>%
      dplyr::select(variable, window_hours, r, p_adj) %>%
      print()
  }
  
  # Plot heatmap
  p_heatmap <- results_all[[site_name]] %>%
    ggplot(aes(x = window_hours, y = reorder(variable, abs(r)), fill = r)) +
    geom_tile() +
    geom_point(
      data = site_best %>% filter(significant),
      aes(x = window_hours, y = variable),
      shape = 8, size = 2, color = "black"
    ) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-0.3, 0.3),
      name = "r"
    ) +
    scale_x_continuous(
      breaks = c(6, 24, 72, 168, 336),
      labels = c("6h", "1d", "3d", "7d", "14d")
    ) +
    labs(
      title = paste0(site_name, " - Residual Correlations"),
      subtitle = paste0("After controlling for ", paste(core_vars, collapse = " × "), " × species"),
      x = "Window", y = NULL
    ) +
    theme_minimal() +
    theme(panel.grid = element_blank())
  
  ggsave(file.path(OUTPUT_DIR, paste0("residual_correlations_", tolower(site_name), ".png")),
         p_heatmap, width = 8, height = 6, dpi = 300)
}

# ============================================================
# COMBINE AND SAVE RESULTS
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("SAVING COMBINED RESULTS")
message(paste(rep("=", 60), collapse = ""))

# Combine all residual correlations
all_residual_cors <- bind_rows(results_all)

# Best windows
best_windows <- all_residual_cors %>%
  group_by(site, variable) %>%
  slice_max(abs(r), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(site) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  mutate(significant = p_adj < FDR_ALPHA)

# Add core driver info
core_info <- tibble(
  site = c("Wetland", "Wetland", "Upland", "Upland"),
  variable = c("TS_Ha2", "bvs_wtd_cm", "TS_Ha2", "SWC_Ha2"),
  window_hours = c(
    base_models$Wetland$core_windows[1], base_models$Wetland$core_windows[2],
    base_models$Upland$core_windows[1], base_models$Upland$core_windows[2]
  ),
  is_core = TRUE,
  r = NA_real_,  # Part of base model, not separately correlated
  p_adj = NA_real_,
  significant = TRUE  # By definition included
)

best_windows <- best_windows %>%
  mutate(is_core = FALSE) %>%
  bind_rows(core_info) %>%
  arrange(site, desc(is_core), p_adj)

write_csv(all_residual_cors, file.path(RESULTS_DIR, "residual_correlations_all.csv"))
write_csv(best_windows, file.path(RESULTS_DIR, "best_windows.csv"))

message("\nSaved:")
message("  ", file.path(RESULTS_DIR, "residual_correlations_all.csv"))
message("  ", file.path(RESULTS_DIR, "best_windows.csv"))

# ============================================================
# SUMMARY
# ============================================================

message("\n", paste(rep("=", 50), collapse = ""))
message("APPROACH B COMPLETE")
message(paste(rep("=", 50), collapse = ""))

cat("\nBase model R²:\n")
cat("  Wetland: ", round(base_models$Wetland$r2, 4), "\n")
cat("  Upland:  ", round(base_models$Upland$r2, 4), "\n")

cat("\nCore driver windows:\n")
cat("  Wetland: TS_Ha2 =", base_models$Wetland$core_windows[1], "h,",
    "bvs_wtd_cm =", base_models$Wetland$core_windows[2], "h\n")
cat("  Upland:  TS_Ha2 =", base_models$Upland$core_windows[1], "h,",
    "SWC_Ha2 =", base_models$Upland$core_windows[2], "h\n")

cat("\nThis approach first accounts for T × moisture × species,\n")
cat("then tests what else explains remaining variance.\n")
cat("Compare with Approach A (02_rolling_correlations.R).\n")

message("\nDone!")
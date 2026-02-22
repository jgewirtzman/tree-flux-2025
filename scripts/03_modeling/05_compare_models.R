# ============================================================
# 106_compare_ems_models.R
# 
# Compare EMS model variants:
#   A: Instantaneous (TS_Ha1 × NEON_SWC_shallow, 1h)
#   B: BGS drivers (TS_Ha2 × bvs_wtd_cm, 282h/132h)
# ============================================================

library(tidyverse)
library(lme4)

# ============================================================
# LOAD MODELS
# ============================================================

m_A <- readRDS("outputs/models/ems_instantaneous/m_final.rds")
m_B <- readRDS("outputs/models/ems_bgs_drivers/m_final.rds")

# ============================================================
# MODEL COMPARISON
# ============================================================

cat("══════════════════════════════════════════════════════════════\n")
cat("           EMS MODEL COMPARISON: A vs B\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# Extract key metrics
get_model_stats <- function(m, name) {
  r2 <- var(predict(m, re.form = NA)) / var(m@frame[[1]])
  tibble(
    model = name,
    R2 = round(r2 * 100, 1),
    AIC = round(AIC(m), 1),
    BIC = round(BIC(m), 1),
    N = nrow(m@frame),
    n_trees = n_distinct(m@frame$Tree),
    residual_SD = round(sigma(m), 4)
  )
}

comparison <- bind_rows(
  get_model_stats(m_A, "A: Instantaneous (TS_Ha1 × SWC)"),
  get_model_stats(m_B, "B: BGS drivers (TS_Ha2 × WTD)")
)

print(comparison, width = 100)

cat("\n──────────────────────────────────────────────────────────────\n")
cat("Model formulas:\n")
cat("──────────────────────────────────────────────────────────────\n")
cat("A:", deparse(formula(m_A)), "\n")
cat("B:", deparse(formula(m_B)), "\n")

# ============================================================
# SPECIES SLOPES COMPARISON
# ============================================================

cat("\n══════════════════════════════════════════════════════════════\n")
cat("           SPECIES-SPECIFIC SLOPES\n")
cat("══════════════════════════════════════════════════════════════\n\n")

slopes_A <- read_csv("outputs/models/ems_instantaneous/species_slopes.csv",
                     show_col_types = FALSE) %>%
  mutate(model = "A")

slopes_B <- read_csv("outputs/models/ems_bgs_drivers/species_slopes.csv",
                     show_col_types = FALSE) %>%
  rename(wtd_slope = any_of(c("wtd_slope", "swc_slope"))) %>%
  mutate(model = "B")

# Standardize column names
if (!"wtd_slope" %in% names(slopes_A)) {
  slopes_A <- slopes_A %>% rename(wtd_slope = swc_slope)
}

slopes_combined <- bind_rows(slopes_A, slopes_B) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

cat("Model A (Instantaneous):\n")
print(slopes_A %>% dplyr::select(-model))

cat("\nModel B (BGS drivers):\n")
print(slopes_B %>% dplyr::select(-model))

# ============================================================
# COEFFICIENT COMPARISON
# ============================================================

cat("\n══════════════════════════════════════════════════════════════\n")
cat("           FIXED EFFECTS SUMMARY\n")
cat("══════════════════════════════════════════════════════════════\n\n")

coef_A <- fixef(m_A)
coef_B <- fixef(m_B)

cat("Model A - significant effects (|t| > 2):\n")
summ_A <- summary(m_A)$coefficients
sig_A <- summ_A[abs(summ_A[, "t value"]) > 2, , drop = FALSE]
if (nrow(sig_A) > 0) {
  print(round(sig_A, 4))
} else {
  cat("  (none)\n")
}

cat("\nModel B - significant effects (|t| > 2):\n")
summ_B <- summary(m_B)$coefficients
sig_B <- summ_B[abs(summ_B[, "t value"]) > 2, , drop = FALSE]
if (nrow(sig_B) > 0) {
  print(round(sig_B, 4))
} else {
  cat("  (none)\n")
}

# ============================================================
# SUMMARY
# ============================================================

cat("\n══════════════════════════════════════════════════════════════\n")
cat("           SUMMARY\n")
cat("══════════════════════════════════════════════════════════════\n\n")

better <- if (comparison$R2[1] > comparison$R2[2]) "A" else "B"
r2_diff <- abs(comparison$R2[1] - comparison$R2[2])

cat(sprintf("Better fit (by R²): Model %s\n", better))
cat(sprintf("R² difference: %.1f percentage points\n", r2_diff))
cat(sprintf("AIC difference: %.1f (lower is better)\n", comparison$AIC[1] - comparison$AIC[2]))

cat("\n──────────────────────────────────────────────────────────────\n")
cat("Interpretation:\n")
cat("──────────────────────────────────────────────────────────────\n")
if (r2_diff < 2) {
  cat("→ Models explain similar variance - neither driver set is clearly better\n")
} else if (better == "A") {
  cat("→ Instantaneous conditions (TS_Ha1 × SWC) explain more variance\n")
  cat("→ Upland CH4 may respond to current rather than integrated conditions\n")
} else {
  cat("→ BGS drivers (TS_Ha2 × WTD) explain more variance\n
")
  cat("→ Same temporal integration as wetland applies to upland\n")
}

cat("\n══════════════════════════════════════════════════════════════\n")






# ============================================================
# UPLAND MODEL RESULTS - BOTH VERSIONS
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("UPLAND MODEL RESULTS - MODEL COMPARISON")
message(paste(rep("=", 60), collapse = ""))

# Determine which model objects exist
model_a_exists <- exists("m_final") && exists("r2_final")
model_b_exists <- exists("m_final_bgs") && exists("r2_final_bgs")

if (!model_a_exists && !model_b_exists) {
  stop("No model objects found. Run either ems_model_A.R or ems_model_B.R first.")
}

# 1. DATA SUMMARY
message("\n--- 1. DATA SUMMARY ---")
message("N observations: ", nrow(model_data_scaled))
message("N trees: ", n_distinct(model_data_scaled$Tree))
message("N species: ", n_distinct(model_data_scaled$species))
message("Date range: ", min(model_data_scaled$datetime), " to ", max(model_data_scaled$datetime))

# 2. MODEL A (INSTANTANEOUS) - if exists
if (model_a_exists) {
  message("\n", paste(rep("=", 60), collapse = ""))
  message("MODEL A: INSTANTANEOUS PREDICTORS")
  message(paste(rep("=", 60), collapse = ""))
  
  message("\nCore predictors:")
  message("  Temperature: ", CORE_TS_VAR, " (1h window)")
  message("  Soil moisture: ", CORE_SWC_VAR, " (1h window)")
  
  message("\nModel formula:")
  print(formula(m_final))
  
  message("\nModel fit:")
  message("  R² (marginal): ", round(r2_final * 100, 1), "%")
  message("  AIC: ", round(AIC(m_final), 1))
  message("  BIC: ", round(BIC(m_final), 1))
  message("  N parameters: ", length(fixef(m_final)))
  
  message("\nFixed effects (standardized):")
  fixef_a <- as.data.frame(summary(m_final)$coefficients)
  fixef_a$term <- rownames(fixef_a)
  fixef_a <- fixef_a %>%
    mutate(
      p_value = 2 * (1 - pt(abs(`t value`), df = nrow(model_data_scaled) - length(fixef(m_final)))),
      sig = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(term, Estimate, `Std. Error`, `t value`, p_value, sig)
  print(fixef_a, digits = 3)
  
  message("\nVariance components:")
  vc_a <- as.data.frame(VarCorr(m_final))
  print(vc_a)
  icc_a <- vc_a$vcov[1] / sum(vc_a$vcov)
  message("ICC: ", round(icc_a, 3))
  
  message("\nVIF:")
  print(vif(m_final))
  
  message("\nSpecies-specific effects:")
  coefs_a <- fixef(m_final)
  
  # Get predictor names dynamically
  ts_pred_a <- CORE_TS_VAR
  swc_pred_a <- CORE_SWC_VAR
  
  message("\nTemperature effects:")
  temp_hem_a <- coefs_a[paste0(ts_pred_a, "_raw_1h")] + 
    coefs_a[paste0(ts_pred_a, "_raw_1h:specieshem")]
  temp_rm_a <- coefs_a[paste0(ts_pred_a, "_raw_1h")] + 
    coefs_a[paste0(ts_pred_a, "_raw_1h:speciesrm")]
  temp_ro_a <- coefs_a[paste0(ts_pred_a, "_raw_1h")]
  
  message("  Q. rubra: ", round(temp_ro_a, 3))
  message("  T. canadensis: ", round(temp_hem_a, 3))
  message("  A. rubrum: ", round(temp_rm_a, 3))
  
  message("\nSoil moisture effects:")
  swc_hem_a <- coefs_a[paste0(swc_pred_a, "_raw_1h")] + 
    coefs_a[paste0(swc_pred_a, "_raw_1h:specieshem")]
  swc_rm_a <- coefs_a[paste0(swc_pred_a, "_raw_1h")] + 
    coefs_a[paste0(swc_pred_a, "_raw_1h:speciesrm")]
  swc_ro_a <- coefs_a[paste0(swc_pred_a, "_raw_1h")]
  
  message("  Q. rubra: ", round(swc_ro_a, 3))
  message("  T. canadensis: ", round(swc_hem_a, 3))
  message("  A. rubrum: ", round(swc_rm_a, 3))
  
  message("\nTemp × SWC interactions:")
  int_hem_a <- coefs_a[paste0(ts_pred_a, "_raw_1h:", swc_pred_a, "_raw_1h")] + 
    coefs_a[paste0(ts_pred_a, "_raw_1h:", swc_pred_a, "_raw_1h:specieshem")]
  int_rm_a <- coefs_a[paste0(ts_pred_a, "_raw_1h:", swc_pred_a, "_raw_1h")] + 
    coefs_a[paste0(ts_pred_a, "_raw_1h:", swc_pred_a, "_raw_1h:speciesrm")]
  int_ro_a <- coefs_a[paste0(ts_pred_a, "_raw_1h:", swc_pred_a, "_raw_1h")]
  
  message("  Q. rubra: ", round(int_ro_a, 3))
  message("  T. canadensis: ", round(int_hem_a, 3))
  message("  A. rubrum: ", round(int_rm_a, 3))
}

# 3. MODEL B (BGS DRIVERS) - if exists  
if (model_b_exists) {
  message("\n", paste(rep("=", 60), collapse = ""))
  message("MODEL B: BGS WETLAND DRIVERS (for comparison)")
  message(paste(rep("=", 60), collapse = ""))
  
  message("\nCore predictors (same as wetland):")
  message("  Temperature: TS_Ha2 (282h window)")
  message("  Water table: bvs_wtd_cm (132h window)")
  
  message("\nModel formula:")
  print(formula(m_final_bgs))
  
  message("\nModel fit:")
  message("  R² (marginal): ", round(r2_final_bgs * 100, 1), "%")
  message("  AIC: ", round(AIC(m_final_bgs), 1))
  message("  BIC: ", round(BIC(m_final_bgs), 1))
  message("  N parameters: ", length(fixef(m_final_bgs)))
  
  message("\nFixed effects (standardized):")
  fixef_b <- as.data.frame(summary(m_final_bgs)$coefficients)
  fixef_b$term <- rownames(fixef_b)
  fixef_b <- fixef_b %>%
    mutate(
      p_value = 2 * (1 - pt(abs(`t value`), df = nrow(model_data_scaled) - length(fixef(m_final_bgs)))),
      sig = case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(term, Estimate, `Std. Error`, `t value`, p_value, sig)
  print(fixef_b, digits = 3)
  
  message("\nVariance components:")
  vc_b <- as.data.frame(VarCorr(m_final_bgs))
  print(vc_b)
  icc_b <- vc_b$vcov[1] / sum(vc_b$vcov)
  message("ICC: ", round(icc_b, 3))
}

# 4. MODEL COMPARISON
if (model_a_exists && model_b_exists) {
  message("\n", paste(rep("=", 60), collapse = ""))
  message("MODEL COMPARISON")
  message(paste(rep("=", 60), collapse = ""))
  
  cat(sprintf("%-40s %6s %8s %8s\n", "Model", "R²", "AIC", "BIC"))
  cat("────────────────────────────────────────────────────\n")
  cat(sprintf("%-40s %5.1f%% %8.1f %8.1f\n", "Model A: Instantaneous (1h windows)", 
              r2_final*100, AIC(m_final), BIC(m_final)))
  cat(sprintf("%-40s %5.1f%% %8.1f %8.1f\n", "Model B: BGS drivers (282h, 132h)", 
              r2_final_bgs*100, AIC(m_final_bgs), BIC(m_final_bgs)))
  cat("────────────────────────────────────────────────────\n")
  cat(sprintf("%-40s %+5.1f%% %+8.1f %+8.1f\n", "Difference (B - A)", 
              (r2_final_bgs - r2_final)*100,
              AIC(m_final_bgs) - AIC(m_final),
              BIC(m_final_bgs) - BIC(m_final)))
}
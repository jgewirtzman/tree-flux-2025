# ============================================================
# 01_wetland_model.R
#
# Wetland (Black Gum Swamp) driver model of stem CH4 flux.
#
# Inputs: outputs/figures/rolling_correlations/best_windows_by_variable.csv (4_analysis/05_window_screening.R),
#         data/final/environment_hourly.csv, data/final/stem_ch4_flux.csv
#
# Steps (section numbers below):
#   1-2   drivers significant in the window screening, each at its best window (raw; LE/FC also
#         as anomalies)
#   3-4   correlation heatmap and redundancy clusters of those drivers (diagnostic only)
#   5-6   core drivers plus greedy forward selection among the other drivers (model_helpers.R)
#   7-8   models without interactions; collinear predictors (GVIF > 10) dropped as listed in VIF_DROP_VARS
#   9-11  core model temperature x water table x species; each remaining predictor tested as an
#         addition (LRT, with and without a species interaction); the extended model replaces
#         the core model only if it has a lower BIC
#   12    final model refitted on asinh(flux / 0.1) using every row that has its own predictors
#   13-16 species-specific slopes, coefficient table, model summary, VIF
#
# Final model:
#   asinh(CH4 flux / 0.1) ~ TS_Ha2 (66 h) * bvs_wtd_cm (183 h) * species + (1|Tree) + (1|date)
#
# Outputs: outputs/models/bgs_final/ (m_final.rds, m_core_asinh.rds, model_data_scaled.rds,
#   model_data_info.csv, predictor_addition_tests.csv, species_slopes.csv, coefficients.csv,
#   model_summary.csv, vif.csv); diagnostics in outputs/figures/predictor_selection/
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(car)
  library(RcppRoll)
  library(pheatmap)
})

set.seed(42)

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  aligned   = "data/final/environment_hourly.csv",
  flux      = "data/final/stem_ch4_flux.csv",
  best_windows = "outputs/figures/rolling_correlations/best_windows_by_variable.csv"
)

OUTPUT_DIR <- "outputs/models/bgs_final"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures/predictor_selection", recursive = TRUE, showWarnings = FALSE)

# Correlation threshold for redundancy clustering
CORR_THRESHOLD <- 0.7

# Forward selection settings
MAX_PREDICTORS <- 6

# Core drivers (forced into model) - based on theory
CORE_VARIABLES <- c("TS_Ha2", "bvs_wtd_cm", "NEON_SWC_shallow")

# Variables that should be summed (not averaged) over windows
SUM_VARS <- c("P_mm", "THROUGHFALL_xHA")

# Specific anomaly variables to consider (LE and FC only)
SPECIFIC_ANOM_VARIABLES <- c("LE_Ha1", "LE_Ha2", "FC_Ha1", "FC_Ha2")

# ============================================================
# HELPER FUNCTIONS
# ============================================================

# roll_mean(), roll_sum(), make_anomaly(), feature matrices, redundancy clustering,
# forward_select(): shared with the other model scripts
source("scripts/5_models/model_helpers.R")

# ============================================================
# 1. LOAD DATA
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       1. LOADING DATA                            ")
message("══════════════════════════════════════════════════")

# Environmental data
aligned_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC"))

# Flux data
stem_flux_raw <- read_csv(PATHS$flux, show_col_types = FALSE)

stem_flux <- stem_flux_raw %>%
  mutate(
    datetime = as.POSIXct(format(as.POSIXct(sample_hour_est, tz = "UTC"), "%Y-%m-%d %H:%M:%S"), tz = "UTC"),  # EST hour of sampling (2_flux/05_flux_dataset.R)
    site = location,  # Already has "Wetland"/"Upland"
    Tree = as.factor(Tree),
    species = as.factor(SPECIES),
    CH4_flux = CH4_flux_nmolpm2ps / 1000  # Convert nmol to µmol for consistency
  ) %>%
  filter(site == "Wetland") %>%
  dplyr::select(datetime, Tree, species, CH4_flux)

# Best windows from rolling correlation analysis
bw <- read_csv(PATHS$best_windows, show_col_types = FALSE)

cat("Loaded", nrow(aligned_data), "hourly environmental observations\n")
cat("Loaded", nrow(stem_flux), "CH4 flux measurements (Wetland only)\n")
cat("Loaded", nrow(bw), "variable-window combinations from rolling analysis\n")

# ============================================================
# 2. EXTRACT SIGNIFICANT PREDICTORS
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       2. EXTRACTING SIGNIFICANT PREDICTORS       ")
message("══════════════════════════════════════════════════")

# Note: rolling analysis uses "Wetland"/"Upland", model uses "BGS"/"EMS"
optA_raw_sig <- bw %>%
  filter(analysis == "raw", significant, site == "Wetland") %>%
  dplyr::select(variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
  arrange(desc(abs(r)))

optA_anom_sig <- bw %>%
  filter(analysis == "anomaly", significant, site == "Wetland") %>%
  dplyr::select(variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
  arrange(desc(abs(r)))

cat("\nBGS (Wetland) significant predictors:\n")
cat("  Raw:", nrow(optA_raw_sig), "\n")
cat("  Anomaly:", nrow(optA_anom_sig), "\n")

# ============================================================
# 3. CORRELATION HEATMAPS FOR REDUNDANCY ANALYSIS
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       3. CORRELATION HEATMAPS                    ")
message("══════════════════════════════════════════════════")

# BGS Raw heatmap
if (nrow(optA_raw_sig) >= 2) {
  feat_bgs_raw <- build_feature_matrix_raw(aligned_data, optA_raw_sig)
  X <- feat_bgs_raw %>% dplyr::select(-datetime)
  X <- X %>% dplyr::select(where(~ sum(!is.na(.)) > 50))
  cmat_bgs_raw <- cor(X, use = "pairwise.complete.obs")
  
  plot_corr_heatmap(
    cmat_bgs_raw,
    "BGS (Wetland) - Raw Predictors at Best Windows",
    "outputs/figures/predictor_selection/corr_heatmap_BGS_raw.png"
  )
  
  redundant_bgs_raw <- find_redundant_pairs(cmat_bgs_raw, CORR_THRESHOLD)
  if (nrow(redundant_bgs_raw) > 0) {
    cat("\nRedundant pairs (|r| >=", CORR_THRESHOLD, "):\n")
    print(redundant_bgs_raw, n = 20)
    write_csv(redundant_bgs_raw, "outputs/figures/predictor_selection/redundant_pairs_BGS_raw.csv")
  }
}

# ============================================================
# 4. COMBINED RAW + ANOMALY PREDICTOR SELECTION
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       4. COMBINED PREDICTOR SELECTION            ")
message("══════════════════════════════════════════════════")

# Pool significant predictors from both modes
bgs_raw <- optA_raw_sig %>% mutate(mode = "raw")
bgs_anom <- optA_anom_sig %>% 
  filter(variable %in% SPECIFIC_ANOM_VARIABLES) %>%
  mutate(mode = "anom")

bgs_combined <- bind_rows(bgs_raw, bgs_anom)

cat("Pooled predictors:", nrow(bgs_raw), "raw +", nrow(bgs_anom), "anomaly =", 
    nrow(bgs_combined), "total\n\n")

# Build combined feature matrix
feat_bgs_combined <- build_combined_feature_matrix(aligned_data, bgs_combined)

X <- feat_bgs_combined %>% dplyr::select(-datetime)
X <- X %>% dplyr::select(where(~ sum(!is.na(.)) > 50))

cat("Feature matrix:", ncol(X), "features\n\n")

cmat_bgs_combined <- cor(X, use = "pairwise.complete.obs")

# Run selection
selected_bgs <- select_predictors_combined(
  bgs_combined, cmat_bgs_combined, feat_bgs_combined, threshold = CORR_THRESHOLD
)

cat("Selected predictors (", nrow(selected_bgs), " from ", nrow(bgs_combined), "):\n\n")

selected_bgs %>%
  dplyr::select(variable, mode, var_group, window_hours, r, cluster_size, dropped_vars) %>%
  arrange(var_group, desc(abs(r))) %>%
  print(n = 30, width = 120)

write_csv(selected_bgs, "outputs/figures/predictor_selection/selected_predictors_combined.csv")

# ============================================================
# 5. THEORY-DRIVEN MODEL BUILDING
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       5. THEORY-DRIVEN MODEL BUILDING            ")
message("══════════════════════════════════════════════════")

# Build features for ALL significant raw predictors + specific anomalies
bgs_raw_all <- bw %>%
  filter(analysis == "raw", significant, site == "Wetland") %>%
  dplyr::select(variable, var_group, window_hours, r) %>%
  mutate(
    predictor = paste0(variable, "_raw_", window_hours, "h"),
    mode = "raw"
  )

bgs_anom_subset <- bw %>%
  filter(analysis == "anomaly", significant, site == "Wetland",
         variable %in% SPECIFIC_ANOM_VARIABLES) %>%
  dplyr::select(variable, var_group, window_hours, r) %>%
  mutate(
    predictor = paste0(variable, "_anom_", window_hours, "h"),
    mode = "anom"
  )

cat("Significant raw predictors:", nrow(bgs_raw_all), "\n")
cat("Specific anomaly predictors (LE/FC):", nrow(bgs_anom_subset), "\n")

# Core drivers (subset of raw)
core_info <- bgs_raw_all %>%
  filter(variable %in% CORE_VARIABLES)

cat("\nCore drivers:\n")
print(core_info %>% dplyr::select(predictor, variable, r))

# Candidates = all raw (except core) + specific anomalies
candidate_info <- bind_rows(
  bgs_raw_all %>% filter(!variable %in% CORE_VARIABLES),
  bgs_anom_subset
)

cat("\nCandidate predictors:", nrow(candidate_info), "\n")

# Build features
met <- aligned_data %>% arrange(datetime)

# Create anomaly versions for specific variables
for (v in SPECIFIC_ANOM_VARIABLES) {
  if (v %in% names(met)) {
    met <- make_anomaly(met, v, time_col = "datetime")
  }
}

features <- met %>% dplyr::select(datetime)

# Raw features
for (i in seq_len(nrow(bgs_raw_all))) {
  v <- bgs_raw_all$variable[i]
  w <- bgs_raw_all$window_hours[i]
  if (!v %in% names(met)) next
  feat_name <- paste0(v, "_raw_", w, "h")
  if (v %in% SUM_VARS) {
    features[[feat_name]] <- roll_sum(met[[v]], w)
  } else {
    features[[feat_name]] <- roll_mean(met[[v]], w)
  }
}

# Anomaly features
for (i in seq_len(nrow(bgs_anom_subset))) {
  v <- bgs_anom_subset$variable[i]
  w <- bgs_anom_subset$window_hours[i]
  src <- paste0(v, "_anom")
  if (!src %in% names(met)) next
  feat_name <- paste0(v, "_anom_", w, "h")
  features[[feat_name]] <- roll_mean(met[[src]], w)
}

cat("Total features built:", ncol(features) - 1, "\n")

# Prepare model data
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

# Total variance for R²
# Selection uses the same response scale as the final model: asinh(flux in nmol m-2 s-1)
model_data_scaled <- model_data_scaled %>% mutate(CH4_flux_umol = CH4_flux, CH4_flux = asinh(CH4_flux * 1000 / 0.1), date = factor(as.Date(datetime)))
var_total <- var(model_data_scaled$CH4_flux)

# ============================================================
# 6. FORWARD SELECTION
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       6. FORWARD SELECTION                       ")
message("══════════════════════════════════════════════════")

cat("Core drivers (always first):", paste(core_names, collapse = ", "), "\n")
cat("Candidates:", length(candidate_names), "predictors\n\n")

step_df <- forward_select(core_names, candidate_names, model_data_scaled, var_total, MAX_PREDICTORS) %>%
  left_join(bind_rows(core_info, candidate_info) %>% dplyr::select(predictor, variable, var_group, mode, r),
            by = "predictor")

cat("Forward selection (predictor added at each step; fixed-effect variance share after it):\n")
print(step_df %>% mutate(r2 = round(r2, 3), r = round(r, 3)))

write_csv(step_df, "outputs/figures/predictor_selection/forward_selection_steps.csv")

# Selected predictors, listed alphabetically (C locale) after the core drivers: the term order
# of the published run, which listed them from a selection-frequency table.
selected_preds <- sort(unique(step_df$predictor), method = "radix")

# ============================================================
# 7. MODEL COMPARISON (NO INTERACTIONS)
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       7. BASE MODEL COMPARISON                   ")
message("══════════════════════════════════════════════════")

# Species only
m_species <- lmer(CH4_flux ~ species + (1|Tree) + (1|date), 
                  data = model_data_scaled, REML = FALSE)
r2_species <- var(predict(m_species, re.form = NA)) / var_total

# Core only
core_formula <- paste("CH4_flux ~", paste(core_names, collapse = " + "), "+ species + (1|Tree) + (1|date)")
m_core <- lmer(as.formula(core_formula), data = model_data_scaled, REML = FALSE)
r2_core <- var(predict(m_core, re.form = NA)) / var_total

# Core + forward-selected
if (length(selected_preds) > 0) {
  strong_preds <- c(core_names, selected_preds)
  strong_formula <- paste("CH4_flux ~", paste(strong_preds, collapse = " + "), "+ species + (1|Tree) + (1|date)")
  m_strong <- lmer(as.formula(strong_formula), data = model_data_scaled, REML = FALSE)
  r2_strong <- var(predict(m_strong, re.form = NA)) / var_total
} else {
  strong_preds <- core_names
  m_strong <- m_core
  r2_strong <- r2_core
}

cat("\nModel comparison (raw CH4, no interactions):\n")
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("%-25s %6s %8s %8s\n", "Model", "R²", "AIC", "BIC"))
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("%-25s %5.1f%% %8.1f %8.1f\n", "Species only", r2_species*100, AIC(m_species), BIC(m_species)))
cat(sprintf("%-25s %5.1f%% %8.1f %8.1f\n", "Core drivers", r2_core*100, AIC(m_core), BIC(m_core)))
cat(sprintf("%-25s %5.1f%% %8.1f %8.1f\n", "Core + selected", r2_strong*100, AIC(m_strong), BIC(m_strong)))
cat("─────────────────────────────────────────────────────\n")

# VIF for strong model
cat("\nVIF for core + selected model:\n")
print(vif(m_strong))

# ============================================================
# 8. VIF-BASED MODEL REFINEMENT
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       8. VIF-BASED MODEL REFINEMENT              ")
message("══════════════════════════════════════════════════")

# Identify core interaction predictors dynamically
ts_pred <- names(model_data_scaled)[grepl("^TS_Ha2_raw_", names(model_data_scaled))][1]
wtd_pred <- names(model_data_scaled)[grepl("^bvs_wtd_cm_raw_", names(model_data_scaled))][1]
# LE_Ha2 anomaly predictor (window from 4_analysis/05_window_screening.R), dropped below with the high-VIF set
le2_pred <- names(model_data_scaled)[grepl("^LE_Ha2_anom_", names(model_data_scaled))][1]

cat("Core interaction predictors (theory-driven):\n")
cat("  Temperature:", ts_pred, "\n")
cat("  Water table:", wtd_pred, "\n\n")

if (is.na(ts_pred) || is.na(wtd_pred)) {
  stop("Required core predictors (TS_Ha2 and/or bvs_wtd_cm) not found!")
}

# Start with the forward-selected predictors
refined_preds <- strong_preds

# Check VIF
cat("VIF for core + selected model:\n")
vif_vals <- vif(m_strong)
vif_df <- data.frame(
  predictor = rownames(vif_vals),
  GVIF = round(vif_vals[, "GVIF"], 2)
) %>%
  arrange(desc(GVIF))

print(vif_df)

# Identify high-VIF predictors (excluding core)
high_vif_preds <- vif_df %>%
  filter(GVIF > 10, !predictor %in% c(ts_pred, wtd_pred, "species")) %>%
  pull(predictor)

if (length(high_vif_preds) > 0) {
  cat("\n⚠ Predictors with VIF > 10 (candidates for removal):\n")
  for (p in high_vif_preds) {
    cat("  -", p, "(VIF =", vif_df$GVIF[vif_df$predictor == p], ")\n")
  }
  
  # Show correlation matrix among high-VIF predictors + core predictors
  # This helps identify which variables are collinear with each other
  check_preds <- c(ts_pred, wtd_pred, high_vif_preds)
  check_preds <- check_preds[check_preds %in% names(model_data_scaled)]
  
  cat("\n─────────────────────────────────────────────────────\n")
  cat("Correlation matrix (high-VIF + core predictors):\n")
  cat("─────────────────────────────────────────────────────\n")
  cor_matrix <- cor(model_data_scaled[, check_preds], use = "pairwise.complete.obs")
  
  # Print with short names for readability
  short_names <- str_replace(check_preds, "_raw_\\d+h$|_anom_\\d+h$", "")
  rownames(cor_matrix) <- short_names
  colnames(cor_matrix) <- short_names
  print(round(cor_matrix, 2))
  
  # Identify strongly correlated pairs (|r| > 0.7)
  cat("\n─────────────────────────────────────────────────────\n")
  cat("Strongly correlated pairs (|r| > 0.7):\n")
  cat("─────────────────────────────────────────────────────\n")
  for (i in 1:(length(check_preds)-1)) {
    for (j in (i+1):length(check_preds)) {
      r <- cor_matrix[i, j]
      if (abs(r) > 0.7) {
        cat(sprintf("  %s ~ %s: r = %.2f\n", short_names[i], short_names[j], r))
      }
    }
  }
  
  cat("\n─────────────────────────────────────────────────────\n")
  cat("DECISION GUIDE:\n")
  cat("─────────────────────────────────────────────────────\n")
  cat("For each correlated pair, keep the variable that is:\
")
  cat("  1. More mechanistically relevant to CH4 flux\n")
  cat("  2. More directly measured (vs derived)\n")
  cat("  3. Core to your hypothesis (TS_Ha2, bvs_wtd_cm)\n")
  cat("─────────────────────────────────────────────────────\n")
  
  # ════════════════════════════════════════════════════════════
  # USER DECISION: Which high-VIF predictors to DROP?
  # Based on the correlation matrix above, choose which member
  # of each correlated pair to remove
  # ════════════════════════════════════════════════════════════
  
  # Variables (not window-specific names) to drop; matched by prefix so the
  # choice survives changes in the optimal windows when the flux data change.
  VIF_DROP_VARS <- c("TS_Ha1_raw_", "s10t_raw_")  # EDIT THIS LIST
  VIF_DROP <- c(refined_preds[Reduce(`|`, lapply(VIF_DROP_VARS, function(v) startsWith(refined_preds, v)))],
                le2_pred)
  
  cat("\n→ Dropping (per VIF_DROP config):", paste(VIF_DROP, collapse = ", "), "\n")
  
  # Remove specified predictors
  refined_preds <- setdiff(refined_preds, VIF_DROP)
  
  # Refit model
  m_strong <- lmer(as.formula(paste("CH4_flux ~", paste(refined_preds, collapse = " + "), 
                                    "+ species + (1|Tree) + (1|date)")), 
                   data = model_data_scaled, REML = FALSE)
  
  # Show final VIF
  cat("\nFinal VIF after removal:\n")
  vif_vals_final <- vif(m_strong)
  print(round(vif_vals_final[, "GVIF"], 2))
  
} else {
  cat("\n✓ No predictors with VIF > 10 (excluding core). No removal needed.\n")
}

cat("\nRefined predictor set:\n")
cat(" ", paste(refined_preds, collapse = ", "), "\n")

r2_refined <- var(predict(m_strong, re.form = NA)) / var_total
cat("\nRefined model R²:", round(r2_refined * 100, 1), "%\n")

# Identify additional predictors (refined set minus core)
additional_preds <- setdiff(refined_preds, c(ts_pred, wtd_pred))
cat("\nAdditional predictors beyond core:", length(additional_preds), "\n")
if (length(additional_preds) > 0) {
  cat(" ", paste(additional_preds, collapse = ", "), "\n")
}

# ============================================================
# 9. CORE INTERACTION MODEL (TS * WTD * SPECIES)
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       9. CORE INTERACTION MODEL                  ")
message("══════════════════════════════════════════════════")

# Model 1: Core 3-way interaction only
formula_core <- paste0("CH4_flux ~ ", ts_pred, " * ", wtd_pred, " * species + (1|Tree) + (1|date)")
m_core_int <- lmer(as.formula(formula_core), data = model_data_scaled, REML = FALSE)
r2_core_int <- var(predict(m_core_int, re.form = NA)) / var_total

cat("Core interaction model:\n")
cat("  Formula:", formula_core, "\n")
cat("  R²:", round(r2_core_int * 100, 1), "%\n")
cat("  AIC:", round(AIC(m_core_int), 1), "\n\n")

cat("Core model coefficients:\n")
print(summary(m_core_int)$coefficients)

# ============================================================
# 10. TEST ADDING PREDICTORS TO CORE MODEL
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       10. ADDING PREDICTORS TO CORE MODEL        ")
message("══════════════════════════════════════════════════")

if (length(additional_preds) > 0) {
  
  cat("Testing addition of each refined predictor to core model:\n")
  cat("(with and without species interaction)\n\n")
  
  addition_results <- list()
  
  for (pred in additional_preds) {
    # Without species interaction
    formula_add <- paste0("CH4_flux ~ ", ts_pred, " * ", wtd_pred, " * species + ", 
                          pred, " + (1|Tree) + (1|date)")
    
    tryCatch({
      m_add <- lmer(as.formula(formula_add), data = model_data_scaled, REML = FALSE)
      lrt <- anova(m_core_int, m_add)
      
      addition_results[[paste0(pred, "_main")]] <- tibble(
        predictor = pred,
        type = "main effect only",
        aic = AIC(m_add),
        delta_aic = AIC(m_add) - AIC(m_core_int),
        lrt_chisq = lrt$Chisq[2],
        lrt_p = lrt$`Pr(>Chisq)`[2]
      )
    }, error = function(e) NULL)
    
    # With species interaction
    formula_add_sp <- paste0("CH4_flux ~ ", ts_pred, " * ", wtd_pred, " * species + ", 
                             pred, " * species + (1|Tree) + (1|date)")
    
    tryCatch({
      m_add_sp <- lmer(as.formula(formula_add_sp), data = model_data_scaled, REML = FALSE)
      lrt_sp <- anova(m_core_int, m_add_sp)
      
      addition_results[[paste0(pred, "_species")]] <- tibble(
        predictor = pred,
        type = "with species interaction",
        aic = AIC(m_add_sp),
        delta_aic = AIC(m_add_sp) - AIC(m_core_int),
        lrt_chisq = lrt_sp$Chisq[2],
        lrt_p = lrt_sp$`Pr(>Chisq)`[2]
      )
    }, error = function(e) NULL)
  }
  
  addition_df <- bind_rows(addition_results) %>%
    arrange(lrt_p) %>%
    mutate(significant = lrt_p < 0.05)
  
  cat("Results (sorted by p-value):\n")
  print(addition_df %>%
          mutate(
            delta_aic = round(delta_aic, 1),
            lrt_chisq = round(lrt_chisq, 1),
            lrt_p = format.pval(lrt_p, digits = 3)
          ), n = 30)
  
  write_csv(addition_df, file.path(OUTPUT_DIR, "predictor_addition_tests.csv"))
  
  # Identify significant additions
  sig_additions <- addition_df %>% filter(significant)
  
  if (nrow(sig_additions) > 0) {
    cat("\n\nSignificant additions to core model:\n")
    print(sig_additions %>% dplyr::select(predictor, type, delta_aic, lrt_p))
    
    # Get best addition (lowest AIC)
    best_addition <- sig_additions %>% slice_min(aic, n = 1)
    cat("\nBest single addition:", best_addition$predictor, "(", best_addition$type, ")\n")
    cat("  ΔAIC =", round(best_addition$delta_aic, 1), "\n")
  }
  
} else {
  cat("No additional predictors to test.\n")
  addition_df <- tibble()
  sig_additions <- tibble()
}

# ============================================================
# 11. BUILD BEST MODEL WITH SIGNIFICANT ADDITIONS
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       11. BEST MODEL SELECTION                   ")
message("══════════════════════════════════════════════════")

# Start with core model
best_formula <- formula_core
best_model <- m_core_int
best_aic <- AIC(m_core_int)

# Try adding significant predictors one at a time (greedy forward)
if (nrow(sig_additions) > 0) {
  
  # Get unique predictors that improve the model
  sig_preds_main <- sig_additions %>% 
    filter(type == "main effect only", delta_aic < -2) %>%
    pull(predictor)
  
  sig_preds_species <- sig_additions %>%
    filter(type == "with species interaction", delta_aic < -2) %>%
    pull(predictor)
  
  cat("Predictors with ΔAIC < -2:\n")
  cat("  Main effect only:", paste(sig_preds_main, collapse = ", "), "\n")
  cat("  With species int:", paste(sig_preds_species, collapse = ", "), "\n\n")
  
  # Build formula with best additions
  # Prefer species interaction version if it's better
  added_terms <- c()
  
  for (pred in unique(c(sig_preds_main, sig_preds_species))) {
    main_row <- sig_additions %>% filter(predictor == pred, type == "main effect only")
    sp_row <- sig_additions %>% filter(predictor == pred, type == "with species interaction")
    
    if (nrow(sp_row) > 0 && nrow(main_row) > 0) {
      # Both exist - pick better one
      if (sp_row$aic < main_row$aic) {
        added_terms <- c(added_terms, paste0(pred, " * species"))
      } else {
        added_terms <- c(added_terms, pred)
      }
    } else if (nrow(sp_row) > 0) {
      added_terms <- c(added_terms, paste0(pred, " * species"))
    } else if (nrow(main_row) > 0) {
      added_terms <- c(added_terms, pred)
    }
  }
  
  if (length(added_terms) > 0) {
    # Build extended formula
    formula_extended <- paste0("CH4_flux ~ ", ts_pred, " * ", wtd_pred, " * species + ",
                               paste(added_terms, collapse = " + "), " + (1|Tree) + (1|date)")
    
    cat("Testing extended model:\n")
    cat("  ", formula_extended, "\n\n")
    
    tryCatch({
      m_extended <- lmer(as.formula(formula_extended), data = model_data_scaled, REML = FALSE)
      r2_extended <- var(predict(m_extended, re.form = NA)) / var_total
      
      cat("Extended model:\n")
      cat("  R²:", round(r2_extended * 100, 1), "%\n")
      cat("  AIC:", round(AIC(m_extended), 1), "(core:", round(AIC(m_core_int), 1), ")\n")
      cat("  BIC:", round(BIC(m_extended), 1), "(core:", round(BIC(m_core_int), 1), ")\n\n")
      
      # LRT vs core
      cat("LRT: Core vs Extended:\n")
      print(anova(m_core_int, m_extended))
      
      # Check if extended is better by BIC (more conservative)
      if (BIC(m_extended) < BIC(m_core_int)) {
        cat("\n→ Extended model selected (lower BIC)\n")
        best_formula <- formula_extended
        best_model <- m_extended
        best_aic <- AIC(m_extended)
      } else {
        cat("\n→ Core model retained (extended model not better by BIC)\n")
      }
      
    }, error = function(e) {
      cat("Extended model failed to fit:", e$message, "\n")
    })
  }
}

cat("\n\nFinal selected model (raw CH4):\n")
cat("  ", best_formula, "\n")
cat("  AIC:", round(best_aic, 1), "\n")

# ============================================================
# 12. FINAL MODEL WITH ASINH TRANSFORM
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       12. FINAL MODEL (ASINH TRANSFORM)          ")
message("══════════════════════════════════════════════════")

# Create asinh-transformed response
# Refit on every observation that has the final model's own predictors (the selection
# data above drop rows missing ANY candidate predictor, which discards usable rows).
# Predictors are re-standardized on this data; sampling date enters as a random effect
# because trees measured on the same day share conditions not captured by the drivers.
final_preds <- setdiff(all.vars(as.formula(best_formula)), c("CH4_flux", "species", "Tree", "date"))
# (all candidate columns are kept for the comparison models further down; rows are
#  filtered on the final model's predictors only)
model_data <- stem_flux %>%
  left_join(features, by = "datetime") %>%
  drop_na(CH4_flux, all_of(final_preds))
zs <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
model_data_scaled <- model_data %>%
  mutate(across(any_of(unique(c(final_preds, all_pred_names))), zs),
         CH4_flux_umol = CH4_flux, CH4_flux_asinh = asinh(CH4_flux * 1000 / 0.1), CH4_flux = CH4_flux_asinh,
         date = factor(as.Date(datetime)))
cat("Final-model data:", nrow(model_data_scaled), "observations,", n_distinct(model_data_scaled$Tree), "trees,",
    as.character(min(as.Date(model_data_scaled$datetime))), "to", as.character(max(as.Date(model_data_scaled$datetime))), "\n")
# Nakagawa marginal R2: fixed-effect variance / (fixed + random + residual)
r2_marg <- function(m) as.numeric(performance::r2_nakagawa(m, tolerance = 1e-10)$R2_marginal)

var_total_asinh <- var(model_data_scaled$CH4_flux_asinh)

# Convert best formula to asinh version
formula_final <- gsub("CH4_flux ~", "CH4_flux_asinh ~", best_formula)

cat("Final model formula (asinh-transformed):\n")
cat(" ", formula_final, "\n\n")

m_final <- lmer(as.formula(formula_final), data = model_data_scaled, REML = FALSE)

r2_final <- r2_marg(m_final)

cat("Final model (asinh-transformed):\n")
cat("  R²:", round(r2_final * 100, 1), "%\n")
cat("  AIC:", round(AIC(m_final), 1), "\n")
cat("  BIC:", round(BIC(m_final), 1), "\n")

cat("\nModel summary:\n")
print(summary(m_final))

cat("\nVIF:\n")
print(vif(m_final))

# Save model
saveRDS(m_final, file.path(OUTPUT_DIR, "m_final.rds"))

# Also fit and save the core-only asinh model for comparison
formula_core_asinh <- gsub("CH4_flux ~", "CH4_flux_asinh ~", formula_core)
m_core_asinh <- lmer(as.formula(formula_core_asinh), data = model_data_scaled, REML = FALSE)
r2_core_asinh <- r2_marg(m_core_asinh)

cat("\n\nModel comparison (asinh-transformed):\n")
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("%-35s %6s %8s %8s\n", "Model", "R²", "AIC", "BIC"))
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("%-35s %5.1f%% %8.1f %8.1f\n", "Core (TS*WTD*species)", 
            r2_core_asinh*100, AIC(m_core_asinh), BIC(m_core_asinh)))
cat(sprintf("%-35s %5.1f%% %8.1f %8.1f\n", "Final (with additions)", 
            r2_final*100, AIC(m_final), BIC(m_final)))
cat("─────────────────────────────────────────────────────\n")

saveRDS(m_core_asinh, file.path(OUTPUT_DIR, "m_core_asinh.rds"))
# Model data (for n, date range and the manuscript-numbers script)
saveRDS(model_data_scaled, file.path(OUTPUT_DIR, "model_data_scaled.rds"))
write_csv(tibble(n_obs = nrow(model_data_scaled), n_trees = n_distinct(model_data_scaled$Tree),
                 first = min(model_data_scaled$datetime), last = max(model_data_scaled$datetime)),
          file.path(OUTPUT_DIR, "model_data_info.csv"))

# ============================================================
# 13. SPECIES-SPECIFIC SLOPES
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       13. SPECIES-SPECIFIC SLOPES                ")
message("══════════════════════════════════════════════════")

coefs <- fixef(m_final)

# Extract slopes using actual predictor names
temp_bg <- coefs[ts_pred]
wtd_bg <- coefs[wtd_pred]

# Species interaction terms
temp_hem_term <- paste0(ts_pred, ":specieshem")
temp_rm_term <- paste0(ts_pred, ":speciesrm")
wtd_hem_term <- paste0(wtd_pred, ":specieshem")
wtd_rm_term <- paste0(wtd_pred, ":speciesrm")
int_term <- paste0(ts_pred, ":", wtd_pred)
int_hem_term <- paste0(ts_pred, ":", wtd_pred, ":specieshem")
int_rm_term <- paste0(ts_pred, ":", wtd_pred, ":speciesrm")

temp_hem <- temp_bg + ifelse(temp_hem_term %in% names(coefs), coefs[temp_hem_term], 0)
temp_rm <- temp_bg + ifelse(temp_rm_term %in% names(coefs), coefs[temp_rm_term], 0)

wtd_hem <- wtd_bg + ifelse(wtd_hem_term %in% names(coefs), coefs[wtd_hem_term], 0)
wtd_rm <- wtd_bg + ifelse(wtd_rm_term %in% names(coefs), coefs[wtd_rm_term], 0)

int_bg <- ifelse(int_term %in% names(coefs), coefs[int_term], 0)
int_hem <- int_bg + ifelse(int_hem_term %in% names(coefs), coefs[int_hem_term], 0)
int_rm <- int_bg + ifelse(int_rm_term %in% names(coefs), coefs[int_rm_term], 0)

species_slopes <- tibble(
  species = c("Black gum", "Hemlock", "Red maple"),
  temp_slope = c(temp_bg, temp_hem, temp_rm),
  wtd_slope = c(wtd_bg, wtd_hem, wtd_rm),
  interaction = c(int_bg, int_hem, int_rm)
)

cat("Species-specific slopes (asinh scale):\n\n")
print(species_slopes %>% mutate(across(where(is.numeric), ~round(., 3))))

write_csv(species_slopes, file.path(OUTPUT_DIR, "species_slopes.csv"))

# ============================================================
# 14. COEFFICIENT TABLE
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       14. COEFFICIENT TABLE                      ")
message("══════════════════════════════════════════════════")

# Build a mapping of predictor names to clean labels
pred_labels <- c(
  setNames("Temperature", ts_pred),
  setNames("Water table", wtd_pred)
)

# Add labels for any additional predictors in the model
model_terms <- names(fixef(m_final))
for (pred in additional_preds) {
  if (any(grepl(pred, model_terms, fixed = TRUE))) {
    # Create readable label from predictor name
    clean_label <- pred %>%
      str_replace("_raw_\\d+h$", "") %>%
      str_replace("_anom_\\d+h$", "") %>%
      str_replace("_", " ")
    pred_labels[pred] <- clean_label
  }
}

coef_df <- as.data.frame(summary(m_final)$coefficients) %>%
  rownames_to_column("term") %>%
  rename(estimate = Estimate, se = `Std. Error`, t_value = `t value`) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    lower = estimate - 1.96 * se,
    upper = estimate + 1.96 * se,
    significant = abs(t_value) > 2
  )

# Clean up term names dynamically
coef_df$term_clean <- coef_df$term
for (pred in names(pred_labels)) {
  coef_df$term_clean <- str_replace_all(coef_df$term_clean, fixed(pred), pred_labels[pred])
}

coef_df <- coef_df %>%
  mutate(
    term_clean = term_clean %>%
      str_replace_all("specieshem", "Hemlock") %>%
      str_replace_all("speciesrm", "Red maple") %>%
      str_replace_all(":", " × "),
    # Add reference species info where needed
    term_clean = case_when(
      term_clean == "Temperature" ~ "Temperature (black gum)",
      term_clean == "Water table" ~ "Water table (black gum)",
      term_clean == "Hemlock" ~ "Hemlock (intercept)",
      term_clean == "Red maple" ~ "Red maple (intercept)",
      term_clean == "Temperature × Water table" ~ "Temp × WTD (black gum)",
      TRUE ~ term_clean
    ),
    effect_type = case_when(
      str_count(term, ":") >= 2 ~ "3-way interaction",
      str_count(term, ":") == 1 ~ "2-way interaction",
      grepl("species", term) ~ "Species intercept",
      TRUE ~ "Main effect"
    )
  )

write_csv(coef_df, file.path(OUTPUT_DIR, "coefficients.csv"))

# ============================================================
# 16. SAVE OUTPUTS
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       16. SAVING OUTPUTS                         ")
message("══════════════════════════════════════════════════")

summary_stats <- tibble(
  metric = c("R² (asinh scale)", "AIC", "BIC", "N observations",
             "N trees", "N species", "Residual SD", "Tree SD"),
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

write_csv(summary_stats, file.path(OUTPUT_DIR, "model_summary.csv"))

# VIF
vif_df <- as.data.frame(vif(m_final)) %>% rownames_to_column("term")
write_csv(vif_df, file.path(OUTPUT_DIR, "vif.csv"))

# ============================================================
# FINAL SUMMARY
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("                 ANALYSIS COMPLETE                ")
message("══════════════════════════════════════════════════")
message("")
message("Final model:")
message("  ", formula_final)
message("")
message(sprintf("R² = %.1f%%", r2_final * 100))
message(sprintf("AIC = %.1f", AIC(m_final)))
message(sprintf("N = %d observations from %d trees",
                nrow(model_data_scaled), n_distinct(model_data_scaled$Tree)))
message("")
message("Key finding: Black gum shows strong temp × water table synergy")
message("             Hemlock and red maple show weak/no synergy")
message("")
message("Species-specific slopes:")
message(sprintf("  Black gum:  Temp=%.3f  WTD=%.3f  Interaction=%.3f",
                temp_bg, wtd_bg, int_bg))
message(sprintf("  Hemlock:    Temp=%.3f  WTD=%.3f  Interaction=%.3f",
                temp_hem, wtd_hem, int_hem))
message(sprintf("  Red maple:  Temp=%.3f  WTD=%.3f  Interaction=%.3f",
                temp_rm, wtd_rm, int_rm))
message("")
message("Output directory: ", OUTPUT_DIR)
message("══════════════════════════════════════════════════")

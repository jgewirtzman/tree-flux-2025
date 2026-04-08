# ============================================================
# 104_bgs_model_complete.R
# 
# Complete BGS CH4 flux modeling pipeline
#
# Requires outputs from:
#   - 103_rolling_windows.R: figures/rolling_correlations/best_windows_by_variable.csv
#   - Data files:
#     - data/processed/aligned_hourly_dataset.csv (environmental)
#     - data/raw/flux_dataset.csv (CH4 flux)
#
# Pipeline:
#   1. Load significant predictors from rolling window analysis
#   2. Build correlation heatmaps for redundancy analysis
#   3. Combined raw + anomaly predictor selection via clustering
#   4. Theory-driven model building with randomized forward selection
#   5. Interaction testing (all 2-way and 3-way)
#   6. Final model with asinh transform
#   7. Visualization (effect sizes + interaction plots)
#
# Final model:
#   asinh(CH4*1000) ~ TS_Ha2 * bvs_wtd_cm * species + SWC_Ha2 + (1|Tree)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(car)
  library(RcppRoll)
  library(cowplot)
  library(pheatmap)
})

set.seed(42)

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  aligned   = "data/processed/aligned_hourly_dataset.csv",
  flux      = "data/processed/flux_with_quality_flags.csv",
  best_windows = "outputs/figures/rolling_correlations/best_windows_by_variable.csv"
)

OUTPUT_DIR <- "outputs/models/bgs_final"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures/predictor_selection", recursive = TRUE, showWarnings = FALSE)

# Correlation threshold for redundancy clustering
CORR_THRESHOLD <- 0.7

# Forward selection settings
N_RANDOM_RUNS <- 100
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

roll_mean <- function(x, n) {
  RcppRoll::roll_mean(x, n = n, align = "right", fill = NA, na.rm = TRUE)
}

roll_sum <- function(x, n) {
  RcppRoll::roll_sum(x, n = n, align = "right", fill = NA, na.rm = TRUE)
}

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

# Build feature matrix (raw mode)
build_feature_matrix_raw <- function(aligned_data, sig_vars_df) {
  met <- aligned_data %>% arrange(datetime)
  result <- met %>% dplyr::select(datetime)
  
  for (i in seq_len(nrow(sig_vars_df))) {
    v <- sig_vars_df$variable[i]
    w <- sig_vars_df$window_hours[i]
    if (!v %in% names(met)) next
    feat_name <- paste0(v, "_", w, "h")
    if (v %in% SUM_VARS) {
      result[[feat_name]] <- roll_sum(met[[v]], w)
    } else {
      result[[feat_name]] <- roll_mean(met[[v]], w)
    }
  }
  result
}

# Build combined feature matrix with both raw and anomaly versions
build_combined_feature_matrix <- function(aligned_data, sig_vars_df) {
  met <- aligned_data %>% arrange(datetime)
  
  # Create anomaly versions for all variables that need them
  anom_vars <- sig_vars_df %>% filter(mode == "anom") %>% pull(variable) %>% unique()
  for (v in anom_vars) {
    if (!v %in% names(met)) next
    met <- make_anomaly(met, v, time_col = "datetime")
  }
  
  result <- met %>% dplyr::select(datetime)
  
  for (i in seq_len(nrow(sig_vars_df))) {
    v <- sig_vars_df$variable[i]
    w <- sig_vars_df$window_hours[i]
    m <- sig_vars_df$mode[i]
    
    src <- if (m == "anom") paste0(v, "_anom") else v
    if (!src %in% names(met)) next
    
    feat_name <- paste0(v, "_", m, "_", w, "h")
    
    if (v %in% SUM_VARS) {
      result[[feat_name]] <- roll_sum(met[[src]], w)
    } else {
      result[[feat_name]] <- roll_mean(met[[src]], w)
    }
  }
  result
}

# Plot correlation heatmap
plot_corr_heatmap <- function(cmat, title, filename) {
  if (nrow(cmat) > 2) {
    hc <- hclust(as.dist(1 - abs(cmat)), method = "complete")
    ord <- hc$order
    cmat <- cmat[ord, ord]
  }
  
  rownames(cmat) <- gsub("_\\d+h$", "", rownames(cmat))
  colnames(cmat) <- gsub("_\\d+h$", "", colnames(cmat))
  
  png(filename, width = 10, height = 9, units = "in", res = 300)
  pheatmap(
    cmat,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    breaks = seq(-1, 1, length.out = 101),
    display_numbers = TRUE,
    number_format = "%.2f",
    number_color = "black",
    fontsize_number = 7,
    fontsize_row = 9,
    fontsize_col = 9,
    main = title,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = "grey90"
  )
  dev.off()
  message("  Saved: ", basename(filename))
}

# Find redundant pairs
find_redundant_pairs <- function(cmat, threshold = 0.7) {
  pairs <- list()
  vars <- colnames(cmat)
  
  for (i in 1:(length(vars) - 1)) {
    for (j in (i + 1):length(vars)) {
      r <- cmat[i, j]
      if (!is.na(r) && abs(r) >= threshold) {
        pairs[[length(pairs) + 1]] <- tibble(
          var1 = vars[i],
          var2 = vars[j],
          r = round(r, 3)
        )
      }
    }
  }
  
  bind_rows(pairs) %>% arrange(desc(abs(r)))
}

# Selection algorithm with clustering
select_predictors_combined <- function(sig_vars_df, cmat, feat_matrix, threshold = 0.7) {
  vars_in_cmat <- colnames(cmat)
  
  sig_vars_df <- sig_vars_df %>%
    mutate(feat_name = paste0(variable, "_", mode, "_", window_hours, "h"))
  
  sig_vars_df <- sig_vars_df %>%
    filter(feat_name %in% vars_in_cmat)
  
  if (nrow(sig_vars_df) == 0) return(tibble())
  
  # Compute missingness
  missingness <- feat_matrix %>%
    dplyr::select(-datetime) %>%
    summarize(across(everything(), ~ mean(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "feat_name", values_to = "missing_frac")
  
  sig_vars_df <- sig_vars_df %>%
    left_join(missingness, by = "feat_name")
  
  # Build adjacency matrix
  adj <- abs(cmat) >= threshold
  diag(adj) <- FALSE
  
  # Find connected components via BFS
  n <- nrow(cmat)
  visited <- rep(FALSE, n)
  cluster_id <- rep(NA_integer_, n)
  current_cluster <- 0
  
  for (start in 1:n) {
    if (visited[start]) next
    current_cluster <- current_cluster + 1
    queue <- start
    
    while (length(queue) > 0) {
      node <- queue[1]
      queue <- queue[-1]
      if (visited[node]) next
      visited[node] <- TRUE
      cluster_id[node] <- current_cluster
      neighbors <- which(adj[node, ] & !visited)
      queue <- c(queue, neighbors)
    }
  }
  
  cluster_map <- tibble(
    feat_name = vars_in_cmat,
    cluster = cluster_id
  )
  
  sig_vars_df <- sig_vars_df %>%
    left_join(cluster_map, by = "feat_name")
  
  # Within each cluster, select best variable
  selected <- sig_vars_df %>%
    group_by(cluster) %>%
    arrange(desc(abs(r)), missing_frac) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(
      selection_reason = "highest |r| in cluster",
      cluster_size = map_int(cluster, ~ sum(sig_vars_df$cluster == .x))
    )
  
  # Track dropped variables
  dropped <- sig_vars_df %>%
    filter(!feat_name %in% selected$feat_name) %>%
    dplyr::select(variable, mode, feat_name, r, cluster)
  
  selected <- selected %>%
    rowwise() %>%
    mutate(
      dropped_vars = paste(
        dropped %>%
          filter(cluster == .data$cluster) %>%
          mutate(label = paste0(variable, "_", mode)) %>%
          pull(label),
        collapse = ", "
      )
    ) %>%
    ungroup()
  
  selected
}

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
    datetime = round_date(as.POSIXct(datetime_posx, tz = "EST"), "hour"),
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

# Store unscaled for plotting later
model_data_unscaled <- model_data

# Total variance for R²
var_total <- var(model_data_scaled$CH4_flux)

# ============================================================
# 6. RANDOMIZED FORWARD SELECTION
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       6. RANDOMIZED FORWARD SELECTION            ")
message("══════════════════════════════════════════════════")

cat("Running", N_RANDOM_RUNS, "iterations...\n")
cat("Core drivers (always first):", paste(core_names, collapse = ", "), "\n")
cat("Candidates:", length(candidate_names), "predictors\n\n")

step_records <- list()

for (run in 1:N_RANDOM_RUNS) {
  if (run %% 20 == 0) cat("  Run", run, "/", N_RANDOM_RUNS, "\n")
  
  selected <- core_names
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
        r2 <- var(predict(fit, re.form = NA)) / var_total
        
        if (r2 > best_r2) {
          best_r2 <- r2
          best_pred <- pred
        }
      }, error = function(e) NULL)
    }
    
    if (is.null(best_pred)) break
    
    step_records[[length(step_records) + 1]] <- tibble(
      run = run,
      step = step,
      predictor = best_pred,
      r2 = best_r2
    )
    
    selected <- c(selected, best_pred)
    remaining <- setdiff(remaining, best_pred)
  }
}

step_df <- bind_rows(step_records)

# Selection frequency
selection_freq <- step_df %>%
  group_by(predictor) %>%
  summarize(
    n_selected = n(),
    pct_selected = n() / N_RANDOM_RUNS * 100,
    mean_step = mean(step),
    .groups = "drop"
  ) %>%
  left_join(
    bind_rows(core_info, candidate_info) %>%
      dplyr::select(predictor, variable, var_group, mode, r),
    by = "predictor"
  ) %>%
  arrange(desc(pct_selected))

cat("\nSelection frequency (top 15):\n")
selection_freq %>%
  head(15) %>%
  mutate(
    pct_selected = round(pct_selected, 1),
    mean_step = round(mean_step, 2),
    r = round(r, 3)
  ) %>%
  print()

write_csv(selection_freq, "outputs/figures/predictor_selection/theory_driven_selection_freq.csv")

# Strong consensus (≥80%)
strong_consensus <- selection_freq %>%
  filter(pct_selected >= 80)

cat("\nStrong consensus predictors (≥80%):", nrow(strong_consensus), "\n")
if (nrow(strong_consensus) > 0) {
  print(strong_consensus %>% dplyr::select(predictor, pct_selected, mean_step))
}

# ============================================================
# 7. MODEL COMPARISON (NO INTERACTIONS)
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       7. BASE MODEL COMPARISON                   ")
message("══════════════════════════════════════════════════")

# Species only
m_species <- lmer(CH4_flux ~ species + (1|Tree), 
                  data = model_data_scaled, REML = FALSE)
r2_species <- var(predict(m_species, re.form = NA)) / var_total

# Core only
core_formula <- paste("CH4_flux ~", paste(core_names, collapse = " + "), "+ species + (1|Tree)")
m_core <- lmer(as.formula(core_formula), data = model_data_scaled, REML = FALSE)
r2_core <- var(predict(m_core, re.form = NA)) / var_total

# Core + consensus
if (nrow(strong_consensus) > 0) {
  strong_preds <- c(core_names, strong_consensus$predictor)
  strong_formula <- paste("CH4_flux ~", paste(strong_preds, collapse = " + "), "+ species + (1|Tree)")
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
cat(sprintf("%-25s %5.1f%% %8.1f %8.1f\n", "Core + consensus", r2_strong*100, AIC(m_strong), BIC(m_strong)))
cat("─────────────────────────────────────────────────────\n")

# VIF for strong model
cat("\nVIF for core + consensus model:\n")
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

cat("Core interaction predictors (theory-driven):\n")
cat("  Temperature:", ts_pred, "\n")
cat("  Water table:", wtd_pred, "\n\n")

if (is.na(ts_pred) || is.na(wtd_pred)) {
  stop("Required core predictors (TS_Ha2 and/or bvs_wtd_cm) not found!")
}

# Start with consensus predictors from forward selection
refined_preds <- strong_preds

# Check VIF
cat("VIF for core + consensus model:\n")
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
  
  VIF_DROP <- c("TS_Ha1_raw_285h", "s10t_raw_264h", "LE_Ha2_anom_9h")  # EDIT THIS LIST
  
  cat("\n→ Dropping (per VIF_DROP config):", paste(VIF_DROP, collapse = ", "), "\n")
  
  # Remove specified predictors
  refined_preds <- setdiff(refined_preds, VIF_DROP)
  
  # Refit model
  m_strong <- lmer(as.formula(paste("CH4_flux ~", paste(refined_preds, collapse = " + "), 
                                    "+ species + (1|Tree)")), 
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
formula_core <- paste0("CH4_flux ~ ", ts_pred, " * ", wtd_pred, " * species + (1|Tree)")
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
                          pred, " + (1|Tree)")
    
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
                             pred, " * species + (1|Tree)")
    
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
                               paste(added_terms, collapse = " + "), " + (1|Tree)")
    
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
model_data_scaled <- model_data_scaled %>%
  mutate(CH4_flux_asinh = asinh(CH4_flux * 1000))

var_total_asinh <- var(model_data_scaled$CH4_flux_asinh)

# Convert best formula to asinh version
formula_final <- gsub("CH4_flux ~", "CH4_flux_asinh ~", best_formula)

cat("Final model formula (asinh-transformed):\n")
cat(" ", formula_final, "\n\n")

m_final <- lmer(as.formula(formula_final), data = model_data_scaled, REML = FALSE)

r2_final <- var(predict(m_final, re.form = NA)) / var_total_asinh

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
r2_core_asinh <- var(predict(m_core_asinh, re.form = NA)) / var_total_asinh

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
# 14. EFFECT SIZE PLOT
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       14. EFFECT SIZE PLOT                       ")
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

effect_colors <- c(
  "Main effect" = "#2166AC",
  "Species intercept" = "#4DAF4A",
  "2-way interaction" = "#FF7F00",
  "3-way interaction" = "#E31A1C"
)

p_effects <- ggplot(coef_df, aes(x = estimate, y = reorder(term_clean, estimate), color = effect_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 0.8) +
  geom_point(aes(shape = significant), size = 3) +
  scale_color_manual(values = effect_colors, name = "Effect type") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                     labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                     name = "Significance") +
  labs(
    x = "Coefficient (± 95% CI)",
    y = NULL,
    title = "Effect sizes for CH₄ flux model",
    subtitle = paste0("R² = ", round(r2_final * 100, 1), "%")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(OUTPUT_DIR, "effect_sizes.png"), p_effects,
       width = 10, height = 7, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "effect_sizes.pdf"), p_effects,
       width = 10, height = 7)

message("Saved: effect_sizes.png/pdf")

write_csv(coef_df, file.path(OUTPUT_DIR, "coefficients.csv"))

# ============================================================
# 15. INTERACTION PLOTS
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       15. INTERACTION PLOTS                      ")
message("══════════════════════════════════════════════════")

# Get scaling parameters using dynamic predictor names
scaling_params <- list(
  TS_mean = mean(model_data_unscaled[[ts_pred]], na.rm = TRUE),
  TS_sd = sd(model_data_unscaled[[ts_pred]], na.rm = TRUE),
  WTD_mean = mean(model_data_unscaled[[wtd_pred]], na.rm = TRUE),
  WTD_sd = sd(model_data_unscaled[[wtd_pred]], na.rm = TRUE)
)

# Prediction grid
temp_range_scaled <- seq(-2, 2, length.out = 100)
wtd_percentiles <- c(10, 25, 50, 75, 90)

wtd_quantiles_unscaled <- quantile(model_data_unscaled[[wtd_pred]],
                                   probs = wtd_percentiles/100, na.rm = TRUE)
wtd_quantiles_scaled <- (wtd_quantiles_unscaled - scaling_params$WTD_mean) / scaling_params$WTD_sd

temp_unscaled <- temp_range_scaled * scaling_params$TS_sd + scaling_params$TS_mean

# Get all predictor names needed for prediction (from the model terms)
model_vars <- all.vars(formula(m_final))
model_vars <- setdiff(model_vars, c("CH4_flux_asinh", "Tree", "species"))

# Generate predictions
pred_list <- list()

for (sp in c("bg", "hem", "rm")) {
  for (i in seq_along(wtd_percentiles)) {
    # Build prediction dataframe - start with the right number of rows
    n_pts <- length(temp_range_scaled)
    
    pred_df <- data.frame(
      species = factor(rep(sp, n_pts), levels = c("bg", "hem", "rm")),
      Tree = factor(rep(model_data_scaled$Tree[1], n_pts))
    )
    
    # Add core predictors
    pred_df[[ts_pred]] <- temp_range_scaled
    pred_df[[wtd_pred]] <- rep(wtd_quantiles_scaled[i], n_pts)
    
    # Add any additional predictors at their mean (0 for scaled data)
    for (v in model_vars) {
      if (!v %in% names(pred_df) && v %in% names(model_data_scaled)) {
        pred_df[[v]] <- rep(0, n_pts)  # Mean for scaled predictors
      }
    }
    
    pred_df$fit_asinh <- predict(m_final, newdata = pred_df, re.form = NA)
    pred_df$fit <- sinh(pred_df$fit_asinh) / 1000  # back-transform
    
    pred_df$temp_unscaled <- temp_unscaled
    pred_df$wtd_percentile <- wtd_percentiles[i]
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

# Colors
wtd_colors <- c("10" = "#8B4513", "25" = "#CD853F", "50" = "#808080",
                "75" = "#6495ED", "90" = "#0000CD")
species_colors <- c("Black gum" = "#1B9E77", "Hemlock" = "#D95F02", "Red maple" = "#7570B3")

# Panel plot by species
p_interaction <- ggplot(pred_all, aes(x = temp_unscaled, y = fit * 1000,
                                      color = wtd_percentile)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ species_name, scales = "free_y") +
  scale_color_manual(values = wtd_colors, labels = paste0(wtd_percentiles, "th"),
                     name = "Water table\npercentile") +
  labs(
    x = "Soil temperature (°C)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Temperature × water table interactions by species"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(OUTPUT_DIR, "interaction_by_species.png"), p_interaction,
       width = 12, height = 5, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "interaction_by_species.pdf"), p_interaction,
       width = 12, height = 5)

message("Saved: interaction_by_species.png/pdf")

# Comparison plot (wet vs dry)
pred_extremes <- pred_all %>%
  filter(wtd_percentile %in% c("10", "90"))

p_comparison <- ggplot(pred_extremes, aes(x = temp_unscaled, y = fit * 1000,
                                          color = species_name, linetype = wtd_percentile)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = species_colors, name = "Species") +
  scale_linetype_manual(values = c("10" = "dashed", "90" = "solid"),
                        labels = c("10" = "Dry (10th)", "90" = "Wet (90th)"),
                        name = "Water table") +
  labs(
    x = "Soil temperature (°C)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Species differ in temperature × water table response"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "interaction_comparison.png"), p_comparison,
       width = 9, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "interaction_comparison.pdf"), p_comparison,
       width = 9, height = 6)

message("Saved: interaction_comparison.png/pdf")

# Combined figure
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








# Model with SWC instead of WTD in the core interaction
formula_swc_core <- "CH4_flux_asinh ~ TS_Ha2_raw_282h * NEON_SWC_shallow_raw_129h * species + (1|Tree)"

m_swc_core <- lmer(as.formula(formula_swc_core), data = model_data_scaled, REML = FALSE)

cat("Model with SWC as core moisture variable:\n")
cat("  R²:", round(var(predict(m_swc_core, re.form = NA)) / var_total_asinh * 100, 1), "%\n")
cat("  AIC:", round(AIC(m_swc_core), 1), "\n")
cat("  BIC:", round(BIC(m_swc_core), 1), "\n\n")

cat("Coefficients:\n")
print(round(summary(m_swc_core)$coefficients, 3))

# Compare to original core model
cat("\n\nComparison:\n")
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("%-30s %6s %8s %8s\n", "Model", "R²", "AIC", "BIC"))
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("%-30s %5.1f%% %8.1f %8.1f\n", "Core (Temp × WTD × species)", 
            r2_core_asinh*100, AIC(m_core_asinh), BIC(m_core_asinh)))
cat(sprintf("%-30s %5.1f%% %8.1f %8.1f\n", "Core (Temp × SWC × species)", 
            var(predict(m_swc_core, re.form = NA)) / var_total_asinh * 100, 
            AIC(m_swc_core), BIC(m_swc_core)))
cat("─────────────────────────────────────────────────────\n")

# Extract species-specific SWC slopes
coefs_swc <- fixef(m_swc_core)
swc_bg <- coefs_swc["NEON_SWC_shallow_raw_129h"]
swc_hem <- swc_bg + coefs_swc["NEON_SWC_shallow_raw_129h:specieshem"]
swc_rm <- swc_bg + coefs_swc["NEON_SWC_shallow_raw_129h:speciesrm"]

cat("\nSWC slopes by species:\n")
cat(sprintf("  Black gum: %.3f\n", swc_bg))
cat(sprintf("  Hemlock:   %.3f\n", swc_hem))
cat(sprintf("  Red maple: %.3f\n", swc_rm))















# Model with LE_Ha1 instead of TEMP in the core interaction
formula_le1_core <- "CH4_flux_asinh ~ LE_Ha1_raw_123h * bvs_wtd_cm_raw_132h * species + (1|Tree)"

m_le1_core <- lmer(as.formula(formula_le1_core), data = model_data_scaled, REML = FALSE)

cat("Model with LE_Ha1 as core energy variable:\n")
cat("  R²:", round(var(predict(m_le1_core, re.form = NA)) / var_total_asinh * 100, 1), "%\n")
cat("  AIC:", round(AIC(m_le1_core), 1), "\n")
cat("  BIC:", round(BIC(m_le1_core), 1), "\n\n")

cat("Coefficients:\n")
print(round(summary(m_le1_core)$coefficients, 3))

# Extract species-specific LE_Ha1 slopes
coefs_le1 <- fixef(m_le1_core)
le1_bg <- coefs_le1["LE_Ha1_raw_123h"]
le1_hem <- le1_bg + coefs_le1["LE_Ha1_raw_123h:specieshem"]
le1_rm <- le1_bg + coefs_le1["LE_Ha1_raw_123h:speciesrm"]

cat("\nLE_Ha1 slopes by species:\n")
cat(sprintf("  Black gum: %.3f\n", le1_bg))
cat(sprintf("  Hemlock:   %.3f\n", le1_hem))
cat(sprintf("  Red maple: %.3f\n", le1_rm))

# Now test LE_Ha2_anom
formula_le2_core <- "CH4_flux_asinh ~ LE_Ha2_anom_9h * bvs_wtd_cm_raw_132h * species + (1|Tree)"

m_le2_core <- lmer(as.formula(formula_le2_core), data = model_data_scaled, REML = FALSE)

cat("\n\n─────────────────────────────────────────────────────\n")
cat("Model with LE_Ha2_anom as core energy variable:\n")
cat("  R²:", round(var(predict(m_le2_core, re.form = NA)) / var_total_asinh * 100, 1), "%\n")
cat("  AIC:", round(AIC(m_le2_core), 1), "\n")
cat("  BIC:", round(BIC(m_le2_core), 1), "\n\n")

cat("Coefficients:\n")
print(round(summary(m_le2_core)$coefficients, 3))

# Extract species-specific LE_Ha2_anom slopes
coefs_le2 <- fixef(m_le2_core)
le2_bg <- coefs_le2["LE_Ha2_anom_9h"]
le2_hem <- le2_bg + coefs_le2["LE_Ha2_anom_9h:specieshem"]
le2_rm <- le2_bg + coefs_le2["LE_Ha2_anom_9h:speciesrm"]

cat("\nLE_Ha2_anom slopes by species:\n")
cat(sprintf("  Black gum: %.3f\n", le2_bg))
cat(sprintf("  Hemlock:   %.3f\n", le2_hem))
cat(sprintf("  Red maple: %.3f\n", le2_rm))

# Summary comparison
cat("\n\n══════════════════════════════════════════════════════\n")
cat("SUMMARY: Core model variants (energy × WTD × species)\n")
cat("══════════════════════════════════════════════════════\n")
cat(sprintf("%-35s %6s %8s %8s\n", "Model", "R²", "AIC", "BIC"))
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("%-35s %5.1f%% %8.1f %8.1f\n", "Temp × WTD × species", 
            r2_core_asinh*100, AIC(m_core_asinh), BIC(m_core_asinh)))
cat(sprintf("%-35s %5.1f%% %8.1f %8.1f\n", "LE_Ha1 × WTD × species", 
            var(predict(m_le1_core, re.form = NA)) / var_total_asinh * 100, 
            AIC(m_le1_core), BIC(m_le1_core)))
cat(sprintf("%-35s %5.1f%% %8.1f %8.1f\n", "LE_Ha2_anom × WTD × species", 
            var(predict(m_le2_core, re.form = NA)) / var_total_asinh * 100, 
            AIC(m_le2_core), BIC(m_le2_core)))
cat("─────────────────────────────────────────────────────\n")











# ============================================================
# BGS MODEL RESULTS SUMMARY FOR MANUSCRIPT
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("WETLAND MODEL RESULTS FOR MANUSCRIPT")
message(paste(rep("=", 60), collapse = ""))

# Final model summary
message("\n--- Final Model Summary ---")
message("Formula: ", formula_final)
message("R² (marginal): ", round(r2_final * 100, 1), "%")
message("AIC: ", round(AIC(m_final), 1))
message("BIC: ", round(BIC(m_final), 1))
message("N obs: ", nrow(model_data_scaled))
message("N trees: ", n_distinct(model_data_scaled$Tree))

# Fixed effects with confidence intervals
message("\n--- Fixed Effects (asinh scale) ---")
fixef_summary <- as.data.frame(summary(m_final)$coefficients)
fixef_summary$term <- rownames(fixef_summary)
fixef_summary <- fixef_summary %>%
  mutate(
    ci_lower = Estimate - 1.96 * `Std. Error`,
    ci_upper = Estimate + 1.96 * `Std. Error`,
    p_value = 2 * (1 - pt(abs(`t value`), df = nrow(model_data_scaled) - length(fixef(m_final))))
  ) %>%
  dplyr::select(term, Estimate, `Std. Error`, ci_lower, ci_upper, `t value`, p_value)
print(fixef_summary, digits = 3)

# Variance components
message("\n--- Variance Components ---")
vc <- as.data.frame(VarCorr(m_final))
print(vc)

# VIF
message("\n--- Variance Inflation Factors ---")
print(vif(m_final))

# Species-specific slopes
message("\n--- Species-Specific Slopes ---")
coefs <- fixef(m_final)
message("\nTemperature effects:")
message("  N. sylvatica: ", round(temp_bg, 3))
message("  T. canadensis: ", round(temp_hem, 3))
message("  A. rubrum: ", round(temp_rm, 3))

message("\nWater table effects:")
message("  N. sylvatica: ", round(wtd_bg, 3))
message("  T. canadensis: ", round(wtd_hem, 3))
message("  A. rubrum: ", round(wtd_rm, 3))

message("\nTemperature × Water table interactions:")
message("  N. sylvatica: ", round(int_bg, 3))
message("  T. canadensis: ", round(int_hem, 3))
message("  A. rubrum: ", round(int_rm, 3))

message("\nSoil water content effect:")
message("  ", round(coefs["NEON_SWC_shallow_raw_129h"], 3))

# Model comparison (alternatives)
message("\n--- Model Comparison ---")
cat(sprintf("%-35s %6s %8s %8s\n", "Model", "R²", "AIC", "BIC"))
cat("─────────────────────────────────────────────────────\n")
cat(sprintf("%-35s %5.1f%% %8.1f %8.1f\n", "Core (Temp × WTD × species)", 
            r2_core_asinh*100, AIC(m_core_asinh), BIC(m_core_asinh)))
cat(sprintf("%-35s %5.1f%% %8.1f %8.1f\n", "Final (+SWC)", 
            r2_final*100, AIC(m_final), BIC(m_final)))
cat(sprintf("%-35s %5.1f%% %8.1f %8.1f\n", "Alt: Temp × SWC × species", 
            var(predict(m_swc_core, re.form = NA)) / var_total_asinh * 100, 
            AIC(m_swc_core), BIC(m_swc_core)))
cat("─────────────────────────────────────────────────────\n")

# Interpretation at extremes
message("\n--- Predicted Flux at Environmental Extremes ---")
# Use the prediction grid to show min/max conditions
pred_summary <- pred_all %>%
  filter(wtd_percentile %in% c("10", "90")) %>%
  group_by(species_name, wtd_percentile) %>%
  summarise(
    temp_range = paste0(round(min(temp_unscaled), 1), "–", round(max(temp_unscaled), 1), "°C"),
    flux_min = round(min(fit * 1000), 2),
    flux_max = round(max(fit * 1000), 2),
    .groups = "drop"
  )
print(pred_summary)









# ============================================================
# COMPLETE BGS MODEL RESULTS FOR MANUSCRIPT
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("WETLAND MODEL - COMPLETE RESULTS")
message(paste(rep("=", 60), collapse = ""))

# 1. DATA SUMMARY
message("\n--- 1. DATA SUMMARY ---")
message("N observations: ", nrow(model_data_scaled))
message("N trees: ", n_distinct(model_data_scaled$Tree))
message("N species: ", n_distinct(model_data_scaled$species))
message("Date range: ", min(model_data_scaled$datetime), " to ", max(model_data_scaled$datetime))

# 2. CORE MODEL (TEMP × WTD × SPECIES)
message("\n--- 2. CORE MODEL: Temp × WTD × Species ---")
message("\nFormula: CH4_flux_asinh ~ TS_Ha2 × bvs_wtd_cm × species + (1|Tree)")
message("All predictors standardized (mean=0, SD=1)")
message("\nR² (marginal): ", round(r2_core_asinh * 100, 1), "%")
message("AIC: ", round(AIC(m_core_asinh), 1))
message("BIC: ", round(BIC(m_core_asinh), 1))
message("N parameters: ", length(fixef(m_core_asinh)))

message("\nFixed effects (standardized scale):")
core_fixef <- as.data.frame(summary(m_core_asinh)$coefficients)
core_fixef$term <- rownames(core_fixef)
core_fixef <- core_fixef %>%
  mutate(
    p_value = 2 * (1 - pt(abs(`t value`), df = nrow(model_data_scaled) - length(fixef(m_core_asinh)))),
    sig = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
   dplyr::select(term, Estimate, `Std. Error`, `t value`, p_value, sig)
print(core_fixef, digits = 3)

message("\nVariance components:")
vc_core <- as.data.frame(VarCorr(m_core_asinh))
print(vc_core)
icc_core <- vc_core$vcov[1] / sum(vc_core$vcov)
message("ICC: ", round(icc_core, 3))

message("\nVIF:")
vif_core <- vif(m_core_asinh)
print(vif_core)

# 3. FULL MODEL (+ LE + SWC)
message("\n--- 3. FULL MODEL: Core + LE + SWC ---")
message("\nFormula: CH4_flux_asinh ~ TS_Ha2 × bvs_wtd_cm × species + LE_Ha1 × species + SWC × species + (1|Tree)")
message("\nR² (marginal): ", round(r2_final * 100, 1), "%")
message("AIC: ", round(AIC(m_final), 1))
message("BIC: ", round(BIC(m_final), 1))
message("N parameters: ", length(fixef(m_final)))
message("\nImprovement over core model:")
message("  ΔR²: +", round((r2_final - r2_core_asinh) * 100, 1), "%")
message("  ΔAIC: ", round(AIC(m_final) - AIC(m_core_asinh), 1))
message("  ΔBIC: ", round(BIC(m_final) - BIC(m_core_asinh), 1))

message("\nFixed effects (standardized scale):")
full_fixef <- as.data.frame(summary(m_final)$coefficients)
full_fixef$term <- rownames(full_fixef)
full_fixef <- full_fixef %>%
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
print(full_fixef, digits = 3)

message("\nVariance components:")
vc_full <- as.data.frame(VarCorr(m_final))
print(vc_full)
icc_full <- vc_full$vcov[1] / sum(vc_full$vcov)
message("ICC: ", round(icc_full, 3))

message("\nVIF (GVIF^(1/(2*Df)) for interactions):")
vif_full <- vif(m_final)
print(vif_full)

# 4. STANDARDIZED EFFECT SIZES - SPECIES-SPECIFIC
message("\n--- 4. STANDARDIZED EFFECT SIZES BY SPECIES ---")

coefs_core <- fixef(m_core_asinh)
coefs_full <- fixef(m_final)

message("\nCORE MODEL:")
message("\nTemperature effects (standardized):")
temp_core_bg <- coefs_core["TS_Ha2_raw_282h"]
temp_core_hem <- coefs_core["TS_Ha2_raw_282h"] + coefs_core["TS_Ha2_raw_282h:specieshem"]
temp_core_rm <- coefs_core["TS_Ha2_raw_282h"] + coefs_core["TS_Ha2_raw_282h:speciesrm"]
message("  N. sylvatica: ", round(temp_core_bg, 3))
message("  T. canadensis: ", round(temp_core_hem, 3))
message("  A. rubrum: ", round(temp_core_rm, 3))

message("\nWater table effects (standardized):")
wtd_core_bg <- coefs_core["bvs_wtd_cm_raw_132h"]
wtd_core_hem <- coefs_core["bvs_wtd_cm_raw_132h"] + coefs_core["bvs_wtd_cm_raw_132h:specieshem"]
wtd_core_rm <- coefs_core["bvs_wtd_cm_raw_132h"] + coefs_core["bvs_wtd_cm_raw_132h:speciesrm"]
message("  N. sylvatica: ", round(wtd_core_bg, 3))
message("  T. canadensis: ", round(wtd_core_hem, 3))
message("  A. rubrum: ", round(wtd_core_rm, 3))

message("\nTemp × WTD interaction (standardized):")
int_core_bg <- coefs_core["TS_Ha2_raw_282h:bvs_wtd_cm_raw_132h"]
int_core_hem <- coefs_core["TS_Ha2_raw_282h:bvs_wtd_cm_raw_132h"] + coefs_core["TS_Ha2_raw_282h:bvs_wtd_cm_raw_132h:specieshem"]
int_core_rm <- coefs_core["TS_Ha2_raw_282h:bvs_wtd_cm_raw_132h"] + coefs_core["TS_Ha2_raw_282h:bvs_wtd_cm_raw_132h:speciesrm"]
message("  N. sylvatica: ", round(int_core_bg, 3))
message("  T. canadensis: ", round(int_core_hem, 3))
message("  A. rubrum: ", round(int_core_rm, 3))

message("\nFULL MODEL:")
message("\nTemperature effects (standardized):")
temp_full_bg <- coefs_full["TS_Ha2_raw_282h"]
temp_full_hem <- coefs_full["TS_Ha2_raw_282h"] + coefs_full["TS_Ha2_raw_282h:specieshem"]
temp_full_rm <- coefs_full["TS_Ha2_raw_282h"] + coefs_full["TS_Ha2_raw_282h:speciesrm"]
message("  N. sylvatica: ", round(temp_full_bg, 3))
message("  T. canadensis: ", round(temp_full_hem, 3))
message("  A. rubrum: ", round(temp_full_rm, 3))

message("\nWater table effects (standardized):")
wtd_full_bg <- coefs_full["bvs_wtd_cm_raw_132h"]
wtd_full_hem <- coefs_full["bvs_wtd_cm_raw_132h"] + coefs_full["bvs_wtd_cm_raw_132h:specieshem"]
wtd_full_rm <- coefs_full["bvs_wtd_cm_raw_132h"] + coefs_full["bvs_wtd_cm_raw_132h:speciesrm"]
message("  N. sylvatica: ", round(wtd_full_bg, 3))
message("  T. canadensis: ", round(wtd_full_hem, 3))
message("  A. rubrum: ", round(wtd_full_rm, 3))

message("\nTemp × WTD interaction (standardized):")
int_full_bg <- coefs_full["TS_Ha2_raw_282h:bvs_wtd_cm_raw_132h"]
int_full_hem <- coefs_full["TS_Ha2_raw_282h:bvs_wtd_cm_raw_132h"] + coefs_full["TS_Ha2_raw_282h:bvs_wtd_cm_raw_132h:specieshem"]
int_full_rm <- coefs_full["TS_Ha2_raw_282h:bvs_wtd_cm_raw_132h"] + coefs_full["TS_Ha2_raw_282h:bvs_wtd_cm_raw_132h:speciesrm"]
message("  N. sylvatica: ", round(int_full_bg, 3))
message("  T. canadensis: ", round(int_full_hem, 3))
message("  A. rubrum: ", round(int_full_rm, 3))

message("\nLatent heat effects (standardized, full model only):")
le_bg <- coefs_full["LE_Ha1_raw_123h"]
le_hem <- coefs_full["LE_Ha1_raw_123h"] + coefs_full["specieshem:LE_Ha1_raw_123h"]
le_rm <- coefs_full["LE_Ha1_raw_123h"] + coefs_full["speciesrm:LE_Ha1_raw_123h"]
message("  N. sylvatica: ", round(le_bg, 3))
message("  T. canadensis: ", round(le_hem, 3))
message("  A. rubrum: ", round(le_rm, 3))

message("\nSoil water content effects (standardized, full model only):")
swc_bg <- coefs_full["NEON_SWC_shallow_raw_129h"]
swc_hem <- coefs_full["NEON_SWC_shallow_raw_129h"] + coefs_full["specieshem:NEON_SWC_shallow_raw_129h"]
swc_rm <- coefs_full["NEON_SWC_shallow_raw_129h"] + coefs_full["speciesrm:NEON_SWC_shallow_raw_129h"]
message("  N. sylvatica: ", round(swc_bg, 3))
message("  T. canadensis: ", round(swc_hem, 3))
message("  A. rubrum: ", round(swc_rm, 3))

# 5. MODEL COMPARISON TABLE
message("\n--- 5. MODEL COMPARISON ---")
cat(sprintf("%-35s %6s %8s %8s %8s\n", "Model", "R²", "AIC", "BIC", "Nparam"))
cat("──────────────────────────────────────────────────────────\n")
cat(sprintf("%-35s %5.1f%% %8.1f %8.1f %8d\n", "Core (Temp × WTD × species)", 
            r2_core_asinh*100, AIC(m_core_asinh), BIC(m_core_asinh), length(fixef(m_core_asinh))))
cat(sprintf("%-35s %5.1f%% %8.1f %8.1f %8d\n", "Full (+LE +SWC interactions)", 
            r2_final*100, AIC(m_final), BIC(m_final), length(fixef(m_final))))
cat("──────────────────────────────────────────────────────────\n")
cat(sprintf("%-35s %+5.1f%% %+8.1f %+8.1f %+8d\n", "Improvement", 
            (r2_final - r2_core_asinh)*100, 
            AIC(m_final) - AIC(m_core_asinh),
            BIC(m_final) - BIC(m_core_asinh),
            length(fixef(m_final)) - length(fixef(m_core_asinh))))

# 6. PREDICTED FLUX RANGE
message("\n--- 6. PREDICTED FLUX AT ENVIRONMENTAL EXTREMES ---")
pred_summary <- pred_all %>%
  filter(wtd_percentile %in% c("10", "90")) %>%
  group_by(species_name, wtd_percentile) %>%
  summarise(
    temp_range = paste0(round(min(temp_unscaled), 1), "–", round(max(temp_unscaled), 1), "°C"),
    flux_min = round(min(fit * 1000), 2),
    flux_max = round(max(fit * 1000), 2),
    flux_range = round(max(fit * 1000) - min(fit * 1000), 2),
    .groups = "drop"
  )
print(pred_summary)



# ============================================================
# ALTERNATIVE MODELS AND ROBUSTNESS CHECKS
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("ALTERNATIVE MODELS - ROBUSTNESS CHECKS")
message(paste(rep("=", 60), collapse = ""))

# 1. SWC AS CORE MOISTURE VARIABLE (instead of WTD)
message("\n--- 1. SWC REPLACING WTD IN CORE INTERACTION ---")
message("\nModel: Temp × SWC × species (instead of Temp × WTD × species)")
message("R²: ", round(var(predict(m_swc_core, re.form = NA)) / var_total_asinh * 100, 1), "%")
message("AIC: ", round(AIC(m_swc_core), 1))
message("BIC: ", round(BIC(m_swc_core), 1))

message("\nSWC slopes by species (standardized):")
coefs_swc <- fixef(m_swc_core)
swc_bg <- coefs_swc["NEON_SWC_shallow_raw_129h"]
swc_hem <- swc_bg + coefs_swc["NEON_SWC_shallow_raw_129h:specieshem"]
swc_rm <- swc_bg + coefs_swc["NEON_SWC_shallow_raw_129h:speciesrm"]
message("  N. sylvatica: ", round(swc_bg, 3))
message("  T. canadensis: ", round(swc_hem, 3))
message("  A. rubrum: ", round(swc_rm, 3))

message("\nComparison to WTD core model:")
message("  WTD core: R² = ", round(r2_core_asinh * 100, 1), "%, AIC = ", round(AIC(m_core_asinh), 1))
message("  SWC core: R² = ", round(var(predict(m_swc_core, re.form = NA)) / var_total_asinh * 100, 1), "%, AIC = ", round(AIC(m_swc_core), 1))
message("  ΔAIC = ", round(AIC(m_swc_core) - AIC(m_core_asinh), 1), " (positive = worse fit)")

# 2. LE_Ha1 AS CORE ENERGY VARIABLE (instead of Temp)
message("\n--- 2. LE_Ha1 REPLACING TEMP IN CORE INTERACTION ---")
message("\nModel: LE_Ha1 × WTD × species (instead of Temp × WTD × species)")
message("R²: ", round(var(predict(m_le1_core, re.form = NA)) / var_total_asinh * 100, 1), "%")
message("AIC: ", round(AIC(m_le1_core), 1))
message("BIC: ", round(BIC(m_le1_core), 1))

message("\nLE_Ha1 slopes by species (standardized):")
coefs_le1 <- fixef(m_le1_core)
le1_bg <- coefs_le1["LE_Ha1_raw_123h"]
le1_hem <- le1_bg + coefs_le1["LE_Ha1_raw_123h:specieshem"]
le1_rm <- le1_bg + coefs_le1["LE_Ha1_raw_123h:speciesrm"]
message("  N. sylvatica: ", round(le1_bg, 3))
message("  T. canadensis: ", round(le1_hem, 3))
message("  A. rubrum: ", round(le1_rm, 3))

message("\nComparison to Temp core model:")
message("  Temp core: R² = ", round(r2_core_asinh * 100, 1), "%, AIC = ", round(AIC(m_core_asinh), 1))
message("  LE_Ha1 core: R² = ", round(var(predict(m_le1_core, re.form = NA)) / var_total_asinh * 100, 1), "%, AIC = ", round(AIC(m_le1_core), 1))
message("  ΔAIC = ", round(AIC(m_le1_core) - AIC(m_core_asinh), 1))

# 3. LE_Ha2_anom AS CORE ENERGY VARIABLE
message("\n--- 3. LE_Ha2_anom REPLACING TEMP IN CORE INTERACTION ---")
message("\nModel: LE_Ha2_anom × WTD × species (instead of Temp × WTD × species)")
message("R²: ", round(var(predict(m_le2_core, re.form = NA)) / var_total_asinh * 100, 1), "%")
message("AIC: ", round(AIC(m_le2_core), 1))
message("BIC: ", round(BIC(m_le2_core), 1))

message("\nLE_Ha2_anom slopes by species (standardized):")
coefs_le2 <- fixef(m_le2_core)
le2_bg <- coefs_le2["LE_Ha2_anom_9h"]
le2_hem <- le2_bg + coefs_le2["LE_Ha2_anom_9h:specieshem"]
le2_rm <- le2_bg + coefs_le2["LE_Ha2_anom_9h:speciesrm"]
message("  N. sylvatica: ", round(le2_bg, 3))
message("  T. canadensis: ", round(le2_hem, 3))
message("  A. rubrum: ", round(le2_rm, 3))

# 4. UNIVARIATE MODELS (to check direction reversals)
message("\n--- 4. UNIVARIATE MODELS (checking effect directions) ---")

# LE_Ha1 alone
m_le_only <- lmer(CH4_flux_asinh ~ LE_Ha1_raw_123h * species + (1|Tree), 
                  data = model_data_scaled, REML = FALSE)
message("\nLE_Ha1 only (no other predictors):")
le_only <- fixef(m_le_only)
message("  N. sylvatica: ", round(le_only["LE_Ha1_raw_123h"], 3))
message("  T. canadensis: ", round(le_only["LE_Ha1_raw_123h"] + le_only["specieshem:LE_Ha1_raw_123h"], 3))
message("  A. rubrum: ", round(le_only["LE_Ha1_raw_123h"] + le_only["speciesrm:LE_Ha1_raw_123h"], 3))
message("  Compare to full model: N. sylvatica = ", round(le_bg, 3), 
        " (direction ", ifelse(sign(le_only["LE_Ha1_raw_123h"]) == sign(le_bg), "same", "REVERSED"), ")")

# SWC alone
m_swc_only <- lmer(CH4_flux_asinh ~ NEON_SWC_shallow_raw_129h * species + (1|Tree), 
                   data = model_data_scaled, REML = FALSE)
message("\nSWC only (no other predictors):")
swc_only <- fixef(m_swc_only)
message("  N. sylvatica: ", round(swc_only["NEON_SWC_shallow_raw_129h"], 3))
message("  T. canadensis: ", round(swc_only["NEON_SWC_shallow_raw_129h"] + swc_only["specieshem:NEON_SWC_shallow_raw_129h"], 3))
message("  A. rubrum: ", round(swc_only["NEON_SWC_shallow_raw_129h"] + swc_only["speciesrm:NEON_SWC_shallow_raw_129h"], 3))
message("  Compare to full model: N. sylvatica = ", round(swc_bg, 3),
        " (direction ", ifelse(sign(swc_only["NEON_SWC_shallow_raw_129h"]) == sign(swc_bg), "same", "REVERSED"), ")")

# Temp alone (for reference)
m_temp_only <- lmer(CH4_flux_asinh ~ TS_Ha2_raw_282h * species + (1|Tree), 
                    data = model_data_scaled, REML = FALSE)
message("\nTemp only (no other predictors):")
temp_only <- fixef(m_temp_only)
message("  N. sylvatica: ", round(temp_only["TS_Ha2_raw_282h"], 3))
message("  T. canadensis: ", round(temp_only["TS_Ha2_raw_282h"] + temp_only["TS_Ha2_raw_282h:specieshem"], 3))
message("  A. rubrum: ", round(temp_only["TS_Ha2_raw_282h"] + temp_only["TS_Ha2_raw_282h:speciesrm"], 3))
message("  Compare to full model: N. sylvatica = ", round(temp_full_bg, 3),
        " (direction ", ifelse(sign(temp_only["TS_Ha2_raw_282h"]) == sign(temp_full_bg), "same", "REVERSED"), ")")

# WTD alone (for reference)
m_wtd_only <- lmer(CH4_flux_asinh ~ bvs_wtd_cm_raw_132h * species + (1|Tree), 
                   data = model_data_scaled, REML = FALSE)
message("\nWTD only (no other predictors):")
wtd_only <- fixef(m_wtd_only)
message("  N. sylvatica: ", round(wtd_only["bvs_wtd_cm_raw_132h"], 3))
message("  T. canadensis: ", round(wtd_only["bvs_wtd_cm_raw_132h"] + wtd_only["bvs_wtd_cm_raw_132h:specieshem"], 3))
message("  A. rubrum: ", round(wtd_only["bvs_wtd_cm_raw_132h"] + wtd_only["bvs_wtd_cm_raw_132h:speciesrm"], 3))
message("  Compare to full model: N. sylvatica = ", round(wtd_full_bg, 3),
        " (direction ", ifelse(sign(wtd_only["bvs_wtd_cm_raw_132h"]) == sign(wtd_full_bg), "same", "REVERSED"), ")")

# 5. SUMMARY TABLE
message("\n--- 5. SUMMARY: MODEL COMPARISON ---")
cat(sprintf("%-40s %6s %8s\n", "Model", "R²", "AIC"))
cat("────────────────────────────────────────────────────\n")
cat(sprintf("%-40s %5.1f%% %8.1f\n", "Core: Temp × WTD × species", 
            r2_core_asinh*100, AIC(m_core_asinh)))
cat(sprintf("%-40s %5.1f%% %8.1f\n", "Full: Core + LE + SWC interactions", 
            r2_final*100, AIC(m_final)))
cat(sprintf("%-40s %5.1f%% %8.1f\n", "Alt: Temp × SWC × species", 
            var(predict(m_swc_core, re.form = NA)) / var_total_asinh * 100, AIC(m_swc_core)))
cat(sprintf("%-40s %5.1f%% %8.1f\n", "Alt: LE_Ha1 × WTD × species", 
            var(predict(m_le1_core, re.form = NA)) / var_total_asinh * 100, AIC(m_le1_core)))
cat(sprintf("%-40s %5.1f%% %8.1f\n", "Alt: LE_Ha2_anom × WTD × species", 
            var(predict(m_le2_core, re.form = NA)) / var_total_asinh * 100, AIC(m_le2_core)))
cat("────────────────────────────────────────────────────\n")
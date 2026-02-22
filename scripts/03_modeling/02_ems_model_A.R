# ============================================================
# 105a_ems_model_instantaneous.R
# 
# EMS (Upland) CH4 flux modeling - INSTANTANEOUS CORE PREDICTORS
#
# Forces TS_Ha1 and NEON_SWC_shallow with 1h window as core,
# regardless of rolling window significance
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
  flux      = "data/input/HF_2023-2025_tree_flux_corrected.csv",
  best_windows = "outputs/figures/rolling_correlations/best_windows_by_variable.csv"
)

OUTPUT_DIR <- "outputs/models/ems_instantaneous"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/figures/predictor_selection_ems", recursive = TRUE, showWarnings = FALSE)

# Site filter
SITE_FILTER <- "Upland"  # Changed from "Wetland" to "Upland"

# Correlation threshold for redundancy clustering
CORR_THRESHOLD <- 0.7

# Forward selection settings
N_RANDOM_RUNS <- 100
MAX_PREDICTORS <- 6

# Core drivers (forced into model) - based on theory
# For upland: temperature and soil moisture (not soil moisture)
# Will be selected based on data completeness below
CORE_VARIABLES <- c("TS_Ha1", "SWC_Ha1", "NEON_SWC_shallow")  # Candidates

# Variables that should be summed (not averaged) over windows
SUM_VARS <- c("P_mm", "THROUGHFALL_xHA")

# Specific anomaly variables to consider (LE and FC only)
SPECIFIC_ANOM_VARIABLES <- c("LE_Ha1", "LE_Ha2", "FC_Ha1", "FC_Ha2")

# ============================================================
# 0. DATA COMPLETENESS CHECK
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       0. DATA COMPLETENESS CHECK                 ")
message("══════════════════════════════════════════════════")

# Load data for completeness check
aligned_data_check <- read_csv(PATHS$aligned, show_col_types = FALSE)
flux_check <- read_csv(PATHS$flux, show_col_types = FALSE)

# Get upland flux period
upland_flux <- flux_check %>% 
  filter(location == SITE_FILTER) %>%
  mutate(datetime = as.POSIXct(datetime_posx))

cat(SITE_FILTER, "flux data range:\n")
cat("  From:", as.character(min(upland_flux$datetime)), "\n")
cat("  To:", as.character(max(upland_flux$datetime)), "\n")
cat("  N observations:", nrow(upland_flux), "\n\n")

# Filter to upland period
upland_period <- aligned_data_check %>%
  mutate(datetime = as.POSIXct(datetime)) %>%
  filter(datetime >= min(upland_flux$datetime),
         datetime <= max(upland_flux$datetime))

# Check temperature variables
cat("Temperature variables:\n")
cat("─────────────────────────────────────────────────────\n")

temp_vars <- c("TS_Ha1", "TS_Ha2", "s10t", "NEON_soilTemp_0cm", "NEON_soilTemp_6cm", "tair_C")
temp_completeness <- list()
for (v in temp_vars) {
  if (v %in% names(upland_period)) {
    n_total <- nrow(upland_period)
    n_valid <- sum(!is.na(upland_period[[v]]))
    pct <- round(n_valid / n_total * 100, 1)
    temp_completeness[[v]] <- pct
    cat(sprintf("  %-25s %6d / %6d (%5.1f%%)\n", v, n_valid, n_total, pct))
  }
}

cat("\nSoil moisture variables:\n")
cat("─────────────────────────────────────────────────────\n")

swc_vars <- c("SWC_Ha1", "SWC_Ha2", "NEON_SWC_shallow", "NEON_SWC_mid", "NEON_SWC_deep")
swc_completeness <- list()
for (v in swc_vars) {
  if (v %in% names(upland_period)) {
    n_total <- nrow(upland_period)
    n_valid <- sum(!is.na(upland_period[[v]]))
    pct <- round(n_valid / n_total * 100, 1)
    swc_completeness[[v]] <- pct
    cat(sprintf("  %-25s %6d / %6d (%5.1f%%)\n", v, n_valid, n_total, pct))
  }
}

cat("\n", SITE_FILTER, "species breakdown:\n")
cat("─────────────────────────────────────────────────────\n")
print(upland_flux %>% count(SPECIES))

# ════════════════════════════════════════════════════════════
# USER DECISION: Select core predictors based on completeness
# Edit these after reviewing the output above
# ════════════════════════════════════════════════════════════

CORE_TS_VAR <- "TS_Ha1"           # EDIT based on completeness
CORE_SWC_VAR <- "NEON_SWC_shallow"         # EDIT based on completeness  

cat("\n─────────────────────────────────────────────────────\n")
cat("SELECTED CORE PREDICTORS:\n")
cat("  Temperature:", CORE_TS_VAR, "\n")
cat("  Soil moisture:", CORE_SWC_VAR, "\n")
cat("─────────────────────────────────────────────────────\n")

# Update CORE_VARIABLES with selected predictors
CORE_VARIABLES <- c(CORE_TS_VAR, CORE_SWC_VAR)

# Clean up
rm(aligned_data_check, flux_check, upland_flux, upland_period)
gc()

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
  filter(site == "Upland") %>%  # All non-BGS plots (E5, etc.)
  dplyr::select(datetime, Tree, species, CH4_flux) %>%
  droplevels()  # Remove unused factor levels (e.g., 'bg' from wetland)

cat("EMS species in data:", paste(levels(stem_flux$species), collapse = ", "), "\n")

# Best windows from rolling correlation analysis
bw <- read_csv(PATHS$best_windows, show_col_types = FALSE)

cat("Loaded", nrow(aligned_data), "hourly environmental observations\n")
cat("Loaded", nrow(stem_flux), "CH4 flux measurements (Upland/EMS - all non-BGS plots)\n")
cat("Loaded", nrow(bw), "variable-window combinations from rolling analysis\n")

# ============================================================
# 2. EXTRACT SIGNIFICANT PREDICTORS
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       2. EXTRACTING SIGNIFICANT PREDICTORS       ")
message("══════════════════════════════════════════════════")

# Note: rolling analysis uses "Wetland"/"Upland"
optA_raw_sig <- bw %>%
  filter(analysis == "raw", significant, site == SITE_FILTER) %>%
  dplyr::select(variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
  arrange(desc(abs(r)))

optA_anom_sig <- bw %>%
  filter(analysis == "anomaly", significant, site == SITE_FILTER) %>%
  dplyr::select(variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
  arrange(desc(abs(r)))

cat("\n", SITE_FILTER, "significant predictors:\n")
cat("  Raw:", nrow(optA_raw_sig), "\n")
cat("  Anomaly:", nrow(optA_anom_sig), "\n")

# ============================================================
# 3. CORRELATION HEATMAPS FOR REDUNDANCY ANALYSIS
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       3. CORRELATION HEATMAPS                    ")
message("══════════════════════════════════════════════════")

# EMS Raw heatmap
if (nrow(optA_raw_sig) >= 2) {
  feat_ems_raw <- build_feature_matrix_raw(aligned_data, optA_raw_sig)
  X <- feat_ems_raw %>% dplyr::select(-datetime)
  X <- X %>% dplyr::select(where(~ sum(!is.na(.)) > 50))
  cmat_ems_raw <- cor(X, use = "pairwise.complete.obs")
  
  plot_corr_heatmap(
    cmat_ems_raw,
    "EMS (Upland) - Raw Predictors at Best Windows",
    "outputs/figures/predictor_selection_ems/corr_heatmap_EMS_raw.png"
  )
  
  redundant_ems_raw <- find_redundant_pairs(cmat_ems_raw, CORR_THRESHOLD)
  if (nrow(redundant_ems_raw) > 0) {
    cat("\nRedundant pairs (|r| >=", CORR_THRESHOLD, "):\n")
    print(redundant_ems_raw, n = 20)
    write_csv(redundant_ems_raw, "outputs/figures/predictor_selection_ems/redundant_pairs_EMS_raw.csv")
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
ems_raw <- optA_raw_sig %>% mutate(mode = "raw")
ems_anom <- optA_anom_sig %>% 
  filter(variable %in% SPECIFIC_ANOM_VARIABLES) %>%
  mutate(mode = "anom")

ems_combined <- bind_rows(ems_raw, ems_anom)

cat("Pooled predictors:", nrow(ems_raw), "raw +", nrow(ems_anom), "anomaly =", 
    nrow(ems_combined), "total\n\n")

# Build combined feature matrix
feat_ems_combined <- build_combined_feature_matrix(aligned_data, ems_combined)

X <- feat_ems_combined %>% dplyr::select(-datetime)
X <- X %>% dplyr::select(where(~ sum(!is.na(.)) > 50))

cat("Feature matrix:", ncol(X), "features\n\n")

cmat_ems_combined <- cor(X, use = "pairwise.complete.obs")

# Run selection
selected_ems <- select_predictors_combined(
  ems_combined, cmat_ems_combined, feat_ems_combined, threshold = CORR_THRESHOLD
)

cat("Selected predictors (", nrow(selected_ems), " from ", nrow(ems_combined), "):\n\n")

selected_ems %>%
  dplyr::select(variable, mode, var_group, window_hours, r, cluster_size, dropped_vars) %>%
  arrange(var_group, desc(abs(r))) %>%
  print(n = 30, width = 120)

write_csv(selected_ems, "outputs/figures/predictor_selection_ems/selected_predictors_combined.csv")

# ============================================================
# 5. THEORY-DRIVEN MODEL BUILDING
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       5. THEORY-DRIVEN MODEL BUILDING            ")
message("══════════════════════════════════════════════════")

# Build features for ALL significant raw predictors + specific anomalies
ems_raw_all <- bw %>%
  filter(analysis == "raw", significant, site == SITE_FILTER) %>%
  dplyr::select(variable, var_group, window_hours, r) %>%
  mutate(
    predictor = paste0(variable, "_raw_", window_hours, "h"),
    mode = "raw"
  )

ems_anom_subset <- bw %>%
  filter(analysis == "anomaly", significant, site == SITE_FILTER,
         variable %in% SPECIFIC_ANOM_VARIABLES) %>%
  dplyr::select(variable, var_group, window_hours, r) %>%
  mutate(
    predictor = paste0(variable, "_anom_", window_hours, "h"),
    mode = "anom"
  )

cat("Significant raw predictors from rolling analysis:", nrow(ems_raw_all), "\n")
cat("Specific anomaly predictors (LE/FC):", nrow(ems_anom_subset), "\n")

# ════════════════════════════════════════════════════════════
# FORCE INSTANTANEOUS CORE PREDICTORS (1h window)
# These are added regardless of rolling window significance
# ════════════════════════════════════════════════════════════

FORCED_CORE_WINDOW <- 1  # Instantaneous (1 hour)

forced_core <- tibble(
  variable = c(CORE_TS_VAR, CORE_SWC_VAR),
  var_group = c("Soil Climate", "Soil Climate"),
  window_hours = FORCED_CORE_WINDOW,
  r = NA_real_,  # Not from rolling analysis
  predictor = paste0(variable, "_raw_", FORCED_CORE_WINDOW, "h"),
  mode = "raw"
)

cat("\n⚠ FORCING instantaneous core predictors (1h window):\n")
print(forced_core %>% dplyr::select(predictor, variable))

# Core info is the forced predictors
core_info <- forced_core

# Candidates = all significant raw (none are core since core wasn't significant) + anomalies
candidate_info <- bind_rows(ems_raw_all, ems_anom_subset)

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

# FORCE: Add core predictors with instantaneous window
for (i in seq_len(nrow(forced_core))) {
  v <- forced_core$variable[i]
  w <- forced_core$window_hours[i]
  if (!v %in% names(met)) {
    cat("  Warning:", v, "not found in met data\n")
    next
  }
  feat_name <- paste0(v, "_raw_", w, "h")
  features[[feat_name]] <- roll_mean(met[[v]], w)
  cat("  Built forced feature:", feat_name, "\n")
}

# Raw features from significant predictors
for (i in seq_len(nrow(ems_raw_all))) {
  v <- ems_raw_all$variable[i]
  w <- ems_raw_all$window_hours[i]
  if (!v %in% names(met)) next
  feat_name <- paste0(v, "_raw_", w, "h")
  if (feat_name %in% names(features)) next  # Skip if already built
  if (v %in% SUM_VARS) {
    features[[feat_name]] <- roll_sum(met[[v]], w)
  } else {
    features[[feat_name]] <- roll_mean(met[[v]], w)
  }
}

# Anomaly features
for (i in seq_len(nrow(ems_anom_subset))) {
  v <- ems_anom_subset$variable[i]
  w <- ems_anom_subset$window_hours[i]
  src <- paste0(v, "_anom")
  if (!src %in% names(met)) next
  feat_name <- paste0(v, "_anom_", w, "h")
  features[[feat_name]] <- roll_mean(met[[src]], w)
}

cat("\nTotal features built:", ncol(features) - 1, "\n")

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

write_csv(selection_freq, "outputs/figures/predictor_selection_ems/theory_driven_selection_freq.csv")

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

# Identify core interaction predictors dynamically based on CORE_TS_VAR and CORE_SWC_VAR
ts_pattern <- paste0("^", CORE_TS_VAR, "_raw_")
swc_pattern <- paste0("^", CORE_SWC_VAR, "_raw_")

ts_pred <- names(model_data_scaled)[grepl(ts_pattern, names(model_data_scaled))][1]
swc_pred <- names(model_data_scaled)[grepl(swc_pattern, names(model_data_scaled))][1]

cat("Core interaction predictors (theory-driven):\n")
cat("  Temperature:", ts_pred, "\n")
cat("  Soil moisture:", swc_pred, "\n\n")

if (is.na(ts_pred) || is.na(swc_pred)) {
  stop("Required core predictors (", CORE_TS_VAR, " and/or ", CORE_SWC_VAR, ") not found!")
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
  filter(GVIF > 10, !predictor %in% c(ts_pred, swc_pred, "species")) %>%
  pull(predictor)

if (length(high_vif_preds) > 0) {
  cat("\n⚠ Predictors with VIF > 10 (candidates for removal):\n")
  for (p in high_vif_preds) {
    cat("  -", p, "(VIF =", vif_df$GVIF[vif_df$predictor == p], ")\n")
  }
  
  # Show correlation matrix among high-VIF predictors + core predictors
  # This helps identify which variables are collinear with each other
  check_preds <- c(ts_pred, swc_pred, high_vif_preds)
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
  cat("For each correlated pair, keep the variable that is:\n")
  cat("  1. More mechanistically relevant to CH4 flux\n")
  cat("  2. More directly measured (vs derived)\n")
  cat("  3. Core to your hypothesis (", CORE_TS_VAR, ",", CORE_SWC_VAR, ")\n")
  cat("─────────────────────────────────────────────────────\n")
  
  # ════════════════════════════════════════════════════════════
  # USER DECISION: Which high-VIF predictors to DROP?
  # Based on the correlation matrix above, choose which member
  # of each correlated pair to remove
  # ════════════════════════════════════════════════════════════
  
  VIF_DROP <- c("USTAR_Ha2_raw_333h")  # EDIT THIS LIST AFTER SEEING OUTPUT
  
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
additional_preds <- setdiff(refined_preds, c(ts_pred, swc_pred))
cat("\nAdditional predictors beyond core:", length(additional_preds), "\n")
if (length(additional_preds) > 0) {
  cat(" ", paste(additional_preds, collapse = ", "), "\n")
}

# ============================================================
# 9. CORE INTERACTION MODEL (TS * SWC * SPECIES)
# ============================================================

message("\n")
message("══════════════════════════════════════════════════")
message("       9. CORE INTERACTION MODEL                  ")
message("══════════════════════════════════════════════════")

# Model 1: Core 3-way interaction only
formula_core <- paste0("CH4_flux ~ ", ts_pred, " * ", swc_pred, " * species + (1|Tree)")
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
    formula_add <- paste0("CH4_flux ~ ", ts_pred, " * ", swc_pred, " * species + ", 
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
    formula_add_sp <- paste0("CH4_flux ~ ", ts_pred, " * ", swc_pred, " * species + ", 
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
    formula_extended <- paste0("CH4_flux ~ ", ts_pred, " * ", swc_pred, " * species + ",
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
cat(sprintf("%-35s %5.1f%% %8.1f %8.1f\n", "Core (TS*SWC*species)", 
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

# Get actual species levels from the data
species_levels <- levels(model_data_scaled$species)
ref_species <- species_levels[1]  # Reference level (alphabetically first)
other_species <- species_levels[-1]

# Species name mapping (used throughout remaining sections)
species_names <- c(
  "bg" = "Black gum",
  "hem" = "Hemlock", 
  "rm" = "Red maple",
  "ro" = "Red oak"
)

cat("Species in model:\n")
cat("  Reference:", ref_species, "(", species_names[ref_species], ")\n")
cat("  Others:", paste(other_species, collapse = ", "), "\n\n")

# Extract slopes for reference species
temp_ref <- coefs[ts_pred]
swc_ref <- coefs[swc_pred]
int_term <- paste0(ts_pred, ":", swc_pred)
int_ref <- ifelse(int_term %in% names(coefs), coefs[int_term], 0)

# Build species slopes dynamically
species_slopes <- tibble(
  species = ref_species,
  temp_slope = temp_ref,
  swc_slope = swc_ref,
  interaction = int_ref
)

for (sp in other_species) {
  temp_sp_term <- paste0(ts_pred, ":species", sp)
  swc_sp_term <- paste0(swc_pred, ":species", sp)
  int_sp_term <- paste0(ts_pred, ":", swc_pred, ":species", sp)
  
  temp_sp <- temp_ref + ifelse(temp_sp_term %in% names(coefs), coefs[temp_sp_term], 0)
  swc_sp <- swc_ref + ifelse(swc_sp_term %in% names(coefs), coefs[swc_sp_term], 0)
  int_sp <- int_ref + ifelse(int_sp_term %in% names(coefs), coefs[int_sp_term], 0)
  
  species_slopes <- bind_rows(species_slopes, tibble(
    species = sp,
    temp_slope = temp_sp,
    swc_slope = swc_sp,
    interaction = int_sp
  ))
}

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
  setNames("Soil moisture", swc_pred)
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

# Get reference species name for labels
ref_species_name <- species_names[ref_species]
if (is.na(ref_species_name)) ref_species_name <- ref_species

coef_df <- coef_df %>%
  mutate(
    term_clean = term_clean %>%
      str_replace_all("specieshem", "Hemlock") %>%
      str_replace_all("speciesrm", "Red maple") %>%
      str_replace_all("speciesro", "Red oak") %>%
      str_replace_all("speciesbg", "Black gum") %>%
      str_replace_all(":", " × "),
    # Add reference species info where needed
    term_clean = case_when(
      term_clean == "Temperature" ~ paste0("Temperature (", ref_species_name, ")"),
      term_clean == "Soil moisture" ~ paste0("Soil moisture (", ref_species_name, ")"),
      term_clean == "Hemlock" ~ "Hemlock (intercept)",
      term_clean == "Red maple" ~ "Red maple (intercept)",
      term_clean == "Red oak" ~ "Red oak (intercept)",
      term_clean == "Black gum" ~ "Black gum (intercept)",
      term_clean == "Temperature × Soil moisture" ~ paste0("Temp × SWC (", ref_species_name, ")"),
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
  SWC_mean = mean(model_data_unscaled[[swc_pred]], na.rm = TRUE),
  SWC_sd = sd(model_data_unscaled[[swc_pred]], na.rm = TRUE)
)

# Prediction grid
temp_range_scaled <- seq(-2, 2, length.out = 100)
swc_percentiles <- c(10, 25, 50, 75, 90)

swc_quantiles_unscaled <- quantile(model_data_unscaled[[swc_pred]],
                                   probs = swc_percentiles/100, na.rm = TRUE)
swc_quantiles_scaled <- (swc_quantiles_unscaled - scaling_params$SWC_mean) / scaling_params$SWC_sd

temp_unscaled <- temp_range_scaled * scaling_params$TS_sd + scaling_params$TS_mean

# Get all predictor names needed for prediction (from the model terms)
model_vars <- all.vars(formula(m_final))
model_vars <- setdiff(model_vars, c("CH4_flux_asinh", "Tree", "species"))

# Get actual species levels from data (species_names already defined in section 13)
species_levels <- levels(model_data_scaled$species)

# Generate predictions
pred_list <- list()

for (sp in species_levels) {
  for (i in seq_along(swc_percentiles)) {
    # Build prediction dataframe - start with the right number of rows
    n_pts <- length(temp_range_scaled)
    
    pred_df <- data.frame(
      species = factor(rep(sp, n_pts), levels = species_levels),
      Tree = factor(rep(model_data_scaled$Tree[1], n_pts), levels = levels(model_data_scaled$Tree))
    )
    
    # Add core predictors
    pred_df[[ts_pred]] <- temp_range_scaled
    pred_df[[swc_pred]] <- rep(swc_quantiles_scaled[i], n_pts)
    
    # Add any additional predictors at their mean (0 for scaled data)
    for (v in model_vars) {
      if (!v %in% names(pred_df) && v %in% names(model_data_scaled)) {
        pred_df[[v]] <- rep(0, n_pts)  # Mean for scaled predictors
      }
    }
    
    pred_df$fit_asinh <- predict(m_final, newdata = pred_df, re.form = NA)
    pred_df$fit <- sinh(pred_df$fit_asinh) / 1000  # back-transform
    
    pred_df$temp_unscaled <- temp_unscaled
    pred_df$swc_percentile <- swc_percentiles[i]
    pred_df$species_name <- ifelse(sp %in% names(species_names), species_names[sp], sp)
    
    pred_list[[length(pred_list) + 1]] <- pred_df
  }
}

pred_all <- bind_rows(pred_list) %>%
  mutate(
    species_name = factor(species_name),
    swc_percentile = factor(swc_percentile)
  )

# Colors
swc_colors <- c("10" = "#8B4513", "25" = "#CD853F", "50" = "#808080",
                "75" = "#6495ED", "90" = "#0000CD")

# Dynamic species colors
unique_species <- unique(pred_all$species_name)
species_color_palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")
species_colors <- setNames(species_color_palette[1:length(unique_species)], unique_species)

# Panel plot by species
p_interaction <- ggplot(pred_all, aes(x = temp_unscaled, y = fit * 1000,
                                      color = swc_percentile)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~ species_name, scales = "free_y") +
  scale_color_manual(values = swc_colors, labels = paste0(swc_percentiles, "th"),
                     name = "Soil moisture\npercentile") +
  labs(
    x = "Soil temperature (°C)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Temperature × soil moisture interactions by species"
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
  filter(swc_percentile %in% c("10", "90"))

p_comparison <- ggplot(pred_extremes, aes(x = temp_unscaled, y = fit * 1000,
                                          color = species_name, linetype = swc_percentile)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = species_colors, name = "Species") +
  scale_linetype_manual(values = c("10" = "dashed", "90" = "solid"),
                        labels = c("10" = "Dry (10th)", "90" = "Wet (90th)"),
                        name = "Soil moisture") +
  labs(
    x = "Soil temperature (°C)",
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
    title = "Species differ in temperature × soil moisture response"
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
message("Site: EMS (Upland)")
message("")
message("Final model:")
message("  ", formula_final)
message("")
message(sprintf("R² = %.1f%%", r2_final * 100))
message(sprintf("AIC = %.1f", AIC(m_final)))
message(sprintf("N = %d observations from %d trees",
                nrow(model_data_scaled), n_distinct(model_data_scaled$Tree)))
message("")
message("Species-specific slopes (see species_slopes.csv for details):")
for (i in 1:nrow(species_slopes)) {
  message(sprintf("  %s:  Temp=%.3f  SWC=%.3f  Interaction=%.3f",
                  species_slopes$species[i], 
                  species_slopes$temp_slope[i], 
                  species_slopes$swc_slope[i], 
                  species_slopes$interaction[i]))
}
message("")
message("Output directory: ", OUTPUT_DIR)
message("══════════════════════════════════════════════════")

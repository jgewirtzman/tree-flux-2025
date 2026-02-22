# ============================================
# Extract Option A significant variables
# Works with output from 08_rolling_window_correlations.R
# ============================================

library(tidyverse)

# Your script saves these objects:
# - best_windows_raw (with 'significant' column based on BH-FDR)
# - best_windows_anom (with 'significant' column based on BH-FDR)
# 
# OR if you only have the combined version:
# - best_windows_all (with 'analysis' column = "raw" or "anomaly")

# If you have best_windows_all (the combined version):
if (exists("best_windows_all") || exists("best_windows")) {
  
  bw <- if (exists("best_windows_all")) best_windows_all else best_windows
  
  cat("\n========== OPTION A: RAW PREDICTORS ==========\n")
  
  optA_raw_sig <- bw %>%
    filter(analysis == "raw", significant) %>%
    dplyr::select(site, variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
    arrange(site, desc(abs(r)))
  
  cat("\nBGS (wetland) - Raw significant:\n")
  print(optA_raw_sig %>% filter(site == "BGS"), n = 50)
  
  cat("\nEMS (upland) - Raw significant:\n")
  print(optA_raw_sig %>% filter(site == "EMS"), n = 50)
  
  cat("\n========== OPTION A: ANOMALY PREDICTORS ==========\n")
  
  optA_anom_sig <- bw %>%
    filter(analysis == "anomaly", significant) %>%
    dplyr::select(site, variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
    arrange(site, desc(abs(r)))
  
  cat("\nBGS (wetland) - Anomaly significant:\n")
  print(optA_anom_sig %>% filter(site == "BGS"), n = 50)
  
  cat("\nEMS (upland) - Anomaly significant:\n")
  print(optA_anom_sig %>% filter(site == "EMS"), n = 50)
  
} else if (exists("best_windows_raw") && exists("best_windows_anom")) {
  
  # If you have the separate objects
  cat("\n========== OPTION A: RAW PREDICTORS ==========\n")
  
  optA_raw_sig <- best_windows_raw %>%
    filter(significant) %>%
    dplyr::select(site, variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
    arrange(site, desc(abs(r)))
  
  cat("\nBGS (wetland) - Raw significant:\n")
  print(optA_raw_sig %>% filter(site == "BGS"), n = 50)
  
  cat("\nEMS (upland) - Raw significant:\n")
  print(optA_raw_sig %>% filter(site == "EMS"), n = 50)
  
  cat("\n========== OPTION A: ANOMALY PREDICTORS ==========\n")
  
  optA_anom_sig <- best_windows_anom %>%
    filter(significant) %>%
    dplyr::select(site, variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
    arrange(site, desc(abs(r)))
  
  cat("\nBGS (wetland) - Anomaly significant:\n")
  print(optA_anom_sig %>% filter(site == "BGS"), n = 50)
  
  cat("\nEMS (upland) - Anomaly significant:\n")
  print(optA_anom_sig %>% filter(site == "EMS"), n = 50)
  
} else {
  cat("ERROR: Cannot find best_windows objects. Run 08_rolling_window_correlations.R first.\n")
  cat("Looking for: best_windows_all, best_windows, best_windows_raw, best_windows_anom\n")
}

# Summary counts
cat("\n========== SUMMARY ==========\n")
cat("BGS Raw:", nrow(optA_raw_sig %>% filter(site == "BGS")), "significant variables\n")
cat("BGS Anom:", nrow(optA_anom_sig %>% filter(site == "BGS")), "significant variables\n")
cat("EMS Raw:", nrow(optA_raw_sig %>% filter(site == "EMS")), "significant variables\n")
cat("EMS Anom:", nrow(optA_anom_sig %>% filter(site == "EMS")), "significant variables\n")

# Overlap analysis
cat("\n========== OVERLAP (sig in both raw AND anomaly) ==========\n")
for (s in c("BGS", "EMS")) {
  raw_vars <- optA_raw_sig %>% filter(site == s) %>% pull(variable)
  anom_vars <- optA_anom_sig %>% filter(site == s) %>% pull(variable)
  both <- intersect(raw_vars, anom_vars)
  raw_only <- setdiff(raw_vars, anom_vars)
  anom_only <- setdiff(anom_vars, raw_vars)
  
  cat("\n", s, ":\n", sep = "")
  cat("  Both (", length(both), "): ", paste(both, collapse = ", "), "\n", sep = "")
  cat("  Raw only (", length(raw_only), "): ", paste(raw_only, collapse = ", "), "\n", sep = "")
  cat("  Anom only (", length(anom_only), "): ", paste(anom_only, collapse = ", "), "\n", sep = "")
}

# OR if you saved to CSV, read from there:
# best_windows_all <- read_csv("figures/rolling_correlations/best_windows_by_variable.csv")














# ============================================
# 15_predictor_correlation_heatmaps.R
#
# Build correlation heatmaps of significant predictors
# at their Option A best windows to identify redundancy
# ============================================

library(tidyverse)
library(RcppRoll)
library(pheatmap)
library(viridis)

# ============================================
# CONFIG
# ============================================

OUTPUT_DIR <- "figures/predictor_select"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Correlation threshold for "redundant"
CORR_THRESHOLD <- 0.7

# Variables that should be summed (not averaged)
SUM_VARS <- c("P_mm", "THROUGHFALL_xHA")

# ============================================
# HELPER FUNCTIONS
# ============================================

# Rolling stats
roll_mean <- function(x, n) RcppRoll::roll_mean(x, n = n, align = "right", fill = NA, na.rm = TRUE)
roll_sum  <- function(x, n) RcppRoll::roll_sum(x, n = n, align = "right", fill = NA, na.rm = TRUE)

# Build feature matrix with each variable at its own best window
build_feature_matrix <- function(aligned_data, sig_vars_df, mode = c("raw", "anomaly")) {
  mode <- match.arg(mode)
  
  met <- aligned_data %>% arrange(datetime)
  
  # If anomaly mode, create anomaly versions
  if (mode == "anomaly") {
    for (v in unique(sig_vars_df$variable)) {
      if (!v %in% names(met)) next
      met <- make_anomaly(met, v, time_col = "datetime")
    }
  }
  
  # Build rolled features
  result <- met %>% dplyr::select(datetime)
  
  for (i in seq_len(nrow(sig_vars_df))) {
    v <- sig_vars_df$variable[i]
    w <- sig_vars_df$window_hours[i]
    
    # Source column
    src <- if (mode == "anomaly") paste0(v, "_anom") else v
    if (!src %in% names(met)) next
    
    # Feature name includes window
    feat_name <- paste0(v, "_", w, "h")
    
    if (v %in% SUM_VARS) {
      result[[feat_name]] <- roll_sum(met[[src]], w)
    } else {
      result[[feat_name]] <- roll_mean(met[[src]], w)
    }
  }
  
  result
}

# Anomaly function (same as in 08)
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

# Plot correlation heatmap
plot_corr_heatmap <- function(cmat, title, filename) {
  # Reorder by hierarchical clustering
  if (nrow(cmat) > 2) {
    hc <- hclust(as.dist(1 - abs(cmat)), method = "complete")
    ord <- hc$order
    cmat <- cmat[ord, ord]
  }
  
  # Clean up names for display (remove _XXh suffix for readability)
  rownames(cmat) <- gsub("_\\d+h$", "", rownames(cmat))
  colnames(cmat) <- gsub("_\\d+h$", "", colnames(cmat))
  
  # Plot
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
    cluster_rows = FALSE,  # Already ordered
    cluster_cols = FALSE,
    border_color = "grey90"
  )
  dev.off()
  
  message("  Saved: ", basename(filename))
  
  invisible(cmat)
}

# Identify redundant pairs
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

# ============================================
# LOAD DATA
# ============================================

message("Loading data...")

# Assumes these exist from running 08_rolling_window_correlations.R:
# - aligned_data (hourly environmental data)
# - optA_raw_sig, optA_anom_sig (significant variables from extraction)

# If not in memory, load from files:
if (!exists("aligned_data")) {
  aligned_data <- read_csv("data/processed/aligned_hourly_dataset.csv", show_col_types = FALSE) %>%
    mutate(datetime = as.POSIXct(datetime, tz = "UTC"))
}

if (!exists("optA_raw_sig") || !exists("optA_anom_sig")) {
  bw <- read_csv("figures/rolling_correlations/best_windows_by_variable.csv", show_col_types = FALSE)
  
  optA_raw_sig <- bw %>%
    filter(analysis == "raw", significant) %>%
    dplyr::select(site, variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
    arrange(site, desc(abs(r)))
  
  optA_anom_sig <- bw %>%
    filter(analysis == "anomaly", significant) %>%
    dplyr::select(site, variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
    arrange(site, desc(abs(r)))
}

# ============================================
# BGS (WETLAND) - RAW
# ============================================

message("\n========== BGS RAW ==========")

bgs_raw_vars <- optA_raw_sig %>% filter(site == "BGS")

if (nrow(bgs_raw_vars) >= 2) {
  # Build feature matrix
  feat_bgs_raw <- build_feature_matrix(aligned_data, bgs_raw_vars, mode = "raw")
  
  # Correlation matrix
  X <- feat_bgs_raw %>% dplyr::select(-datetime)
  X <- X %>% dplyr::select(where(~ sum(!is.na(.)) > 50))
  
  cmat_bgs_raw <- cor(X, use = "pairwise.complete.obs")
  
  # Plot
  plot_corr_heatmap(
    cmat_bgs_raw,
    "BGS (Wetland) - Raw Predictors at Best Windows",
    file.path(OUTPUT_DIR, "corr_heatmap_BGS_raw.png")
  )
  
  # Find redundant pairs
  redundant_bgs_raw <- find_redundant_pairs(cmat_bgs_raw, CORR_THRESHOLD)
  
  if (nrow(redundant_bgs_raw) > 0) {
    cat("\nRedundant pairs (|r| >=", CORR_THRESHOLD, "):\n")
    print(redundant_bgs_raw, n = 30)
  }
}

# ============================================
# BGS (WETLAND) - ANOMALY
# ============================================

message("\n========== BGS ANOMALY ==========")

bgs_anom_vars <- optA_anom_sig %>% filter(site == "BGS")

if (nrow(bgs_anom_vars) >= 2) {
  feat_bgs_anom <- build_feature_matrix(aligned_data, bgs_anom_vars, mode = "anomaly")
  
  X <- feat_bgs_anom %>% dplyr::select(-datetime)
  X <- X %>% dplyr::select(where(~ sum(!is.na(.)) > 50))
  
  cmat_bgs_anom <- cor(X, use = "pairwise.complete.obs")
  
  plot_corr_heatmap(
    cmat_bgs_anom,
    "BGS (Wetland) - Anomaly Predictors at Best Windows",
    file.path(OUTPUT_DIR, "corr_heatmap_BGS_anomaly.png")
  )
  
  redundant_bgs_anom <- find_redundant_pairs(cmat_bgs_anom, CORR_THRESHOLD)
  
  if (nrow(redundant_bgs_anom) > 0) {
    cat("\nRedundant pairs (|r| >=", CORR_THRESHOLD, "):\n")
    print(redundant_bgs_anom, n = 30)
  }
}

# ============================================
# EMS (UPLAND) - ANOMALY ONLY
# ============================================

message("\n========== EMS ANOMALY ==========")

ems_anom_vars <- optA_anom_sig %>% filter(site == "EMS")

if (nrow(ems_anom_vars) >= 2) {
  feat_ems_anom <- build_feature_matrix(aligned_data, ems_anom_vars, mode = "anomaly")
  
  X <- feat_ems_anom %>% dplyr::select(-datetime)
  X <- X %>% dplyr::select(where(~ sum(!is.na(.)) > 50))
  
  cmat_ems_anom <- cor(X, use = "pairwise.complete.obs")
  
  plot_corr_heatmap(
    cmat_ems_anom,
    "EMS (Upland) - Anomaly Predictors at Best Windows",
    file.path(OUTPUT_DIR, "corr_heatmap_EMS_anomaly.png")
  )
  
  redundant_ems_anom <- find_redundant_pairs(cmat_ems_anom, CORR_THRESHOLD)
  
  if (nrow(redundant_ems_anom) > 0) {
    cat("\nRedundant pairs (|r| >=", CORR_THRESHOLD, "):\n")
    print(redundant_ems_anom, n = 30)
  }
} else {
  message("  Only ", nrow(ems_anom_vars), " significant variables - skipping heatmap")
}

# ============================================
# SUMMARY: GROUPED BY PROCESS
# ============================================

message("\n========== PROCESS GROUPINGS FOR FINAL MODEL ==========\n")

cat("Based on Option A results, here are the candidate predictors by process:\n\n")

# Define process groups
process_groups <- list(
  "Water Table" = c("bvs_wtd_cm", "bgs_wtd_cm"),
  "Soil Temperature" = c("TS_Ha1", "TS_Ha2", "s10t", "T_CANOPY_xHA", "tair_C"),
  "Soil Moisture" = c("NEON_SWC_shallow", "NEON_SWC_mid", "NEON_SWC_deep", "SWC_Ha1", "SWC_Ha2"),
  "Radiation" = c("PAR", "slrr", "rnet"),
  "Atmospheric" = c("VPD_kPa", "RH", "p_kPa", "USTAR_Ha1", "USTAR_Ha2"),
  "Ecosystem Fluxes" = c("LE_Ha1", "LE_Ha2", "H_Ha1", "H_Ha2", "FC_Ha1", "FC_Ha2"),
  "Gases" = c("CO2_MR_xHA", "CO2_MR_Ha2", "CH4_MR_xHA"),
  "Precipitation" = c("P_mm", "THROUGHFALL_xHA"),
  "Phenology" = c("gcc", "ndvi")
)

for (grp in names(process_groups)) {
  grp_vars <- process_groups[[grp]]
  
  # BGS raw
  bgs_raw_in_grp <- bgs_raw_vars %>% filter(variable %in% grp_vars)
  # BGS anom
  bgs_anom_in_grp <- bgs_anom_vars %>% filter(variable %in% grp_vars)
  
  if (nrow(bgs_raw_in_grp) > 0 || nrow(bgs_anom_in_grp) > 0) {
    cat("---", grp, "---\n")
    
    if (nrow(bgs_raw_in_grp) > 0) {
      cat("  BGS Raw:\n")
      for (i in seq_len(nrow(bgs_raw_in_grp))) {
        cat(sprintf("    %s: r=%.3f @ %dh (%.1fd)\n",
                    bgs_raw_in_grp$variable[i],
                    bgs_raw_in_grp$r[i],
                    bgs_raw_in_grp$window_hours[i],
                    bgs_raw_in_grp$window_days[i]))
      }
    }
    
    if (nrow(bgs_anom_in_grp) > 0) {
      cat("  BGS Anomaly:\n")
      for (i in seq_len(nrow(bgs_anom_in_grp))) {
        cat(sprintf("    %s: r=%.3f @ %dh (%.1fd)\n",
                    bgs_anom_in_grp$variable[i],
                    bgs_anom_in_grp$r[i],
                    bgs_anom_in_grp$window_hours[i],
                    bgs_anom_in_grp$window_days[i]))
      }
    }
    cat("\n")
  }
}

# ============================================
# SAVE RESULTS
# ============================================

# Save redundant pairs
if (exists("redundant_bgs_raw") && nrow(redundant_bgs_raw) > 0) {
  write_csv(redundant_bgs_raw, file.path(OUTPUT_DIR, "redundant_pairs_BGS_raw.csv"))
}
if (exists("redundant_bgs_anom") && nrow(redundant_bgs_anom) > 0) {
  write_csv(redundant_bgs_anom, file.path(OUTPUT_DIR, "redundant_pairs_BGS_anomaly.csv"))
}
if (exists("redundant_ems_anom") && nrow(redundant_ems_anom) > 0) {
  write_csv(redundant_ems_anom, file.path(OUTPUT_DIR, "redundant_pairs_EMS_anomaly.csv"))
}

message("\nOutputs saved to: ", OUTPUT_DIR)
message("Done!")





















# ============================================
# 16b_predictor_selection_combined.R
#
# Combined raw + anomaly predictor selection
# 
# Pools all significant predictors from both raw and anomaly
# analyses, then runs correlation clustering on the combined set.
# This allows the algorithm to choose the best representation
# (raw vs anomaly) for each process, or keep both if uncorrelated.
# ============================================

library(tidyverse)
library(RcppRoll)

# ============================================
# CONFIG
# ============================================

OUTPUT_DIR <- "figures/predictor_selection"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Selection parameters (DOCUMENT THESE IN METHODS)
CORR_THRESHOLD <- 0.7    # |r| >= this = "redundant"
SUM_VARS <- c("P_mm", "THROUGHFALL_xHA")

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

# Build combined feature matrix with both raw and anomaly versions
build_combined_feature_matrix <- function(aligned_data, sig_vars_df) {
  # sig_vars_df must have columns: variable, window_hours, mode (raw/anomaly)
  
  met <- aligned_data %>% arrange(datetime)
  
  # Create anomaly versions for all variables that need them
  anom_vars <- sig_vars_df %>% filter(mode == "anomaly") %>% pull(variable) %>% unique()
  for (v in anom_vars) {
    if (!v %in% names(met)) next
    met <- make_anomaly(met, v, time_col = "datetime")
  }
  
  # Build rolled features
  result <- met %>% dplyr::select(datetime)
  
  for (i in seq_len(nrow(sig_vars_df))) {
    v <- sig_vars_df$variable[i]
    w <- sig_vars_df$window_hours[i]
    m <- sig_vars_df$mode[i]
    
    # Source column depends on mode
    src <- if (m == "anomaly") paste0(v, "_anom") else v
    if (!src %in% names(met)) next
    
    # Feature name includes mode and window
    feat_name <- paste0(v, "_", m, "_", w, "h")
    
    if (v %in% SUM_VARS) {
      result[[feat_name]] <- roll_sum(met[[src]], w)
    } else {
      result[[feat_name]] <- roll_mean(met[[src]], w)
    }
  }
  
  result
}

# Selection algorithm (same as before but handles combined naming)
select_predictors_combined <- function(sig_vars_df, cmat, feat_matrix, threshold = 0.7) {
  
  vars_in_cmat <- colnames(cmat)
  
  # Map back to sig_vars_df
  sig_vars_df <- sig_vars_df %>%
    mutate(feat_name = paste0(variable, "_", mode, "_", window_hours, "h"))
  
  # Only keep variables that are in both
  sig_vars_df <- sig_vars_df %>%
    filter(feat_name %in% vars_in_cmat)
  
  if (nrow(sig_vars_df) == 0) {
    return(tibble())
  }
  
  # Compute missingness for each feature
  missingness <- feat_matrix %>%
    dplyr::select(-datetime) %>%
    summarize(across(everything(), ~ mean(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "feat_name", values_to = "missing_frac")
  
  sig_vars_df <- sig_vars_df %>%
    left_join(missingness, by = "feat_name")
  
  # Build adjacency matrix (|r| >= threshold)
  adj <- abs(cmat) >= threshold
  diag(adj) <- FALSE
  
  # Find connected components using simple BFS
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
  
  # Map cluster IDs back to variables
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
  
  # Add info about what was dropped
  dropped <- sig_vars_df %>%
    filter(!feat_name %in% selected$feat_name) %>%
    dplyr::select(variable, mode, feat_name, r, cluster)
  
  # For each selected variable, note what it beat
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

# ============================================
# LOAD DATA
# ============================================

message("Loading data...")

if (!exists("aligned_data")) {
  aligned_data <- read_csv("data/processed/aligned_hourly_dataset.csv", show_col_types = FALSE) %>%
    mutate(datetime = as.POSIXct(datetime, tz = "UTC"))
}

if (!exists("optA_raw_sig") || !exists("optA_anom_sig")) {
  bw <- read_csv("figures/rolling_correlations/best_windows_by_variable.csv", show_col_types = FALSE)
  
  optA_raw_sig <- bw %>%
    filter(analysis == "raw", significant) %>%
    dplyr::select(site, variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
    arrange(site, desc(abs(r)))
  
  optA_anom_sig <- bw %>%
    filter(analysis == "anomaly", significant) %>%
    dplyr::select(site, variable, var_label, var_group, window_hours, window_days, r, p_adj) %>%
    arrange(site, desc(abs(r)))
}

# ============================================
# BGS: COMBINE RAW + ANOMALY
# ============================================

message("\n========== BGS COMBINED (RAW + ANOMALY) ==========\n")

# Pool significant predictors from both modes
bgs_raw <- optA_raw_sig %>% 
  filter(site == "BGS") %>%
  mutate(mode = "raw")

bgs_anom <- optA_anom_sig %>% 
  filter(site == "BGS") %>%
  mutate(mode = "anom")

bgs_combined <- bind_rows(bgs_raw, bgs_anom)

cat("Pooled predictors:", nrow(bgs_raw), "raw +", nrow(bgs_anom), "anomaly =", nrow(bgs_combined), "total\n\n")

# Check for variables significant in both modes
both_modes <- bgs_raw %>%
  inner_join(bgs_anom %>% dplyr::select(variable), by = "variable") %>%
  pull(variable) %>%
  unique()

cat("Variables significant in BOTH modes:", length(both_modes), "\n")
cat(" ", paste(both_modes, collapse = ", "), "\n\n")

# Build combined feature matrix
feat_bgs_combined <- build_combined_feature_matrix(aligned_data, bgs_combined)

X <- feat_bgs_combined %>% dplyr::select(-datetime)
X <- X %>% dplyr::select(where(~ sum(!is.na(.)) > 50))

cat("Feature matrix:", ncol(X), "features\n\n")

cmat_bgs_combined <- cor(X, use = "pairwise.complete.obs")

# Run selection
selected_bgs_combined <- select_predictors_combined(
  bgs_combined, cmat_bgs_combined, feat_bgs_combined, threshold = CORR_THRESHOLD
) %>%
  mutate(site = "BGS")

cat("Selected predictors (", nrow(selected_bgs_combined), " from ", nrow(bgs_combined), "):\n\n")

selected_bgs_combined %>%
  dplyr::select(variable, mode, var_group, window_hours, r, cluster_size, dropped_vars) %>%
  arrange(var_group, desc(abs(r))) %>%
  print(n = 40, width = 120)

# ============================================
# EMS: COMBINE RAW + ANOMALY
# ============================================

message("\n========== EMS COMBINED (RAW + ANOMALY) ==========\n")

ems_raw <- optA_raw_sig %>% 
  filter(site == "EMS") %>%
  mutate(mode = "raw")

ems_anom <- optA_anom_sig %>% 
  filter(site == "EMS") %>%
  mutate(mode = "anom")

ems_combined <- bind_rows(ems_raw, ems_anom)

cat("Pooled predictors:", nrow(ems_raw), "raw +", nrow(ems_anom), "anomaly =", nrow(ems_combined), "total\n\n")

if (nrow(ems_combined) >= 2) {
  feat_ems_combined <- build_combined_feature_matrix(aligned_data, ems_combined)
  
  X <- feat_ems_combined %>% dplyr::select(-datetime)
  X <- X %>% dplyr::select(where(~ sum(!is.na(.)) > 50))
  
  cmat_ems_combined <- cor(X, use = "pairwise.complete.obs")
  
  selected_ems_combined <- select_predictors_combined(
    ems_combined, cmat_ems_combined, feat_ems_combined, threshold = CORR_THRESHOLD
  ) %>%
    mutate(site = "EMS")
  
  cat("Selected predictors (", nrow(selected_ems_combined), " from ", nrow(ems_combined), "):\n\n")
  
  selected_ems_combined %>%
    dplyr::select(variable, mode, var_group, window_hours, r, cluster_size, dropped_vars) %>%
    arrange(var_group, desc(abs(r))) %>%
    print(n = 20, width = 120)
} else {
  selected_ems_combined <- tibble()
  message("  Insufficient predictors")
}

# ============================================
# SUMMARY
# ============================================

cat("\n\n========== FINAL PREDICTOR SETS ==========\n\n")

cat("--- BGS (Wetland) Combined ---\n")
selected_bgs_combined %>%
  dplyr::select(variable, mode, var_group, window_hours, r) %>%
  mutate(
    window = paste0(round(window_hours/24, 1), "d"),
    r = round(r, 3)
  ) %>%
  dplyr::select(Process = var_group, Variable = variable, Mode = mode, Window = window, r) %>%
  arrange(desc(abs(r))) %>%
  print(n = 30)

# Mode breakdown
cat("\nMode breakdown:\n")
selected_bgs_combined %>%
  count(mode) %>%
  print()

cat("\n--- EMS (Upland) Combined ---\n")
if (nrow(selected_ems_combined) > 0) {
  selected_ems_combined %>%
    dplyr::select(variable, mode, var_group, window_hours, r) %>%
    mutate(
      window = paste0(round(window_hours/24, 1), "d"),
      r = round(r, 3)
    ) %>%
    dplyr::select(Process = var_group, Variable = variable, Mode = mode, Window = window, r) %>%
    arrange(desc(abs(r))) %>%
    print(n = 20)
  
  cat("\nMode breakdown:\n")
  selected_ems_combined %>%
    count(mode) %>%
    print()
}

# ============================================
# SAVE RESULTS
# ============================================

all_selected <- bind_rows(
  selected_bgs_combined,
  selected_ems_combined
)

write_csv(all_selected, file.path(OUTPUT_DIR, "selected_predictors_combined.csv"))

# ============================================
# METHODS TEXT
# ============================================

cat("\n\n========== METHODS TEXT ==========\n\n")

cat("We selected environmental predictors using a two-stage procedure applied to\n")
cat("both raw and diurnal/seasonal anomaly-filtered representations. First, we\n")
cat("identified predictors with significant associations with CH4 flux at any\n")
cat("integration window (1h-14d) using variable-level Benjamini-Hochberg FDR\n")
cat("correction (α = 0.05), separately for raw and anomaly representations.\n")
cat("We then pooled all significant predictors from both representations and\n")
cat("computed pairwise correlations among the windowed features. We clustered\n")
cat("variables with |r| ≥", CORR_THRESHOLD, "using connected components and selected\n")
cat("the predictor with the highest absolute correlation with CH4 flux from each\n")
cat("cluster. This approach allows the algorithm to choose the optimal representation\n")
cat("(raw vs anomaly) for each process, or retain both if they capture distinct signals.\n")

cat("\n\nThis procedure reduced the predictor set from:\n")
cat("  BGS:", nrow(bgs_combined), "→", nrow(selected_bgs_combined), 
    "(", sum(selected_bgs_combined$mode == "raw"), "raw,", 
    sum(selected_bgs_combined$mode == "anom"), "anomaly)\n")
if (nrow(ems_combined) > 0) {
  cat("  EMS:", nrow(ems_combined), "→", nrow(selected_ems_combined),
      "(", sum(selected_ems_combined$mode == "raw"), "raw,",
      sum(selected_ems_combined$mode == "anom"), "anomaly)\n")
}

message("\n\nOutputs saved to: ", OUTPUT_DIR)
message("Done!")
# ============================================================
# model_selection.R
#
# Extracts significant predictors from rolling correlation analysis
# and builds feature matrices for model fitting.
#
# Run AFTER: rolling_correlations.R
#
# Outputs:
#   - figures/model_selection/predictor_correlations_{site}.png
#   - figures/model_selection/selected_predictors.csv
#   - data/processed/model_data_{site}.csv (features + flux for each site)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(RcppRoll)
  library(pheatmap)
})

set.seed(42)

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  flux = "data/HF_2023-2025_tree_flux_corrected.csv",
  aligned = "data/processed/aligned_hourly_dataset.csv",
  best_windows = "figures/rolling_correlations/best_windows_by_variable.csv"
)

OUTPUT_DIR <- "figures/model_selection"
DATA_OUTPUT_DIR <- "data/processed"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Correlation threshold for redundancy clustering
CORR_THRESHOLD <- 0.7

# Variables that should be summed (not averaged) over windows
SUM_VARS <- c("P_mm", "THROUGHFALL_xHA")

# Core variables - always use RAW mode for these (theory-driven)
# These are fundamental drivers where raw values are most interpretable
CORE_VARIABLES <- c(
  "TS_Ha2",           # Soil temperature
  "s10t",             # Alternative soil temp (Fisher)
  "bvs_wtd_cm",       # Water table depth
  "SWC_Ha2",          # Soil water content (Harvard)
  "NEON_SWC_shallow", # Soil water content (NEON shallow)
  "NEON_SWC_mid"      # Soil water content (NEON mid)
)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

roll_mean <- function(x, n) RcppRoll::roll_mean(x, n = n, align = "right", fill = NA, na.rm = TRUE)
roll_sum  <- function(x, n) RcppRoll::roll_sum(x, n = n, align = "right", fill = NA, na.rm = TRUE)

#' Create anomaly version of a variable (remove DOY + hour effects)
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

#' Build feature matrix with rolling windows
build_features <- function(met, predictors_df) {
  
  features <- met %>% dplyr::select(datetime)
  
  for (i in seq_len(nrow(predictors_df))) {
    v <- predictors_df$variable[i]
    w <- predictors_df$window_hours[i]
    mode <- predictors_df$mode[i]
    
    # Source column
    src <- if (mode == "anom") paste0(v, "_anom") else v
    if (!src %in% names(met)) {
      message("  Warning: ", src, " not found, skipping")
      next
    }
    
    # Feature name
    feat_name <- paste0(v, "_", mode, "_", w, "h")
    
    if (v %in% SUM_VARS) {
      features[[feat_name]] <- roll_sum(met[[src]], w)
    } else {
      features[[feat_name]] <- roll_mean(met[[src]], w)
    }
  }
  
  features
}

#' Select non-redundant predictors using correlation clustering
#' For CORE_VARIABLES, always prefer raw mode regardless of |r|
select_nonredundant <- function(predictors_df, cor_matrix, threshold = 0.7) {
  
  vars_in_mat <- intersect(predictors_df$feat_name, rownames(cor_matrix))
  if (length(vars_in_mat) < 2) return(predictors_df)
  
  cmat <- cor_matrix[vars_in_mat, vars_in_mat]
  
  # Build adjacency (|r| >= threshold)
  adj <- abs(cmat) >= threshold
  diag(adj) <- FALSE
  
  # Connected components via BFS
  visited <- rep(FALSE, length(vars_in_mat))
  names(visited) <- vars_in_mat
  cluster_id <- rep(NA_integer_, length(vars_in_mat))
  names(cluster_id) <- vars_in_mat
  current_cluster <- 0
  
  for (start in vars_in_mat) {
    if (visited[start]) next
    current_cluster <- current_cluster + 1
    queue <- start
    
    while (length(queue) > 0) {
      node <- queue[1]
      queue <- queue[-1]
      if (visited[node]) next
      visited[node] <- TRUE
      cluster_id[node] <- current_cluster
      neighbors <- vars_in_mat[adj[node, ] & !visited[vars_in_mat]]
      queue <- c(queue, neighbors)
    }
  }
  
  # Map clusters back
  predictors_df <- predictors_df %>%
    left_join(tibble(feat_name = names(cluster_id), cluster = cluster_id), by = "feat_name")
  
  # Select best from each cluster
  # KEY CHANGE: For core variables, prefer raw mode even if anomaly has higher |r|
  selected <- predictors_df %>%
    filter(feat_name %in% vars_in_mat) %>%
    group_by(cluster) %>%
    arrange(
      # First: prefer raw for core variables
      !(variable %in% CORE_VARIABLES & mode == "raw"),
      # Second: highest |r|
      desc(abs(r))
    ) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(cluster_size = map_int(cluster, ~ sum(predictors_df$cluster == .x, na.rm = TRUE)))
  
  selected
}

# ============================================================
# LOAD DATA
# ============================================================

message("Loading data...")

# Aligned environmental data
aligned_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC"))

message("  Aligned data: ", nrow(aligned_data), " rows, ", 
        min(aligned_data$datetime), " to ", max(aligned_data$datetime))

# Flux data
flux_raw <- read_csv(PATHS$flux, show_col_types = FALSE)

flux_data <- flux_raw %>%
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  filter(!is.na(PLOT)) %>%
  mutate(
    datetime = as.POSIXct(datetime_posx, tz = "UTC"),
    datetime = round_date(datetime, "hour"),
    site = ifelse(PLOT == "BGS", "Wetland", "Upland"),
    Tree = as.factor(Tree),
    species = as.factor(SPECIES),
    CH4_flux_asinh = asinh(CH4_flux_nmolpm2ps)
  ) %>%
  dplyr::select(datetime, site, Tree, species, CH4_flux_nmolpm2ps, CH4_flux_asinh)

message("  Flux data: ", nrow(flux_data), " observations")
message("    Wetland: ", sum(flux_data$site == "Wetland"))
message("    Upland: ", sum(flux_data$site == "Upland"))

# Best windows from rolling correlations
best_windows <- read_csv(PATHS$best_windows, show_col_types = FALSE)

# Update site names if needed
if ("BGS" %in% best_windows$site) {
  best_windows <- best_windows %>%
    mutate(site = case_when(
      site == "BGS" ~ "Wetland",
      site == "EMS" ~ "Upland",
      TRUE ~ site
    ))
}

message("  Best windows: ", nrow(best_windows), " variable-window combinations")

# ============================================================
# EXTRACT SIGNIFICANT PREDICTORS
# ============================================================

message("\n========== EXTRACTING SIGNIFICANT PREDICTORS ==========\n")

# Separate raw and anomaly, add mode column
sig_raw <- best_windows %>%
  filter(analysis == "raw", significant) %>%
  mutate(mode = "raw")

sig_anom <- best_windows %>%
  filter(analysis == "anomaly", significant) %>%
  mutate(mode = "anom")

# Combine
sig_all <- bind_rows(sig_raw, sig_anom) %>%
  mutate(feat_name = paste0(variable, "_", mode, "_", window_hours, "h"))

for (s in c("Wetland", "Upland")) {
  cat("\n---", s, "---\n")
  
  raw_count <- sum(sig_all$site == s & sig_all$mode == "raw")
  anom_count <- sum(sig_all$site == s & sig_all$mode == "anom")
  
  cat("  Raw significant: ", raw_count, "\n")
  cat("  Anomaly significant: ", anom_count, "\n")
  cat("  Total: ", raw_count + anom_count, "\n")
  
  # Show top predictors
  cat("\n  Top 10 by |r|:\n")
  sig_all %>%
    filter(site == s) %>%
    arrange(desc(abs(r))) %>%
    head(10) %>%
    dplyr::select(variable, mode, window_hours, r, p_adj) %>%
    print()
}

# ============================================================
# BUILD FEATURES FOR EACH SITE
# ============================================================

message("\n========== BUILDING FEATURES ==========\n")

# Prepare met data with anomalies
met <- aligned_data %>% arrange(datetime)

# Get all unique anomaly variables needed
anom_vars <- sig_anom %>% pull(variable) %>% unique()

message("Creating anomaly versions for ", length(anom_vars), " variables...")

for (v in anom_vars) {
  if (v %in% names(met)) {
    met <- make_anomaly(met, v, time_col = "datetime")
  }
}

# Process each site
model_data_list <- list()

for (s in c("Wetland", "Upland")) {
  
  message("\n--- ", s, " ---")
  
  site_predictors <- sig_all %>% filter(site == s)
  
  if (nrow(site_predictors) == 0) {
    message("  No significant predictors, skipping")
    next
  }
  
  # Build features
  message("  Building ", nrow(site_predictors), " features...")
  features <- build_features(met, site_predictors)
  
  # Join with flux data
  site_flux <- flux_data %>% filter(site == s)
  
  model_data <- site_flux %>%
    left_join(features, by = "datetime") %>%
    drop_na()
  
  message("  Model data: ", nrow(model_data), " observations (", 
          n_distinct(model_data$Tree), " trees)")
  
  # Compute correlation matrix among predictors
  pred_cols <- site_predictors$feat_name
  pred_cols <- pred_cols[pred_cols %in% names(model_data)]
  
  X <- model_data %>% dplyr::select(all_of(pred_cols))
  X <- X %>% dplyr::select(where(~ sum(!is.na(.)) > 50))
  
  cor_matrix <- cor(X, use = "pairwise.complete.obs")
  
  # Plot correlation heatmap
  if (ncol(X) > 1) {
    message("  Creating correlation heatmap...")
    
    # Clean names for display
    display_names <- gsub("_raw_", " (raw ", colnames(cor_matrix))
    display_names <- gsub("_anom_", " (anom ", display_names)
    display_names <- gsub("h$", "h)", display_names)
    
    rownames(cor_matrix) <- display_names
    colnames(cor_matrix) <- display_names
    
    png(file.path(OUTPUT_DIR, paste0("predictor_correlations_", tolower(s), ".png")),
        width = 12, height = 11, units = "in", res = 300)
    
    pheatmap(
      cor_matrix,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      breaks = seq(-1, 1, length.out = 101),
      display_numbers = TRUE,
      number_format = "%.2f",
      number_color = "black",
      fontsize_number = 6,
      fontsize_row = 8,
      fontsize_col = 8,
      main = paste0(s, " - Predictor Correlations"),
      clustering_method = "complete"
    )
    
    dev.off()
    
    # Reset names for selection
    rownames(cor_matrix) <- pred_cols[pred_cols %in% names(X)]
    colnames(cor_matrix) <- pred_cols[pred_cols %in% names(X)]
  }
  
  # Select non-redundant predictors
  message("  Selecting non-redundant predictors (threshold = ", CORR_THRESHOLD, ")...")
  
  site_predictors <- site_predictors %>%
    filter(feat_name %in% colnames(cor_matrix))
  
  selected <- select_nonredundant(site_predictors, cor_matrix, threshold = CORR_THRESHOLD)
  
  message("  Selected ", nrow(selected), " from ", nrow(site_predictors), " predictors")
  
  cat("\n  Selected predictors:\n")
  selected %>%
    dplyr::select(variable, mode, window_hours, var_group, r, cluster_size) %>%
    arrange(var_group, desc(abs(r))) %>%
    print(n = 30)
  
  # Save model data
  model_data_list[[s]] <- list(
    data = model_data,
    all_predictors = site_predictors,
    selected_predictors = selected,
    cor_matrix = cor_matrix
  )
  
  # Save to CSV
  write_csv(model_data, file.path(DATA_OUTPUT_DIR, paste0("model_data_", tolower(s), ".csv")))
  message("  Saved: model_data_", tolower(s), ".csv")
}

# ============================================================
# SAVE SELECTED PREDICTORS
# ============================================================

message("\n========== SAVING RESULTS ==========\n")

# Combine selected predictors from both sites
all_selected <- bind_rows(
  model_data_list$Wetland$selected_predictors %>% mutate(site = "Wetland"),
  model_data_list$Upland$selected_predictors %>% mutate(site = "Upland")
)

write_csv(all_selected, file.path(OUTPUT_DIR, "selected_predictors.csv"))
message("Saved: selected_predictors.csv")

# Summary
cat("\n========== SUMMARY ==========\n\n")

for (s in c("Wetland", "Upland")) {
  if (!s %in% names(model_data_list)) next
  
  sel <- model_data_list[[s]]$selected_predictors
  
  cat(s, ":\n")
  cat("  Observations: ", nrow(model_data_list[[s]]$data), "\n")
  cat("  Trees: ", n_distinct(model_data_list[[s]]$data$Tree), "\n")
  cat("  Selected predictors: ", nrow(sel), "\n")
  cat("    Raw: ", sum(sel$mode == "raw"), "\n")
  cat("    Anomaly: ", sum(sel$mode == "anom"), "\n")
  cat("\n")
}

# Print predictor table for methods
cat("\n========== PREDICTOR TABLE FOR METHODS ==========\n\n")

all_selected %>%
  dplyr::select(site, variable, mode, window_hours, var_group, r) %>%
  mutate(
    window = paste0(round(window_hours/24, 1), "d"),
    r = round(r, 3)
  ) %>%
  dplyr::select(Site = site, Variable = variable, Mode = mode, Window = window, 
                Process = var_group, r) %>%
  arrange(Site, Process, desc(abs(r))) %>%
  print(n = 50)

message("\n\nOutputs saved to: ", OUTPUT_DIR, " and ", DATA_OUTPUT_DIR)
message("Done!")
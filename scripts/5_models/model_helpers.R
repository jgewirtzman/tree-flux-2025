# ============================================================
# model_helpers.R
# Functions shared by the three driver-model scripts (sourced, not run on its own):
#   5_models/01_wetland_model.R, 02_upland_model_A.R, 03_upland_model_B.R
#
#   roll_mean(), roll_sum(), make_anomaly()      rolling windows, anomalies (scripts/helpers/windows.R)
#   build_feature_matrix_raw()                   drivers at their best windows (raw)
#   build_combined_feature_matrix()              drivers at their best windows (raw + anomaly)
#   plot_corr_heatmap(), find_redundant_pairs()  redundancy among the screened drivers
#   select_predictors_combined()                 one driver per cluster of |r| >= threshold
#   forward_select()                             greedy forward selection after the core drivers
#
# The calling script defines SUM_VARS (drivers summed, not averaged, over a window).
# ============================================================
source("scripts/helpers/windows.R")

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

# Greedy forward selection: starting from the core drivers, add at each step the candidate
# that most increases the fixed-effect variance share var(fixed prediction) / var_total of
#   CH4_flux ~ <selected> + species + (1|Tree) + (1|date)   (ML fit),
# for up to max_predictors steps. Deterministic: the candidate order matters only for exact ties.
# Returns one row per step: the predictor added and the variance share after adding it.
forward_select <- function(core_names, candidate_names, data, var_total, max_predictors) {
  selected <- core_names
  remaining <- candidate_names
  steps <- list()
  for (step in 1:max_predictors) {
    best_r2 <- -Inf
    best_pred <- NULL
    for (pred in remaining) {
      formula_str <- paste("CH4_flux ~", paste(c(selected, pred), collapse = " + "),
                           "+ species + (1|Tree) + (1|date)")
      tryCatch({
        fit <- lmer(as.formula(formula_str), data = data, REML = FALSE)
        r2 <- var(predict(fit, re.form = NA)) / var_total
        if (r2 > best_r2) {
          best_r2 <- r2
          best_pred <- pred
        }
      }, error = function(e) NULL)
    }
    if (is.null(best_pred)) break
    steps[[length(steps) + 1]] <- tibble(step = step, predictor = best_pred, r2 = best_r2)
    selected <- c(selected, best_pred)
    remaining <- setdiff(remaining, best_pred)
  }
  bind_rows(steps)
}

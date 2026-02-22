# ============================================================
# 09_tree_predictor_plots.R
#
# Simple descriptive plots: each predictor vs CH4 flux
# faceted by location × species
#
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  stem_flux   = "/Users/jongewirtzman/Google Drive/Research/tree-flux-2025/data/HF_2023-2025_tree_flux.csv",
  tomography  = "/Users/jongewirtzman/Downloads/tomography_results_compiled.csv",
  hummock     = "/Users/jongewirtzman/Downloads/hummock_hollow.csv"
)

OUTPUT_DIR <- "figures/tree_predictor_plots"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD DATA
# ============================================================

message("\nLoading data...")

# --- Stem flux data ---
stem_flux_raw <- read_csv(PATHS$stem_flux, show_col_types = FALSE)

stem_flux <- stem_flux_raw %>%
  mutate(CH4_flux_nmolpm2ps = ifelse(year < 2025, CH4_flux_nmolpm2ps * 1000, CH4_flux_nmolpm2ps)) %>%
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  filter(CH4_flux_nmolpm2ps >= -1, !is.na(PLOT)) %>%
  mutate(
    tree = Tree,
    date = as.Date(datetime_posx),
    location = ifelse(PLOT == "BGS", "wetland", "upland"),
    species = SPECIES,
    CH4_flux = CH4_flux_nmolpm2ps,
    species_full = case_when(
      species == "bg"  ~ "Nyssa sylvatica",
      species == "hem" ~ "Tsuga canadensis",
      species == "rm"  ~ "Acer rubrum",
      species == "ro"  ~ "Quercus rubra",
      TRUE ~ species
    )
  )

message("  Stem flux: ", nrow(stem_flux), " measurements from ", 
        n_distinct(stem_flux$tree), " trees")

# --- Tomography data ---
tomography <- read_csv(PATHS$tomography, show_col_types = FALSE) %>%
  dplyr::select(-filename) %>%
  rename(
    ert_radial_gradient = ert_radialgradiant,
    pct_solid = sot_solid,
    pct_damaged = sot_damaged
  )

message("  Tomography: ", nrow(tomography), " trees")

# --- Hummock/hollow data ---
hummock <- read_csv(PATHS$hummock, show_col_types = FALSE) %>%
  rename(tree = Tag, microsite = Classification) %>%
  mutate(
    microsite_simple = case_when(
      microsite == "hummock" ~ "hummock",
      microsite == "hollow" ~ "hollow",
      TRUE ~ "intermediate"
    )
  ) %>%
  dplyr::select(tree, microsite, microsite_simple)

message("  Hummock/hollow: ", nrow(hummock), " trees")

# ============================================================
# COMPUTE TREE-LEVEL MEAN FLUX
# ============================================================

message("\nComputing tree-level summaries...")

tree_flux <- stem_flux %>%
  group_by(tree, location, species, species_full) %>%
  summarize(
    n_obs = n(),
    CH4_mean = mean(CH4_flux, na.rm = TRUE),
    .groups = "drop"
  )

message("  Tree-level summaries: ", nrow(tree_flux), " trees")

# ============================================================
# MERGE DATASETS
# ============================================================

tree_data <- tree_flux %>%
  left_join(tomography, by = "tree") %>%
  left_join(hummock, by = "tree")

message("  Trees with tomography: ", sum(!is.na(tree_data$ert_mean)), " / ", nrow(tree_data))
message("  Trees with microsite: ", sum(!is.na(tree_data$microsite)), " / ", nrow(tree_data))

# Save merged data
write_csv(tree_data, file.path(OUTPUT_DIR, "tree_level_data.csv"))

# ============================================================
# PREDICTOR LABELS
# ============================================================

PRED_LABELS <- c(
  ert_mean = "ERT Mean",
  ert_median = "ERT Median", 
  ert_sd = "ERT Std Dev",
  ert_cv = "ERT CV",
  ert_gini = "ERT Gini",
  ert_entropy = "ERT Entropy",
  ert_cma = "ERT CMA",
  ert_radial_gradient = "ERT Radial Gradient",
  pct_solid = "% Solid Wood",
  pct_damaged = "% Damaged Wood",
  microsite_simple = "Microsite"
)

# ============================================================
# CONTINUOUS PREDICTORS - SCATTER PLOTS
# ============================================================

message("\nGenerating scatter plots for continuous predictors...")

cont_vars <- c("ert_mean", "ert_median", "ert_sd", "ert_cv", "ert_gini", 
               "ert_entropy", "ert_cma", "ert_radial_gradient", 
               "pct_solid", "pct_damaged")

# Store correlation results
cor_results <- tibble()

for (var in cont_vars) {
  if (!var %in% names(tree_data)) next
  
  plot_data <- tree_data %>%
    filter(!is.na(.data[[var]]), !is.na(CH4_mean))
  
  if (nrow(plot_data) < 3) next
  
  # Create facet label
  plot_data <- plot_data %>%
    mutate(facet_label = paste(location, "-", species_full))
  
  # Calculate correlation and p-value for each facet
  facet_stats <- plot_data %>%
    group_by(facet_label, location, species_full) %>%
    summarize(
      n = n(),
      r = if (n() >= 3) cor(.data[[var]], CH4_mean, use = "complete.obs") else NA_real_,
      p = if (n() >= 3) {
        tryCatch(cor.test(.data[[var]], CH4_mean)$p.value, error = function(e) NA_real_)
      } else NA_real_,
      .groups = "drop"
    ) %>%
    mutate(
      sig = p < 0.05 & !is.na(p),
      label = paste0("r = ", round(r, 2), "\np = ", formatC(p, format = "f", digits = 3)),
      line_color = ifelse(sig, "red", "steelblue")
    )
  
  # Store results
  cor_results <- bind_rows(cor_results, 
                           facet_stats %>% mutate(variable = var, var_label = PRED_LABELS[var]))
  
  # Join stats back to plot data
  plot_data <- plot_data %>%
    left_join(facet_stats %>% dplyr::select(facet_label, sig, line_color, label), by = "facet_label")
  
  # Get x position for labels (90% of range)
  label_pos <- plot_data %>%
    group_by(facet_label) %>%
    summarize(
      x_pos = min(.data[[var]], na.rm = TRUE) + 0.05 * (max(.data[[var]], na.rm = TRUE) - min(.data[[var]], na.rm = TRUE)),
      y_pos = max(CH4_mean, na.rm = TRUE) - 0.05 * (max(CH4_mean, na.rm = TRUE) - min(CH4_mean, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    left_join(facet_stats %>% dplyr::select(facet_label, label), by = "facet_label")
  
  p <- ggplot(plot_data, aes(x = .data[[var]], y = CH4_mean)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, aes(color = sig), linewidth = 0.8, show.legend = FALSE) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "steelblue")) +
    geom_text_repel(aes(label = tree), size = 2.5, max.overlaps = 15, color = "gray40") +
    geom_text(data = label_pos, aes(x = x_pos, y = y_pos, label = label), 
              hjust = 0, vjust = 1, size = 3, color = "gray20") +
    facet_wrap(~ facet_label, scales = "free") +
    labs(
      title = paste("CH4 Flux vs", PRED_LABELS[var]),
      subtitle = "Red line = p < 0.05",
      x = PRED_LABELS[var],
      y = expression("Mean CH"[4]*" Flux (nmol m"^-2*" s"^-1*")")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "gray95", color = NA)
    )
  
  # Determine dimensions based on number of facets
  n_facets <- n_distinct(plot_data$facet_label)
  ncol <- min(3, n_facets)
  nrow <- ceiling(n_facets / ncol)
  
  ggsave(file.path(OUTPUT_DIR, paste0("scatter_", var, ".png")), 
         p, width = 4 * ncol, height = 3.5 * nrow, dpi = 300)
  
  message("  Saved: scatter_", var, ".png")
}

# Save correlation results
write_csv(cor_results, file.path(OUTPUT_DIR, "correlation_by_group.csv"))
message("  Saved: correlation_by_group.csv")

# ============================================================
# MICROSITE - BOX PLOTS
# ============================================================

message("\nGenerating microsite box plots...")

microsite_data <- tree_data %>%
  filter(!is.na(microsite_simple), !is.na(CH4_mean))

if (nrow(microsite_data) >= 3) {
  
  microsite_data <- microsite_data %>%
    mutate(facet_label = paste(location, "-", species_full))
  
  p_micro <- ggplot(microsite_data, aes(x = microsite_simple, y = CH4_mean, fill = microsite_simple)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 2.5, alpha = 0.7) +
    geom_text_repel(aes(label = tree), size = 2.5, max.overlaps = 10, color = "gray40") +
    facet_wrap(~ facet_label, scales = "free_y") +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
      title = "CH4 Flux by Microsite",
      x = "Microsite",
      y = expression("Mean CH"[4]*" Flux (nmol m"^-2*" s"^-1*")")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_rect(fill = "gray95", color = NA),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  n_facets <- n_distinct(microsite_data$facet_label)
  ncol <- min(3, n_facets)
  nrow <- ceiling(n_facets / ncol)
  
  ggsave(file.path(OUTPUT_DIR, "boxplot_microsite.png"), 
         p_micro, width = 4 * ncol, height = 3.5 * nrow, dpi = 300)
  
  message("  Saved: boxplot_microsite.png")
}

# ============================================================
# COMBINED PANEL - ALL PREDICTORS
# ============================================================

message("\nGenerating combined panel plot...")

# Pivot all continuous predictors to long format
long_data <- tree_data %>%
  filter(!is.na(CH4_mean)) %>%
  dplyr::select(tree, location, species_full, CH4_mean, all_of(cont_vars)) %>%
  pivot_longer(cols = all_of(cont_vars), names_to = "predictor", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(pred_label = PRED_LABELS[predictor])

# Calculate overall correlation per predictor
overall_stats <- long_data %>%
  group_by(predictor, pred_label) %>%
  summarize(
    n = n(),
    r = cor(value, CH4_mean, use = "complete.obs"),
    p = tryCatch(cor.test(value, CH4_mean)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    sig = p < 0.05 & !is.na(p),
    label = paste0("r = ", round(r, 2), ", p = ", formatC(p, format = "f", digits = 3))
  )

# Join significance to long_data
long_data <- long_data %>%
  left_join(overall_stats %>% dplyr::select(predictor, sig), by = "predictor")

p_all <- ggplot(long_data, aes(x = value, y = CH4_mean)) +
  geom_point(aes(color = species_full, shape = location), size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, aes(linetype = sig), 
              color = "gray30", linewidth = 0.6, show.legend = FALSE) +
  geom_smooth(data = long_data %>% filter(sig), method = "lm", se = TRUE,
              color = "red", linewidth = 0.8, show.legend = FALSE) +
  geom_smooth(data = long_data %>% filter(!sig), method = "lm", se = TRUE,
              color = "gray50", linewidth = 0.6, show.legend = FALSE) +
  facet_wrap(~ pred_label, scales = "free_x", ncol = 4) +
  scale_color_brewer(palette = "Dark2", name = "Species") +
  scale_shape_manual(values = c("wetland" = 16, "upland" = 17), name = "Location") +
  labs(
    title = "CH4 Flux vs All Tomography Predictors",
    subtitle = "Red line = p < 0.05; Color = species, shape = location",
    x = "Predictor Value",
    y = expression("Mean CH"[4]*" Flux (nmol m"^-2*" s"^-1*")")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    strip.text = element_text(face = "bold", size = 9),
    strip.background = element_rect(fill = "gray95", color = NA),
    legend.position = "bottom"
  )

ggsave(file.path(OUTPUT_DIR, "panel_all_predictors.png"), 
       p_all, width = 14, height = 10, dpi = 300)

message("  Saved: panel_all_predictors.png")

# Print overall stats
message("\nOverall correlations (ignoring location × species):")
print(overall_stats %>% arrange(p))

# ============================================================
# SUMMARY TABLE
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("SUMMARY")
message(paste(rep("=", 60), collapse = ""))

# Sample sizes by location × species
sample_sizes <- tree_data %>%
  group_by(location, species_full) %>%
  summarize(
    n_trees = n(),
    n_with_tomo = sum(!is.na(ert_mean)),
    n_with_microsite = sum(!is.na(microsite)),
    .groups = "drop"
  )

message("\nSample sizes by location × species:")
print(sample_sizes)

message("\nOutputs saved to: ", OUTPUT_DIR)
message("Done!")
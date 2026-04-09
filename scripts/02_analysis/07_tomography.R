# ============================================================
# tomography.R
#
# ERT/Sonic tomography visualization with CH4 flux bars,
# plus PCA-based decay classification and comprehensive
# ERT metric vs flux correlation analysis.
#
# ERT sign convention:
#   - Resistivity (Ohm-m): HIGHER = DRIER wood
#   - ert_mean, ert_median: raw resistivity, higher = drier
#   - ert_cv, ert_gini: variability metrics, higher = more heterogeneous
#   - ert_cma: Central Moisture Accumulation, positive = wet center (decay)
#   - ert_entropy: higher = more uniform moisture distribution
#   - PC1 (species-normalized): higher = WETTER / more anomalous
#     (mean/median load negatively, cv/gini load positively)
#
# Required external data:
#   - tomography_results_compiled.csv
#   - hummock_hollow.csv (for Nyssa identification)
#   - ERT image directory
#   - Sonic image directory
#
# Outputs:
#   - tomography_specialists.png/pdf (Nyssa + Oak)
#   - tomography_generalists.png/pdf (Maple + Hemlock by location)
#   - pca_biplot_flux.png/pdf
#   - ert_metric_vs_flux_*.png/pdf (per-metric scatter grids)
#   - ert_correlation_summary.csv
#   - ert_metric_correlations.csv
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(magick)
  library(patchwork)
  library(grid)
  library(scales)
})

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  tomography = "data/input/tomography_results_compiled.csv",
  hummock = "data/input/hummock_hollow.csv",
  flux = "data/processed/flux_with_quality_flags.csv",
  ert_images = "data/input/tomography/ERTs_Absolute",
  sonic_images = "data/input/tomography/CH4_PITs"
)

OUTPUT_DIR <- "outputs/figures/tomography"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Project-wide color palettes (consistent across all figures)
species_colors <- c("N. sylvatica" = "#2A7F7A", "Q. rubra" = "#6E8B3D",
                     "A. rubrum" = "#A7DAD1", "T. canadensis" = "#C9D6A4")
site_colors <- c("Wetland" = "#2A7F7A", "Upland" = "#6E8B3D")
spp_shapes <- c("N. sylvatica" = 15, "A. rubrum" = 17,
                "T. canadensis" = 16, "Q. rubra" = 18)

# ============================================================
# LOAD AND PREPARE DATA
# ============================================================

message("Loading data...")

if (!file.exists(PATHS$tomography)) stop("Tomography data not found: ", PATHS$tomography)
if (!file.exists(PATHS$flux)) stop("Flux data not found: ", PATHS$flux)

tomography <- read_csv(PATHS$tomography, show_col_types = FALSE)

if (file.exists(PATHS$hummock)) {
  hummock <- read_csv(PATHS$hummock, show_col_types = FALSE) %>% rename(tree = Tag)
} else {
  message("  Warning: hummock_hollow.csv not found")
  hummock <- tibble(tree = integer(), Species = character())
}

flux_data <- read_csv(PATHS$flux, show_col_types = FALSE) %>%
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  filter(!is.na(PLOT)) %>%
  mutate(location = factor(ifelse(PLOT == "BGS", "Wetland", "Upland"),
                           levels = c("Wetland", "Upland")))

tree_flux_means <- flux_data %>%
  group_by(Tree, SPECIES, PLOT, location) %>%
  summarize(CH4_mean = mean(CH4_flux_nmolpm2ps, na.rm = TRUE), .groups = "drop")

message("  Flux data: ", nrow(tree_flux_means), " trees")

# Image file lists
ert_images <- if (dir.exists(PATHS$ert_images)) {
  list.files(PATHS$ert_images, pattern = "\\.jpg$", full.names = TRUE)
} else { character(0) }
sonic_images <- if (dir.exists(PATHS$sonic_images)) {
  list.files(PATHS$sonic_images, pattern = "\\.jpg$", full.names = TRUE)
} else { character(0) }
message("  ERT images: ", length(ert_images), "  Sonic images: ", length(sonic_images))

# ============================================================
# MERGE TOMOGRAPHY + FLUX FOR FULL ANALYSIS DATASET
# ============================================================

species_labels_full <- c(bg = "N. sylvatica", rm = "A. rubrum",
                         hem = "T. canadensis", ro = "Q. rubra")

tomo_flux <- tomography %>%
  left_join(tree_flux_means, by = c("tree" = "Tree")) %>%
  filter(!is.na(CH4_mean), !is.na(SPECIES)) %>%
  mutate(
    species_full = species_labels_full[SPECIES],
    # Derived metrics (absolute values for PCA)
    abs_cma = abs(ert_cma),
    abs_radgrad = abs(ert_radialgradiant)
  )

message("  Merged tomo+flux: ", nrow(tomo_flux), " trees")

# ============================================================
# PCA ON 8 ERT METRICS (species-normalized)
# ============================================================

message("\nComputing species-normalized PCA on ERT metrics...")

# ERT metrics for PCA
# NOTE: These are all in resistivity space (Ohm-m).
#   ert_mean, ert_median: higher = drier
#   ert_sd: absolute variability
#   ert_cv, ert_gini: relative variability, higher = more heterogeneous
#   ert_entropy: higher = more uniform distribution
#   abs_cma, abs_radgrad: spatial pattern magnitude
pca_metrics <- c("ert_mean", "ert_median", "ert_sd", "ert_cv",
                 "ert_gini", "ert_entropy", "abs_cma", "abs_radgrad")

# Species z-score normalization
spp_stats <- tomo_flux %>%
  group_by(species_full) %>%
  summarise(across(all_of(pca_metrics),
                   list(mu = ~ mean(., na.rm = TRUE),
                        sigma = ~ sd(., na.rm = TRUE))),
            .groups = "drop")

pca_normed <- tomo_flux %>%
  select(tree, species_full, all_of(pca_metrics)) %>%
  left_join(spp_stats, by = "species_full")

for (m in pca_metrics) {
  mu_col <- paste0(m, "_mu")
  sd_col <- paste0(m, "_sigma")
  pca_normed[[m]] <- (pca_normed[[m]] - pca_normed[[mu_col]]) / pca_normed[[sd_col]]
}

pca_input <- pca_normed %>% select(all_of(pca_metrics)) %>% as.matrix()
pca_input[is.nan(pca_input)] <- 0

pca_fit <- prcomp(pca_input, center = FALSE, scale. = FALSE)

ve <- summary(pca_fit)$importance[2, 1:3] * 100
message("  PC1: ", round(ve[1], 1), "% | PC2: ", round(ve[2], 1), "% | PC3: ", round(ve[3], 1), "%")

tomo_flux$pc1 <- pca_fit$x[, 1]
tomo_flux$pc2 <- pca_fit$x[, 2]

# Ensure PC1 direction: high PC1 = low resistivity = WETTER / more anomalous
# (mean/median should load negatively on final PC1)
if (pca_fit$rotation["ert_mean", 1] > 0) {
  tomo_flux$pc1 <- -tomo_flux$pc1
  pca_fit$rotation[, 1] <- -pca_fit$rotation[, 1]
  message("  Flipped PC1: high PC1 = wetter/more anomalous (lower resistivity)")
}

message("\n  PC1 loadings (high PC1 = wetter/more anomalous):")
pc1_load <- round(pca_fit$rotation[, 1], 3)
for (nm in names(pc1_load)) message("    ", nm, ": ", pc1_load[nm])

# --- Decay classification (CV-based) ---
# SoT threshold: 1% structural damage
# ERT threshold: species-normalized CV z-score = 0
#   Higher CV = more moisture heterogeneity = more decay
sot_threshold <- 1

cv_spp_stats <- tomo_flux %>%
  group_by(species_full) %>%
  summarise(cv_mu = mean(ert_cv, na.rm = TRUE),
            cv_sd = sd(ert_cv, na.rm = TRUE), .groups = "drop")

tomo_flux <- tomo_flux %>%
  left_join(cv_spp_stats, by = "species_full") %>%
  mutate(
    cv_z = (ert_cv - cv_mu) / cv_sd,
    decay_phase = case_when(
      sot_damaged <= sot_threshold & cv_z <= 0 ~ "I: Sound",
      sot_damaged <= sot_threshold & cv_z >  0 ~ "II: Incipient",
      sot_damaged >  sot_threshold & cv_z >  0 ~ "III: Active",
      sot_damaged >  sot_threshold & cv_z <= 0 ~ "IV: Cavity"
    ),
    decay_phase_short = str_extract(decay_phase, "^[IV]+")
  ) %>%
  select(-cv_mu, -cv_sd)

message("\n  Decay phase distribution:")
print(table(tomo_flux$decay_phase))

# Merge PC1 and decay phase back into tomography for image panels
tomography <- tomography %>%
  left_join(tomo_flux %>% select(tree, pc1, decay_phase, decay_phase_short),
            by = "tree")

# ============================================================
# SPECIES/LOCATION GROUPS FOR IMAGE PANELS
# ============================================================

nyssa_trees <- hummock %>% filter(Species == "bg") %>% pull(tree)
oak_trees <- flux_data %>% filter(SPECIES == "ro") %>% distinct(Tree) %>% pull(Tree)
wetland_rm_trees <- flux_data %>% filter(location == "Wetland", SPECIES == "rm") %>% distinct(Tree) %>% pull(Tree)
wetland_hem_trees <- flux_data %>% filter(location == "Wetland", SPECIES == "hem") %>% distinct(Tree) %>% pull(Tree)
upland_rm_trees <- flux_data %>% filter(location == "Upland", SPECIES == "rm") %>% distinct(Tree) %>% pull(Tree)
upland_hem_trees <- flux_data %>% filter(location == "Upland", SPECIES == "hem") %>% distinct(Tree) %>% pull(Tree)

# ============================================================
# HELPER FUNCTIONS (image panels)
# ============================================================

prepare_species_data <- function(tree_ids, species_name, sort_metric = "ert_cv") {
  df <- tomography %>%
    filter(tree %in% tree_ids) %>%
    arrange(.data[[sort_metric]]) %>%
    mutate(order = row_number()) %>%
    left_join(tree_flux_means %>% select(Tree, CH4_mean), by = c("tree" = "Tree")) %>%
    mutate(
      ert_path = map_chr(tree, function(t) {
        matches <- ert_images[grepl(paste0("^", t, "_"), basename(ert_images))]
        if (length(matches) > 0) matches[1] else NA_character_
      }),
      sonic_path = map_chr(tree, function(t) {
        matches <- sonic_images[grepl(paste0("^", t, "_"), basename(sonic_images))]
        if (length(matches) > 0) matches[1] else NA_character_
      }),
      species = species_name
    ) %>%
    filter(!is.na(ert_path), !is.na(sonic_path), !is.na(CH4_mean))

  if (nrow(df) > 0) {
    df <- df %>% mutate(
      flux_transformed = asinh(CH4_mean),
      flux_norm = (flux_transformed - min(flux_transformed)) /
        (max(flux_transformed) - min(flux_transformed) + 1e-10)
    )
  }
  df
}

read_image_as_raster <- function(path, target_size = 200) {
  img <- image_read(path)
  img <- image_resize(img, paste0(target_size, "x", target_size, "!"))
  as.raster(img)
}

create_species_panel <- function(data, species_label, bad_sonic_indices = c(),
                                  annotation_pos = "topright",
                                  metric = "ert_cv", metric_label = "ERT CV",
                                  metric_x_label = "ERT CV") {
  n <- nrow(data)
  if (n == 0) return(NULL)

  message("  Reading images for ", species_label, "...")
  ert_rasters <- map(data$ert_path, read_image_as_raster)
  sonic_rasters <- map(1:n, function(i) {
    if (i %in% bad_sonic_indices) NULL else read_image_as_raster(data$sonic_path[i])
  })

  p_images <- ggplot() +
    coord_fixed(ratio = 1, xlim = c(0, n), ylim = c(0, 3.6), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(5, 10, 5, 40),
          plot.title = element_text(hjust = 0.5, size = 14, face = "italic")) +
    ggtitle(species_label)

  # ERT metric axis guide
  arrow_y <- 3.4
  p_images <- p_images +
    annotate("segment", x = 0, xend = n, y = arrow_y, yend = arrow_y,
             arrow = arrow(ends = "both", length = unit(0.08, "inches")), linewidth = 0.4) +
    annotate("text", x = 0, y = arrow_y + 0.15, label = "lower", size = 3, hjust = 0) +
    annotate("text", x = n / 2, y = arrow_y + 0.15, label = metric_label, size = 3.5,
             hjust = 0.5, fontface = "bold") +
    annotate("text", x = n, y = arrow_y + 0.15, label = "higher", size = 3, hjust = 1)

  # Row labels
  p_images <- p_images +
    annotate("text", x = -0.4, y = 2.6, label = "ERT", size = 5, fontface = "bold", hjust = 0.5) +
    annotate("text", x = -0.4, y = 1.55, label = "SoT", size = 5, fontface = "bold", hjust = 0.5) +
    annotate("text", x = -0.4, y = 0.7, label = "Phase", size = 5, fontface = "bold", hjust = 0.5) +
    annotate("text", x = -0.4, y = 0.2, label = "CH[4]", size = 5, fontface = "bold",
             hjust = 0.5, parse = TRUE)

  # ERT images (top row)
  for (i in 1:n) {
    p_images <- p_images + annotation_raster(ert_rasters[[i]],
      xmin = i - 1, xmax = i, ymin = 2.1, ymax = 3.1)
  }

  # Sonic images
  for (i in 1:n) {
    if (!is.null(sonic_rasters[[i]])) {
      p_images <- p_images + annotation_raster(sonic_rasters[[i]],
        xmin = i - 1, xmax = i, ymin = 1.05, ymax = 2.05)
    }
  }

  # Decay phase labels (PCA-based quadrant classification)
  phase_colors <- c("I" = "#4E79A7", "II" = "#E5C460",
                     "III" = "#D4873F", "IV" = "#C4524E")
  for (i in 1:n) {
    if (!(i %in% bad_sonic_indices) && !is.na(data$decay_phase_short[i])) {
      p_images <- p_images + annotate(
        "text", x = i - 0.5, y = 0.7,
        label = data$decay_phase_short[i],
        size = 5, fontface = "bold",
        color = phase_colors[data$decay_phase_short[i]]
      )
    }
  }

  # Flux bars
  cor_colors <- c("#4393C3", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#D6604D")
  flux_bar_data <- data %>%
    mutate(xmin = order - 1, xmax = order, ymin = 0, ymax = 0.4,
           fill_color = cor_colors[pmin(pmax(round(flux_norm * 4) + 1, 1), 5)],
           text_color = ifelse(flux_norm > 0.25 & flux_norm < 0.75, "black", "white"),
           flux_label = sprintf("%.2f", CH4_mean))

  for (i in 1:n) {
    row <- flux_bar_data[i, ]
    p_images <- p_images +
      annotate("rect", xmin = row$xmin, xmax = row$xmax,
               ymin = row$ymin, ymax = row$ymax,
               fill = row$fill_color, color = NA) +
      annotate("text", x = (row$xmin + row$xmax) / 2, y = (row$ymin + row$ymax) / 2,
               label = row$flux_label, color = row$text_color, size = 5)
  }

  # Scatter: ERT metric vs CH4 flux
  cor_test <- cor.test(data[[metric]], data$CH4_mean)
  r_val <- cor_test$estimate; p_val <- cor_test$p.value

  if (annotation_pos == "bottomright") {
    ann_x <- Inf; ann_y <- -Inf; ann_hjust <- 1.1; ann_vjust <- -0.5
  } else {
    ann_x <- Inf; ann_y <- Inf; ann_hjust <- 1.1; ann_vjust <- 1.5
  }

  p_scatter <- ggplot(data, aes(x = .data[[metric]], y = CH4_mean)) +
    geom_point(size = 2.5, shape = 16)
  if (p_val < 0.05) {
    p_scatter <- p_scatter + geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8)
  }
  p_scatter <- p_scatter +
    labs(x = metric_x_label, y = expression(CH[4]~flux)) +
    annotate("text", x = ann_x, y = ann_y,
             label = sprintf("r = %.2f\np = %.3f", r_val, p_val),
             hjust = ann_hjust, vjust = ann_vjust, size = 3.5) +
    theme_classic(base_size = 11) +
    theme(axis.line = element_line(linewidth = 0.5),
          axis.ticks = element_line(linewidth = 0.5),
          aspect.ratio = 1,
          plot.margin = margin(15, 15, 15, 15))

  p_images + p_scatter + plot_layout(widths = c(3, 1.2))
}

# ============================================================
# BUILD IMAGE PANEL FIGURES
# ============================================================

message("\nPreparing species data...")

nyssa_data <- prepare_species_data(nyssa_trees, "Nyssa sylvatica")
oak_data <- prepare_species_data(oak_trees, "Quercus rubra")
wetland_rm_data <- prepare_species_data(wetland_rm_trees, "Acer rubrum (Wetland)")
wetland_hem_data <- prepare_species_data(wetland_hem_trees, "Tsuga canadensis (Wetland)")
upland_rm_data <- prepare_species_data(upland_rm_trees, "Acer rubrum (Upland)")
upland_hem_data <- prepare_species_data(upland_hem_trees, "Tsuga canadensis (Upland)")

message("  Nyssa: ", nrow(nyssa_data), " | Oak: ", nrow(oak_data),
        " | Wet rm: ", nrow(wetland_rm_data), " | Wet hem: ", nrow(wetland_hem_data),
        " | Up rm: ", nrow(upland_rm_data), " | Up hem: ", nrow(upland_hem_data))

# --- Site-level scatter plot helper ---
# Creates a scatter of ERT metric vs CH4 flux for one site,
# with per-species regression lines + pooled line + pooled r,p
site_scatter <- function(site_name, metric = "ert_cv", x_label = "ERT CV",
                         annotation_pos = "topright") {
  site_data <- tomo_flux %>% filter(location == site_name)

  # Pooled stats
  ct <- cor.test(site_data[[metric]], site_data$CH4_mean)
  pooled_p <- ct$p.value
  pooled_label <- sprintf("r = %.2f, p = %.3f", ct$estimate, pooled_p)

  # Per-species significance (only show regression line if p < 0.05)
  sig_species <- site_data %>%
    group_by(species_full) %>%
    filter(sum(!is.na(.data[[metric]]) & !is.na(CH4_mean)) >= 3) %>%
    summarise(p = cor.test(.data[[metric]], CH4_mean)$p.value, .groups = "drop") %>%
    filter(p < 0.05) %>%
    pull(species_full)

  sig_data <- site_data %>% filter(species_full %in% sig_species)

  p <- ggplot(site_data, aes(x = .data[[metric]], y = CH4_mean))

  # Per-species regression lines (dashed, only significant)
  if (nrow(sig_data) > 0) {
    p <- p + geom_smooth(data = sig_data, aes(color = species_full),
                         method = "lm", se = FALSE, linetype = "dashed",
                         linewidth = 0.6, alpha = 0.6)
  }

  # Pooled regression line (solid, only if significant)
  if (pooled_p < 0.05) {
    p <- p + geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.9,
                         alpha = 0.15)
  }

  # Points
  p <- p +
    geom_point(aes(color = species_full, shape = species_full), size = 3, alpha = 0.8) +
    scale_color_manual(values = species_colors, name = "Species") +
    scale_shape_manual(values = spp_shapes, name = "Species") +
    annotate("text",
             x = if (annotation_pos == "topleft") -Inf else Inf,
             y = Inf, label = pooled_label,
             hjust = if (annotation_pos == "topleft") -0.1 else 1.1,
             vjust = 1.5, size = 4, color = "grey30") +
    labs(
      x = x_label,
      y = expression(Mean~CH[4]~flux~(nmol~m^{-2}~s^{-1})),
      title = site_name
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
      aspect.ratio = 1,
      legend.position = "bottom"
    )
  p
}

# --- Specialists (main figure: CV) ---
message("\nBuilding specialist panels...")
if (nrow(nyssa_data) > 0 && nrow(oak_data) > 0) {
  p_nyssa <- create_species_panel(nyssa_data, "Nyssa sylvatica",
                                  bad_sonic_indices = c(7), annotation_pos = "topright",
                                  metric = "ert_cv", metric_label = "ERT CV",
                                  metric_x_label = "ERT CV")
  p_oak <- create_species_panel(oak_data, "Quercus rubra",
                                bad_sonic_indices = c(6), annotation_pos = "bottomright",
                                metric = "ert_cv", metric_label = "ERT CV",
                                metric_x_label = "ERT CV")

  # Site-level ERT CV vs CH4 flux scatter panels
  p_wet_cv <- site_scatter("Wetland", metric = "ert_cv",
                           x_label = "ERT CV")
  p_up_cv  <- site_scatter("Upland",  metric = "ert_cv",
                           x_label = "ERT CV", annotation_pos = "topleft")

  p_specialists <- p_nyssa / plot_spacer() / p_oak / plot_spacer() /
    (p_wet_cv + plot_spacer() + p_up_cv + plot_layout(widths = c(1, 0.1, 1))) +
    plot_layout(heights = c(1, 0.05, 1, 0.05, 1.5))
  print(p_specialists)
  ggsave(file.path(OUTPUT_DIR, "tomography_specialists.png"), p_specialists,
         width = 14, height = 14, dpi = 300, bg = "white")
  ggsave(file.path(OUTPUT_DIR, "tomography_specialists.pdf"), p_specialists,
         width = 14, height = 14, bg = "white")
  message("  Saved: tomography_specialists.png/pdf")

  # --- SI figure: ERT mean throughout (re-sort images + scatter by mean) ---
  nyssa_data_mean <- prepare_species_data(nyssa_trees, "Nyssa sylvatica", sort_metric = "ert_mean")
  oak_data_mean   <- prepare_species_data(oak_trees, "Quercus rubra", sort_metric = "ert_mean")

  p_nyssa_mean <- create_species_panel(nyssa_data_mean, "Nyssa sylvatica",
                                       bad_sonic_indices = c(), annotation_pos = "topright",
                                       metric = "ert_mean", metric_label = "ERT Mean",
                                       metric_x_label = "Mean resistivity (Ohm-m)")
  p_oak_mean <- create_species_panel(oak_data_mean, "Quercus rubra",
                                     bad_sonic_indices = c(), annotation_pos = "bottomright",
                                     metric = "ert_mean", metric_label = "ERT Mean",
                                     metric_x_label = "Mean resistivity (Ohm-m)")

  p_wet_mean <- site_scatter("Wetland", metric = "ert_mean",
                             x_label = "Mean resistivity (Ohm-m)")
  p_up_mean  <- site_scatter("Upland",  metric = "ert_mean",
                             x_label = "Mean resistivity (Ohm-m)")

  p_specialists_si <- p_nyssa_mean / plot_spacer() / p_oak_mean / plot_spacer() /
    (p_wet_mean + plot_spacer() + p_up_mean + plot_layout(widths = c(1, 0.1, 1))) +
    plot_layout(heights = c(1, 0.05, 1, 0.05, 1.5))
  print(p_specialists_si)
  ggsave(file.path(OUTPUT_DIR, "tomography_specialists_SI_mean.png"), p_specialists_si,
         width = 14, height = 14, dpi = 300, bg = "white")
  ggsave(file.path(OUTPUT_DIR, "tomography_specialists_SI_mean.pdf"), p_specialists_si,
         width = 14, height = 14, bg = "white")
  message("  Saved: tomography_specialists_SI_mean.png/pdf")
}

# --- Generalists ---
message("\nBuilding generalist panels...")
panel_list <- list(); height_list <- c()
for (d in list(
  list(data = wetland_rm_data, label = "Acer rubrum (Wetland)"),
  list(data = wetland_hem_data, label = "Tsuga canadensis (Wetland)"),
  list(data = upland_rm_data, label = "Acer rubrum (Upland)"),
  list(data = upland_hem_data, label = "Tsuga canadensis (Upland)")
)) {
  if (nrow(d$data) > 0) {
    if (length(panel_list) > 0) {
      panel_list <- c(panel_list, list(plot_spacer())); height_list <- c(height_list, 0.05)
    }
    p <- create_species_panel(d$data, d$label, annotation_pos = "topright")
    panel_list <- c(panel_list, list(p)); height_list <- c(height_list, 1)
  }
}
if (length(panel_list) > 0) {
  p_generalists <- wrap_plots(panel_list, ncol = 1) + plot_layout(heights = height_list)
  n_sp <- sum(sapply(list(wetland_rm_data, wetland_hem_data, upland_rm_data, upland_hem_data), nrow) > 0)
  print(p_generalists)
  ggsave(file.path(OUTPUT_DIR, "tomography_generalists.png"), p_generalists,
         width = 14, height = 4.5 * n_sp, dpi = 300, bg = "white")
  ggsave(file.path(OUTPUT_DIR, "tomography_generalists.pdf"), p_generalists,
         width = 14, height = 4.5 * n_sp, bg = "white")
  message("  Saved: tomography_generalists.png/pdf")
}

# ============================================================
# COMPREHENSIVE ERT METRIC vs CH4 FLUX ANALYSIS
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("ERT METRIC vs CH4 FLUX CORRELATIONS")
message(paste(rep("=", 60), collapse = ""))

# Metrics to test, with sign convention labels
# NOTE: All raw ERT values are resistivity (Ohm-m): HIGHER = DRIER
metric_info <- tribble(
  ~metric,              ~label,                         ~direction,
  "ert_mean",           "Mean resistivity",             "higher = drier",
  "ert_median",         "Median resistivity",           "higher = drier",
  "ert_sd",             "SD resistivity",               "higher = more variable",
  "ert_cv",             "CV",                           "higher = more heterogeneous",
  "ert_gini",           "Gini coefficient",             "higher = more unequal moisture",
  "ert_entropy",        "Shannon entropy",              "higher = more uniform",
  "ert_cma",            "CMA (signed)",                 "positive = wet center",
  "ert_radialgradiant", "Radial gradient (signed)",     "positive = drier edges",
  "pc1",                "PC1 (species-normalized)",     "higher = wetter/more anomalous"
)

# --- Compute correlations for all groupings ---
# Groups: Overall, by location, by species, by species x location
compute_cors <- function(data, metric_name) {
  results <- tibble()

  # Helper: safe cor.test
  safe_cor <- function(x, y, grp_label) {
    complete <- complete.cases(x, y)
    x <- x[complete]; y <- y[complete]
    if (length(x) < 4) return(tibble(
      metric = metric_name, group = grp_label, n = length(x),
      r = NA_real_, r2 = NA_real_, p = NA_real_
    ))
    ct <- cor.test(x, y)
    tibble(metric = metric_name, group = grp_label, n = length(x),
           r = round(ct$estimate, 3), r2 = round(ct$estimate^2, 3),
           p = round(ct$p.value, 4))
  }

  # Overall
  results <- bind_rows(results, safe_cor(data[[metric_name]], data$CH4_mean, "Overall"))

  # By location
  for (loc in unique(data$location)) {
    sub <- data %>% filter(location == loc)
    results <- bind_rows(results, safe_cor(sub[[metric_name]], sub$CH4_mean, as.character(loc)))
  }

  # By species
  for (sp in unique(data$species_full)) {
    sub <- data %>% filter(species_full == sp)
    results <- bind_rows(results, safe_cor(sub[[metric_name]], sub$CH4_mean, sp))
  }

  # By species x location
  for (sp in unique(data$species_full)) {
    for (loc in unique(data$location)) {
      sub <- data %>% filter(species_full == sp, location == loc)
      if (nrow(sub) >= 3) {
        results <- bind_rows(results,
          safe_cor(sub[[metric_name]], sub$CH4_mean, paste0(sp, " (", loc, ")")))
      }
    }
  }

  results
}

all_cors <- tibble()
for (metric in metric_info$metric) {
  if (metric %in% names(tomo_flux)) {
    all_cors <- bind_rows(all_cors, compute_cors(tomo_flux, metric))
  }
}

all_cors <- all_cors %>%
  mutate(sig = case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    p < 0.1 ~ ".",
    TRUE ~ ""
  ))

# Print summary: Overall, Wetland, Upland for each metric
message("\n--- CORRELATIONS: Overall / Wetland / Upland ---")
summary_cors <- all_cors %>%
  filter(group %in% c("Overall", "Wetland", "Upland")) %>%
  left_join(metric_info %>% select(metric, direction), by = "metric") %>%
  select(metric, direction, group, n, r, r2, p, sig) %>%
  pivot_wider(names_from = group, values_from = c(n, r, r2, p, sig),
              names_glue = "{group}_{.value}")
print(summary_cors, width = 140)

# Print species-level
message("\n--- CORRELATIONS: By Species ---")
species_cors <- all_cors %>%
  filter(group %in% unique(tomo_flux$species_full)) %>%
  select(metric, group, n, r, r2, p, sig)
print(species_cors, width = 120)

# Print species x location
message("\n--- CORRELATIONS: By Species x Location ---")
sxl_cors <- all_cors %>%
  filter(str_detect(group, "\\(")) %>%
  select(metric, group, n, r, r2, p, sig)
print(sxl_cors, width = 120)

# Save full table
write.csv(all_cors, file.path(OUTPUT_DIR, "ert_correlation_summary.csv"), row.names = FALSE)
message("\nSaved: ert_correlation_summary.csv")

# ============================================================
# SCATTER PLOT GRIDS: Each metric vs CH4 flux
# ============================================================

message("\nGenerating per-metric scatter plots...")

for (i in seq_len(nrow(metric_info))) {
  m <- metric_info$metric[i]
  m_label <- metric_info$label[i]
  m_dir <- metric_info$direction[i]

  if (!m %in% names(tomo_flux)) next

  # Compute per-facet stats
  facet_stats <- tomo_flux %>%
    group_by(species_full, location) %>%
    filter(sum(!is.na(.data[[m]]) & !is.na(CH4_mean)) >= 3) %>%
    summarise(
      n = sum(!is.na(.data[[m]]) & !is.na(CH4_mean)),
      r = cor(.data[[m]], CH4_mean, use = "complete.obs"),
      r2 = r^2,
      p = cor.test(.data[[m]], CH4_mean)$p.value,
      label = sprintf("r=%.2f  R\u00b2=%.2f\np=%.3f  n=%d", r, r2, p, n),
      # Position for annotation
      x_pos = max(.data[[m]], na.rm = TRUE),
      y_pos = max(CH4_mean, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot(tomo_flux, aes(x = .data[[m]], y = CH4_mean)) +
    geom_point(aes(color = species_full), size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6,
                linetype = "dashed", alpha = 0.15) +
    geom_text(data = facet_stats,
              aes(x = x_pos, y = y_pos, label = label),
              hjust = 1, vjust = 1, size = 3, color = "grey30", inherit.aes = FALSE) +
    facet_grid(location ~ species_full, scales = "free") +
    labs(
      x = m_label,
      y = expression(Mean~CH[4]~flux~(nmol~m^{-2}~s^{-1})),
      title = paste0(m_label, " vs CH4 flux"),
      color = "Species"
    ) +
    theme_classic(base_size = 11) +
    theme(
      strip.text = element_text(face = "bold", size = 10),
      strip.background = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 13)
    )

  fname <- paste0("ert_metric_vs_flux_", gsub("[^a-z0-9]", "_", m), ".png")
  ggsave(file.path(OUTPUT_DIR, fname), p,
         width = 12, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(OUTPUT_DIR, gsub(".png$", ".pdf", fname)), p,
         width = 12, height = 6, bg = "white")
  message("  Saved: ", fname)

  # --- Pooled by location ---
  pooled_stats <- tomo_flux %>%
    group_by(location) %>%
    filter(sum(!is.na(.data[[m]]) & !is.na(CH4_mean)) >= 3) %>%
    summarise(
      n = sum(!is.na(.data[[m]]) & !is.na(CH4_mean)),
      r = cor(.data[[m]], CH4_mean, use = "complete.obs"),
      r2 = r^2,
      p = cor.test(.data[[m]], CH4_mean)$p.value,
      label = sprintf("r=%.2f  R\u00b2=%.2f\np=%.3f  n=%d", r, r2, p, n),
      x_pos = max(.data[[m]], na.rm = TRUE),
      y_pos = max(CH4_mean, na.rm = TRUE),
      .groups = "drop"
    )

  p_pooled <- ggplot(tomo_flux, aes(x = .data[[m]], y = CH4_mean)) +
    geom_point(aes(color = species_full, shape = species_full), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7,
                linetype = "dashed", alpha = 0.15) +
    geom_text(data = pooled_stats,
              aes(x = x_pos, y = y_pos, label = label),
              hjust = 1, vjust = 1, size = 3.5, color = "grey30", inherit.aes = FALSE) +
    facet_wrap(~ location, scales = "free") +
    labs(
      x = m_label,
      y = expression(Mean~CH[4]~flux~(nmol~m^{-2}~s^{-1})),
      title = paste0(m_label, " vs CH4 flux — pooled by location"),
      color = "Species", shape = "Species"
    ) +
    theme_classic(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 12),
      strip.background = element_blank(),
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 13)
    )

  fname_pooled <- paste0("ert_pooled_vs_flux_", gsub("[^a-z0-9]", "_", m), ".png")
  ggsave(file.path(OUTPUT_DIR, fname_pooled), p_pooled,
         width = 10, height = 5, dpi = 300, bg = "white")
  ggsave(file.path(OUTPUT_DIR, gsub(".png$", ".pdf", fname_pooled)), p_pooled,
         width = 10, height = 5, bg = "white")
  message("  Saved: ", fname_pooled)

  # --- Overall (all trees, colored by location) ---
  overall_ct <- cor.test(tomo_flux[[m]], tomo_flux$CH4_mean)
  overall_label <- sprintf("r=%.2f  R\u00b2=%.2f\np=%.4f  n=%d",
                           overall_ct$estimate, overall_ct$estimate^2,
                           overall_ct$p.value, nrow(tomo_flux))

  p_overall <- ggplot(tomo_flux, aes(x = .data[[m]], y = CH4_mean)) +
    geom_point(aes(color = location, shape = species_full), size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7,
                linetype = "dashed", alpha = 0.15) +
    scale_color_manual(values = site_colors) +
    annotate("text", x = Inf, y = Inf, label = overall_label,
             hjust = 1.1, vjust = 1.3, size = 4, color = "grey30") +
    labs(
      x = m_label,
      y = expression(Mean~CH[4]~flux~(nmol~m^{-2}~s^{-1})),
      title = paste0(m_label, " vs CH4 flux — all trees"),
      color = "Site", shape = "Species"
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold", size = 13))

  fname_overall <- paste0("ert_overall_vs_flux_", gsub("[^a-z0-9]", "_", m), ".png")
  ggsave(file.path(OUTPUT_DIR, fname_overall), p_overall,
         width = 7, height = 6, dpi = 300, bg = "white")
  ggsave(file.path(OUTPUT_DIR, gsub(".png$", ".pdf", fname_overall)), p_overall,
         width = 7, height = 6, bg = "white")
  message("  Saved: ", fname_overall)
}

# ============================================================
# PCA BIPLOT with flux coloring
# ============================================================

message("\nGenerating PCA biplot...")

loadings <- as.data.frame(pca_fit$rotation[, 1:2])
loadings$metric <- rownames(loadings)
names(loadings)[1:2] <- c("PC1", "PC2")

arrow_scale <- min(max(abs(tomo_flux$pc1)), max(abs(tomo_flux$pc2))) * 0.8 /
  max(sqrt(loadings$PC1^2 + loadings$PC2^2))
loadings$PC1 <- loadings$PC1 * arrow_scale
loadings$PC2 <- loadings$PC2 * arrow_scale

metric_labels_pca <- c(ert_mean = "Mean", ert_median = "Median", ert_sd = "SD",
                       ert_cv = "CV", ert_gini = "Gini", ert_entropy = "Entropy",
                       abs_cma = "|CMA|", abs_radgrad = "|RadGrad|")
loadings$label <- metric_labels_pca[loadings$metric]

p_biplot <- ggplot(tomo_flux, aes(x = pc1, y = pc2)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(aes(shape = species_full, color = location), size = 3, alpha = 0.7) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "#B22222", linewidth = 0.7, inherit.aes = FALSE) +
  geom_text(data = loadings,
            aes(x = PC1 * 1.15, y = PC2 * 1.15, label = label),
            color = "#B22222", size = 3.5, fontface = "bold", inherit.aes = FALSE) +
  scale_shape_manual(name = "Species", values = spp_shapes) +
  scale_color_manual(name = "Site", values = site_colors) +
  labs(x = paste0("PC1 (", round(ve[1], 1), "% var) — high = wetter/more anomalous"),
       y = paste0("PC2 (", round(ve[2], 1), "% var)")) +
  theme_classic(base_size = 13) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        legend.position = "right")

ggsave(file.path(OUTPUT_DIR, "pca_biplot_flux.png"), p_biplot,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "pca_biplot_flux.pdf"), p_biplot,
       width = 10, height = 8, bg = "white")
message("Saved: pca_biplot_flux.png/pdf")

# ============================================================
# SUMMARY STATISTICS
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("TOMOGRAPHY ANALYSIS SUMMARY")
message(paste(rep("=", 60), collapse = ""))

message("\nTotal trees: ", nrow(tomo_flux))
message("  Wetland: ", sum(tomo_flux$location == "Wetland"),
        "  Upland: ", sum(tomo_flux$location == "Upland"))

message("\n--- DECAY PHASE DISTRIBUTION ---")
print(table(tomo_flux$decay_phase))

message("\nBy location:")
print(table(tomo_flux$location, tomo_flux$decay_phase))

message("\n--- CH4 FLUX BY DECAY PHASE ---")
flux_by_phase <- tomo_flux %>%
  group_by(decay_phase) %>%
  summarise(n = n(), mean = round(mean(CH4_mean), 3),
            se = round(sd(CH4_mean) / sqrt(n()), 3),
            median = round(median(CH4_mean), 3), .groups = "drop")
print(flux_by_phase)

message("\n--- ERT METRICS BY LOCATION (medians) ---")
ert_loc <- tomo_flux %>%
  group_by(location) %>%
  summarise(across(c(ert_mean, ert_median, ert_cv, ert_gini, ert_cma, pc1),
                   ~ round(median(., na.rm = TRUE), 3)),
            .groups = "drop")
print(ert_loc)

message("\n", paste(rep("=", 60), collapse = ""))
message("Done!")
# ============================================================
# tomography.R
#
# ERT/Sonic tomography visualization with CH4 flux bars.
# Creates composite figures showing tree internal structure
# alongside flux measurements.
#
# Run AFTER: timeseries.R (which creates the corrected flux file)
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
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(magick)
  library(patchwork)
  library(grid)
})

# ============================================================
# CONFIGURATION - UPDATE PATHS AS NEEDED
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

# ============================================================
# LOAD AND PREPARE DATA
# ============================================================

message("Loading data...")

# Check if files exist
if (!file.exists(PATHS$tomography)) {
  
  stop("Tomography data not found: ", PATHS$tomography)
}
if (!file.exists(PATHS$flux)) {
  stop("Flux data not found: ", PATHS$flux)
}

# Tomography results
tomography <- read_csv(PATHS$tomography, show_col_types = FALSE)

# Hummock data (for Nyssa identification)
if (file.exists(PATHS$hummock)) {
  hummock <- read_csv(PATHS$hummock, show_col_types = FALSE) %>%
    rename(tree = Tag)
} else {
  message("  Warning: hummock_hollow.csv not found, Nyssa identification may be incomplete")
  hummock <- tibble(tree = integer(), Species = character())
}

# Flux data
flux_data <- read_csv(PATHS$flux, show_col_types = FALSE) %>%
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  filter(!is.na(PLOT)) %>%
  mutate(
    location = factor(
      ifelse(PLOT == "BGS", "Wetland", "Upland"),
      levels = c("Wetland", "Upland")
    )
  )

tree_flux_means <- flux_data %>%
  group_by(Tree, SPECIES, location) %>%
  summarize(CH4_mean = mean(CH4_flux_nmolpm2ps, na.rm = TRUE), .groups = "drop")

message("  Flux data: ", nrow(tree_flux_means), " trees")

# Image file lists
if (dir.exists(PATHS$ert_images)) {
  ert_images <- list.files(PATHS$ert_images, pattern = "\\.jpg$", full.names = TRUE)
  message("  ERT images: ", length(ert_images))
} else {
  ert_images <- character(0)
  message("  Warning: ERT image directory not found")
}

if (dir.exists(PATHS$sonic_images)) {
  sonic_images <- list.files(PATHS$sonic_images, pattern = "\\.jpg$", full.names = TRUE)
  message("  Sonic images: ", length(sonic_images))
} else {
  sonic_images <- character(0)
  message("  Warning: Sonic image directory not found")
}

# ============================================================
# IDENTIFY SPECIES/LOCATION GROUPS
# ============================================================

# Nyssa from hummock data
nyssa_trees <- hummock %>% filter(Species == "bg") %>% pull(tree)

# Oak from flux data
oak_trees <- flux_data %>% filter(SPECIES == "ro") %>% distinct(Tree) %>% pull(Tree)

# Generalists by location
wetland_rm_trees <- flux_data %>% filter(location == "Wetland", SPECIES == "rm") %>% distinct(Tree) %>% pull(Tree)
wetland_hem_trees <- flux_data %>% filter(location == "Wetland", SPECIES == "hem") %>% distinct(Tree) %>% pull(Tree)
upland_rm_trees <- flux_data %>% filter(location == "Upland", SPECIES == "rm") %>% distinct(Tree) %>% pull(Tree)
upland_hem_trees <- flux_data %>% filter(location == "Upland", SPECIES == "hem") %>% distinct(Tree) %>% pull(Tree)

# ============================================================
# HELPER FUNCTIONS
# ============================================================

#' Prepare species data for plotting
prepare_species_data <- function(tree_ids, species_name) {
  
  df <- tomography %>%
    filter(tree %in% tree_ids) %>%
    arrange(ert_cv) %>%
    mutate(order = row_number()) %>%
    left_join(tree_flux_means %>% dplyr::select(Tree, CH4_mean), by = c("tree" = "Tree")) %>%
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
  
  # Normalize flux within species using asinh
  if (nrow(df) > 0) {
    df <- df %>%
      mutate(
        flux_transformed = asinh(CH4_mean),
        flux_norm = (flux_transformed - min(flux_transformed)) / 
          (max(flux_transformed) - min(flux_transformed) + 1e-10)
      )
  }
  
  df
}

#' Read image as raster for ggplot
read_image_as_raster <- function(path, target_size = 200) {
  img <- image_read(path)
  img <- image_resize(img, paste0(target_size, "x", target_size, "!"))
  as.raster(img)
}

#' Create a single species panel
create_species_panel <- function(data, species_label, bad_sonic_indices = c(), annotation_pos = "topright") {
  
  n <- nrow(data)
  if (n == 0) return(NULL)
  
  # Read all images
  message("  Reading images for ", species_label, "...")
  ert_rasters <- map(data$ert_path, read_image_as_raster)
  
  # For sonic, replace bad indices with NULL
  sonic_rasters <- map(1:n, function(i) {
    if (i %in% bad_sonic_indices) NULL else read_image_as_raster(data$sonic_path[i])
  })
  
  # Create base plot
  p_images <- ggplot() +
    xlim(-0.8, n) +
    ylim(0, 3.6) +
    coord_fixed(ratio = 1) +
    theme_void() +
    theme(
      plot.margin = margin(5, 5, 5, 5),
      plot.title = element_text(hjust = 0.5, size = 14, face = "italic")
    ) +
    ggtitle(species_label)
  
  # Add ERT CV axis guide at top
  arrow_y <- 3.4
  p_images <- p_images + 
    annotate("segment", x = 0, xend = n, y = arrow_y, yend = arrow_y,
             arrow = arrow(ends = "both", length = unit(0.08, "inches")), 
             linewidth = 0.4) +
    annotate("text", x = 0, y = arrow_y + 0.15, label = "lower", size = 3, hjust = 0) +
    annotate("text", x = n / 2, y = arrow_y + 0.15, label = "ERT CV", size = 3.5, hjust = 0.5, fontface = "bold") +
    annotate("text", x = n, y = arrow_y + 0.15, label = "higher", size = 3, hjust = 1)
  
  # Add row labels on left
  p_images <- p_images + 
    annotate("text", x = -0.4, y = 2.6, label = "ERT", size = 5, fontface = "bold", hjust = 0.5) +
    annotate("text", x = -0.4, y = 1.55, label = "SoT", size = 5, fontface = "bold", hjust = 0.5) +
    annotate("text", x = -0.4, y = 0.7, label = "Class", size = 5, fontface = "bold", hjust = 0.5) +
    annotate("text", x = -0.4, y = 0.2, label = "CH[4]", size = 5, fontface = "bold", hjust = 0.5, parse = TRUE)
  
  # Add ERT images (top row)
  for (i in 1:n) {
    p_images <- p_images + annotation_raster(
      ert_rasters[[i]],
      xmin = i - 1, xmax = i,
      ymin = 2.1, ymax = 3.1
    )
  }
  
  # Add Sonic images (skip bad indices)
  for (i in 1:n) {
    if (!is.null(sonic_rasters[[i]])) {
      p_images <- p_images + annotation_raster(
        sonic_rasters[[i]],
        xmin = i - 1, xmax = i,
        ymin = 1.05, ymax = 2.05
      )
    }
  }
  
  # Add decay category row
  decay_data <- data %>%
    mutate(
      low_sonic = sot_damaged > 1,
      low_resistance = ert_cv > 0.5,
      decay_cat = case_when(
        !low_sonic & !low_resistance ~ "I",
        !low_sonic & low_resistance ~ "II",
        low_sonic & low_resistance ~ "III",
        low_sonic & !low_resistance ~ "IV"
      )
    )
  
  for (i in 1:n) {
    if (!(i %in% bad_sonic_indices)) {
      p_images <- p_images + annotate(
        "text", x = i - 0.5, y = 0.7,
        label = decay_data$decay_cat[i],
        size = 5, fontface = "bold"
      )
    }
  }
  
  # Add flux bars (bottom row)
  cor_colors <- c("#4393C3", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#D6604D")
  
  flux_bar_data <- data %>%
    mutate(
      xmin = order - 1,
      xmax = order,
      ymin = 0,
      ymax = 0.4,
      fill_color = cor_colors[pmin(pmax(round(flux_norm * 4) + 1, 1), 5)],
      text_color = ifelse(flux_norm > 0.25 & flux_norm < 0.75, "black", "white"),
      flux_label = sprintf("%.2f", CH4_mean)
    )
  
  for (i in 1:n) {
    row <- flux_bar_data[i, ]
    p_images <- p_images + 
      annotate("rect", xmin = row$xmin, xmax = row$xmax,
               ymin = row$ymin, ymax = row$ymax,
               fill = row$fill_color, color = NA) +
      annotate("text", x = (row$xmin + row$xmax) / 2, y = (row$ymin + row$ymax) / 2,
               label = row$flux_label, color = row$text_color, size = 5)
  }
  
  # Create scatterplot: ERT CV vs CH4 flux
  cor_test <- cor.test(data$ert_cv, data$CH4_mean)
  r_val <- cor_test$estimate
  p_val <- cor_test$p.value
  
  if (annotation_pos == "bottomright") {
    ann_x <- Inf; ann_y <- -Inf; ann_hjust <- 1.1; ann_vjust <- -0.5
  } else {
    ann_x <- Inf; ann_y <- Inf; ann_hjust <- 1.1; ann_vjust <- 1.5
  }
  
  p_scatter <- ggplot(data, aes(x = ert_cv, y = CH4_mean)) +
    geom_point(size = 2.5, shape = 16)
  
  if (p_val < 0.05) {
    p_scatter <- p_scatter + geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8)
  }
  
  p_scatter <- p_scatter +
    labs(x = "ERT CV", y = expression(CH[4]~flux)) +
    annotate("text", x = ann_x, y = ann_y, 
             label = sprintf("r = %.2f\np = %.3f", r_val, p_val),
             hjust = ann_hjust, vjust = ann_vjust, size = 3.5) +
    theme_classic(base_size = 11) +
    theme(
      axis.line = element_line(linewidth = 0.5),
      axis.ticks = element_line(linewidth = 0.5),
      plot.margin = margin(15, 15, 15, 15)
    )
  
  # Combine images and scatterplot
  p_combined <- p_images + p_scatter + plot_layout(widths = c(4, 1))
  
  p_combined
}

# ============================================================
# BUILD FIGURES
# ============================================================

message("\nPreparing species data...")

nyssa_data <- prepare_species_data(nyssa_trees, "Nyssa sylvatica")
oak_data <- prepare_species_data(oak_trees, "Quercus rubra")
wetland_rm_data <- prepare_species_data(wetland_rm_trees, "Acer rubrum (Wetland)")
wetland_hem_data <- prepare_species_data(wetland_hem_trees, "Tsuga canadensis (Wetland)")
upland_rm_data <- prepare_species_data(upland_rm_trees, "Acer rubrum (Upland)")
upland_hem_data <- prepare_species_data(upland_hem_trees, "Tsuga canadensis (Upland)")

message("  Nyssa: ", nrow(nyssa_data), " trees")
message("  Oak: ", nrow(oak_data), " trees")
message("  Wetland red maple: ", nrow(wetland_rm_data), " trees")
message("  Wetland hemlock: ", nrow(wetland_hem_data), " trees")
message("  Upland red maple: ", nrow(upland_rm_data), " trees")
message("  Upland hemlock: ", nrow(upland_hem_data), " trees")

# ============================================================
# FIGURE 1: Specialists (Nyssa + Oak)
# ============================================================

message("\nBuilding specialist panels...")

if (nrow(nyssa_data) > 0 && nrow(oak_data) > 0) {
  p_nyssa <- create_species_panel(nyssa_data, "Nyssa sylvatica", 
                                  bad_sonic_indices = c(7), annotation_pos = "topright")
  p_oak <- create_species_panel(oak_data, "Quercus rubra", 
                                bad_sonic_indices = c(6), annotation_pos = "bottomright")
  
  p_specialists <- p_nyssa / plot_spacer() / p_oak +
    plot_layout(heights = c(1, 0.05, 1))
  
  print(p_specialists)
  
  ggsave(file.path(OUTPUT_DIR, "tomography_specialists.png"), p_specialists, 
         width = 14, height = 10, dpi = 300, bg = "white")
  ggsave(file.path(OUTPUT_DIR, "tomography_specialists.pdf"), p_specialists, 
         width = 14, height = 10, bg = "white")
  
  message("  Saved: tomography_specialists.png/pdf")
} else {
  message("  Skipping specialists figure (insufficient data)")
}

# ============================================================
# FIGURE 2: Generalists (Maple + Hemlock by location)
# ============================================================

message("\nBuilding generalist panels...")

panel_list <- list()
height_list <- c()

# Wetland maples
if (nrow(wetland_rm_data) > 0) {
  p <- create_species_panel(wetland_rm_data, "Acer rubrum (Wetland)", annotation_pos = "topright")
  panel_list <- c(panel_list, list(p))
  height_list <- c(height_list, 1)
}

# Wetland hemlocks
if (nrow(wetland_hem_data) > 0) {
  if (length(panel_list) > 0) {
    panel_list <- c(panel_list, list(plot_spacer()))
    height_list <- c(height_list, 0.05)
  }
  p <- create_species_panel(wetland_hem_data, "Tsuga canadensis (Wetland)", annotation_pos = "topright")
  panel_list <- c(panel_list, list(p))
  height_list <- c(height_list, 1)
}

# Upland maples
if (nrow(upland_rm_data) > 0) {
  if (length(panel_list) > 0) {
    panel_list <- c(panel_list, list(plot_spacer()))
    height_list <- c(height_list, 0.05)
  }
  p <- create_species_panel(upland_rm_data, "Acer rubrum (Upland)", annotation_pos = "topright")
  panel_list <- c(panel_list, list(p))
  height_list <- c(height_list, 1)
}

# Upland hemlocks
if (nrow(upland_hem_data) > 0) {
  if (length(panel_list) > 0) {
    panel_list <- c(panel_list, list(plot_spacer()))
    height_list <- c(height_list, 0.05)
  }
  p <- create_species_panel(upland_hem_data, "Tsuga canadensis (Upland)", annotation_pos = "topright")
  panel_list <- c(panel_list, list(p))
  height_list <- c(height_list, 1)
}

if (length(panel_list) > 0) {
  p_generalists <- wrap_plots(panel_list, ncol = 1) + plot_layout(heights = height_list)
  
  n_species <- sum(nrow(wetland_rm_data) > 0, nrow(wetland_hem_data) > 0, 
                   nrow(upland_rm_data) > 0, nrow(upland_hem_data) > 0)
  fig_height <- 4.5 * n_species
  
  print(p_generalists)
  
  ggsave(file.path(OUTPUT_DIR, "tomography_generalists.png"), p_generalists, 
         width = 14, height = fig_height, dpi = 300, bg = "white")
  ggsave(file.path(OUTPUT_DIR, "tomography_generalists.pdf"), p_generalists, 
         width = 14, height = fig_height, bg = "white")
  
  message("  Saved: tomography_generalists.png/pdf")
} else {
  message("  Skipping generalists figure (no data)")
}

message("\nOutputs saved to: ", OUTPUT_DIR)
message("Done!")


# ============================================================
# TOMOGRAPHY RESULTS SUMMARY
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("TOMOGRAPHY ANALYSIS SUMMARY")
message(paste(rep("=", 60), collapse = ""))

# Combine all species data
all_tomo_data <- bind_rows(
  nyssa_data %>% mutate(species_group = "N. sylvatica", location = "Wetland"),
  oak_data %>% mutate(species_group = "Q. rubra", location = "Upland"),
  wetland_rm_data %>% mutate(species_group = "A. rubrum", location = "Wetland"),
  wetland_hem_data %>% mutate(species_group = "T. canadensis", location = "Wetland"),
  upland_rm_data %>% mutate(species_group = "A. rubrum", location = "Upland"),
  upland_hem_data %>% mutate(species_group = "T. canadensis", location = "Upland")
) %>%
  mutate(
    low_sonic = sot_damaged > 1,
    low_resistance = ert_cv > 0.5,
    decay_cat = case_when(
      !low_sonic & !low_resistance ~ "I",
      !low_sonic & low_resistance ~ "II",
      low_sonic & low_resistance ~ "III",
      low_sonic & !low_resistance ~ "IV"
    )
  )

message("\n--- OVERALL SUMMARY ---")
message("Total trees with tomography: ", nrow(all_tomo_data))
message("  Wetland: ", sum(all_tomo_data$location == "Wetland"))
message("  Upland: ", sum(all_tomo_data$location == "Upland"))

# Decay category distribution
message("\n--- DECAY CATEGORY DISTRIBUTION ---")
decay_summary <- all_tomo_data %>%
  group_by(decay_cat) %>%
  summarise(
    n = n(),
    pct = round(100 * n() / nrow(all_tomo_data), 1),
    .groups = "drop"
  ) %>%
  arrange(decay_cat)
print(decay_summary)

message("\nBy location:")
decay_by_loc <- all_tomo_data %>%
  group_by(location, decay_cat) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = decay_cat, values_from = n, values_fill = 0)
print(decay_by_loc)

# ERT CV summary
message("\n--- ERT CV SUMMARY ---")
ert_summary <- all_tomo_data %>%
  group_by(location) %>%
  summarise(
    n = n(),
    mean = round(mean(ert_cv), 3),
    sd = round(sd(ert_cv), 3),
    min = round(min(ert_cv), 3),
    max = round(max(ert_cv), 3),
    .groups = "drop"
  )
print(ert_summary)

# Sonic damage summary  
message("\n--- SONIC DAMAGE (% > 1) SUMMARY ---")
sonic_summary <- all_tomo_data %>%
  group_by(location) %>%
  summarise(
    n = n(),
    mean = round(mean(sot_damaged), 1),
    sd = round(sd(sot_damaged), 1),
    min = round(min(sot_damaged), 1),
    max = round(max(sot_damaged), 1),
    pct_damaged = round(100 * mean(low_sonic), 1),
    .groups = "drop"
  )
print(sonic_summary)

# CH4 flux by decay category
message("\n--- CH4 FLUX BY DECAY CATEGORY ---")
flux_by_decay <- all_tomo_data %>%
  group_by(decay_cat) %>%
  summarise(
    n = n(),
    mean_flux = round(mean(CH4_mean), 2),
    se_flux = round(sd(CH4_mean) / sqrt(n()), 2),
    median_flux = round(median(CH4_mean), 2),
    .groups = "drop"
  ) %>%
  arrange(decay_cat)
print(flux_by_decay)

# Correlations by species/location
message("\n--- CORRELATIONS: ERT CV vs CH4 FLUX ---")

correlation_results <- all_tomo_data %>%
  group_by(species_group, location) %>%
  summarise(
    n = n(),
    r = round(cor(ert_cv, CH4_mean), 3),
    .groups = "drop"
  ) %>%
  arrange(location, species_group)

# Add p-values
for (i in 1:nrow(correlation_results)) {
  sp <- correlation_results$species_group[i]
  loc <- correlation_results$location[i]
  data_subset <- all_tomo_data %>% filter(species_group == sp, location == loc)
  
  if (nrow(data_subset) >= 3) {
    cor_test <- cor.test(data_subset$ert_cv, data_subset$CH4_mean)
    correlation_results$p_value[i] <- round(cor_test$p.value, 4)
    correlation_results$sig[i] <- ifelse(cor_test$p.value < 0.05, "*", "")
  } else {
    correlation_results$p_value[i] <- NA
    correlation_results$sig[i] <- ""
  }
}

print(correlation_results)

# Overall correlation
message("\n--- OVERALL CORRELATION (all trees) ---")
overall_cor <- cor.test(all_tomo_data$ert_cv, all_tomo_data$CH4_mean)
message("r = ", round(overall_cor$estimate, 3))
message("p = ", round(overall_cor$p.value, 4))

# By location only
message("\n--- CORRELATION BY LOCATION ---")
wetland_cor <- cor.test(
  all_tomo_data %>% filter(location == "Wetland") %>% pull(ert_cv),
  all_tomo_data %>% filter(location == "Wetland") %>% pull(CH4_mean)
)
upland_cor <- cor.test(
  all_tomo_data %>% filter(location == "Upland") %>% pull(ert_cv),
  all_tomo_data %>% filter(location == "Upland") %>% pull(CH4_mean)
)

message("Wetland: r = ", round(wetland_cor$estimate, 3), ", p = ", round(wetland_cor$p.value, 4))
message("Upland: r = ", round(upland_cor$estimate, 3), ", p = ", round(upland_cor$p.value, 4))

message("\n", paste(rep("=", 60), collapse = ""))

# ============================================================
# PCA-BASED ERT ANALYSIS
# (following approach from Tree-Tomography project)
# ============================================================

library(scales)

message("\n", paste(rep("=", 60), collapse = ""))
message("PCA-BASED ERT ANALYSIS")
message(paste(rep("=", 60), collapse = ""))

# --- Prepare metrics ---
pca_metrics <- c("ert_mean", "ert_median", "ert_sd", "ert_cv",
                 "ert_gini", "ert_entropy", "abs_cma", "abs_radgrad")

all_tomo_data <- all_tomo_data %>%
  mutate(abs_cma = abs(ert_cma),
         abs_radgrad = abs(ert_radialgradiant))

# --- Species-normalize using within-species z-scores ---
spp_stats <- all_tomo_data %>%
  group_by(species_group) %>%
  summarise(across(all_of(pca_metrics),
                   list(mu = ~ mean(., na.rm = TRUE),
                        sigma = ~ sd(., na.rm = TRUE))),
            .groups = "drop")

pca_normed <- all_tomo_data %>%
  select(tree, species_group, all_of(pca_metrics)) %>%
  left_join(spp_stats, by = "species_group")

for (m in pca_metrics) {
  mu_col <- paste0(m, "_mu")
  sd_col <- paste0(m, "_sigma")
  pca_normed[[m]] <- (pca_normed[[m]] - pca_normed[[mu_col]]) / pca_normed[[sd_col]]
}

pca_input <- pca_normed %>% select(all_of(pca_metrics)) %>% as.matrix()
pca_input[is.nan(pca_input)] <- 0

# --- Fit PCA ---
pca_fit <- prcomp(pca_input, center = FALSE, scale. = FALSE)

ve <- summary(pca_fit)$importance[2, 1:3] * 100
message("\nPCA variance explained:")
message("  PC1: ", round(ve[1], 1), "%  PC2: ", round(ve[2], 1), "%  PC3: ", round(ve[3], 1), "%")

message("\nPC1 loadings:")
pc1_load <- round(pca_fit$rotation[, 1], 3)
for (nm in names(pc1_load)) message("  ", nm, ": ", pc1_load[nm])

all_tomo_data$pc1 <- pca_fit$x[, 1]
all_tomo_data$pc2 <- pca_fit$x[, 2]

# Flip so high PC1 = wet/anomalous (mean loading should be negative)
if (pca_fit$rotation["ert_mean", 1] > 0) {
  all_tomo_data$pc1 <- -all_tomo_data$pc1
  pca_fit$rotation[, 1] <- -pca_fit$rotation[, 1]
  message("Flipped PC1 so high = wet/anomalous")
}

# --- Phase diagram: PC1 x SoT with flux as color ---
sot_threshold <- 1
ert_threshold <- mean(all_tomo_data$pc1, na.rm = TRUE)

all_tomo_data <- all_tomo_data %>%
  mutate(
    pca_quadrant = case_when(
      sot_damaged <= sot_threshold & pc1 <= ert_threshold ~ "I: Sound",
      sot_damaged <= sot_threshold & pc1 >  ert_threshold ~ "II: Incipient",
      sot_damaged >  sot_threshold & pc1 >  ert_threshold ~ "III: Active",
      sot_damaged >  sot_threshold & pc1 <= ert_threshold ~ "IV: Cavity"
    )
  )

message("\n--- PCA PHASE DIAGRAM QUADRANTS ---")
print(table(all_tomo_data$pca_quadrant))

# Quadrant colors (severity ramp)
quad_fill <- c("I: Sound" = "#4E79A7", "II: Incipient" = "#E5C460",
               "III: Active" = "#D4873F", "IV: Cavity" = "#C4524E")
quad_text <- c("I: Sound" = "#355570", "II: Incipient" = "#A89030",
               "III: Active" = "#9A5F28", "IV: Cavity" = "#8C3535")

# Signed sqrt transform
signed_sqrt_trans <- trans_new(
  name = "signed_sqrt",
  transform = function(x) sign(x) * sqrt(abs(x)),
  inverse = function(x) sign(x) * x^2
)

pc1_range <- range(all_tomo_data$pc1, na.rm = TRUE)
pc1_max_dev <- max(abs(pc1_range - ert_threshold))
pc1_lim <- ert_threshold + c(-1, 1) * pc1_max_dev * 1.15
sot_max <- max(all_tomo_data$sot_damaged, na.rm = TRUE) * 1.1
sot_lim <- c(-1.5, sot_max)

spp_shapes <- c("N. sylvatica" = 15, "A. rubrum" = 17,
                "T. canadensis" = 16, "Q. rubra" = 18)

# Phase diagram with flux as color
p_phase <- ggplot(all_tomo_data, aes(x = pc1, y = sot_damaged)) +
  # Quadrant shading
  annotate("rect", xmin = pc1_lim[1], xmax = ert_threshold,
           ymin = -2, ymax = sot_threshold, fill = quad_fill[1], alpha = 0.07) +
  annotate("rect", xmin = ert_threshold, xmax = pc1_lim[2],
           ymin = -2, ymax = sot_threshold, fill = quad_fill[2], alpha = 0.07) +
  annotate("rect", xmin = ert_threshold, xmax = pc1_lim[2],
           ymin = sot_threshold, ymax = sot_max, fill = quad_fill[3], alpha = 0.07) +
  annotate("rect", xmin = pc1_lim[1], xmax = ert_threshold,
           ymin = sot_threshold, ymax = sot_max, fill = quad_fill[4], alpha = 0.07) +
  # Threshold lines
  geom_hline(yintercept = sot_threshold, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = ert_threshold, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  # Quadrant labels
  annotate("text", x = pc1_lim[1], y = sot_lim[1], label = "I: Sound",
           hjust = 0, vjust = 0, color = quad_text[1], fontface = "bold", size = 4) +
  annotate("text", x = pc1_lim[2], y = sot_lim[1], label = "II: Incipient",
           hjust = 1, vjust = 0, color = quad_text[2], fontface = "bold", size = 4) +
  annotate("text", x = pc1_lim[2], y = sot_max, label = "III: Active",
           hjust = 1, vjust = 1, color = quad_text[3], fontface = "bold", size = 4) +
  annotate("text", x = pc1_lim[1], y = sot_max, label = "IV: Cavity",
           hjust = 0, vjust = 1, color = quad_text[4], fontface = "bold", size = 4) +
  # Points: shape=species, color=flux
  geom_point(aes(shape = species_group, fill = asinh(CH4_mean)),
             size = 3.5, alpha = 0.85, color = "grey30", stroke = 0.3,
             position = position_jitter(width = 0, height = 0.1, seed = 42)) +
  scale_shape_manual(name = "Species", values = c(21, 22, 23, 24),
                     labels = names(spp_shapes)) +
  scale_fill_viridis_c(option = "inferno", name = expression(asinh(CH[4]~flux)),
                       direction = 1) +
  scale_x_continuous(trans = signed_sqrt_trans) +
  scale_y_continuous(trans = signed_sqrt_trans,
                     breaks = c(0, 1, 5, 10, 20, 30, 40, 50)) +
  coord_cartesian(xlim = pc1_lim, ylim = sot_lim) +
  labs(x = "ERT PC1 (species-normalized)\nhigh = wet / anomalous",
       y = "Structural loss (%, SoT)") +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "right"
  ) +
  guides(shape = guide_legend(override.aes = list(size = 3, fill = "grey50")))

ggsave(file.path(OUTPUT_DIR, "phase_diagram_flux.png"), p_phase,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "phase_diagram_flux.pdf"), p_phase,
       width = 10, height = 8, bg = "white")
message("\nSaved: phase_diagram_flux.png/pdf")

# --- CH4 flux by PCA-derived quadrant ---
message("\n--- CH4 FLUX BY PCA QUADRANT ---")
flux_by_pca_quad <- all_tomo_data %>%
  group_by(pca_quadrant) %>%
  summarise(
    n = n(),
    mean_flux = round(mean(CH4_mean), 3),
    se_flux = round(sd(CH4_mean) / sqrt(n()), 3),
    median_flux = round(median(CH4_mean), 3),
    .groups = "drop"
  ) %>%
  arrange(pca_quadrant)
print(flux_by_pca_quad)

# --- Multi-metric correlations: all ERT metrics + PC1 vs CH4 flux ---
message("\n--- MULTI-METRIC CORRELATIONS: ERT METRICS vs CH4 FLUX ---")

ert_metrics_to_test <- c("ert_mean", "ert_median", "ert_sd", "ert_cv",
                          "ert_gini", "ert_entropy", "ert_cma",
                          "ert_radialgradiant", "pc1")

metric_cor_results <- tibble()

for (metric in ert_metrics_to_test) {
  if (!metric %in% names(all_tomo_data)) next

  # Overall
  ct <- cor.test(all_tomo_data[[metric]], all_tomo_data$CH4_mean)
  metric_cor_results <- bind_rows(metric_cor_results, tibble(
    metric = metric, group = "Overall",
    n = nrow(all_tomo_data),
    r = round(ct$estimate, 3),
    p = round(ct$p.value, 4)
  ))

  # By location
  for (loc in c("Wetland", "Upland")) {
    sub <- all_tomo_data %>% filter(location == loc)
    if (nrow(sub) >= 3) {
      ct <- cor.test(sub[[metric]], sub$CH4_mean)
      metric_cor_results <- bind_rows(metric_cor_results, tibble(
        metric = metric, group = loc,
        n = nrow(sub),
        r = round(ct$estimate, 3),
        p = round(ct$p.value, 4)
      ))
    }
  }
}

metric_cor_results <- metric_cor_results %>%
  mutate(sig = ifelse(p < 0.05, "*", ""))

# Print as wide table
metric_cor_wide <- metric_cor_results %>%
  select(metric, group, r, p, sig) %>%
  pivot_wider(names_from = group,
              values_from = c(r, p, sig),
              names_glue = "{group}_{.value}")
print(metric_cor_wide, width = 120)

# Save correlation table
write.csv(metric_cor_results,
          file.path(OUTPUT_DIR, "ert_metric_correlations.csv"),
          row.names = FALSE)
message("Saved: ert_metric_correlations.csv")

# --- PCA biplot ---
loadings <- as.data.frame(pca_fit$rotation[, 1:2])
loadings$metric <- rownames(loadings)
names(loadings)[1:2] <- c("PC1", "PC2")

# Scale arrows to score range
arrow_scale <- min(max(abs(all_tomo_data$pc1)), max(abs(all_tomo_data$pc2))) * 0.8 /
  max(sqrt(loadings$PC1^2 + loadings$PC2^2))
loadings$PC1 <- loadings$PC1 * arrow_scale
loadings$PC2 <- loadings$PC2 * arrow_scale

metric_labels <- c(ert_mean = "Mean", ert_median = "Median", ert_sd = "SD",
                   ert_cv = "CV", ert_gini = "Gini", ert_entropy = "Entropy",
                   abs_cma = "|CMA|", abs_radgrad = "|RadGrad|")
loadings$label <- metric_labels[loadings$metric]

p_biplot <- ggplot(all_tomo_data, aes(x = pc1, y = pc2)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(aes(shape = species_group, fill = asinh(CH4_mean)),
             size = 3, alpha = 0.7, color = "grey30", stroke = 0.3) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")),
               color = "#B22222", linewidth = 0.7, inherit.aes = FALSE) +
  geom_text(data = loadings,
            aes(x = PC1 * 1.15, y = PC2 * 1.15, label = label),
            color = "#B22222", size = 3.5, fontface = "bold", inherit.aes = FALSE) +
  scale_shape_manual(name = "Species", values = c(21, 22, 23, 24)) +
  scale_fill_viridis_c(option = "inferno", name = expression(asinh(CH[4]~flux))) +
  labs(x = paste0("PC1 (", round(ve[1], 1), "% variance)\nhigh = wet / anomalous"),
       y = paste0("PC2 (", round(ve[2], 1), "% variance)")) +
  theme_classic(base_size = 13) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    legend.position = "right"
  ) +
  guides(shape = guide_legend(override.aes = list(size = 3, fill = "grey50")))

ggsave(file.path(OUTPUT_DIR, "pca_biplot_flux.png"), p_biplot,
       width = 10, height = 8, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "pca_biplot_flux.pdf"), p_biplot,
       width = 10, height = 8, bg = "white")
message("Saved: pca_biplot_flux.png/pdf")

# --- PC1 vs CH4 flux scatter by species ---
p_pc1_flux <- ggplot(all_tomo_data, aes(x = pc1, y = CH4_mean, color = species_group)) +
  geom_point(size = 2.5, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, linetype = "dashed") +
  facet_wrap(~ location, scales = "free") +
  labs(x = "ERT PC1 (species-normalized)\nhigh = wet / anomalous",
       y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})),
       color = "Species") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "pc1_vs_flux.png"), p_pc1_flux,
       width = 10, height = 5, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "pc1_vs_flux.pdf"), p_pc1_flux,
       width = 10, height = 5, bg = "white")
message("Saved: pc1_vs_flux.png/pdf")

message("\n--- PC1 vs CH4 FLUX CORRELATIONS ---")
for (loc in c("Overall", "Wetland", "Upland")) {
  sub <- if (loc == "Overall") all_tomo_data else all_tomo_data %>% filter(location == loc)
  ct <- cor.test(sub$pc1, sub$CH4_mean)
  message(sprintf("  %s: r = %.3f, p = %.4f (n = %d)", loc, ct$estimate, ct$p.value, nrow(sub)))
}

message("\n", paste(rep("=", 60), collapse = ""))
message("Done!")
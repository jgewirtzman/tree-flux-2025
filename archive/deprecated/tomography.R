# ============================================================
# Clean ggplot layout of ERT/Sonic tomography with flux bars
# ============================================================

library(tidyverse)
library(magick)
library(patchwork)
library(grid)

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  tomography = "/Users/jongewirtzman/Downloads/tomography_results_compiled.csv",
  hummock = "/Users/jongewirtzman/Downloads/hummock_hollow.csv",
  flux = "/Users/jongewirtzman/Google Drive/Research/tree-flux-2025/data/HF_2023-2025_tree_flux.csv",
  ert_images = "/Users/jongewirtzman/My Drive/Research/Tomography/HarvardForest_Tomography/ERTs_Absolute",
  sonic_images = "/Users/jongewirtzman/My Drive/Research/Tomography/HarvardForest_Tomography/CH4 PITs"
)

OUTPUT_DIR <- "figures/ert_layouts"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD AND PREPARE DATA
# ============================================================

# Tomography
tomography <- read_csv(PATHS$tomography, show_col_types = FALSE)

# Hummock (for Nyssa identification)
hummock <- read_csv(PATHS$hummock, show_col_types = FALSE) %>%
  rename(tree = Tag)

# Flux data
flux_data <- read_csv(PATHS$flux, show_col_types = FALSE) %>%
  mutate(CH4_flux_nmolpm2ps = ifelse(year < 2025, CH4_flux_nmolpm2ps * 1000, CH4_flux_nmolpm2ps)) %>%
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  filter(CH4_flux_nmolpm2ps >= -1, !is.na(PLOT))

tree_flux_means <- flux_data %>%
  group_by(Tree, SPECIES) %>%
  summarize(CH4_mean = mean(CH4_flux_nmolpm2ps, na.rm = TRUE), .groups = "drop")

# Image file lists
ert_images <- list.files(PATHS$ert_images, pattern = "\\.jpg$", full.names = TRUE)
sonic_images <- list.files(PATHS$sonic_images, pattern = "\\.jpg$", full.names = TRUE)

# ============================================================
# IDENTIFY SPECIES/LOCATION GROUPS
# ============================================================

# Nyssa from hummock
nyssa_trees <- hummock %>% filter(Species == "bg") %>% pull(tree)

# Oak from flux
oak_trees <- flux_data %>% filter(SPECIES == "ro") %>% distinct(Tree) %>% pull(Tree)

# Red maple and hemlock by location (BGS = wetland, everything else = upland/EMS)
upland_rm_trees <- flux_data %>% filter(PLOT != "BGS", SPECIES == "rm") %>% distinct(Tree) %>% pull(Tree)
upland_hem_trees <- flux_data %>% filter(PLOT != "BGS", SPECIES == "hem") %>% distinct(Tree) %>% pull(Tree)
wetland_rm_trees <- flux_data %>% filter(PLOT == "BGS", SPECIES == "rm") %>% distinct(Tree) %>% pull(Tree)
wetland_hem_trees <- flux_data %>% filter(PLOT == "BGS", SPECIES == "hem") %>% distinct(Tree) %>% pull(Tree)

# ============================================================
# HELPER FUNCTION: Prepare species data
# ============================================================

prepare_species_data <- function(tree_ids, species_name) {
  
  df <- tomography %>%
    filter(tree %in% tree_ids) %>%
    arrange(ert_cv) %>%
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
  
  # Normalize flux within species using asinh
  df <- df %>%
    mutate(
      flux_transformed = asinh(CH4_mean),
      flux_norm = (flux_transformed - min(flux_transformed)) / 
        (max(flux_transformed) - min(flux_transformed))
    )
  
  df
}

nyssa_data <- prepare_species_data(nyssa_trees, "Nyssa sylvatica")
oak_data <- prepare_species_data(oak_trees, "Quercus rubra")
upland_rm_data <- prepare_species_data(upland_rm_trees, "Acer rubrum (upland)")
upland_hem_data <- prepare_species_data(upland_hem_trees, "Tsuga canadensis (upland)")
wetland_rm_data <- prepare_species_data(wetland_rm_trees, "Acer rubrum (wetland)")
wetland_hem_data <- prepare_species_data(wetland_hem_trees, "Tsuga canadensis (wetland)")

message("Nyssa trees: ", nrow(nyssa_data))
message("Oak trees: ", nrow(oak_data))
message("Upland red maple trees: ", nrow(upland_rm_data))
message("Upland hemlock trees: ", nrow(upland_hem_data))
message("Wetland red maple trees: ", nrow(wetland_rm_data))
message("Wetland hemlock trees: ", nrow(wetland_hem_data))

# Check which datasets have data
has_upland_rm <- nrow(upland_rm_data) > 0
has_upland_hem <- nrow(upland_hem_data) > 0
has_wetland_rm <- nrow(wetland_rm_data) > 0
has_wetland_hem <- nrow(wetland_hem_data) > 0

# ============================================================
# HELPER: Read image as raster for ggplot (keep original colors)
# ============================================================

read_image_as_raster <- function(path, target_size = 200) {
  img <- image_read(path)
  img <- image_resize(img, paste0(target_size, "x", target_size, "!"))
  as.raster(img)
}

# ============================================================
# BUILD GGPLOT
# ============================================================

# Function to create a single species panel
create_species_panel <- function(data, species_label, bad_sonic_indices = c(), annotation_pos = "topright") {
  
  n <- nrow(data)
  img_size <- 1  # Each image takes 1 unit
  
  # Read all images
  message("  Reading images for ", species_label, "...")
  ert_rasters <- map(data$ert_path, read_image_as_raster)
  
  # For sonic, replace bad indices with NULL (will be blank space)
  sonic_rasters <- map(1:n, function(i) {
    if (i %in% bad_sonic_indices) {
      NULL  # Will be blank space
    } else {
      read_image_as_raster(data$sonic_path[i])
    }
  })
  
  # Create base plot for images
  # Layout (bottom to top):
  #   y = 0 to 0.4: CH4 flux bars
  #   y = 0.5 to 0.9: Decay class
  #   y = 1.05 to 2.05: Sonic images
  #   y = 2.1 to 3.1: ERT images
  #   y = 3.4: ERT CV axis guide
  p_images <- ggplot() +
    xlim(-0.8, n) +
    ylim(0, 3.6) +
    coord_fixed(ratio = 1) +
    theme_void() +
    theme(
      plot.margin = margin(5, 5, 5, 5),
      plot.title = element_text(hjust = 0.5, size = 14, face = "italic")  # Smaller species name
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
  p_images <- p_images + annotate("text", x = -0.4, y = 2.6, label = "ERT", size = 5, fontface = "bold", hjust = 0.5)
  p_images <- p_images + annotate("text", x = -0.4, y = 1.55, label = "SoT", size = 5, fontface = "bold", hjust = 0.5)
  p_images <- p_images + annotate("text", x = -0.4, y = 0.7, label = "Class", size = 5, fontface = "bold", hjust = 0.5)
  p_images <- p_images + annotate("text", x = -0.4, y = 0.2, label = "CH[4]", size = 5, fontface = "bold", hjust = 0.5, parse = TRUE)
  
  # Add ERT images (top row, y = 2.1 to 3.1)
  for (i in 1:n) {
    p_images <- p_images + annotation_raster(
      ert_rasters[[i]],
      xmin = i - 1, xmax = i,
      ymin = 2.1, ymax = 3.1
    )
  }
  
  # Add Sonic images (y = 1.05 to 2.05) - skip bad indices (blank space)
  for (i in 1:n) {
    if (!is.null(sonic_rasters[[i]])) {
      p_images <- p_images + annotation_raster(
        sonic_rasters[[i]],
        xmin = i - 1, xmax = i,
        ymin = 1.05, ymax = 2.05
      )
    }
  }
  
  # Add decay category row (y = 0.5 to 0.9)
  # Categories based on Brazee et al. (2011), Marra et al. (2018):
  # I (Sound): high sonic velocity (sot_damaged <= 1) + high resistance (ert_cv <= 0.5)
  # II (Incipient): high sonic velocity (sot_damaged <= 1) + low resistance (ert_cv > 0.5)
  # III (Advanced): low sonic velocity (sot_damaged > 1) + low resistance (ert_cv > 0.5)
  # IV (Cavity): low sonic velocity (sot_damaged > 1) + high resistance (ert_cv <= 0.5)
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
    # Skip class label if this tree has a missing sonic image
    if (i %in% bad_sonic_indices) next
    
    p_images <- p_images + annotate(
      "text",
      x = i - 0.5,
      y = 0.7,
      label = decay_data$decay_cat[i],
      size = 5,
      fontface = "bold"
    )
  }
  
  # Add flux bars (bottom row, y = 0 to 0.4)
  cor_colors <- c("#4393C3", "#D1E5F0", "#F7F7F7", "#FDDBC7", "#D6604D")
  
  flux_bar_data <- data %>%
    mutate(
      xmin = order - 1,
      xmax = order,
      ymin = 0,
      ymax = 0.4,
      fill_color = cor_colors[round(flux_norm * 4) + 1],
      text_color = ifelse(flux_norm > 0.25 & flux_norm < 0.75, "black", "white"),
      flux_label = sprintf("%.2f", CH4_mean)
    )
  
  # Add colored rectangles for flux
  for (i in 1:n) {
    row <- flux_bar_data[i, ]
    p_images <- p_images + annotate(
      "rect",
      xmin = row$xmin, xmax = row$xmax,
      ymin = row$ymin, ymax = row$ymax,
      fill = row$fill_color,
      color = NA
    )
    p_images <- p_images + annotate(
      "text",
      x = (row$xmin + row$xmax) / 2,
      y = (row$ymin + row$ymax) / 2,
      label = row$flux_label,
      color = row$text_color,
      size = 5
    )
  }
  
  # Create scatterplot: ERT CV vs CH4 flux
  cor_test <- cor.test(data$ert_cv, data$CH4_mean)
  r_val <- cor_test$estimate
  p_val <- cor_test$p.value
  
  if (annotation_pos == "bottomright") {
    ann_x <- Inf
    ann_y <- -Inf
    ann_hjust <- 1.1
    ann_vjust <- -0.5
  } else {
    ann_x <- Inf
    ann_y <- Inf
    ann_hjust <- 1.1
    ann_vjust <- 1.5
  }
  
  p_scatter <- ggplot(data, aes(x = ert_cv, y = CH4_mean)) +
    geom_point(size = 2.5, shape = 16)
  
  # Only add regression line if significant (p < 0.05)
  if (p_val < 0.05) {
    p_scatter <- p_scatter + geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8)
  }
  
  p_scatter <- p_scatter +
    labs(
      x = "ERT CV",
      y = expression(CH[4]~flux)
    ) +
    annotate("text", x = ann_x, y = ann_y, 
             label = sprintf("r = %.2f\np = %.3f", r_val, p_val),
             hjust = ann_hjust, vjust = ann_vjust, size = 3.5) +
    coord_fixed(ratio = diff(range(data$ert_cv)) / diff(range(data$CH4_mean))) +
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

message("\nBuilding panels...")

p_nyssa <- create_species_panel(nyssa_data, "Nyssa sylvatica", bad_sonic_indices = c(7), annotation_pos = "topright")
p_oak <- create_species_panel(oak_data, "Quercus rubra", bad_sonic_indices = c(6), annotation_pos = "bottomright")

# Combine wetland species (Nyssa + Oak)
p_wetland_main <- p_nyssa / plot_spacer() / p_oak +
  plot_layout(heights = c(1, 0.05, 1))

# Save wetland main figure
output_path <- file.path(OUTPUT_DIR, "tomography_ggplot.png")
ggsave(output_path, p_wetland_main, width = 14, height = 10, dpi = 300, bg = "white")
message("\nSaved: ", output_path)

output_path_pdf <- file.path(OUTPUT_DIR, "tomography_ggplot.pdf")
ggsave(output_path_pdf, p_wetland_main, width = 14, height = 10, bg = "white")
message("Saved: ", output_path_pdf)

# Build maple/hemlock panels only if data exists
panel_list <- list()
height_list <- c()

# Maples first
if (has_upland_rm) {
  p_upland_rm <- create_species_panel(upland_rm_data, "Acer rubrum (upland)", bad_sonic_indices = c(), annotation_pos = "topright")
  panel_list <- c(panel_list, list(p_upland_rm))
  height_list <- c(height_list, 1)
}

if (has_wetland_rm) {
  if (length(panel_list) > 0) {
    panel_list <- c(panel_list, list(plot_spacer()))
    height_list <- c(height_list, 0.05)
  }
  p_wetland_rm <- create_species_panel(wetland_rm_data, "Acer rubrum (wetland)", bad_sonic_indices = c(), annotation_pos = "topright")
  panel_list <- c(panel_list, list(p_wetland_rm))
  height_list <- c(height_list, 1)
}

# Then hemlocks
if (has_upland_hem) {
  if (length(panel_list) > 0) {
    panel_list <- c(panel_list, list(plot_spacer()))
    height_list <- c(height_list, 0.05)
  }
  p_upland_hem <- create_species_panel(upland_hem_data, "Tsuga canadensis (upland)", bad_sonic_indices = c(), annotation_pos = "topright")
  panel_list <- c(panel_list, list(p_upland_hem))
  height_list <- c(height_list, 1)
}

if (has_wetland_hem) {
  if (length(panel_list) > 0) {
    panel_list <- c(panel_list, list(plot_spacer()))
    height_list <- c(height_list, 0.05)
  }
  p_wetland_hem <- create_species_panel(wetland_hem_data, "Tsuga canadensis (wetland)", bad_sonic_indices = c(), annotation_pos = "topright")
  panel_list <- c(panel_list, list(p_wetland_hem))
  height_list <- c(height_list, 1)
}

# Save maple/hemlock figure if any panels exist
if (length(panel_list) > 0) {
  p_other <- wrap_plots(panel_list, ncol = 1) + plot_layout(heights = height_list)
  
  # Adjust height based on number of species panels
  n_species <- sum(has_wetland_rm, has_wetland_hem, has_upland_rm, has_upland_hem)
  fig_height <- 4.5 * n_species
  
  output_path2 <- file.path(OUTPUT_DIR, "tomography_ggplot_maple_hemlock.png")
  ggsave(output_path2, p_other, width = 14, height = fig_height, dpi = 300, bg = "white")
  message("Saved: ", output_path2)
  
  output_path2_pdf <- file.path(OUTPUT_DIR, "tomography_ggplot_maple_hemlock.pdf")
  ggsave(output_path2_pdf, p_other, width = 14, height = fig_height, bg = "white")
  message("Saved: ", output_path2_pdf)
} else {
  message("No maple/hemlock data available for second figure")
}

message("\nDone!")
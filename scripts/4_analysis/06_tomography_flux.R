# ============================================================
# 06_tomography_flux.R
#
# Internal wood condition (tomography) and stem CH4 flux: Figure 5 and Figure S5.
# Each panel shows the ERT and SoT tomograms of every tree of one species x site, ordered by
# ERT CV, with its decay class (3_trees/01_tomography_classes.R) and mean CH4 flux, next to
# a scatter of flux against ERT CV. Tests: every measurement enters
# flux ~ ERT CV (per SD) + (1 | tree) + (1 | date), Kenward-Roger df (obs_test below).
#
# ERT sign convention: resistivity (Ohm-m) higher = drier wood; ERT CV higher = more
# heterogeneous moisture within the cross-section.
#
# Inputs: data/package/tomography/ERT_application_results.csv (ERT metrics per tomogram),
#         data/package/tomography/images/ERT_absolute/ and SoT/ (tomograms),
#         data/package/trees.csv, data/final/stem_ch4_flux.csv, data/final/tomography_classes.csv
# Outputs (outputs/figures/tomography/):
#   tomography_specialists.png/pdf  Figure 5 (N. sylvatica, Q. rubra; site panels)
#   tomography_generalists.png/pdf  Figure S5 (A. rubrum and T. canadensis at each site)
#   outputs/tables/ert_allobs_tests.csv  every test of Figures 5 and S5 (slope per SD, SE, KR df, p)
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
  tomography = "data/package/tomography/ERT_application_results.csv",
  hummock = "data/package/trees.csv",
  flux = "data/final/stem_ch4_flux.csv",
  ert_images = "data/package/tomography/images/ERT_absolute",
  sonic_images = "data/package/tomography/images/SoT"
)

OUTPUT_DIR <- "outputs/figures/tomography"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Project-wide color palettes (consistent across all figures)
species_colors <- c("N. sylvatica" = "#2A7F7A", "Q. rubra" = "#6E8B3D",
                     "A. rubrum" = "#A7DAD1", "T. canadensis" = "#C9D6A4")
spp_shapes <- c("N. sylvatica" = 15, "A. rubrum" = 17,
                "T. canadensis" = 16, "Q. rubra" = 18)

# ============================================================
# LOAD AND PREPARE DATA
# ============================================================

message("Loading data...")

if (!file.exists(PATHS$tomography)) stop("Tomography data not found: ", PATHS$tomography)
if (!file.exists(PATHS$flux)) stop("Flux data not found: ", PATHS$flux)

tomography <- read_csv(PATHS$tomography, show_col_types = FALSE) %>%
  rename_with(~ paste0("ert_", sub("radialgradient", "radialgradiant", .x)),
              c(mean, median, sd, cv, gini, entropy, cma, radialgradient))

if (file.exists(PATHS$hummock)) {
  hummock <- read_csv(PATHS$hummock, show_col_types = FALSE) %>% transmute(tree = Tree, Species = species)
} else {
  message("  Warning: trees.csv not found")
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

# ---- Tree-level predictor against repeated flux measurements ----
# Wood condition is measured once per tree and flux ~30 times, so the unit of replication is
# the tree. Every measurement enters a mixed model with random intercepts for tree and sampling
# date, flux ~ metric (per SD) + (1 | tree) + (1 | date), on untransformed flux so that the
# effect reads as nmol m-2 s-1 per SD of the metric; the metric is tested with Kenward-Roger
# degrees of freedom (about 8 for 10 trees), because likelihood-ratio tests are anticonservative
# with so few groups. Figures show flux on linear axes, scaled to each panel.
suppressPackageStartupMessages(library(lmerTest))
obs_for <- function(trees, metric = "ert_cv") {
  flux_data %>% filter(!is.na(CH4_flux_nmolpm2ps)) %>%
    inner_join(trees %>% transmute(Tree = as.numeric(tree), m = .data[[metric]]) %>% filter(!is.na(m)), by = "Tree") %>%
    mutate(flux = CH4_flux_nmolpm2ps, y = flux)
}
obs_test <- function(o, adjust_species = FALSE) {
  mu <- mean(unique(o[c("Tree", "m")])$m); sdv <- sd(unique(o[c("Tree", "m")])$m)
  o$x <- (o$m - mu) / sdv
  f <- if (adjust_species) flux ~ x + SPECIES + (1 | Tree) + (1 | date) else flux ~ x + (1 | Tree) + (1 | date)
  fit <- suppressMessages(lmer(f, data = o, REML = TRUE))
  a <- anova(fit, ddf = "Kenward-Roger")["x", ]
  grid <- data.frame(m = seq(min(o$m), max(o$m), length.out = 60)); X <- cbind(1, (grid$m - mu) / sdv)
  b <- fixef(fit)[c("(Intercept)", "x")]; V <- as.matrix(vcov(fit))[c("(Intercept)", "x"), c("(Intercept)", "x")]
  if (adjust_species) b[1] <- mean(predict(fit, re.form = NA)) - b[["x"]] * mean(o$x)   # line at the average species
  fr <- drop(X %*% b); se <- sqrt(rowSums((X %*% V) * X))
  grid$fit <- fr; grid$lo <- fr - 1.96 * se; grid$hi <- fr + 1.96 * se
  list(p = a$`Pr(>F)`, df = a$DenDF, beta = fixef(fit)[["x"]], se_beta = sqrt(V[2, 2]), metric_sd = sdv,
       n_obs = nrow(o), n_trees = n_distinct(o$Tree), line = grid)
}
fmt_p <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("= %.3f", p))
# Linear y axis zoomed to the tree means (the tests use every measurement). Measurements beyond
# the axis are marked at the panel edge, one small grey triangle per tree (explained in the caption).
zoom_edge <- function(o, tm) {
  r <- range(c(tm$y, 0)); top <- r[2] + 0.6 * diff(r); bot <- r[1] - 0.4 * diff(r); pad <- 0 * (top - bot)
  hi <- o %>% filter(flux > top) %>% distinct(Tree, m); lo <- o %>% filter(flux < bot) %>% distinct(Tree, m)
  list(coord_cartesian(ylim = c(bot, top)),   # default 5% margin keeps edge points and markers whole
       if (nrow(hi)) geom_point(data = hi, aes(m, top - pad), inherit.aes = FALSE, shape = 24, size = 1.3,
                                colour = "grey50", fill = "grey50", alpha = 0.6, stroke = 0.2),
       if (nrow(lo)) geom_point(data = lo, aes(m, bot + pad), inherit.aes = FALSE, shape = 25, size = 1.3,
                                colour = "grey50", fill = "grey50", alpha = 0.6, stroke = 0.2))
}

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
  mutate(species_full = species_labels_full[SPECIES])

message("  Merged tomo+flux: ", nrow(tomo_flux), " trees")

# --- Decay classification: taken from the companion tomography paper ---
# scripts/3_trees/01_tomography_classes.R reproduces that paper's classes
# exactly: SoT structural loss > 1 % crossed with species-normalized ERT PC1
# (8 metrics, PCA on the 57 trees with SoT) above the study-set mean.
tomo_classes <- read_csv("data/final/tomography_classes.csv", show_col_types = FALSE) %>%
  transmute(tree = as.numeric(tree), pc1_paper = ert_pc1, sot_structural_loss,
            decay_phase = decay_class)
tomo_flux <- tomo_flux %>%
  left_join(tomo_classes, by = "tree") %>%
  mutate(decay_phase_short = str_extract(decay_phase, "^[IV]+"))

message("\n  Decay phase distribution:")
print(table(tomo_flux$decay_phase))

# Merge decay class back into tomography for image panels
tomography <- tomography %>%
  left_join(tomo_flux %>% select(tree, pc1_paper, decay_phase, decay_phase_short),
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
    annotate("text", x = -0.4, y = 0.7, label = "Class", size = 5, fontface = "bold", hjust = 0.5) +
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

  # Scatter: ERT metric vs every CH4 measurement (grey) and tree means; mixed-model test (see obs_test)
  o <- obs_for(data, metric); tt <- obs_test(o)
  tm <- o %>% group_by(Tree) %>% summarise(m = first(m), y = mean(flux), .groups = "drop")
  p_scatter <- ggplot(tm, aes(m, y)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_point(data = o, colour = "grey55", size = 0.7, alpha = 0.35, stroke = 0)
  if (tt$p < 0.10) {   # alpha = 0.10 for tree-level tests (10 trees)
    p_scatter <- p_scatter +
      geom_ribbon(data = tt$line, aes(m, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = "grey30", alpha = 0.15) +
      geom_line(data = tt$line, aes(m, fit), inherit.aes = FALSE, colour = "black", linewidth = 0.8)
  }
  p_scatter <- p_scatter + geom_point(size = 2.5, shape = 16) + zoom_edge(o, tm) +
    labs(x = metric_x_label, y = expression(CH[4]~flux), subtitle = sprintf("p %s", fmt_p(tt$p))) +
    theme_classic(base_size = 11) +
    theme(axis.line = element_line(linewidth = 0.5),
          axis.ticks = element_line(linewidth = 0.5),
          aspect.ratio = 1,
          plot.margin = margin(15, 15, 15, 15))

  list(images = p_images, scatter = p_scatter)
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

# --- Specialists: Figure 5 ---
message("\nBuilding specialist panels...")
if (nrow(nyssa_data) > 0 && nrow(oak_data) > 0) {
  # ------------------------------------------------------------------
  # Figure 5: trees ordered by ERT CV (moisture heterogeneity; the bioRxiv v1 metric; the other
  # metrics are in Table S2, 07_decay_definitions.R); decay classes from the tomography paper.
  # ------------------------------------------------------------------
  SIG_P <- 0.1   # Figure 5: fit lines and * where p < 0.1 (stated in the caption)
  class_cols <- c("I" = "#4E79A7", "II" = "#E0B83C", "III" = "#D4873F", "IV" = "#B03A2E")
  flux_pal <- colorRampPalette(c("#F4F4F4", "#9ECAE1", "#2171B5", "#08306B"))(101)

  strip_panel <- function(data, title, sp_col) {
    n <- nrow(data)
    er <- map(data$ert_path, read_image_as_raster)
    so <- map(seq_len(n), function(i) if (is.na(data$decay_phase[i])) NULL else read_image_as_raster(data$sonic_path[i]))
    g <- ggplot() + coord_fixed(ratio = 1, xlim = c(0, n), ylim = c(-0.05, 3.55), clip = "off") +
      theme_void() + theme(plot.margin = margin(0, 4, 0, 34),
                           plot.title = element_text(hjust = 0.5, size = 13, face = "italic", margin = margin(b = 2))) +
      ggtitle(title) +
      annotate("segment", x = 0, xend = n, y = 3.32, yend = 3.32,
               arrow = arrow(length = unit(0.07, "inches")), linewidth = 0.35, colour = "grey30") +
      annotate("text", x = n / 2, y = 3.47, label = "ERT CV (heterogeneity of wood moisture)",
               size = 3, colour = "grey20") +
      annotate("text", x = -0.15, y = c(2.6, 1.55, 0.72), label = c("ERT", "SoT", "Class"),
               size = 4, fontface = "bold", hjust = 1) +
      annotate("text", x = -0.15, y = 0.2, label = "CH[4]", parse = TRUE, size = 4, fontface = "bold", hjust = 1)
    for (i in seq_len(n)) {
      g <- g + annotation_raster(er[[i]], xmin = i - 0.97, xmax = i - 0.03, ymin = 2.12, ymax = 3.08)
      if (!is.null(so[[i]])) g <- g + annotation_raster(so[[i]], xmin = i - 0.97, xmax = i - 0.03, ymin = 1.07, ymax = 2.03)
      else g <- g + annotate("text", x = i - 0.5, y = 1.55, label = "no scan", size = 2.6, colour = "grey50")
      cl <- data$decay_phase_short[i]
      if (!is.na(cl)) g <- g + annotate("text", x = i - 0.5, y = 0.72, label = cl, size = 4.2, fontface = "bold", colour = class_cols[cl])
      fn <- data$flux_norm[i]; fc <- flux_pal[round(fn * 100) + 1]
      g <- g + annotate("rect", xmin = i - 0.98, xmax = i - 0.02, ymin = 0.02, ymax = 0.38, fill = fc, colour = NA) +
        annotate("text", x = i - 0.5, y = 0.2, label = formatC(data$CH4_mean[i], format = "fg", digits = 2),
                 size = 3.1, colour = ifelse(fn > 0.55, "white", "black"))
    }
    g
  }
  class_scatter <- function(data, sp_col, thr) {
    o <- obs_for(data); tt <- obs_test(o)
    lab <- sprintf("p %s\n%d trees, %d measurements", fmt_p(tt$p), tt$n_trees, tt$n_obs)
    tm <- o %>% group_by(Tree) %>% summarise(m = first(m), y = mean(flux), .groups = "drop") %>%
      left_join(data %>% transmute(Tree = as.numeric(tree), cls = ifelse(is.na(decay_phase_short), "no SoT", decay_phase_short)), by = "Tree")
    ggplot(tm, aes(m, y)) +
      geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
      geom_point(data = o, colour = "grey55", size = 0.7, alpha = 0.35, stroke = 0) +
      { if (tt$p < SIG_P) list(geom_ribbon(data = tt$line, aes(m, ymin = lo, ymax = hi), inherit.aes = FALSE, fill = sp_col, alpha = 0.15),
                               geom_line(data = tt$line, aes(m, fit), inherit.aes = FALSE, colour = sp_col, linewidth = 0.8)) } +
      geom_point(aes(fill = cls), shape = 21, size = 3.2, colour = "grey20", stroke = 0.3) +
      scale_fill_manual(values = c(class_cols, "no SoT" = "white"), name = "Class", drop = TRUE) +
      zoom_edge(o, tm) +
      labs(x = "ERT CV", y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})), subtitle = lab) +
      theme_classic(base_size = 12) +
      theme(aspect.ratio = 1, legend.position = "none", plot.margin = margin(4, 8, 4, 4),
            plot.title = element_text(face = "italic", hjust = 0.5, size = 12),
            plot.subtitle = element_text(size = 9.5, colour = "grey20", lineheight = 1.05, hjust = 0.5))
  }
  thr <- read_csv("data/final/tomography_classes.csv", show_col_types = FALSE)$ert_pc1_threshold[1]
  ny2 <- prepare_species_data(nyssa_trees, "Nyssa sylvatica", sort_metric = "ert_cv")
  oa2 <- prepare_species_data(oak_trees, "Quercus rubra", sort_metric = "ert_cv")
  tg <- theme(plot.tag = element_text(size = 14, face = "bold"))
  pa <- strip_panel(ny2, "Nyssa sylvatica (wetland)", species_colors[["N. sylvatica"]]) + labs(tag = "a") + tg
  pb <- class_scatter(ny2, species_colors[["N. sylvatica"]], thr) + labs(tag = "b", title = "N. sylvatica") + tg
  pc <- strip_panel(oa2, "Quercus rubra (upland)", species_colors[["Q. rubra"]]) + labs(tag = "c") + tg
  pd <- class_scatter(oa2, species_colors[["Q. rubra"]], thr) + labs(tag = "d", title = "Q. rubra") + tg

  # site-level scatters (all trees), species coloured; each species' p in the legend;
  # a line only where that species' test has p < SIG_P; the species-adjusted test in the subtitle
  site_pc1 <- function(site_name, tag) {
    sd0 <- tomo_flux %>% filter(location == site_name, !is.na(ert_cv))
    o_all <- obs_for(sd0) %>% left_join(sd0 %>% transmute(Tree = as.numeric(tree), species_full), by = "Tree")
    tests <- lapply(split(o_all, o_all$species_full), obs_test)
    st <- tibble(species_full = names(tests), beta = sapply(tests, `[[`, "beta"), p = sapply(tests, `[[`, "p")) %>%
      mutate(lab = sprintf("%s  p %s%s", species_full, fmt_p(p), ifelse(p < SIG_P, " *", "")))
    labs_v <- setNames(st$lab, st$species_full)
    adj <- obs_test(o_all, adjust_species = TRUE)
    sub <- sprintf("All species, species-adjusted: p %s", fmt_p(adj$p))
    lines <- bind_rows(lapply(names(tests)[st$p < SIG_P], function(k) tests[[k]]$line %>% mutate(species_full = k)))
    tm <- o_all %>% group_by(Tree, species_full) %>% summarise(m = first(m), y = mean(flux), .groups = "drop")
    g <- ggplot(tm, aes(m, y, colour = species_full, shape = species_full)) +
      geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
      geom_point(data = o_all, size = 0.6, alpha = 0.25, stroke = 0, show.legend = FALSE)
    if (nrow(lines)) g <- g + geom_ribbon(data = lines, aes(m, ymin = lo, ymax = hi, fill = species_full),
                                          inherit.aes = FALSE, alpha = 0.12, show.legend = FALSE) +
      geom_line(data = lines, aes(m, fit, colour = species_full), inherit.aes = FALSE, linewidth = 0.8, show.legend = FALSE)
    g + geom_point(size = 2.6, alpha = 0.95) +
      scale_colour_manual(values = species_colors, labels = labs_v, name = NULL) +
      scale_fill_manual(values = species_colors, guide = "none") +
      scale_shape_manual(values = spp_shapes, labels = labs_v, name = NULL) +
      zoom_edge(o_all, tm) +
      labs(x = "ERT CV",
           y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1})), title = site_name, subtitle = sub, tag = tag) +
      theme_classic(base_size = 12) + tg +
      theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
            plot.subtitle = element_text(size = 8.5, colour = "grey30", hjust = 0.5, lineheight = 1.05),
            legend.position = "bottom", legend.direction = "vertical", legend.text = element_text(size = 9),
            legend.key.size = unit(0.35, "cm"), aspect.ratio = 1, plot.margin = margin(4, 8, 4, 4))
  }
  # every test of Figure 5 and Figure S5, for the text (outputs/tables/ert_allobs_tests.csv)
  allobs_tests <- bind_rows(lapply(c("Wetland", "Upland"), function(s) {
    sd0 <- tomo_flux %>% filter(location == s, !is.na(ert_cv))
    o_all <- obs_for(sd0) %>% left_join(sd0 %>% transmute(Tree = as.numeric(tree), species_full), by = "Tree")
    bind_rows(lapply(c(sort(unique(o_all$species_full)), "species-adjusted"), function(k) {
      t <- if (k == "species-adjusted") obs_test(o_all, adjust_species = TRUE) else obs_test(o_all %>% filter(species_full == k))
      tibble(site = s, group = k, metric = "ert_cv", flux_change_per_sd = t$beta, se = t$se_beta, kr_df = t$df, p = t$p,
             ert_cv_sd = t$metric_sd, n_trees = t$n_trees, n_obs = t$n_obs)
    }))
  }))
  write_csv(allobs_tests, "outputs/tables/ert_allobs_tests.csv")
  print(allobs_tests)
  pe <- site_pc1("Wetland", "e"); pf <- site_pc1("Upland", "f")
  p_fig5 <- (pa + pb + plot_layout(widths = c(3.1, 1))) / (pc + pd + plot_layout(widths = c(3.1, 1))) /
    (pe + pf) + plot_layout(heights = c(1, 1, 1.35))
  ggsave(file.path(OUTPUT_DIR, "tomography_specialists.png"), p_fig5, width = 11, height = 11.5, dpi = 300, bg = "white")
  ggsave(file.path(OUTPUT_DIR, "tomography_specialists.pdf"), p_fig5, width = 11, height = 11.5, bg = "white")
  message("  Saved: tomography_specialists.png/pdf (Figure 5, ERT CV ordering)")
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
    panels <- create_species_panel(d$data, d$label, annotation_pos = "topright")
    p <- panels$images + panels$scatter + plot_layout(widths = c(3, 1.2))
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

message("Done!")
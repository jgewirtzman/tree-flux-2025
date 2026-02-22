# ============================================================
# flux_figures.R
#
# Main figure for CH4 tree flux by species and location.
# Also runs statistical models (mixed effects, emmeans).
#
# Run AFTER: timeseries.R (which creates the corrected flux file)
#
# Outputs:
#   - fig2_final.png/pdf (main boxplot with inset)
#   - fig2_main_only.png/pdf (without inset)
#   - Statistical model summaries printed to console
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lme4)
  library(emmeans)
  library(patchwork)
})

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  flux = "data/input/HF_2023-2025_tree_flux_corrected.csv"
)

OUTPUT_DIR <- "outputs/figures"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# DATA LOADING AND CLEANING
# ============================================================

message("Loading data...")

fluxes <- read_csv(PATHS$flux, show_col_types = FALSE) %>%
  # Fill missing species from same tree
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  # Filter bad data
  filter(CH4_flux_nmolpm2ps >= -1,
         !is.na(PLOT)) %>%
  # Add derived columns
  mutate(
    date = as.Date(datetime_posx),
    location = factor(
      ifelse(PLOT == "BGS", "Wetland", "Upland"),
      levels = c("Wetland", "Upland")
    ),
    month = month(date),
    doy = yday(date),
    species_full = case_when(
      SPECIES == "bg"  ~ "Nyssa sylvatica",
      SPECIES == "hem" ~ "Tsuga canadensis",
      SPECIES == "rm"  ~ "Acer rubrum",
      SPECIES == "ro"  ~ "Quercus rubra",
      TRUE ~ SPECIES
    ),
    strategy = case_when(
      SPECIES == "bg" ~ "Wetland specialist",
      SPECIES == "ro" ~ "Upland specialist",
      SPECIES %in% c("hem", "rm") ~ "Generalist",
      TRUE ~ NA_character_
    )
  )

message("  Loaded ", nrow(fluxes), " observations")

# ============================================================
# TREE-LEVEL MEANS
# ============================================================

tree_means <- fluxes %>%
  group_by(Tree, SPECIES, location) %>%
  summarise(
    n_obs = n(),
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    species_full = case_when(
      SPECIES == "bg"  ~ "Nyssa sylvatica",
      SPECIES == "hem" ~ "Tsuga canadensis",
      SPECIES == "rm"  ~ "Acer rubrum",
      SPECIES == "ro"  ~ "Quercus rubra",
      TRUE ~ SPECIES
    ),
    strategy = case_when(
      SPECIES == "bg" ~ "Wetland specialist",
      SPECIES == "ro" ~ "Upland specialist",
      SPECIES %in% c("hem", "rm") ~ "Generalist",
      TRUE ~ NA_character_
    )
  )

# Set up factors
species_order <- c("Nyssa sylvatica", "Acer rubrum", "Tsuga canadensis", "Quercus rubra")

tree_means <- tree_means %>%
  mutate(
    species_full = factor(species_full, levels = species_order),
    strategy = factor(strategy, levels = c("Wetland specialist", "Generalist", "Upland specialist"))
  )

# ============================================================
# STATISTICAL MODELS
# ============================================================

message("\nRunning statistical models...")

# Mixed model: species × location with tree as random effect
model_full <- lmer(CH4_flux_nmolpm2ps ~ SPECIES * location + (1 | Tree), 
                   data = fluxes, REML = TRUE)

message("\n--- Full Model Summary ---")
print(summary(model_full))

# Marginal means
emm <- emmeans(model_full, ~ SPECIES | location)
message("\n--- Marginal Means ---")
print(emm)

# Pairwise comparisons
message("\n--- Pairwise Comparisons ---")
print(pairs(emm))

# --------------------------------------------
# Models reflecting ecological design
# --------------------------------------------

# 1. Wetland: specialist (bg) vs generalists (hem, rm)
message("\n--- Wetland Model (specialist vs generalists) ---")
model_wetland <- lmer(CH4_flux_nmolpm2ps ~ SPECIES + (1 | Tree), 
                      data = filter(fluxes, location == "Wetland"))
print(summary(model_wetland))
emm_wetland <- emmeans(model_wetland, ~ SPECIES)
print(pairs(emm_wetland))

# 2. Upland: specialist (ro) vs generalists (hem, rm)
message("\n--- Upland Model (specialist vs generalists) ---")
model_upland <- lmer(CH4_flux_nmolpm2ps ~ SPECIES + (1 | Tree), 
                     data = filter(fluxes, location == "Upland"))
print(summary(model_upland))
emm_upland <- emmeans(model_upland, ~ SPECIES)
print(pairs(emm_upland))

# 3. Generalists only: location effect for hem and rm
message("\n--- Generalists Model (location effect) ---")
model_generalists <- lmer(CH4_flux_nmolpm2ps ~ SPECIES * location + (1 | Tree), 
                          data = filter(fluxes, SPECIES %in% c("hem", "rm")))
print(summary(model_generalists))
emm_gen <- emmeans(model_generalists, ~ location | SPECIES)
print(pairs(emm_gen))

# ============================================================
# FIGURE 2: MAIN BOXPLOT WITH INSET
# ============================================================

message("\nGenerating figures...")

# --------------------------------------------
# Theme and colors
# --------------------------------------------

theme_clean <- theme_classic(base_size = 16) +
  theme(
    text = element_text(color = "gray20"),
    axis.line = element_line(color = "gray40", linewidth = 0.3),
    axis.ticks = element_line(color = "gray40", linewidth = 0.3),
    axis.text = element_text(color = "gray30"),
    legend.key.size = unit(1.2, "lines"),
    plot.margin = margin(5, 5, 5, 5)
  )

strategy_colors <- c(
  "Wetland specialist" = "#3D7C9C",
  "Generalist" = "#888888",
  "Upland specialist" = "#C2703D"
)

# --------------------------------------------
# Main plot (faceted by species)
# --------------------------------------------

fig2_facet <- ggplot(tree_means, 
                     aes(x = location, y = mean_CH4, fill = strategy)) +
  geom_boxplot(outlier.shape = NA, 
               width = 0.6,
               linewidth = 0.4,
               color = "gray30") +
  geom_jitter(width = 0.12, 
              size = 2,
              alpha = 0.5, 
              color = "gray20",
              shape = 16) +
  scale_fill_manual(values = strategy_colors, name = "Ecological strategy") +
  scale_x_discrete(drop = TRUE) +
  facet_grid(~ species_full, 
             scales = "free_x", 
             space = "free_x",
             switch = "x") +
  labs(
    x = NULL,
    y = expression("Tree mean CH"[4]~"flux (nmol m"^-2~"s"^-1*")")
  ) +
  theme_clean +
  theme(
    axis.text.x = element_text(size = 14),
    legend.position = "bottom",
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.background = element_rect(color = "gray40", fill = "white", linewidth = 0.3),
    strip.text = element_text(face = "italic", size = 14, margin = margin(t = 5)),
    panel.spacing = unit(0.3, "lines")
  )

# --------------------------------------------
# Inset panels with custom ticks
# --------------------------------------------

short_labels <- c(
  "Nyssa sylvatica" = expression(italic("N. sylvatica")),
  "Acer rubrum" = expression(italic("A. rubrum")),
  "Tsuga canadensis" = expression(italic("T. canadensis")),
  "Quercus rubra" = expression(italic("Q. rubra"))
)

theme_inset <- theme_clean +
  theme(
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 9),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.title = element_blank(),
    strip.text = element_text(size = 10, face = "bold", margin = margin(2, 0, 2, 0)),
    strip.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(2, 2, 2, 2)
  )

# Wetland inset
inset_wetland <- tree_means %>%
  filter(location == "Wetland") %>%
  ggplot(aes(x = species_full, y = mean_CH4, fill = strategy)) +
  geom_boxplot(outlier.shape = NA, 
               width = 0.7,
               linewidth = 0.25,
               color = "gray30") +
  geom_jitter(width = 0.12, alpha = 0.4, size = 0.9, color = "gray30") +
  scale_fill_manual(values = strategy_colors) +
  scale_y_continuous(trans = "pseudo_log",
                     breaks = c(0, 2, 10),
                     labels = c("0", "2", "10")) +
  scale_x_discrete(labels = short_labels, drop = TRUE) +
  ggtitle("Wetland") +
  theme_inset +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 2)))

# Upland inset
inset_upland <- tree_means %>%
  filter(location == "Upland") %>%
  ggplot(aes(x = species_full, y = mean_CH4, fill = strategy)) +
  geom_boxplot(outlier.shape = NA, 
               width = 0.7,
               linewidth = 0.25,
               color = "gray30") +
  geom_jitter(width = 0.12, alpha = 0.4, size = 0.9, color = "gray30") +
  scale_fill_manual(values = strategy_colors) +
  scale_y_continuous(trans = "pseudo_log",
                     breaks = c(0, 0.1, 0.2),
                     labels = c("0", "0.1", "0.2")) +
  scale_x_discrete(labels = short_labels, drop = TRUE) +
  ggtitle("Upland") +
  theme_inset +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5, margin = margin(b = 2)))

# Combine inset panels
p_eco_inset <- (inset_wetland | inset_upland) + plot_layout(widths = c(1, 1))

# --------------------------------------------
# Final combined figure
# --------------------------------------------

fig2_final <- fig2_facet +
  inset_element(
    p_eco_inset,
    left = 0.50,
    bottom = 0.50,
    right = 1,
    top = 1
  )

print(fig2_final)

# ============================================================
# SAVE
# ============================================================

ggsave(file.path(OUTPUT_DIR, "fig2_final.pdf"), fig2_final, 
       width = 9, height = 6, units = "in")
ggsave(file.path(OUTPUT_DIR, "fig2_final.png"), fig2_final, 
       width = 9, height = 6, units = "in", dpi = 300)

# Also save without inset
ggsave(file.path(OUTPUT_DIR, "fig2_main_only.pdf"), fig2_facet, 
       width = 8, height = 5.5, units = "in")
ggsave(file.path(OUTPUT_DIR, "fig2_main_only.png"), fig2_facet, 
       width = 8, height = 5.5, units = "in", dpi = 300)

message("\nSaved figures to: ", OUTPUT_DIR)
message("Done!")








# ============================================================
# SUMMARY STATISTICS FOR RESULTS
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("STATISTICAL RESULTS SUMMARY")
message(paste(rep("=", 60), collapse = ""))

# Tree-level mean summary by site × species
message("\n--- Tree-Level Mean CH4 by Site × Species ---")
tree_means %>%
  group_by(location, species_full) %>%
  summarise(
    n_trees = n(),
    mean = round(mean(mean_CH4), 3),
    se = round(sd(mean_CH4) / sqrt(n()), 3),
    median = round(median(mean_CH4), 3),
    .groups = "drop"
  ) %>%
  arrange(location, species_full) %>%
  print()

# Full model: extract fixed effects with CIs
message("\n--- Full Model Fixed Effects ---")
full_fixef <- as.data.frame(summary(model_full)$coefficients)
full_fixef$term <- rownames(full_fixef)
full_fixef <- full_fixef %>%
  mutate(
    ci_lower = Estimate - 1.96 * `Std. Error`,
    ci_upper = Estimate + 1.96 * `Std. Error`
  ) %>%
  dplyr::select(term, Estimate, `Std. Error`, ci_lower, ci_upper, `t value`)
print(full_fixef)

# Variance components
message("\n--- Variance Components (Full Model) ---")
vc <- as.data.frame(VarCorr(model_full))
print(vc)
icc_full <- vc$vcov[1] / sum(vc$vcov)
message("ICC (tree): ", round(icc_full, 3))

# Marginal means with CIs
message("\n--- Marginal Means (emmeans) ---")
emm_df <- as.data.frame(emm)
print(emm_df)

# Pairwise contrasts with effect sizes
message("\n--- Pairwise Contrasts (with 95% CI) ---")
pairs_df <- as.data.frame(pairs(emm))
print(pairs_df)

# Wetland model contrasts
message("\n--- Wetland Pairwise Contrasts ---")
pairs_wetland <- as.data.frame(pairs(emm_wetland))
print(pairs_wetland)

# Upland model contrasts
message("\n--- Upland Pairwise Contrasts ---")
pairs_upland <- as.data.frame(pairs(emm_upland))
print(pairs_upland)

# Generalists: location effect
message("\n--- Generalists Location Effect ---")
pairs_gen <- as.data.frame(pairs(emm_gen))
print(pairs_gen)

# Effect size: wetland vs upland (for generalists)
message("\n--- Location Effect Size for Generalists ---")
emm_gen_df <- as.data.frame(emmeans(model_generalists, ~ location | SPECIES))
emm_gen_df %>%
  dplyr::select(SPECIES, location, emmean, SE) %>%
  pivot_wider(names_from = location, values_from = c(emmean, SE)) %>%
  mutate(
    diff = emmean_Wetland - emmean_Upland,
    fold_change = round(emmean_Wetland / pmax(emmean_Upland, 0.001), 1)
  ) %>%
  print()

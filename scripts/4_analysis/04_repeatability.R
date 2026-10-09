# ============================================================
# repeatability.R
#
# Tree-level repeatability and consistency analyses.
# ICC, likelihood ratio tests, Spearman correlations, z-score tracks.
#
# Run AFTER: timeseries.R (which creates the corrected flux file)
#
# Outputs:
#   All statistics on asinh(flux), the scale of the models.
#   - repeatability_summary.csv: ICC, Spearman rho by species/location
#   - fig_period_cor.png/pdf: Early vs late period scatter
#   - fig_blups.png/pdf: Tree random effects (exploratory)
#   - zscore_tracks.png/pdf: Z-score time series (Figure 3)
#   - outputs/tables/synchrony_by_group.csv: date share of within-tree variance
#   - fig_tree_effects.png/pdf: Tree random effects with 95% intervals (Figure S3)
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(lme4)
  library(performance)
  library(ggtext)
  library(scales)
})

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  flux = "data/final/stem_ch4_flux.csv"
)

OUTPUT_DIR <- "outputs/figures/repeatability"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# DATA LOADING AND CLEANING
# ============================================================

message("Loading data...")

fluxes <- read_csv(PATHS$flux, show_col_types = FALSE) %>%
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  filter(!is.na(PLOT)) %>%
  mutate(
    date = as.Date(datetime_posx),
    location = factor(
      ifelse(PLOT == "BGS", "Wetland", "Upland"),
      levels = c("Wetland", "Upland")
    ),
    species_full = case_when(
      SPECIES == "bg"  ~ "Nyssa sylvatica",
      SPECIES == "hem" ~ "Tsuga canadensis",
      SPECIES == "rm"  ~ "Acer rubrum",
      SPECIES == "ro"  ~ "Quercus rubra",
      TRUE ~ SPECIES
    ),
    facet_label = paste(location, "-", species_full),
    flux_asinh = asinh(CH4_flux_nmolpm2ps)   # all repeatability statistics use the model scale
  )

message("  Loaded ", nrow(fluxes), " observations")

# ============================================================
# 1. ICC (Intraclass Correlation Coefficient)
# ============================================================
# ICC = proportion of variance due to between-tree differences
# High ICC = trees are consistently different from each other

message("\nCalculating ICC by group...")

icc_by_group <- fluxes %>%
  filter(!is.na(facet_label)) %>%
  group_by(SPECIES, location) %>%
  summarise(
    icc = {
      mod <- lmer(flux_asinh ~ 1 + (1 | Tree), data = cur_data())
      var_tree <- as.numeric(VarCorr(mod)$Tree)
      var_resid <- sigma(mod)^2
      var_tree / (var_tree + var_resid)
    },
    n_trees = n_distinct(Tree),
    n_obs = n(),
    .groups = "drop"
  )

message("\n--- ICC by Species × Location ---")
print(icc_by_group)
# ICC > 0.5 = strong consistency; 0.2-0.5 = moderate; < 0.2 = weak

# ============================================================
# 2. Likelihood Ratio Test: does Tree matter?
# ============================================================

message("\nRunning likelihood ratio tests...")

lrt_by_group <- fluxes %>%
  filter(!is.na(facet_label)) %>%
  group_by(SPECIES, location) %>%
  reframe(
    chisq = {
      mod_tree <- lmer(flux_asinh ~ 1 + (1 | Tree), data = cur_data(), REML = FALSE)
      mod_null <- lm(flux_asinh ~ 1, data = cur_data())
      as.numeric(2 * (logLik(mod_tree) - logLik(mod_null)))
    },
    p_value = pchisq(chisq, df = 1, lower.tail = FALSE)
  )

message("\n--- Likelihood Ratio Tests ---")
print(lrt_by_group)

# ============================================================
# 3. Spearman Correlation Across Time Periods
# ============================================================
# Do tree ranks correlate between early and late periods?

message("\nCalculating temporal consistency (Spearman correlation)...")

tree_by_period <- fluxes %>%
  filter(!is.na(facet_label)) %>%
  mutate(period = ifelse(date < as.Date("2024-06-01"), "early", "late")) %>%
  group_by(Tree, SPECIES, location, period) %>%
  summarise(mean_CH4 = mean(flux_asinh, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_CH4) %>%
  filter(!is.na(early), !is.na(late))

cor_by_group <- tree_by_period %>%
  group_by(SPECIES, location) %>%
  summarise(
    spearman_rho = cor(early, late, method = "spearman"),
    p_value = cor.test(early, late, method = "spearman")$p.value,
    n_trees = n(),
    .groups = "drop"
  )

message("\n--- Spearman Correlations (Early vs Late) ---")
print(cor_by_group)

# ============================================================
# 4. Summary Table for Publication
# ============================================================

summary_table <- icc_by_group %>%
  left_join(cor_by_group %>% dplyr::select(SPECIES, location, spearman_rho, p_value), 
            by = c("SPECIES", "location")) %>%
  rename(cor_p = p_value) %>%
  mutate(
    species_full = case_when(
      SPECIES == "bg"  ~ "N. sylvatica",
      SPECIES == "hem" ~ "T. canadensis",
      SPECIES == "rm"  ~ "A. rubrum",
      SPECIES == "ro"  ~ "Q. rubra"
    ),
    sig = ifelse(cor_p < 0.05, "*", "")
  ) %>%
  dplyr::select(species_full, location, n_trees, n_obs, icc, spearman_rho, cor_p, sig)

message("\n", paste(rep("=", 60), collapse = ""))
message("SUMMARY TABLE FOR PUBLICATION")
message(paste(rep("=", 60), collapse = ""))
print(summary_table)

write_csv(summary_table, file.path(OUTPUT_DIR, "repeatability_summary.csv"))

# ============================================================
# FIGURE 1: Early vs Late Period Correlation
# ============================================================

message("\nGenerating figures...")

# Add species labels
tree_by_period <- tree_by_period %>%
  mutate(
    species_full = case_when(
      SPECIES == "bg"  ~ "N. sylvatica",
      SPECIES == "hem" ~ "T. canadensis",
      SPECIES == "rm"  ~ "A. rubrum",
      SPECIES == "ro"  ~ "Q. rubra"
    ),
    species_full = factor(species_full, 
                          levels = c("A. rubrum", "T. canadensis", "Q. rubra", "N. sylvatica")),
    location = factor(location, levels = c("Wetland", "Upland"))
  )

# Correlation labels for plot
cor_labels <- cor_by_group %>%
  mutate(
    species_full = case_when(
      SPECIES == "bg"  ~ "N. sylvatica",
      SPECIES == "hem" ~ "T. canadensis",
      SPECIES == "rm"  ~ "A. rubrum",
      SPECIES == "ro"  ~ "Q. rubra"
    ),
    species_full = factor(species_full, 
                          levels = c("A. rubrum", "T. canadensis", "Q. rubra", "N. sylvatica")),
    location = factor(location, levels = c("Wetland", "Upland")),
    label = paste0("ρ = ", round(spearman_rho, 2), 
                   ifelse(p_value < 0.05, "*", ""),
                   ifelse(p_value < 0.01, "*", ""))
  )

# one panel per species x site (the six study groups), free scales; Spearman rho as in the text
grp_levels <- c("Wetland · N. sylvatica", "Wetland · A. rubrum", "Wetland · T. canadensis",
                "Upland · Q. rubra", "Upland · A. rubrum", "Upland · T. canadensis")
tbp <- tree_by_period %>% mutate(group = factor(paste(location, "·", species_full), levels = grp_levels))
lab <- cor_labels %>% mutate(group = factor(paste(location, "·", species_full), levels = grp_levels),
                             label = sprintf("ρ = %.2f, p %s", spearman_rho, ifelse(p_value < 0.001, "< 0.001", sprintf("= %.3f", p_value))))
fig_period_cor <- ggplot(tbp, aes(x = early, y = late)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_point(size = 2.6, alpha = 0.8) +
  geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.4, size = 3.3) +
  facet_wrap(~ group, nrow = 2, scales = "free") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.2))) +
  labs(x = expression("Tree mean of asinh(CH"[4]~"flux), before June 2024"),
       y = expression("Tree mean of asinh(CH"[4]~"flux), after June 2024")) +
  theme_classic(base_size = 11) +
  theme(strip.text = element_text(face = "italic"), strip.background = element_blank())

print(fig_period_cor)

ggsave(file.path(OUTPUT_DIR, "fig_period_cor.png"), fig_period_cor, 
       width = 9, height = 5.6, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "fig_period_cor.pdf"), fig_period_cor, 
       width = 9, height = 6)

message("  Saved: fig_period_cor.png/pdf")

# ============================================================
# FIGURE 2: Tree Random Effects (BLUPs) - for SI
# ============================================================

model <- lmer(flux_asinh ~ SPECIES + location + (1 | Tree), data = fluxes)

tree_effects <- ranef(model)$Tree %>%
  as.data.frame() %>%
  rownames_to_column("Tree") %>%
  rename(tree_effect = `(Intercept)`) %>%
  mutate(Tree = as.integer(Tree)) %>%
  left_join(fluxes %>% distinct(Tree, SPECIES, location), by = "Tree") %>%
  arrange(desc(tree_effect))

fig_blups <- ggplot(tree_effects, aes(x = reorder(factor(Tree), tree_effect), 
                                      y = tree_effect, fill = SPECIES)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~location, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Tree", y = "Random effect (deviation from mean)",
       title = "Tree-level CH4 flux effects (BLUPs)") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 6))

print(fig_blups)

ggsave(file.path(OUTPUT_DIR, "fig_blups.png"), fig_blups, 
       width = 12, height = 8, dpi = 300)
ggsave(file.path(OUTPUT_DIR, "fig_blups.pdf"), fig_blups, 
       width = 12, height = 8)

message("  Saved: fig_blups.png/pdf")

# ============================================================
# FIGURE 3: Z-Score Tracks (Main Figure)
# ============================================================

message("  Creating z-score tracks figure...")

# Z-score within year × location × species
fluxes_zscore <- fluxes %>%
  mutate(
    date = as.Date(date),
    year = year(date)
  ) %>%
  group_by(year, location, SPECIES) %>%
  mutate(
    group_mean = mean(flux_asinh, na.rm = TRUE),
    group_sd   = sd(flux_asinh, na.rm = TRUE),
    z_score    = (flux_asinh - group_mean) / group_sd
  ) %>%
  ungroup() %>%
  mutate(
    z_score = ifelse(is.finite(z_score), z_score, 0)
  )

# Mean z-score per tree (for ordering)
tree_mean_zscore <- fluxes_zscore %>%
  group_by(Tree, location, SPECIES) %>%
  summarise(mean_zscore = mean(z_score, na.rm = TRUE), .groups = "drop")

fluxes_zscore <- fluxes_zscore %>%
  left_join(tree_mean_zscore, by = c("Tree", "location", "SPECIES"))

# Facet labels - Wetland first
fluxes_zscore <- fluxes_zscore %>%
  mutate(
    species_full = factor(
      species_full,
      levels = c("Acer rubrum", "Tsuga canadensis", "Nyssa sylvatica", "Quercus rubra")
    ),
    facet_key = paste(as.character(location), as.character(species_full), sep = " - "),
    facet_label_z = case_when(
      facet_key == "Upland - Acer rubrum"        ~ "<i>Acer rubrum</i><br>Upland",
      facet_key == "Upland - Tsuga canadensis"   ~ "<i>Tsuga canadensis</i><br>Upland",
      facet_key == "Upland - Quercus rubra"      ~ "<i>Quercus rubra</i><br>Upland",
      facet_key == "Wetland - Acer rubrum"       ~ "<i>Acer rubrum</i><br>Wetland",
      facet_key == "Wetland - Tsuga canadensis"  ~ "<i>Tsuga canadensis</i><br>Wetland",
      facet_key == "Wetland - Nyssa sylvatica"   ~ "<i>Nyssa sylvatica</i><br>Wetland",
      TRUE ~ NA_character_
    ),
    facet_label_z = factor(
      facet_label_z,
      levels = c(
        "<i>Acer rubrum</i><br>Wetland",
        "<i>Tsuga canadensis</i><br>Wetland",
        "<i>Nyssa sylvatica</i><br>Wetland",
        "<i>Acer rubrum</i><br>Upland",
        "<i>Tsuga canadensis</i><br>Upland",
        "<i>Quercus rubra</i><br>Upland"
      )
    )
  ) %>%
  filter(!is.na(facet_label_z))

# Global measurement index (even spacing across all measurement dates)
all_dates <- sort(unique(fluxes_zscore$date))
fluxes_zscore <- fluxes_zscore %>%
  mutate(
    date_idx = match(date, all_dates),
    z_score_abs = pmin(abs(z_score), 4)
  ) %>%
  group_by(facet_label_z) %>%
  mutate(
    tree_ordered = reorder(factor(Tree), mean_zscore)
  ) %>%
  ungroup()

# X-axis labels
n_breaks <- 7
break_idx <- unique(round(seq(1, length(all_dates), length.out = n_breaks)))
break_labs <- format(all_dates[break_idx], "%b '%y")

# Plot
fig_zscore_tracks <- ggplot(
  fluxes_zscore,
  aes(
    x = date_idx,
    y = tree_ordered,
    color = z_score,
    size = z_score_abs
  )
) +
  # Subtle connecting lines
  geom_line(
    aes(group = factor(Tree), color = mean_zscore), 
    alpha = 0.25, 
    linewidth = 0.3, 
    show.legend = FALSE
  ) +
  # Points
  geom_point(alpha = 0.8, shape = 16) +
  # Color scale
  scale_color_gradient2(
    low = "#2A7F7A",
    mid = "grey92",
    high = "#C46A2F",
    midpoint = 0,
    limits = c(-2, 4),
    oob = squish,
    name = "Z-score"
  ) +
  scale_size_continuous(range = c(0.5, 4.5), guide = "none") +
  scale_x_continuous(
    breaks = break_idx,
    labels = break_labs,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  facet_wrap(~ facet_label_z, ncol = 3, scales = "free_y") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8, color = "grey30", angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    strip.background = element_blank(),
    strip.text = element_markdown(size = 10, lineheight = 1.1, color = "grey20"),
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.6, "lines"),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "plain", color = "grey30"),
    legend.text = element_text(size = 8, color = "grey40"),
    legend.key.height = unit(35, "pt"),
    legend.key.width = unit(8, "pt"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "grey98", color = NA)
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      frame.colour = "grey80",
      frame.linewidth = 0.3,
      ticks = FALSE
    )
  )

print(fig_zscore_tracks)

ggsave(file.path(OUTPUT_DIR, "zscore_tracks.png"), fig_zscore_tracks, 
       width = 7, height = 5, dpi = 600, bg = "white")
ggsave(file.path(OUTPUT_DIR, "zscore_tracks.pdf"), fig_zscore_tracks, 
       width = 10, height = 7)

  message("  Saved: zscore_tracks.png/pdf")

 message("\nOutputs saved to: ", OUTPUT_DIR)
message("Done!")








# ============================================================
# REPEATABILITY RESULTS SUMMARY
# ============================================================

message("\n", paste(rep("=", 60), collapse = ""))
message("REPEATABILITY RESULTS FOR MANUSCRIPT")
message(paste(rep("=", 60), collapse = ""))

# Combined ICC + LRT + Spearman table
message("\n--- Combined Repeatability Table ---")
repeatability_full <- icc_by_group %>%
  left_join(lrt_by_group, by = c("SPECIES", "location")) %>%
  left_join(cor_by_group, by = c("SPECIES", "location", "n_trees")) %>%
  mutate(
    species_full = case_when(
      SPECIES == "bg"  ~ "N. sylvatica",
      SPECIES == "hem" ~ "T. canadensis",
      SPECIES == "rm"  ~ "A. rubrum",
      SPECIES == "ro"  ~ "Q. rubra"
    ),
    icc = round(icc, 3),
    spearman_rho = round(spearman_rho, 2),
    lrt_p = ifelse(p_value.x < 0.001, "<0.001", round(p_value.x, 3)),
    cor_p = ifelse(p_value.y < 0.001, "<0.001", round(p_value.y, 3))
  ) %>%
  dplyr::select(species_full, location, n_trees, n_obs, icc, lrt_p, spearman_rho, cor_p)
print(repeatability_full)

# Overall ICC summary
message("\n--- ICC Summary Across Groups ---")
message("Range: ", round(min(icc_by_group$icc), 3), " - ", round(max(icc_by_group$icc), 3))
message("Mean: ", round(mean(icc_by_group$icc), 3))
message("Groups with ICC > 0.5 (strong): ", sum(icc_by_group$icc > 0.5))
message("Groups with ICC > 0.2 (moderate+): ", sum(icc_by_group$icc > 0.2))

# Spearman summary
message("\n--- Temporal Consistency Summary ---")
message("Spearman rho range: ", round(min(cor_by_group$spearman_rho), 2), 
        " - ", round(max(cor_by_group$spearman_rho), 2))
message("Significant correlations (p < 0.05): ", sum(cor_by_group$p_value < 0.05), 
        " of ", nrow(cor_by_group))

# Which groups show strongest repeatability?
message("\n--- Strongest Repeatability ---")
repeatability_full %>%
  arrange(desc(icc)) %>%
  head(3) %>%
  print()

# Tree random effects range (from BLUPs)
message("\n--- Tree Random Effects (BLUPs) ---")
message("Range: ", round(min(tree_effects$tree_effect), 2), 
        " to ", round(max(tree_effects$tree_effect), 2))
tree_effects %>%
  group_by(location) %>%
  summarise(
    min_effect = round(min(tree_effect), 2),
    max_effect = round(max(tree_effect), 2),
    range = round(max(tree_effect) - min(tree_effect), 2),
    .groups = "drop"
  ) %>%
  print()
# ============================================================
# SYNCHRONY: how much of the variation within trees is shared across trees on the same date
# asinh(flux) ~ 1 + (1|Tree) + (1|date) per species x site; date share of within-tree variance
# = var(date) / (var(date) + var(residual)); LRT for the date effect (ML fits).
# Output: outputs/tables/synchrony_by_group.csv (read by 5_models/07_manuscript_numbers.R)
# ============================================================
sync_one <- function(d) {
  d$y <- asinh(d$CH4_flux_nmolpm2ps)
  m <- lmer(y ~ 1 + (1 | Tree) + (1 | date), data = d)
  v <- as.data.frame(VarCorr(m)); s <- setNames(v$vcov, v$grp)
  m0 <- lmer(y ~ 1 + (1 | Tree), data = d, REML = FALSE)
  tibble(n_obs = nrow(d), n_dates = n_distinct(d$date),
         tree_share = s[["Tree"]] / sum(s), date_share = s[["date"]] / sum(s), resid_share = s[["Residual"]] / sum(s),
         date_share_within_tree = s[["date"]] / (s[["date"]] + s[["Residual"]]),
         date_lrt_p = anova(m0, update(m, REML = FALSE))$`Pr(>Chisq)`[2])
}
synchrony <- fluxes %>% filter(!is.na(facet_label)) %>% group_by(location, SPECIES) %>%
  group_modify(~ sync_one(.x)) %>% ungroup() %>%
  mutate(species_full = c(rm = "A. rubrum", hem = "T. canadensis", ro = "Q. rubra", bg = "N. sylvatica")[SPECIES])
write_csv(synchrony, "outputs/tables/synchrony_by_group.csv")
message("\n--- Synchrony (date share of within-tree variance, asinh scale) ---"); print(synchrony)

# ============================================================
# SI FIGURE: tree random effects with 95% intervals
# ============================================================
suppressPackageStartupMessages({ library(patchwork); library(lme4) })
grp_levels <- c("A. rubrum|Upland", "T. canadensis|Upland", "Q. rubra|Upland",
                "A. rubrum|Wetland", "T. canadensis|Wetland", "N. sylvatica|Wetland")
sp_lab <- c(rm = "A. rubrum", hem = "T. canadensis", ro = "Q. rubra", bg = "N. sylvatica")
traj <- fluxes %>% mutate(sp = sp_lab[SPECIES], grp = factor(paste(sp, location, sep = "|"), levels = grp_levels),
                          date = as.Date(date)) %>%
  group_by(grp) %>% mutate(rel_mean = { m <- ave(CH4_flux_nmolpm2ps, Tree); r <- rank(m, ties.method = "average"); (r - min(r)) / max(1, max(r) - min(r)) }) %>%
  ungroup()
blup <- traj %>% group_by(grp) %>% group_modify(~ {
  m <- lmer(flux_asinh ~ 1 + (1 | Tree), data = .x)
  re <- ranef(m, condVar = TRUE)$Tree; pv <- attr(re, "postVar")
  tibble(Tree = rownames(re), effect = re[, 1], se = sqrt(pv[1, 1, ]))
}) %>% ungroup() %>% mutate(lo = effect - 1.96 * se, hi = effect + 1.96 * se, excl0 = lo > 0 | hi < 0,
                            grp_lab = sub("\\|", " - ", grp)) %>%
  group_by(grp) %>% mutate(Tree = reorder(Tree, effect)) %>% ungroup()
pB <- ggplot(blup, aes(effect, Tree)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", colour = "grey55") +
  geom_point(aes(colour = excl0), size = 1.8) +
  scale_colour_manual(values = c(`FALSE` = "grey40", `TRUE` = "#B03A2E"), guide = "none") +
  facet_wrap(~ factor(grp_lab, levels = sub("\\|", " - ", grp_levels)), scales = "free", nrow = 2) +
  labs(x = expression(Tree~random~effect~(deviation~from~group~mean*","~asinh~CH[4]~flux)), y = "Tree", tag = "B") +
  theme_classic(base_size = 9) + theme(strip.background = element_blank(), axis.text.y = element_text(size = 6))
pB <- pB + labs(tag = NULL)
ggsave(file.path(OUTPUT_DIR, "fig_tree_effects.png"), pB, width = 11, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(OUTPUT_DIR, "fig_tree_effects.pdf"), pB, width = 11, height = 5.5, bg = "white")
message("  Saved: fig_tree_effects.png/pdf")

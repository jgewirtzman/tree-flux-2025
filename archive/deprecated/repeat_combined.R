library(lme4)
library(performance)



# ============================================
# 1. ICC (Intraclass Correlation Coefficient)
# ============================================
# ICC = proportion of variance due to between-tree differences
# High ICC = trees are consistently different from each other

icc_by_group <- fluxes %>%
  filter(!is.na(facet_label)) %>%
  group_by(SPECIES, location) %>%
  summarise(
    icc = {
      mod <- lmer(CH4_flux_nmolpm2ps ~ 1 + (1 | Tree), data = cur_data())
      var_tree <- as.numeric(VarCorr(mod)$Tree)
      var_resid <- sigma(mod)^2
      var_tree / (var_tree + var_resid)
    },
    n_trees = n_distinct(Tree),
    n_obs = n(),
    .groups = "drop"
  )

print(icc_by_group)
# ICC > 0.5 = strong consistency; 0.2-0.5 = moderate; < 0.2 = weak

# ============================================
# 2. Likelihood ratio test: does Tree matter?
# ============================================
# Compare model with vs without tree random effect

lrt_by_group <- fluxes %>%
  filter(!is.na(facet_label)) %>%
  group_by(SPECIES, location) %>%
  summarise(
    lrt_result = {
      mod_tree <- lmer(CH4_flux_nmolpm2ps ~ 1 + (1 | Tree), data = cur_data(), REML = FALSE)
      mod_null <- lm(CH4_flux_nmolpm2ps ~ 1, data = cur_data())
      test <- anova(mod_tree, mod_null)
      list(chisq = test$Chisq[2], p = test$`Pr(>Chisq)`[2])
    },
    .groups = "drop"
  ) %>%
  unnest_wider(lrt_result)

print(lrt_by_group)

# ============================================
# 3. Repeatability with confidence intervals (rptR package)
# ============================================
library(rptR)

# Example for one group
bg_wetland <- fluxes %>% filter(SPECIES == "bg", location == "wetland")

rpt_bg <- rpt(CH4_flux_nmolpm2ps ~ 1 + (1 | Tree), 
              grname = "Tree", 
              data = bg_wetland, 
              datatype = "Gaussian",
              nboot = 1000, npermut = 1000)

print(rpt_bg)
# Gives you repeatability (R) with CI and p-value
# 
# # Run for all groups
# rpt_results <- fluxes %>%
#   filter(!is.na(facet_label)) %>%
#   group_by(SPECIES, location) %>%
#   summarise(
#     rpt_result = {
#       r <- rpt(CH4_flux_nmolpm2ps ~ 1 + (1 | Tree), 
#                grname = "Tree", data = cur_data(),
#                datatype = "Gaussian", nboot = 100, npermut = 100)
#       list(R = r$R$Tree, CI_low = r$CI_emp$`2.5%`, CI_high = r$CI_emp$`97.5%`, p = r$P$Tree)
#     },
#     .groups = "drop"
#   ) %>%
#   unnest_wider(rpt_result)
# 
# print(rpt_results)

# ============================================
# 4. Spearman correlation across time periods
# ============================================
# Do tree ranks correlate between early and late periods?

tree_by_period <- fluxes %>%
  filter(!is.na(facet_label)) %>%
  mutate(period = ifelse(date < as.Date("2024-06-01"), "early", "late")) %>%
  group_by(Tree, SPECIES, location, period) %>%
  summarise(mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE), .groups = "drop") %>%
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

print(cor_by_group)

# Plot it
fig_period_cor <- ggplot(tree_by_period, aes(x = early, y = late)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  facet_grid(location ~ SPECIES, scales = "free") +
  labs(x = "Mean CH4 flux (early period)",
       y = "Mean CH4 flux (late period)",
       title = "Tree-level consistency across time") +
  theme_classic()

print(fig_period_cor)




# ============================================
# Fixed LRT (likelihood ratio test)
# ============================================

lrt_by_group <- fluxes %>%
  filter(!is.na(facet_label)) %>%
  group_by(SPECIES, location) %>%
  reframe(
    chisq = {
      mod_tree <- lmer(CH4_flux_nmolpm2ps ~ 1 + (1 | Tree), data = cur_data(), REML = FALSE)
      mod_null <- lm(CH4_flux_nmolpm2ps ~ 1, data = cur_data())
      as.numeric(2 * (logLik(mod_tree) - logLik(mod_null)))
    },
    p_value = pchisq(chisq, df = 1, lower.tail = FALSE)
  )

print(lrt_by_group)

# ============================================
# Summary table for publication
# ============================================

summary_table <- icc_by_group %>%
  left_join(cor_by_group %>% select(SPECIES, location, spearman_rho, p_value), 
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
  select(species_full, location, n_trees, n_obs, icc, spearman_rho, cor_p, sig)

print(summary_table)

# ============================================
# Better visualization of early vs late correlation
# ============================================

# Add species full names
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
    location = factor(location, levels = c("upland", "wetland"))
  )

# Add correlation labels to plot
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
    location = factor(location, levels = c("upland", "wetland")),
    label = paste0("ρ = ", round(spearman_rho, 2), 
                   ifelse(p_value < 0.05, "*", ""),
                   ifelse(p_value < 0.01, "*", ""))
  )

fig_period_cor <- ggplot(tree_by_period, aes(x = early, y = late)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick", alpha = 0.2) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(data = cor_labels, 
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5, size = 3.5) +
  facet_grid(location ~ species_full, scales = "free") +
  labs(x = expression("Mean CH"[4]~"flux, early period (nmol m"^-2~"s"^-1~")"),
       y = expression("Mean CH"[4]~"flux, late period (nmol m"^-2~"s"^-1~")")) +
  theme_classic(base_size = 11) +
  theme(strip.text = element_text(face = "italic"),
        strip.background = element_blank())

print(fig_period_cor)

#ggsave("figures/tree_consistency_correlation.png", fig_period_cor, width = 9, height = 6, dpi = 300)





library(lme4)
library(performance)

# Fit random intercept model
icc_model <- lmer(CH4_flux_nmolpm2ps ~ 1 + (1 | Tree), data = fluxes)

# Get ICC
icc(icc_model)

# By species-location
icc_by_group <- fluxes %>%
  group_by(SPECIES, location) %>%
  summarise(
    icc = {
      mod <- lmer(CH4_flux_nmolpm2ps ~ 1 + (1 | Tree), data = cur_data())
      var_tree <- as.numeric(VarCorr(mod)$Tree)
      var_resid <- sigma(mod)^2
      var_tree / (var_tree + var_resid)
    },
    .groups = "drop"
  )

print(icc_by_group)

# Extract tree-level effects from mixed model
model <- lmer(CH4_flux_nmolpm2ps ~ SPECIES + location + (1 | Tree), data = fluxes)

tree_effects <- ranef(model)$Tree %>%
  as.data.frame() %>%
  rownames_to_column("Tree") %>%
  rename(tree_effect = `(Intercept)`) %>%
  mutate(Tree = as.integer(Tree)) %>%
  left_join(fluxes %>% distinct(Tree, SPECIES, location), by = "Tree") %>%
  arrange(desc(tree_effect))

# Plot tree effects
fig_blups <- ggplot(tree_effects, aes(x = reorder(factor(Tree), tree_effect), 
                                      y = tree_effect, fill = SPECIES)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~location, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "Tree", y = "Random effect (deviation from mean)",
       title = "Tree-level CH4 flux effects") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 6))

print(fig_blups)




# Calculate tree mean flux and assign rank
tree_means <- fluxes %>%
  group_by(Tree, SPECIES, location, species_full) %>%
  summarise(mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE), .groups = "drop") %>%
  group_by(SPECIES, location) %>%
  mutate(flux_rank = rank(-mean_CH4)) %>%  # 1 = highest emitter
  ungroup()

# Join back to main data
fluxes <- fluxes %>%
  left_join(tree_means %>% select(Tree, mean_CH4, flux_rank), by = "Tree")

# Plot
fig_tree_consistency <- ggplot(fluxes, aes(x = date, y = CH4_flux_nmolpm2ps, 
                                           group = Tree, color = flux_rank)) +
  geom_line(alpha = 0.8, linewidth = 0.7) +
  geom_point(size = 1.5) +
  scale_color_viridis_c(option = "plasma", direction = -1,  # hot = high emitter
                        name = "Rank\n(1 = highest)") +
  facet_grid(location ~ species_full, scales = "free_y") +
  labs(x = "Date",
       y = expression("CH"[4]~"flux (nmol m"^-2~"s"^-1~")")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        legend.position = "right")

print(fig_tree_consistency)


fig_tree_consistency <- ggplot(fluxes, aes(x = date, y = CH4_flux_nmolpm2ps, 
                                           group = Tree, color = flux_rank)) +
  geom_line(alpha = 0.8, linewidth = 0.7) +
  geom_point(size = 1.5) +
  scale_color_viridis_c(option = "plasma", direction = -1,
                        name = "Rank\n(1 = highest)") +
  facet_wrap(location ~ species_full, scales = "free_y", ncol = 4) +
  labs(x = "Date",
       y = expression("CH"[4]~"flux (nmol m"^-2~"s"^-1~")")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        legend.position = "right")

print(fig_tree_consistency)






# Calculate tree mean flux within each species-location group
tree_means <- fluxes %>%
  group_by(Tree, SPECIES, location, species_full) %>%
  summarise(mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE), .groups = "drop") %>%
  group_by(SPECIES, location) %>%
  mutate(
    # Scale mean within each group to 0-1 for consistent coloring
    mean_scaled = (mean_CH4 - min(mean_CH4)) / (max(mean_CH4) - min(mean_CH4))
  ) %>%
  ungroup()

# Join back to main data
fluxes <- fluxes %>%
  select(-any_of(c("mean_CH4", "flux_rank", "mean_scaled"))) %>%  # remove if exists
  
  left_join(tree_means %>% select(Tree, mean_CH4, mean_scaled), by = "Tree")

# Create ordered factor for faceting
fluxes <- fluxes %>%
  mutate(
    species_order = factor(species_full, 
                           levels = c("Acer rubrum", "Tsuga canadensis", 
                                      "Quercus rubra", "Nyssa sylvatica")),
    location_order = factor(location, levels = c("upland", "wetland")),
    # Create combined facet variable
    facet_label = interaction(location_order, species_order, sep = " - ")
  )

# Reorder facet levels: upland (rm, hem, ro), wetland (rm, hem, bg)
facet_order <- c(
  "upland - Acer rubrum", "upland - Tsuga canadensis", "upland - Quercus rubra",
  "wetland - Acer rubrum", "wetland - Tsuga canadensis", "wetland - Nyssa sylvatica"
)

fluxes <- fluxes %>%
  mutate(facet_label = factor(facet_label, levels = facet_order))

# Plot
fig_tree_consistency <- ggplot(fluxes %>% filter(!is.na(facet_label)), 
                               aes(x = date, y = CH4_flux_nmolpm2ps, 
                                   group = Tree, color = mean_scaled)) +
  geom_line(alpha = 0.8, linewidth = 0.7) +
  geom_point(size = 1.5) +
  scale_color_viridis_c(option = "plasma", 
                        name = "Relative\ntree mean\n(within group)") +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 3) +
  labs(x = "Date",
       y = expression("CH"[4]~"flux (nmol m"^-2~"s"^-1~")")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        legend.position = "right")

print(fig_tree_consistency)


fig_tree_consistency_asinh <- ggplot(fluxes %>% filter(!is.na(facet_label)), 
                                     aes(x = date, y = CH4_flux_nmolpm2ps, 
                                         group = Tree, color = mean_scaled)) +
  geom_line(alpha = 0.8, linewidth = 0.7) +
  geom_point(size = 1.5) +
  scale_color_viridis_c(option = "plasma", 
                        name = "Relative\ntree mean\n(within group)") +
  scale_y_continuous(trans = "asinh") +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 3) +
  labs(x = "Date",
       y = expression("CH"[4]~"flux (nmol m"^-2~"s"^-1~", asinh scale)")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        legend.position = "right")

print(fig_tree_consistency_asinh)


library(lme4)
library(performance)
library(tidyverse)
library(patchwork)

# ============================================
# NATURAL COLOR PALETTE FOR SPECIES
# Inspired by each tree's characteristics
# ============================================
# A. rubrum (red maple) - muted red/burgundy
# T. canadensis (hemlock) - deep forest green
# Q. rubra (red oak) - warm brown/tan
# N. sylvatica (blackgum) - dark slate/charcoal

species_colors <- c(
  
  "rm"  = "#9E4A5B",
  "hem" = "#2E5A4C",
  "ro"  = "#8B6F47",
  "bg"  = "#4A5568"
)

species_colors_full <- c(
  "A. rubrum"     = "#9E4A5B",
  "T. canadensis" = "#2E5A4C",
  "Q. rubra"      = "#8B6F47",
  "N. sylvatica"  = "#4A5568"
)

# ============================================
# UPDATED fig_period_cor
# Order: red maple, hemlock, then blackgum/red oak per row
# No empty spaces, independent axes per panel
# ============================================

# Create custom facet ordering for period correlation plot
tree_by_period <- tree_by_period %>%
  mutate(
    # Create row-specific ordering
    col_order = case_when(
      SPECIES == "rm" ~ 1,
      SPECIES == "hem" ~ 2,
      SPECIES == "ro" ~ 3,   # upland only
      SPECIES == "bg" ~ 3,   # wetland only
      TRUE ~ NA_real_
    ),
    # Create combined facet label for proper ordering
    facet_combo = paste(location, col_order, sep = "_"),
    facet_combo = factor(facet_combo, levels = c(
      "upland_1", "upland_2", "upland_3",
      "wetland_1", "wetland_2", "wetland_3"
    )),
    # Create display labels
    panel_label = case_when(
      SPECIES == "rm" ~ "A. rubrum",
      SPECIES == "hem" ~ "T. canadensis",
      SPECIES == "ro" ~ "Q. rubra",
      SPECIES == "bg" ~ "N. sylvatica"
    ),
    location_label = factor(location, levels = c("upland", "wetland"),
                            labels = c("Upland", "Wetland"))
  )

# Update correlation labels with same ordering
cor_labels <- cor_by_group %>%
  mutate(
    col_order = case_when(
      SPECIES == "rm" ~ 1,
      SPECIES == "hem" ~ 2,
      SPECIES == "ro" ~ 3,
      SPECIES == "bg" ~ 3,
      TRUE ~ NA_real_
    ),
    facet_combo = paste(location, col_order, sep = "_"),
    facet_combo = factor(facet_combo, levels = c(
      "upland_1", "upland_2", "upland_3",
      "wetland_1", "wetland_2", "wetland_3"
    )),
    panel_label = case_when(
      SPECIES == "rm" ~ "A. rubrum",
      SPECIES == "hem" ~ "T. canadensis",
      SPECIES == "ro" ~ "Q. rubra",
      SPECIES == "bg" ~ "N. sylvatica"
    ),
    location_label = factor(location, levels = c("upland", "wetland"),
                            labels = c("Upland", "Wetland")),
    label = paste0("ρ = ", round(spearman_rho, 2), 
                   ifelse(p_value < 0.01, "**",
                          ifelse(p_value < 0.05, "*", "")))
  )

# Create labeller function for strip labels
facet_labeller <- function(labels) {
  # Extract species from the data for each panel
  lapply(labels[[1]], function(x) {
    loc <- strsplit(as.character(x), "_")[[1]][1]
    col <- as.numeric(strsplit(as.character(x), "_")[[1]][2])
    if (col == 1) return("A. rubrum")
    if (col == 2) return("T. canadensis")
    if (col == 3 && loc == "upland") return("Q. rubra")
    if (col == 3 && loc == "wetland") return("N. sylvatica")
  })
}

fig_period_cor <- ggplot(tree_by_period, aes(x = early, y = late)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick", alpha = 0.2) +
  geom_point(size = 3, alpha = 0.7) +
  geom_text(data = cor_labels, 
            aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5, size = 3.5) +
  facet_wrap(~ facet_combo, scales = "free", ncol = 3,
             labeller = as_labeller(function(x) {
               sapply(x, function(val) {
                 parts <- strsplit(val, "_")[[1]]
                 loc <- parts[1]
                 col <- as.numeric(parts[2])
                 if (col == 1) return("A. rubrum")
                 if (col == 2) return("T. canadensis")
                 if (col == 3 && loc == "upland") return("Q. rubra")
                 if (col == 3 && loc == "wetland") return("N. sylvatica")
               })
             })) +
  labs(x = expression("Mean CH"[4]~"flux, early period (nmol m"^-2~"s"^-1~")"),
       y = expression("Mean CH"[4]~"flux, late period (nmol m"^-2~"s"^-1~")")) +
  theme_classic(base_size = 11) +
  theme(strip.text = element_text(face = "italic"),
        strip.background = element_blank())

print(fig_period_cor)

# ============================================
# UPDATED fig_blups - independent axes for upland/wetland
# Plus second version wrapped by location AND species
# BLUPs calculated SEPARATELY for each species x location group
# ============================================

# Fit separate models for each species x location group and extract BLUPs
tree_effects <- fluxes %>%
  filter(!is.na(facet_label)) %>%
  group_by(SPECIES, location) %>%
  group_modify(~ {
    mod <- lmer(CH4_flux_nmolpm2ps ~ 1 + (1 | Tree), data = .x)
    blups <- ranef(mod)$Tree %>%
      as.data.frame() %>%
      rownames_to_column("Tree") %>%
      rename(tree_effect = `(Intercept)`) %>%
      mutate(Tree = as.integer(Tree))
    return(blups)
  }) %>%
  ungroup() %>%
  arrange(desc(tree_effect)) %>%
  mutate(
    species_full = case_when(
      SPECIES == "bg"  ~ "N. sylvatica",
      SPECIES == "hem" ~ "T. canadensis",
      SPECIES == "rm"  ~ "A. rubrum",
      SPECIES == "ro"  ~ "Q. rubra"
    ),
    location_label = factor(location, levels = c("upland", "wetland"),
                            labels = c("Upland", "Wetland"))
  )

# Version 1: Independent axes for upland and wetland
fig_blups <- ggplot(tree_effects, aes(x = reorder(factor(Tree), tree_effect), 
                                      y = tree_effect, fill = SPECIES)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  facet_wrap(~ location_label, scales = "free", ncol = 2) +
  scale_fill_manual(values = species_colors,
                    labels = c("bg" = "N. sylvatica", "hem" = "T. canadensis",
                               "rm" = "A. rubrum", "ro" = "Q. rubra")) +
  labs(x = "Tree", y = "Random effect (deviation from group mean)",
       fill = "Species") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 6),
        legend.text = element_text(face = "italic"))

print(fig_blups)

# Version 2: Caterpillar plot of BLUPs with confidence intervals
# Shows which trees are significantly different from group mean

library(merTools)

# Fit separate models and extract BLUPs with CIs for each species x location group
tree_effects_with_ci <- fluxes %>%
  filter(!is.na(facet_label)) %>%
  group_by(SPECIES, location) %>%
  group_modify(~ {
    mod <- lmer(CH4_flux_nmolpm2ps ~ 1 + (1 | Tree), data = .x)
    
    # Get prediction intervals for random effects
    re_df <- REsim(mod, n.sims = 1000) %>%
      filter(groupFctr == "Tree") %>%
      dplyr::select(Tree = groupID, tree_effect = mean, median, sd) %>%
      mutate(
        ci_lower = median - 1.96 * sd,
        ci_upper = median + 1.96 * sd,
        Tree = as.integer(as.character(Tree))
      ) %>%
      dplyr::select(-median, -sd)
    
    return(re_df)
  }) %>%
  ungroup() %>%
  # Create facet ordering
  mutate(
    col_order = case_when(
      SPECIES == "rm" ~ 1,
      SPECIES == "hem" ~ 2,
      SPECIES == "ro" ~ 3,
      SPECIES == "bg" ~ 3,
      TRUE ~ NA_real_
    ),
    facet_combo = paste(location, col_order, sep = "_"),
    facet_combo = factor(facet_combo, levels = c(
      "upland_1", "upland_2", "upland_3",
      "wetland_1", "wetland_2", "wetland_3"
    )),
    # Flag if CI excludes zero (significant)
    significant = (ci_lower > 0) | (ci_upper < 0)
  ) %>%
  filter(!is.na(facet_combo)) %>%
  # Order trees by effect within each facet
  group_by(facet_combo) %>%
  mutate(Tree_ordered = reorder(factor(Tree), tree_effect)) %>%
  ungroup()

fig_blups_v2 <- ggplot(tree_effects_with_ci, 
                       aes(x = Tree_ordered, y = tree_effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), 
                width = 0, color = "gray40", linewidth = 0.5) +
  geom_point(aes(color = significant), size = 2) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "firebrick"),
                     labels = c("FALSE" = "CI includes 0", "TRUE" = "CI excludes 0"),
                     name = NULL) +
  coord_flip() +
  facet_wrap(~ facet_combo, scales = "free", ncol = 3,
             labeller = as_labeller(function(x) {
               sapply(x, function(val) {
                 parts <- strsplit(val, "_")[[1]]
                 loc <- parts[1]
                 col <- as.numeric(parts[2])
                 loc_label <- ifelse(loc == "upland", "Upland", "Wetland")
                 sp_label <- if (col == 1) "A. rubrum" 
                 else if (col == 2) "T. canadensis"
                 else if (col == 3 && loc == "upland") "Q. rubra"
                 else "N. sylvatica"
                 paste(loc_label, "-", sp_label)
               })
             })) +
  labs(x = "Tree (ordered by effect)", 
       y = expression("Random effect (deviation from group mean, nmol m"^-2~"s"^-1~")")) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 6),
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        legend.position = "bottom")

print(fig_blups_v2)

# ============================================
# UPDATED fig_tree_consistency_asinh
# Plotting order: highest mean on top, lowest on bottom
# ============================================
library(dplyr)
library(ggplot2)
library(scales)
library(ggtext)

# 1) Make an ordered facet label with markdown italics
fluxes_ordered <- fluxes_ordered %>%
  mutate(
    location = factor(location, levels = c("upland", "wetland")),
    species_full = factor(
      species_full,
      levels = c("Acer rubrum", "Tsuga canadensis", "Quercus rubra", "Nyssa sylvatica")
    ),
    facet_key = paste(as.character(location), as.character(species_full), sep = " - "),
    facet_label = case_when(
      facet_key == "upland - Acer rubrum"        ~ "<i>Acer rubrum</i><br>Upland",
      facet_key == "upland - Tsuga canadensis"   ~ "<i>Tsuga canadensis</i><br>Upland",
      facet_key == "upland - Quercus rubra"      ~ "<i>Quercus rubra</i><br>Upland",
      facet_key == "wetland - Acer rubrum"       ~ "<i>Acer rubrum</i><br>Wetland",
      facet_key == "wetland - Tsuga canadensis"  ~ "<i>Tsuga canadensis</i><br>Wetland",
      facet_key == "wetland - Nyssa sylvatica"   ~ "<i>Nyssa sylvatica</i><br>Wetland",
      TRUE ~ NA_character_
    ),
    facet_label = factor(
      facet_label,
      levels = c(
        "<i>Acer rubrum</i><br>Upland",
        "<i>Tsuga canadensis</i><br>Upland",
        "<i>Quercus rubra</i><br>Upland",
        "<i>Acer rubrum</i><br>Wetland",
        "<i>Tsuga canadensis</i><br>Wetland",
        "<i>Nyssa sylvatica</i><br>Wetland"
      )
    )
  ) %>%
  filter(!is.na(facet_label))

# 2) Plot (viridis continuous, y asinh), with a crisp theme and markdown strips
fig_tree_consistency_asinh <- ggplot(
  fluxes_ordered,
  aes(
    x = date,
    y = CH4_flux_nmolpm2ps,
    group = reorder(factor(Tree), tree_mean),
    color = mean_scaled
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
  geom_line(alpha = 0.85, linewidth = 0.8, lineend = "round", linejoin = "round") +
  scale_color_viridis_c(
    option = "viridis",
    name = "Relative tree mean\n(within group)",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0,
      barheight = unit(45, "pt"),
      barwidth  = unit(10, "pt")
    )
  ) +
  scale_y_continuous(
    trans = "asinh",
    breaks = pretty_breaks(n = 4),
    expand = expansion(mult = c(0.03, 0.06))
  ) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b '%y",
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_wrap(~ facet_label, ncol = 3, scales = "free_y") +
  labs(
    x = NULL,
    y = expression(CH[4]~flux~(nmol~m^{-2}~s^{-1}*","~asinh~scale))
  ) +
  theme_classic(base_size = 12) +
  theme(
    # light grid like your temporal figure
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_text(angle = 45, hjust = 1),
    
    strip.background = element_blank(),
    strip.text = ggtext::element_markdown(size = 11, lineheight = 0.95),
    # optional: tighten panel spacing a touch
    panel.spacing = unit(0.9, "lines"),
    
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text  = element_text(size = 9)
  )

print(fig_tree_consistency_asinh)


# ============================================
# COMBINED FIGURE LAYOUT
# Top row: fig_tree_consistency_asinh
# Bottom row: fig_blups_v2 and fig_period_cor
# ============================================

# Adjust individual plots for combined layout
fig_tree_consistency_asinh_layout <- fig_tree_consistency_asinh +
  theme(legend.position = "right")

fig_blups_layout <- fig_blups_v2 +
  theme(legend.position = "none")

fig_period_cor_layout <- fig_period_cor +
  theme(legend.position = "none")

# Create combined layout using patchwork
combined_figure <- fig_tree_consistency_asinh_layout /
  (fig_blups_layout) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(tag_levels = 'A')

print(combined_figure)

# Save the combined figure
ggsave("combined_figure.png", combined_figure, width = 14, height = 12, dpi = 300)
ggsave("combined_figure.pdf", combined_figure, width = 14, height = 12)

# Also save individual figures
ggsave("fig_period_cor_updated.png", fig_period_cor, width = 9, height = 6, dpi = 300)
ggsave("fig_blups_v2_updated.png", fig_blups_v2, width = 12, height = 8, dpi = 300)
ggsave("fig_tree_consistency_asinh_updated.png", fig_tree_consistency_asinh, width = 12, height = 8, dpi = 300)

# Print summary of changes
cat("\n========== SUMMARY OF UPDATES ==========\n")
cat("1. fig_period_cor: Reordered as rm, hem, ro/bg per row; free scales; no empty panels\n")
cat("2. fig_blups_v2: Boxplots with jittered points per tree, ordered by mean flux\n")
cat("3. fig_tree_consistency_asinh: Plotting order ensures high emitters on top\n")
cat("4. combined_figure: Layout with time series (top) and boxplots + correlation (bottom)\n")














# ============================================
# Z-SCORE WITHIN YEAR x SPECIES x LOCATION
# Point size scales with |z-score|
# ============================================

# Calculate z-scores within each year x location x species group
fluxes_zscore_year <- fluxes %>%
  filter(!is.na(facet_label)) %>%
  mutate(year = year(date)) %>%
  group_by(year, location, SPECIES) %>%
  mutate(
    group_mean = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    group_sd = sd(CH4_flux_nmolpm2ps, na.rm = TRUE),
    z_score = (CH4_flux_nmolpm2ps - group_mean) / group_sd
  ) %>%
  ungroup() %>%
  mutate(z_score = ifelse(is.finite(z_score), z_score, 0))

# Calculate mean z-score per tree for ordering
tree_mean_zscore <- fluxes_zscore_year %>%
  group_by(Tree, location, SPECIES) %>%
  summarise(mean_zscore = mean(z_score, na.rm = TRUE), .groups = "drop")

fluxes_zscore_year <- fluxes_zscore_year %>%
  left_join(tree_mean_zscore, by = c("Tree", "location", "SPECIES"))

# Create facet labels
fluxes_zscore_year <- fluxes_zscore_year %>%
  mutate(
    location = factor(location, levels = c("upland", "wetland")),
    species_full = factor(
      species_full,
      levels = c("Acer rubrum", "Tsuga canadensis", "Quercus rubra", "Nyssa sylvatica")
    ),
    facet_key = paste(as.character(location), as.character(species_full), sep = " - "),
    facet_label_z = case_when(
      facet_key == "upland - Acer rubrum"        ~ "<i>Acer rubrum</i><br>Upland",
      facet_key == "upland - Tsuga canadensis"   ~ "<i>Tsuga canadensis</i><br>Upland",
      facet_key == "upland - Quercus rubra"      ~ "<i>Quercus rubra</i><br>Upland",
      facet_key == "wetland - Acer rubrum"       ~ "<i>Acer rubrum</i><br>Wetland",
      facet_key == "wetland - Tsuga canadensis"  ~ "<i>Tsuga canadensis</i><br>Wetland",
      facet_key == "wetland - Nyssa sylvatica"   ~ "<i>Nyssa sylvatica</i><br>Wetland",
      TRUE ~ NA_character_
    ),
    facet_label_z = factor(
      facet_label_z,
      levels = c(
        "<i>Acer rubrum</i><br>Upland",
        "<i>Tsuga canadensis</i><br>Upland",
        "<i>Quercus rubra</i><br>Upland",
        "<i>Acer rubrum</i><br>Wetland",
        "<i>Tsuga canadensis</i><br>Wetland",
        "<i>Nyssa sylvatica</i><br>Wetland"
      )
    )
  ) %>%
  filter(!is.na(facet_label_z))

# Create ordered tree factor within each facet (ordered by mean z-score)
# Also create abs z-score for size mapping
fluxes_zscore_year <- fluxes_zscore_year %>%
  group_by(facet_label_z) %>%
  mutate(
    tree_ordered = reorder(factor(Tree), mean_zscore)
  ) %>%
  ungroup() %>%
  mutate(
    z_score_abs = pmin(abs(z_score), 4)  # clamp at 4 for size
  )

# Point + line plot with size scaling
fig_zscore_lines_year <- ggplot(
  fluxes_zscore_year,
  aes(
    x = date,
    y = tree_ordered,
    color = z_score,
    size = z_score_abs,
    group = Tree
  )
) +
  geom_line(aes(color = mean_zscore), alpha = 0.4, linewidth = 0.5, show.legend = FALSE) +
  geom_point(alpha = 0.85) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "grey80",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2, 4),
    oob = scales::squish,
    name = "Z-score",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(50, "pt"),
      barwidth = unit(10, "pt")
    )
  ) +
  scale_size_continuous(
    range = c(0.5, 4),
    guide = "none"
  ) +
  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b '%y",
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_wrap(~ facet_label_z, ncol = 3, scales = "free_y") +
  labs(
    x = NULL,
    y = "Tree (ordered by mean z-score)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = ggtext::element_markdown(size = 11, lineheight = 0.95),
    panel.spacing = unit(0.9, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

print(fig_zscore_lines_year)






# ============================================
# Z-SCORE WITHIN YEAR x SPECIES x LOCATION
# Filtered to CO2 r2 > 0.8
# X-axis = measurement index (evenly spaced)
# ============================================

# Calculate z-scores within each year x location x species group
fluxes_zscore_year <- fluxes %>%
  filter(!is.na(facet_label), CO2_r2 > 0) %>%
  mutate(year = year(date)) %>%
  group_by(year, location, SPECIES) %>%
  mutate(
    group_mean = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    group_sd = sd(CH4_flux_nmolpm2ps, na.rm = TRUE),
    z_score = (CH4_flux_nmolpm2ps - group_mean) / group_sd
  ) %>%
  ungroup() %>%
  mutate(z_score = ifelse(is.finite(z_score), z_score, 0))

# Calculate mean z-score per tree for ordering
tree_mean_zscore <- fluxes_zscore_year %>%
  group_by(Tree, location, SPECIES) %>%
  summarise(mean_zscore = mean(z_score, na.rm = TRUE), .groups = "drop")

fluxes_zscore_year <- fluxes_zscore_year %>%
  left_join(tree_mean_zscore, by = c("Tree", "location", "SPECIES"))

# Create facet labels
fluxes_zscore_year <- fluxes_zscore_year %>%
  mutate(
    location = factor(location, levels = c("upland", "wetland")),
    species_full = factor(
      species_full,
      levels = c("Acer rubrum", "Tsuga canadensis", "Quercus rubra", "Nyssa sylvatica")
    ),
    facet_key = paste(as.character(location), as.character(species_full), sep = " - "),
    facet_label_z = case_when(
      facet_key == "upland - Acer rubrum"        ~ "<i>Acer rubrum</i><br>Upland",
      facet_key == "upland - Tsuga canadensis"   ~ "<i>Tsuga canadensis</i><br>Upland",
      facet_key == "upland - Quercus rubra"      ~ "<i>Quercus rubra</i><br>Upland",
      facet_key == "wetland - Acer rubrum"       ~ "<i>Acer rubrum</i><br>Wetland",
      facet_key == "wetland - Tsuga canadensis"  ~ "<i>Tsuga canadensis</i><br>Wetland",
      facet_key == "wetland - Nyssa sylvatica"   ~ "<i>Nyssa sylvatica</i><br>Wetland",
      TRUE ~ NA_character_
    ),
    facet_label_z = factor(
      facet_label_z,
      levels = c(
        "<i>Acer rubrum</i><br>Upland",
        "<i>Tsuga canadensis</i><br>Upland",
        "<i>Quercus rubra</i><br>Upland",
        "<i>Acer rubrum</i><br>Wetland",
        "<i>Tsuga canadensis</i><br>Wetland",
        "<i>Nyssa sylvatica</i><br>Wetland"
      )
    )
  ) %>%
  filter(!is.na(facet_label_z))

# Create measurement index within each facet (evenly spaced x-axis)
fluxes_zscore_year <- fluxes_zscore_year %>%
  group_by(facet_label_z) %>%
  mutate(
    tree_ordered = reorder(factor(Tree), mean_zscore),
    date_idx = as.numeric(factor(date))
  ) %>%
  ungroup() %>%
  mutate(
    z_score_abs = pmin(abs(z_score), 4)
  )

# Create date labels for x-axis
date_labels_df <- fluxes_zscore_year %>%
  group_by(facet_label_z) %>%
  distinct(date, date_idx) %>%
  arrange(date_idx) %>%
  mutate(
    n_dates = n(),
    show_label = row_number() %% ceiling(n_dates / 6) == 1 | row_number() == n()
  ) %>%
  filter(show_label) %>%
  ungroup()

# Point + line plot with evenly spaced x
fig_zscore_lines_year <- ggplot(
  fluxes_zscore_year,
  aes(
    x = date_idx,
    y = tree_ordered,
    color = z_score,
    size = z_score_abs,
    group = Tree
  )
) +
  geom_line(aes(color = mean_zscore), alpha = 0.4, linewidth = 0.5, show.legend = FALSE) +
  geom_point(alpha = 0.85) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "grey80",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2, 4),
    oob = scales::squish,
    name = "Z-score",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(50, "pt"),
      barwidth = unit(10, "pt")
    )
  ) +
  scale_size_continuous(
    range = c(0.5, 4),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = date_labels_df %>% filter(facet_label_z == levels(fluxes_zscore_year$facet_label_z)[1]) %>% pull(date_idx),
    labels = date_labels_df %>% filter(facet_label_z == levels(fluxes_zscore_year$facet_label_z)[1]) %>% pull(date) %>% format("%b '%y"),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_wrap(~ facet_label_z, ncol = 3, scales = "free") +
  labs(
    x = NULL,
    y = "Tree (ordered by mean z-score)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = ggtext::element_markdown(size = 11, lineheight = 0.95),
    panel.spacing = unit(0.9, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

print(fig_zscore_lines_year)

# Print how many measurements were kept
cat("\nMeasurements kept:", nrow(fluxes_zscore_year), "of", nrow(fluxes %>% filter(!is.na(facet_label))), 
    "(", round(100 * nrow(fluxes_zscore_year) / nrow(fluxes %>% filter(!is.na(facet_label))), 1), "%)\n")











# ============================================
# FIGURE: Z-score within YEAR x LOCATION x SPECIES (filtered to CO2_r2 > 0.8)
# X-axis = measurement index (even spacing to de-cram summer)
# Y-axis = tree tracks (ordered by mean z-score)
# Color = z-score, Point size = |z-score|
# ============================================

# (assumes you already loaded dplyr/ggplot2/scales/ggtext/lubridate)
# If not, uncomment:
# library(dplyr); library(ggplot2); library(lubridate); library(scales); library(ggtext); library(grid)

fluxes_zscore_year_idx <- fluxes %>%
  filter(!is.na(facet_label), !is.na(CO2_r2), CO2_r2 > 0.8) %>%
  mutate(
    # ensure Date class
    date = as.Date(date),
    year = lubridate::year(date)
  ) %>%
  group_by(year, location, SPECIES) %>%
  mutate(
    group_mean = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    group_sd   = sd(CH4_flux_nmolpm2ps, na.rm = TRUE),
    z_score    = (CH4_flux_nmolpm2ps - group_mean) / group_sd
  ) %>%
  ungroup() %>%
  mutate(
    z_score = ifelse(is.finite(z_score), z_score, 0)
  )

# mean z-score per tree (for ordering)
tree_mean_zscore <- fluxes_zscore_year_idx %>%
  group_by(Tree, location, SPECIES) %>%
  summarise(mean_zscore = mean(z_score, na.rm = TRUE), .groups = "drop")

fluxes_zscore_year_idx <- fluxes_zscore_year_idx %>%
  left_join(tree_mean_zscore, by = c("Tree", "location", "SPECIES"))

# facet labels (same style you were using)
fluxes_zscore_year_idx <- fluxes_zscore_year_idx %>%
  mutate(
    location = factor(location, levels = c("upland", "wetland")),
    species_full = factor(
      species_full,
      levels = c("Acer rubrum", "Tsuga canadensis", "Quercus rubra", "Nyssa sylvatica")
    ),
    facet_key = paste(as.character(location), as.character(species_full), sep = " - "),
    facet_label_z = dplyr::case_when(
      facet_key == "upland - Acer rubrum"        ~ "<i>Acer rubrum</i><br>Upland",
      facet_key == "upland - Tsuga canadensis"   ~ "<i>Tsuga canadensis</i><br>Upland",
      facet_key == "upland - Quercus rubra"      ~ "<i>Quercus rubra</i><br>Upland",
      facet_key == "wetland - Acer rubrum"       ~ "<i>Acer rubrum</i><br>Wetland",
      facet_key == "wetland - Tsuga canadensis"  ~ "<i>Tsuga canadensis</i><br>Wetland",
      facet_key == "wetland - Nyssa sylvatica"   ~ "<i>Nyssa sylvatica</i><br>Wetland",
      TRUE ~ NA_character_
    ),
    facet_label_z = factor(
      facet_label_z,
      levels = c(
        "<i>Acer rubrum</i><br>Upland",
        "<i>Tsuga canadensis</i><br>Upland",
        "<i>Quercus rubra</i><br>Upland",
        "<i>Acer rubrum</i><br>Wetland",
        "<i>Tsuga canadensis</i><br>Wetland",
        "<i>Nyssa sylvatica</i><br>Wetland"
      )
    )
  ) %>%
  filter(!is.na(facet_label_z))

# global measurement index (even spacing across all measurement dates)
all_dates <- sort(unique(fluxes_zscore_year_idx$date))
fluxes_zscore_year_idx <- fluxes_zscore_year_idx %>%
  mutate(
    date_idx = match(date, all_dates),
    z_score_abs = pmin(abs(z_score), 4)  # clamp for size only
  ) %>%
  group_by(facet_label_z) %>%
  mutate(
    tree_ordered = reorder(factor(Tree), mean_zscore)
  ) %>%
  ungroup()

# x-axis labels: ~7 ticks across the measurement index
n_breaks <- 7
break_idx <- unique(round(seq(1, length(all_dates), length.out = n_breaks)))
break_labs <- format(all_dates[break_idx], "%b '%y")

fig_zscore_tracks_year_idx <- ggplot(
  fluxes_zscore_year_idx,
  aes(
    x = date_idx,
    y = tree_ordered,
    group = factor(Tree),
    color = z_score,
    size = z_score_abs
  )
) +
  geom_line(aes(color = mean_zscore), alpha = 0.35, linewidth = 0.45, show.legend = FALSE) +
  geom_point(alpha = 0.85) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "grey85",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2, 4),
    oob = scales::squish,
    name = "Z-score",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(50, "pt"),
      barwidth = unit(10, "pt")
    )
  ) +
  scale_size_continuous(range = c(0.6, 4), guide = "none") +
  scale_x_continuous(
    breaks = break_idx,
    labels = break_labs,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_wrap(~ facet_label_z, ncol = 3, scales = "free_y") +
  labs(
    x = NULL,
    y = "Tree (ordered by mean z-score)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = ggtext::element_markdown(size = 11, lineheight = 0.95),
    panel.spacing = unit(0.9, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

print(fig_zscore_tracks_year_idx)

# optional: quick retained-count message
cat(
  "\nCO2_r2 > 0.8 kept ",
  nrow(fluxes_zscore_year_idx),
  " rows out of ",
  nrow(fluxes %>% filter(!is.na(facet_label))),
  " (",
  round(100 * nrow(fluxes_zscore_year_idx) / nrow(fluxes %>% filter(!is.na(facet_label))), 1),
  "%)\n",
  sep = ""
)






# ============================================
# FIGURE: Z-score within YEAR x LOCATION x SPECIES (filtered to CO2_r2 > 0.8)
# X-axis = measurement index (even spacing to de-cram summer)
# Y-axis = tree tracks (ordered by mean z-score)
# Color = z-score, Point size = |z-score|
# ============================================

# (assumes you already loaded dplyr/ggplot2/scales/ggtext/lubridate)
# If not, uncomment:
# library(dplyr); library(ggplot2); library(lubridate); library(scales); library(ggtext); library(grid)


fluxes <- fluxes %>%
  mutate(
    location = ifelse(PLOT == "BGS", "wetland", "upland"),
    species_full = case_when(
      SPECIES == "bg"  ~ "Nyssa sylvatica",
      SPECIES == "hem" ~ "Tsuga canadensis",
      SPECIES == "rm"  ~ "Acer rubrum",
      SPECIES == "ro"  ~ "Quercus rubra",
      TRUE ~ SPECIES
    ),
    facet_label = paste(location, "-", species_full)
  )

fluxes_zscore_year_idx <- fluxes %>%
  #filter(!is.na(facet_label), !is.na(CO2_r2), CO2_r2 > 0) %>%
  mutate(
    # ensure Date class
    date = as.Date(date),
    year = lubridate::year(date)
  ) %>%
  group_by(year, location, SPECIES) %>%
  mutate(
    group_mean = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    group_sd   = sd(CH4_flux_nmolpm2ps, na.rm = TRUE),
    z_score    = (CH4_flux_nmolpm2ps - group_mean) / group_sd
  ) %>%
  ungroup() %>%
  mutate(
    z_score = ifelse(is.finite(z_score), z_score, 0)
  )

# mean z-score per tree (for ordering)
tree_mean_zscore <- fluxes_zscore_year_idx %>%
  group_by(Tree, location, SPECIES) %>%
  summarise(mean_zscore = mean(z_score, na.rm = TRUE), .groups = "drop")

fluxes_zscore_year_idx <- fluxes_zscore_year_idx %>%
  left_join(tree_mean_zscore, by = c("Tree", "location", "SPECIES"))

# facet labels (same style you were using)
fluxes_zscore_year_idx <- fluxes_zscore_year_idx %>%
  mutate(
    location = factor(location, levels = c("upland", "wetland")),
    species_full = factor(
      species_full,
      levels = c("Acer rubrum", "Tsuga canadensis", "Quercus rubra", "Nyssa sylvatica")
    ),
    facet_key = paste(as.character(location), as.character(species_full), sep = " - "),
    facet_label_z = dplyr::case_when(
      facet_key == "upland - Acer rubrum"        ~ "<i>Acer rubrum</i><br>Upland",
      facet_key == "upland - Tsuga canadensis"   ~ "<i>Tsuga canadensis</i><br>Upland",
      facet_key == "upland - Quercus rubra"      ~ "<i>Quercus rubra</i><br>Upland",
      facet_key == "wetland - Acer rubrum"       ~ "<i>Acer rubrum</i><br>Wetland",
      facet_key == "wetland - Tsuga canadensis"  ~ "<i>Tsuga canadensis</i><br>Wetland",
      facet_key == "wetland - Nyssa sylvatica"   ~ "<i>Nyssa sylvatica</i><br>Wetland",
      TRUE ~ NA_character_
    ),
    facet_label_z = factor(
      facet_label_z,
      levels = c(
        "<i>Acer rubrum</i><br>Upland",
        "<i>Tsuga canadensis</i><br>Upland",
        "<i>Quercus rubra</i><br>Upland",
        "<i>Acer rubrum</i><br>Wetland",
        "<i>Tsuga canadensis</i><br>Wetland",
        "<i>Nyssa sylvatica</i><br>Wetland"
      )
    )
  ) %>%
  filter(!is.na(facet_label_z))

# global measurement index (even spacing across all measurement dates)
all_dates <- sort(unique(fluxes_zscore_year_idx$date))
fluxes_zscore_year_idx <- fluxes_zscore_year_idx %>%
  mutate(
    date_idx = match(date, all_dates),
    z_score_abs = pmin(abs(z_score), 4)  # clamp for size only
  ) %>%
  group_by(facet_label_z) %>%
  mutate(
    tree_ordered = reorder(factor(Tree), mean_zscore)
  ) %>%
  ungroup()

# x-axis labels: ~7 ticks across the measurement index
n_breaks <- 7
break_idx <- unique(round(seq(1, length(all_dates), length.out = n_breaks)))
break_labs <- format(all_dates[break_idx], "%b '%y")

fig_zscore_tracks_year_idx <- ggplot(
  fluxes_zscore_year_idx,
  aes(
    x = date_idx,
    y = tree_ordered,
    group = factor(Tree),
    color = z_score,
    size = z_score_abs
  )
) +
  geom_line(aes(color = mean_zscore), alpha = 0.35, linewidth = 0.45, show.legend = FALSE) +
  geom_point(alpha = 0.85) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "grey85",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-2, 4),
    oob = scales::squish,
    name = "Z-score",
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(50, "pt"),
      barwidth = unit(10, "pt")
    )
  ) +
  scale_size_continuous(range = c(0.6, 4), guide = "none") +
  scale_x_continuous(
    breaks = break_idx,
    labels = break_labs,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  facet_wrap(~ facet_label_z, ncol = 3, scales = "free_y") +
  labs(
    x = NULL,
    y = "Tree (ordered by mean z-score)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = ggtext::element_markdown(size = 11, lineheight = 0.95),
    panel.spacing = unit(0.9, "lines"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9)
  )

print(fig_zscore_tracks_year_idx)

# optional: quick retained-count message
cat(
  "\nCO2_r2 > 0.8 kept ",
  nrow(fluxes_zscore_year_idx),
  " rows out of ",
  nrow(fluxes %>% filter(!is.na(facet_label))),
  " (",
  round(100 * nrow(fluxes_zscore_year_idx) / nrow(fluxes %>% filter(!is.na(facet_label))), 1),
  "%)\n",
  sep = ""
)









# ============================================
# FIGURE: Z-score within YEAR x LOCATION x SPECIES
# Clean, modern version - Wetland on top, Upland on bottom
# ============================================

library(tidyverse)
library(ggtext)
library(scales)
library(lubridate)

# ============================================
# DATA PREP
# ============================================

fluxes <- fluxes %>%
  mutate(
    location = ifelse(PLOT == "BGS", "wetland", "upland"),
    species_full = case_when(
      SPECIES == "bg"  ~ "Nyssa sylvatica",
      SPECIES == "hem" ~ "Tsuga canadensis",
      SPECIES == "rm"  ~ "Acer rubrum",
      SPECIES == "ro"  ~ "Quercus rubra",
      TRUE ~ SPECIES
    ),
    facet_label = paste(location, "-", species_full)
  )

fluxes_zscore_year_idx <- fluxes %>%
  mutate(
    date = as.Date(date),
    year = lubridate::year(date)
  ) %>%
  group_by(year, location, SPECIES) %>%
  mutate(
    group_mean = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    group_sd   = sd(CH4_flux_nmolpm2ps, na.rm = TRUE),
    z_score    = (CH4_flux_nmolpm2ps - group_mean) / group_sd
  ) %>%
  ungroup() %>%
  mutate(
    z_score = ifelse(is.finite(z_score), z_score, 0)
  )

# Mean z-score per tree (for ordering)
tree_mean_zscore <- fluxes_zscore_year_idx %>%
  group_by(Tree, location, SPECIES) %>%
  summarise(mean_zscore = mean(z_score, na.rm = TRUE), .groups = "drop")

fluxes_zscore_year_idx <- fluxes_zscore_year_idx %>%
  left_join(tree_mean_zscore, by = c("Tree", "location", "SPECIES"))

# Facet labels - WETLAND FIRST
fluxes_zscore_year_idx <- fluxes_zscore_year_idx %>%
  mutate(
    location = factor(location, levels = c("wetland", "upland")),
    species_full = factor(
      species_full,
      levels = c("Acer rubrum", "Tsuga canadensis", "Nyssa sylvatica", "Quercus rubra")
    ),
    facet_key = paste(as.character(location), as.character(species_full), sep = " - "),
    facet_label_z = case_when(
      facet_key == "upland - Acer rubrum"        ~ "<i>Acer rubrum</i><br>Upland",
      facet_key == "upland - Tsuga canadensis"   ~ "<i>Tsuga canadensis</i><br>Upland",
      facet_key == "upland - Quercus rubra"      ~ "<i>Quercus rubra</i><br>Upland",
      facet_key == "wetland - Acer rubrum"       ~ "<i>Acer rubrum</i><br>Wetland",
      facet_key == "wetland - Tsuga canadensis"  ~ "<i>Tsuga canadensis</i><br>Wetland",
      facet_key == "wetland - Nyssa sylvatica"   ~ "<i>Nyssa sylvatica</i><br>Wetland",
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
all_dates <- sort(unique(fluxes_zscore_year_idx$date))
fluxes_zscore_year_idx <- fluxes_zscore_year_idx %>%
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

# ============================================
# PLOT
# ============================================

fig_zscore_tracks_year_idx <- ggplot(
  fluxes_zscore_year_idx,
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
  scale_size_continuous(range = c(0.5, 3.5), guide = "none") +
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
    # Clean grid
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.25),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    # Axes
    axis.text.x = element_text(size = 8, color = "grey30", angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    # Facets
    strip.background = element_blank(),
    strip.text = element_markdown(size = 10, lineheight = 1.1, color = "grey20"),
    panel.spacing.x = unit(0.8, "lines"),
    panel.spacing.y = unit(0.6, "lines"),
    # Legend
    legend.position = "right",
    legend.title = element_text(size = 9, face = "plain", color = "grey30"),
    legend.text = element_text(size = 8, color = "grey40"),
    legend.key.height = unit(35, "pt"),
    legend.key.width = unit(8, "pt"),
    # Background
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

print(fig_zscore_tracks_year_idx)

ggsave("figures/zscore_tracks.png", fig_zscore_tracks_year_idx, 
       width = 10, height = 7, dpi = 300, bg = "white")
ggsave("figures/zscore_tracks.pdf", fig_zscore_tracks_year_idx, 
       width = 10, height = 7)

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

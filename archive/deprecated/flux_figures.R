# ============================================
# CH4 Tree Flux Analysis - Harvard Forest
# ============================================

library(tidyverse)
library(lubridate)
library(lme4)
library(emmeans)

# --------------------------------------------
# 1. DATA LOADING AND CLEANING
# --------------------------------------------

fluxes <- read.csv('/Users/jongewirtzman/Google Drive/Research/tree-flux-2025/data/HF_2023-2025_tree_flux.csv')

fluxes <- fluxes %>%
  # Fix units for pre-2025 data
  
  mutate(CH4_flux_nmolpm2ps = ifelse(year < 2025, CH4_flux_nmolpm2ps * 1000, CH4_flux_nmolpm2ps)) %>%
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
    location = ifelse(PLOT == "BGS", "wetland", "upland"),
    month = month(date),
    doy = yday(date),
    species_full = case_when(
      SPECIES == "bg"  ~ "Nyssa sylvatica",
      SPECIES == "hem" ~ "Tsuga canadensis",
      SPECIES == "rm"  ~ "Acer rubrum",
      SPECIES == "ro"  ~ "Quercus rubra",
      TRUE ~ SPECIES
    )
  )

# --------------------------------------------
# 2. SUMMARY STATISTICS
# --------------------------------------------

# By species and location
flux_summary_species <- fluxes %>%
  group_by(SPECIES, location) %>%
  summarise(
    n_obs = n(),
    n_trees = n_distinct(Tree),
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    sd_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd_CH4 / sqrt(n_obs),
    median_CH4 = median(CH4_flux_nmolpm2ps, na.rm = TRUE),
    .groups = "drop"
  )

print(flux_summary_species)

# Tree-level means (for understanding individual variation)
tree_means <- fluxes %>%
  group_by(Tree, SPECIES, location) %>%
  summarise(
    n_obs = n(),
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    .groups = "drop"
  )

# --------------------------------------------
# 3. PLOTS
# --------------------------------------------

# Plot 1: Time series by species, faceted by location
p1 <- ggplot(fluxes, aes(x = date, y = CH4_flux_nmolpm2ps, color = SPECIES)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~location, ncol = 1, scales = "free_y") +
  labs(x = "Date", y = expression("CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       title = "CH4 flux over time by species and location") +
  theme_classic()

# Plot 2: Boxplot of tree means by species and location
p2 <- ggplot(tree_means, aes(x = SPECIES, y = mean_CH4, fill = location)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha = 0.5) +
  labs(x = "Species", y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       title = "Tree-level mean CH4 flux by species and location") +
  theme_classic()

# Plot 3: Distribution of fluxes
p3 <- ggplot(fluxes, aes(x = CH4_flux_nmolpm2ps, fill = location)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~SPECIES, scales = "free") +
  labs(x = expression("CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       title = "Distribution of CH4 flux by species") +
  theme_classic()

# Plot 4: CO2 vs CH4 relationship
p4 <- ggplot(fluxes, aes(x = CO2_flux_umolpm2ps, y = CH4_flux_nmolpm2ps, color = location)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~SPECIES) +
  labs(x = expression("CO"[2]~"flux (µmol m"^-2~"s"^-1~")"),
       y = expression("CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       title = "CO2 vs CH4 flux") +
  theme_classic()

print(p1)
print(p2)
print(p3)
print(p4)

# --------------------------------------------
# 4. STATISTICAL MODELS
# --------------------------------------------

# Mixed model: species × location with tree as random effect
model <- lmer(CH4_flux_nmolpm2ps ~ SPECIES * location + (1 | Tree), 
              data = fluxes, REML = TRUE)

summary(model)

# Marginal means
emm <- emmeans(model, ~ SPECIES | location)
print(emm)

# Pairwise comparisons
pairs(emm)

# Plot marginal means
emm_df <- as.data.frame(emm)

p5 <- ggplot(emm_df, aes(x = SPECIES, y = emmean, fill = location)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.8), width = 0.2) +
  labs(x = "Species", y = expression("Estimated CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       title = "Marginal mean CH4 flux by species and location") +
  theme_classic()

print(p5)


# --------------------------------------------
# MODELS REFLECTING ECOLOGICAL DESIGN
# --------------------------------------------

# 1. Wetland: specialist (bg) vs generalists (hem, rm)
model_wetland <- lmer(CH4_flux_nmolpm2ps ~ SPECIES + (1 | Tree), 
                      data = filter(fluxes, location == "wetland"))
summary(model_wetland)
emm_wetland <- emmeans(model_wetland, ~ SPECIES)
pairs(emm_wetland)

# 2. Upland: specialist (ro) vs generalists (hem, rm)
model_upland <- lmer(CH4_flux_nmolpm2ps ~ SPECIES + (1 | Tree), 
                     data = filter(fluxes, location == "upland"))
summary(model_upland)
emm_upland <- emmeans(model_upland, ~ SPECIES)
pairs(emm_upland)

# 3. Generalists only: location effect for hem and rm
model_generalists <- lmer(CH4_flux_nmolpm2ps ~ SPECIES * location + (1 | Tree), 
                          data = filter(fluxes, SPECIES %in% c("hem", "rm")))
summary(model_generalists)
emm_gen <- emmeans(model_generalists, ~ location | SPECIES)
pairs(emm_gen)

# --------------------------------------------
# CLEANER PLOT FOR YOUR DESIGN
# --------------------------------------------

# Add specialist/generalist designation
fluxes <- fluxes %>%
  mutate(strategy = case_when(
    SPECIES == "bg" ~ "wetland specialist",
    SPECIES == "ro" ~ "upland specialist",
    SPECIES %in% c("hem", "rm") ~ "generalist"
  ))

tree_means <- tree_means %>%
  mutate(strategy = case_when(
    SPECIES == "bg" ~ "wetland specialist",
    SPECIES == "ro" ~ "upland specialist",
    SPECIES %in% c("hem", "rm") ~ "generalist"
  ))

# Plot emphasizing the ecological structure
p_eco <- ggplot(tree_means, aes(x = SPECIES, y = mean_CH4, fill = strategy)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~location, scales = "free") +
  scale_fill_manual(values = c("generalist" = "gray70", 
                               "wetland specialist" = "steelblue", 
                               "upland specialist" = "tan")) +
  scale_y_continuous(trans = "pseudo_log") +
  labs(x = "Species", y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       title = "CH4 flux: specialists vs generalists") +
  theme_classic()

print(p_eco)



library(patchwork)

# Make sure species_full is in tree_means
tree_means <- tree_means %>%
  mutate(species_full = case_when(
    SPECIES == "bg"  ~ "Nyssa sylvatica",
    SPECIES == "hem" ~ "Tsuga canadensis",
    SPECIES == "rm"  ~ "Acer rubrum",
    SPECIES == "ro"  ~ "Quercus rubra"
  ))

# --------------------------------------------
# FIGURE 1: Main comparison (log scale)
# --------------------------------------------

fig1 <- ggplot(tree_means, aes(x = species_full, y = mean_CH4, fill = location)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("upland" = "#D2691E", "wetland" = "#4682B4"),
                    labels = c("Upland", "Wetland")) +
  scale_y_continuous(trans = "pseudo_log", breaks = c(0, 0.1, 1, 10)) +
  labs(x = NULL, 
       y = expression("Tree mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       fill = "Location") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
        legend.position = "top")

print(fig1)

# --------------------------------------------
# FIGURE 2: Linear scale showing bg dominance
# --------------------------------------------

fig2 <- ggplot(tree_means, aes(x = species_full, y = mean_CH4, fill = location)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("upland" = "#D2691E", "wetland" = "#4682B4"),
                    labels = c("Upland", "Wetland")) +
  labs(x = NULL, 
       y = expression("Tree mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       fill = "Location") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
        legend.position = "top")

print(fig2)

# --------------------------------------------
# FIGURE 3: Temporal trends by location
# --------------------------------------------

# Daily/weekly means by location
temporal_location <- fluxes %>%
  group_by(date, location) %>%
  summarise(
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

fig3 <- ggplot(temporal_location, aes(x = date, y = mean_CH4, color = location)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, alpha = 0.2) +
  scale_color_manual(values = c("upland" = "#D2691E", "wetland" = "#4682B4"),
                     labels = c("Upland", "Wetland")) +
  labs(x = "Date", 
       y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       color = "Location") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

print(fig3)

# --------------------------------------------
# FIGURE 4: Temporal trends by location × species
# --------------------------------------------

temporal_species <- fluxes %>%
  mutate(species_full = case_when(
    SPECIES == "bg"  ~ "Nyssa sylvatica",
    SPECIES == "hem" ~ "Tsuga canadensis",
    SPECIES == "rm"  ~ "Acer rubrum",
    SPECIES == "ro"  ~ "Quercus rubra"
  )) %>%
  group_by(date, location, species_full) %>%
  summarise(
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

fig4 <- ggplot(temporal_species, aes(x = date, y = mean_CH4, color = species_full)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~location, ncol = 1, scales = "free_y") +
  labs(x = "Date", 
       y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       color = "Species") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        legend.text = element_text(face = "italic"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12))

print(fig4)

# --------------------------------------------
# FIGURE 5: Combined panel figure
# --------------------------------------------

# Inset showing non-bg species on linear scale
fig2_inset <- tree_means %>%
  filter(SPECIES != "bg") %>%
  ggplot(aes(x = species_full, y = mean_CH4, fill = location)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = c("upland" = "#D2691E", "wetland" = "#4682B4")) +
  labs(x = NULL, y = NULL, title = "Excluding N. sylvatica") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1),
        legend.position = "none",
        plot.title = element_text(size = 10))

# Combine main + inset
fig_combined <- fig2 + 
  inset_element(fig2_inset, left = 0.5, bottom = 0.5, right = 0.98, top = 0.98)

print(fig_combined)


# --------------------------------------------
# FIGURE 3: Temporal trends by location (lines + SE)
# --------------------------------------------

temporal_location <- fluxes %>%
  group_by(date, location) %>%
  summarise(
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

fig3 <- ggplot(temporal_location, aes(x = date, y = mean_CH4, color = location)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 5) +
  scale_color_manual(values = c("upland" = "#D2691E", "wetland" = "#4682B4"),
                     labels = c("Upland", "Wetland")) +
  labs(x = "Date", 
       y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       color = "Location") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

print(fig3)

# --------------------------------------------
# FIGURE 4: Temporal trends by location × species (lines + SE)
# --------------------------------------------

temporal_species <- fluxes %>%
  mutate(species_full = case_when(
    SPECIES == "bg"  ~ "Nyssa sylvatica",
    SPECIES == "hem" ~ "Tsuga canadensis",
    SPECIES == "rm"  ~ "Acer rubrum",
    SPECIES == "ro"  ~ "Quercus rubra"
  )) %>%
  group_by(date, location, species_full) %>%
  summarise(
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

fig4 <- ggplot(temporal_species, aes(x = date, y = mean_CH4, color = species_full)) +
  geom_line() +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 5) +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~location, ncol = 1, scales = "free_y") +
  labs(x = "Date", 
       y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       color = "Species") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        legend.text = element_text(face = "italic"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12))

print(fig4)


fluxes <- fluxes %>%
  arrange(date) %>%
  mutate(
    days_since_last = as.numeric(date - lag(date), units = "days"),
    days_since_last = replace_na(days_since_last, 0),
    new_round = days_since_last > 2,
    sampling_round = cumsum(new_round) + 1
  )

# Check it worked
fluxes %>% 
  group_by(sampling_round) %>% 
  summarise(start = min(date), end = max(date), n_days = n_distinct(date), n_obs = n())

# Aggregate by sampling round and location
temporal_round <- fluxes %>%
  group_by(sampling_round, location) %>%
  summarise(
    date = mean(date),  # midpoint of sampling round
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

fig3 <- ggplot(temporal_round, aes(x = date, y = mean_CH4, color = location)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 5) +
  scale_color_manual(values = c("upland" = "#D2691E", "wetland" = "#4682B4"),
                     labels = c("Upland", "Wetland")) +
  labs(x = "Date", 
       y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       color = "Location") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

print(fig3)

# By species
temporal_round_species <- fluxes %>%
  mutate(species_full = case_when(
    SPECIES == "bg"  ~ "Nyssa sylvatica",
    SPECIES == "hem" ~ "Tsuga canadensis",
    SPECIES == "rm"  ~ "Acer rubrum",
    SPECIES == "ro"  ~ "Quercus rubra"
  )) %>%
  group_by(sampling_round, location, species_full) %>%
  summarise(
    date = mean(date),
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

fig4 <- ggplot(temporal_round_species, aes(x = date, y = mean_CH4, color = species_full)) +
  geom_line() +
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 5) +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~location, ncol = 1, scales = "free_y") +
  labs(x = "Date", 
       y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       color = "Species") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right",
        legend.text = element_text(face = "italic"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 12))

print(fig4)
























library(tidyverse)
library(viridis)
library(patchwork)

# ============================================
# HARVARD FOREST - INDIVIDUAL TREES OVER TIME
# ============================================

# Assuming fluxes is already loaded and cleaned from earlier
# If not, reload and clean:

# fluxes <- read.csv('/Users/jongewirtzman/Google Drive/Research/tree-flux-2025/data/HF_2023-2025_tree_flux.csv')
# 
# fluxes <- fluxes %>%
#   mutate(CH4_flux_nmolpm2ps = ifelse(year < 2025, CH4_flux_nmolpm2ps * 1000, CH4_flux_nmolpm2ps)) %>%
#   group_by(Tree) %>%
#   mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
#   ungroup() %>%
#   filter(CH4_flux_nmolpm2ps >= -1, !is.na(PLOT)) %>%
#   mutate(
#     date = as.Date(datetime_posx),
#     location = ifelse(PLOT == "BGS", "wetland", "upland"),
#     species_full = case_when(
#       SPECIES == "bg"  ~ "Nyssa sylvatica",
#       SPECIES == "hem" ~ "Tsuga canadensis",
#       SPECIES == "rm"  ~ "Acer rubrum",
#       SPECIES == "ro"  ~ "Quercus rubra"
#     ),
#     Tree = factor(Tree)
#   )

# Add sampling rounds (>7 day gap = new round)
fluxes <- fluxes %>%
  arrange(date) %>%
  mutate(
    days_since_last = as.numeric(date - lag(date), units = "days"),
    days_since_last = replace_na(days_since_last, 0),
    new_round = days_since_last > 7,
    sampling_round = cumsum(new_round) + 1
  )

# --------------------------------------------
# Figure 1: Individual tree time series (faceted by species × location)
# --------------------------------------------

fig_trees <- ggplot(fluxes, aes(x = date, y = CH4_flux_nmolpm2ps, 
                                group = Tree, color = factor(Tree))) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.5) +
  scale_color_viridis_d() +
  facet_grid(location ~ species_full, scales = "free_y") +
  labs(x = "Date",
       y = expression("CH"[4]~"flux (nmol m"^-2~"s"^-1~")")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        legend.position = "none")

print(fig_trees)

# --------------------------------------------
# Figure 2: Rank over time (are high emitters consistent?)
# --------------------------------------------

fluxes_ranked <- fluxes %>%
  group_by(date, SPECIES, location) %>%
  mutate(rank = rank(-CH4_flux_nmolpm2ps, ties.method = "average")) %>%
  ungroup()

fig_ranks <- ggplot(fluxes_ranked, aes(x = date, y = rank, 
                                       group = Tree, color = factor(Tree))) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.5) +
  scale_color_viridis_d() +
  scale_y_reverse() +
  facet_grid(location ~ species_full) +
  labs(x = "Date",
       y = "Rank (1 = highest flux)") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        legend.position = "none")

print(fig_ranks)

# --------------------------------------------
# Figure 3: Highlight top emitters vs others
# --------------------------------------------

top_emitters <- fluxes %>%
  group_by(Tree, SPECIES, location) %>%
  summarise(mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE), .groups = "drop") %>%
  group_by(SPECIES, location) %>%
  slice_max(mean_CH4, n = 2) %>%
  pull(Tree)

fluxes <- fluxes %>%
  mutate(is_top = Tree %in% top_emitters)

fig_top <- ggplot(fluxes, aes(x = date, y = CH4_flux_nmolpm2ps, group = Tree)) +
  geom_line(aes(alpha = is_top, color = is_top)) +
  geom_point(aes(alpha = is_top, color = is_top), size = 1.5) +
  scale_alpha_manual(values = c("FALSE" = 0.15, "TRUE" = 0.9)) +
  scale_color_manual(values = c("FALSE" = "gray50", "TRUE" = "firebrick")) +
  facet_grid(location ~ species_full, scales = "free_y") +
  labs(x = "Date",
       y = expression("CH"[4]~"flux (nmol m"^-2~"s"^-1~")")) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        legend.position = "none")

print(fig_top)

# --------------------------------------------
# Figure 4: Mean ± SE by sampling round (by location)
# --------------------------------------------

temporal_round <- fluxes %>%
  group_by(sampling_round, location) %>%
  summarise(
    date = mean(date),
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    se_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

fig_temporal <- ggplot(temporal_round, aes(x = date, y = mean_CH4, color = location)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 5) +
  scale_color_manual(values = c("upland" = "#C2703D", "wetland" = "#3D7C9C"),
                     labels = c("Upland", "Wetland")) +
  labs(x = "Date",
       y = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       color = "Location") +
  theme_classic(base_size = 12) +
  theme(legend.position = "top")

print(fig_temporal)



# Split into time periods and correlate
tree_by_period <- fluxes %>%
  mutate(period = case_when(
    date < as.Date("2024-01-01") ~ "2023",
    date < as.Date("2025-01-01") ~ "2024",
    TRUE ~ "2025"
  )) %>%
  group_by(Tree, SPECIES, location, period) %>%
  summarise(mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_CH4)

# Plot year-to-year consistency
fig_consistency_23_24 <- ggplot(tree_by_period, aes(x = `2023`, y = `2024`, color = SPECIES)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_brewer(palette = "Dark2") +
  facet_wrap(~location, scales = "free") +
  labs(x = "Mean CH4 flux 2023", y = "Mean CH4 flux 2024",
       title = "Tree-level consistency: 2023 vs 2024") +
  theme_classic()

print(fig_consistency_23_24)

# Correlation test
cor.test(tree_by_period$`2023`, tree_by_period$`2024`, use = "complete.obs")

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


variance_decomp <- fluxes %>%
  group_by(SPECIES, location) %>%
  summarise(
    # Between-tree variance (variance of tree means)
    var_between = var(tapply(CH4_flux_nmolpm2ps, Tree, mean, na.rm = TRUE), na.rm = TRUE),
    # Within-tree variance (mean of tree variances)
    var_within = mean(tapply(CH4_flux_nmolpm2ps, Tree, var, na.rm = TRUE), na.rm = TRUE),
    # Ratio
    var_ratio = var_between / var_within,
    .groups = "drop"
  )

print(variance_decomp)
# var_ratio > 1 means trees are more different from each other than from themselves over time

# Calculate tree means and rank once
tree_means <- fluxes %>%
  group_by(Tree, SPECIES, location) %>%
  summarise(
    mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE),
    sd_CH4 = sd(CH4_flux_nmolpm2ps, na.rm = TRUE),
    cv = sd_CH4 / mean_CH4,
    n = n(),
    .groups = "drop"
  )

# Plot: mean vs CV (are high emitters more or less variable?)
fig_mean_cv <- ggplot(tree_means, aes(x = mean_CH4, y = cv, color = SPECIES)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_x_continuous(trans = "pseudo_log") +
  facet_wrap(~location, scales = "free") +
  scale_color_brewer(palette = "Dark2") +
  labs(x = expression("Mean CH"[4]~"flux (nmol m"^-2~"s"^-1~")"),
       y = "Coefficient of variation",
       title = "Are high emitters more consistent?") +
  theme_classic()

print(fig_mean_cv)







# Calculate tree mean flux and assign rank
tree_means <- fluxes %>%
  group_by(Tree, SPECIES, location, species_full) %>%
  summarise(mean_CH4 = mean(CH4_flux_nmolpm2ps, na.rm = TRUE), .groups = "drop") %>%
  group_by(SPECIES, location) %>%
  mutate(flux_rank = rank(-mean_CH4)) %>%  # 1 = highest emitter
  ungroup()

# Join back to main data
fluxes <- fluxes %>%
  left_join(tree_means %>% dplyr::select(Tree, mean_CH4, flux_rank), by = "Tree")

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
  dplyr::select(-any_of(c("mean_CH4", "flux_rank", "mean_scaled"))) %>%  # remove if exists
  
  left_join(tree_means %>% dplyr::select(Tree, mean_CH4, mean_scaled), by = "Tree")

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























# ============================================
# Figure 2: CH4 Flux - Final Version
# Faceted main plot with inset, custom tick marks
# DOUBLED TEXT SIZES
# ============================================

library(tidyverse)
library(patchwork)

# --------------------------------------------
# 1. DATA PREPARATION
# --------------------------------------------

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
    location = factor(location, levels = c("wetland", "upland")),
    strategy = factor(strategy, levels = c("Wetland specialist", "Generalist", "Upland specialist"))
  )

# --------------------------------------------
# 2. THEME AND COLORS (1.5x base_size: 11 -> 16)
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
# 3. MAIN PLOT (faceted by species)
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
  scale_x_discrete(labels = c("wetland" = "Wetland", "upland" = "Upland"),
                   drop = TRUE) +
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
    #legend.margin = margin(t = 5),
    strip.placement = "outside",
    strip.background = element_blank(),
    legend.background = element_rect(color = "gray40", fill = "white", linewidth = 0.3),
    strip.text = element_text(face = "italic", size = 14, margin = margin(t = 5)),
    panel.spacing = unit(0.3, "lines")
  )

# --------------------------------------------
# 4. INSET - SEPARATE PANELS WITH CUSTOM TICKS
# --------------------------------------------

short_labels <- c(
  "Nyssa sylvatica" = expression(italic("N. sylvatica")),
  "Acer rubrum" = expression(italic("A. rubrum")),
  "Tsuga canadensis" = expression(italic("T. canadensis")),
  "Quercus rubra" = expression(italic("Q. rubra"))
)

# Common theme for inset panels (1.5x sizes)
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
  filter(location == "wetland") %>%
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
  filter(location == "upland") %>%
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

# Combine inset panels side by side
p_eco_inset <- (inset_wetland | inset_upland) + plot_layout(widths = c(1, 1))

# --------------------------------------------
# 5. COMBINE MAIN + INSET
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

# --------------------------------------------
# 6. SAVE
# --------------------------------------------

ggsave("fig2_final.pdf", fig2_final, width = 9, height = 6, units = "in")
ggsave("fig2_final.png", fig2_final, width = 9, height = 6, units = "in", dpi = 300)

# Also save without inset for flexibility
ggsave("fig2_main_only.pdf", fig2_facet, width = 8, height = 5.5, units = "in")
ggsave("fig2_main_only.png", fig2_facet, width = 8, height = 5.5, units = "in", dpi = 300)


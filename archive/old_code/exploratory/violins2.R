ymf <- read.csv("data/raw/ymf_dataset.csv")
library(ggplot2)
library(dplyr)
library(viridis)
library(gghalves)

# Create a mapping from species code to common name
species_names <- c(
  "PIST" = "White pine",
  "ACRU" = "Red maple",
  "QURU" = "Red oak",
  "BELE" = "Yellow birch",
  "TSCA" = "Eastern hemlock",
  "BEAL" = "Black birch",
  "FRAM" = "White ash",
  "BEPA" = "Paper birch",
  "ACSA" = "Sugar maple",
  "CALA" = "American chestnut",
  "QUAL" = "White oak",
  "QUVE" = "Black oak",
  "BEPO" = "Gray birch",
  "PRSE" = "Black cherry"
)

# Modify the dataset: merge BEPO into BEPA, remove CALA
ymf_cleaned <- ymf %>%
  mutate(Species.Code = ifelse(Species.Code %in% c("BEPA", "BEPO"), "BEPA", Species.Code)) %>%
  filter(Species.Code != "CALA") %>%
  # Add a new column with common names
  mutate(Species.Common = species_names[Species.Code])

# Compute mean, SE, and sample size for each species
summary_stats <- ymf_cleaned %>%
  group_by(Species.Code, Species.Common) %>%
  summarise(
    mean_flux = mean(CH4_flux, na.rm = TRUE),
    se_flux = sd(CH4_flux, na.rm = TRUE) / sqrt(n()),  # Standard Error
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(
    Species.Label = paste0(Species.Common, " (n=", n, ")"),  # Append n to species name
    Species.Label = factor(Species.Label, levels = Species.Label[order(-mean_flux)])  # Order by mean_flux descending
  )

# Add the Species.Label to ymf_cleaned to ensure alignment
ymf_cleaned <- ymf_cleaned %>%
  left_join(summary_stats %>% select(Species.Code, Species.Label), by = "Species.Code")

# Calculate y-axis limits to avoid clipping points
flux_range <- range(ymf_cleaned$CH4_flux * 1000, na.rm = TRUE)
buffer <- (flux_range[2] - flux_range[1]) * 0.1
y_limits <- c(flux_range[1] - buffer, flux_range[2] + buffer)

# 1. Create the bar plot with individual points
mean_plot <- ggplot() +
  # Add bars for mean values
  geom_bar(data = summary_stats, 
           aes(x = Species.Label, y = mean_flux * 1000, fill = mean_flux),
           stat = "identity", color = "black") +
  # Add error bars
  geom_errorbar(data = summary_stats,
                aes(x = Species.Label, 
                    ymin = (mean_flux - se_flux) * 1000, 
                    ymax = (mean_flux + se_flux) * 1000),
                width = 0.2) +
  # Add individual data points
  geom_jitter(data = ymf_cleaned, 
              aes(x = Species.Label, y = CH4_flux * 1000), 
              width = 0.2, height = 0, 
              size = 1.5, color = "darkblue") +
  scale_fill_viridis(option = "D", direction = -1) +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(y = "CH4 Flux (nmol CH4/m²/s)", 
       x = "Species", 
       title = "CH4 Flux by Tree Species with Mean, SE, and Individual Measurements") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")

print(mean_plot)

# 2. Create a half-violin plot with jittered points
half_violin_plot <- ggplot() +
  # Half violin on the left
  geom_half_violin(data = ymf_cleaned, 
                   aes(x = Species.Label, y = CH4_flux * 1000, fill = Species.Code), 
                   side = "l", alpha = 0.7) +
  # Add jittered points on the right
  geom_jitter(data = ymf_cleaned,
              aes(x = Species.Label, y = CH4_flux * 1000),
              width = 0.1, height = 0, 
              alpha = 0.7, size = 1.5, color = "darkblue") +
  # Add mean and SE
  geom_errorbar(data = summary_stats,
                aes(x = Species.Label, 
                    ymin = (mean_flux - se_flux) * 1000, 
                    ymax = (mean_flux + se_flux) * 1000),
                width = 0.1) +
  geom_point(data = summary_stats,
             aes(x = Species.Label, y = mean_flux * 1000),
             size = 2.5, shape = 18, color = "black") +
  # Horizontal line at y=0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Use viridis color palette for discrete values
  scale_fill_viridis_d(option = "D") +
  coord_cartesian(ylim = y_limits) +
  theme_minimal() +
  labs(y = "CH4 Flux (nmol CH4/m²/s)", 
       x = "Species", 
       title = "CH4 Flux by Tree Species (Half-Violin with Points)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

print(half_violin_plot)

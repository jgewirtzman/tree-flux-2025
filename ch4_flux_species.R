library(ggplot2)
library(dplyr)
library(viridis)

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

# Tree families with a relatedness order
# Order based on APG IV phylogeny (rough approximation):
# Magnoliids -> Rosids -> Asterids, and within them by family
tree_families <- c(
  "White pine" = "Pinaceae",        # Gymnosperms
  "Eastern hemlock" = "Pinaceae",   # Gymnosperms
  "Black cherry" = "Rosaceae",      # Rosids - earlier diverging
  "Yellow birch" = "Betulaceae",    # Rosids
  "Black birch" = "Betulaceae",     # Rosids
  "Paper birch" = "Betulaceae",     # Rosids
  "Gray birch" = "Betulaceae",      # Rosids
  "American chestnut" = "Fagaceae", # Rosids
  "Red oak" = "Fagaceae",           # Rosids
  "White oak" = "Fagaceae",         # Rosids
  "Black oak" = "Fagaceae",         # Rosids
  "Red maple" = "Sapindaceae",      # Rosids - was Aceraceae, now in Sapindaceae
  "Sugar maple" = "Sapindaceae",    # Rosids - was Aceraceae, now in Sapindaceae
  "White ash" = "Oleaceae"          # Asterids - most derived
)

# Numerical relatedness index for viridis scale
# We'll create a position value for each family based on their phylogenetic position
family_order <- c(
  "Pinaceae" = 1,      # Gymnosperms - earliest diverging
  "Rosaceae" = 2,      # Early diverging angiosperms
  "Betulaceae" = 3,    # Fagales order
  "Fagaceae" = 4,      # Fagales order
  "Sapindaceae" = 5,   # Malvids clade
  "Oleaceae" = 6       # Asterids - most derived
)

# Modify the dataset: merge BEPO into BEPA, remove CALA
ymf_cleaned <- ymf %>%
  mutate(Species.Code = ifelse(Species.Code %in% c("BEPA", "BEPO"), "BEPA", Species.Code)) %>%
  filter(Species.Code != "CALA")

# Compute mean, SE, and sample size for each species
summary_stats <- ymf_cleaned %>%
  group_by(Species.Code) %>%
  summarise(
    mean_flux = mean(CH4_flux, na.rm = TRUE),
    se_flux = sd(CH4_flux, na.rm = TRUE) / sqrt(n()),  # Standard Error
    n = n()
  ) %>%
  mutate(
    # Create labels with common names
    Species.Label = paste0(species_names[Species.Code], " (n=", n, ")"),
    # Keep original ordering from the code you provided
    Species.Label = factor(Species.Label, levels = Species.Label[order(-mean_flux)])
  )

# Add family to summary_stats
summary_stats$common_name <- gsub(" \\(n=.*\\)$", "", summary_stats$Species.Label)
summary_stats$family <- tree_families[summary_stats$common_name]
summary_stats$family_order <- family_order[summary_stats$family]

# Color by family with viridis scale based on relatedness
mean_plot_family <- ggplot(summary_stats, aes(x = Species.Label, y = mean_flux * 1000, 
                                              fill = family_order)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = (mean_flux - se_flux) * 1000, ymax = (mean_flux + se_flux) * 1000),
                width = 0.2) +
  # Use viridis scale to color by family relatedness
  scale_fill_viridis_c(
    option = "viridis",  # You can also try "magma", "plasma", "cividis", or "inferno"
    name = "Plant Family",
    labels = function(x) names(family_order)[match(x, family_order)],
    breaks = family_order,
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = 10
    )
  ) +
  theme_minimal() +
  labs(y = expression("Mean CH"[4]~"Flux (nmol CH"[4]*"/m²/s)"), 
       x = "Species", 
       title = expression("CH"[4]~"Flux per Species Colored by Phylogenetic Relatedness")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Display the plot
print(mean_plot_family)

# Alternative: Discrete colors by family with a viridis palette
# This approach is more visually distinct for different families
mean_plot_family_discrete <- ggplot(summary_stats, aes(x = Species.Label, y = mean_flux * 1000, 
                                                       fill = family)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = (mean_flux - se_flux) * 1000, ymax = (mean_flux + se_flux) * 1000),
                width = 0.2) +
  # Use a discrete viridis palette
  scale_fill_viridis_d(
    option = "viridis",
    name = "Plant Family",
    # Reorder the family legend by phylogenetic position
    breaks = names(sort(family_order))
  ) +
  theme_minimal() +
  labs(y = expression("Mean CH"[4]~"Flux (nmol CH"[4]*"/m²/s)"), 
       x = "Species", 
       title = expression("CH"[4]~"Flux per Species by Plant Family")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

# Display the discrete version
print(mean_plot_family_discrete)

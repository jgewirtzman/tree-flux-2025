


library(lubridate)
library(tidyverse)

# Load data
fluxes <- read.csv("data/raw/flux_dataset.csv")
fluxes<-fluxes%>%filter(CH4_flux>=-.02)

# Convert jday and year into a date format
fluxes$date <- as.Date(fluxes$jday, origin = paste0(fluxes$year, "-01-01"))

# Add the location column based on the plot values
fluxes <- fluxes %>%
  filter(!is.na(Plot) & !is.na(species) & !is.na(Tree) & !is.na(CH4_flux)) %>%
  mutate(location = ifelse(Plot == "BGS", "wetland", "upland"))

# Get unique species-location combinations
species_location_combos <- fluxes %>%
  distinct(species, location) %>%
  arrange(species, location)

# Generate separate plots for each species-location combination
plots <- list()

for (i in 1:nrow(species_location_combos)) {
  sp <- species_location_combos$species[i]
  loc <- species_location_combos$location[i]
  
  # Filter data for current species and location
  plot_data <- fluxes %>% filter(species == sp, location == loc)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = date, y = CH4_flux*1000, group = Tree, color = Tree)) +
    geom_line(alpha = 0.8, linewidth = 0.7) +  # Connect the same tree across time
    geom_point(size = 2) +  # Points at each time
    labs(
      title = paste("CH4 Flux Rate Over Time for", sp, "in", loc),
      x = "Date",
      y = "CH4 Flux Rate",
      color = "Tree ID"
    ) +
    scale_color_viridis_d() +  # Unique colors for each tree
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Store plot in list
  plots[[paste(sp, loc, sep = "_")]] <- p
}

# Display all plots
#plots
patchwork::wrap_plots(plots)



library(lubridate)
library(tidyverse)

# Load data
fluxes <- read.csv("data/raw/flux_dataset.csv")

# Convert jday and year into a date format
fluxes$date <- as.Date(fluxes$jday, origin = paste0(fluxes$year, "-01-01"))

# Add the location column based on the plot values
fluxes <- fluxes %>%
  filter(!is.na(Plot) & !is.na(species) & !is.na(Tree)) %>%
  mutate(location = ifelse(Plot == "BGS", "wetland", "upland"))

# Rank CH4 flux within each date, species, and location
fluxes <- fluxes %>%
  group_by(date, species, location) %>%
  mutate(rank = rank(-CH4_flux, na.last = "keep")) %>%  # Descending rank (- sign)
  ungroup()

# Get unique species-location combinations
species_location_combos <- fluxes %>%
  distinct(species, location) %>%
  arrange(species, location)

# Generate separate plots for each species-location combination
plots <- list()

for (i in 1:nrow(species_location_combos)) {
  sp <- species_location_combos$species[i]
  loc <- species_location_combos$location[i]
  
  # Filter data for current species and location
  plot_data <- fluxes %>% filter(species == sp, location == loc)
  
  # Create plot
  p <- ggplot(plot_data, aes(x = date, y = rank, group = Tree, color = Tree)) +
    geom_line(alpha = 0.8, linewidth = 0.7) +  # Connect the same tree across time
    geom_point(size = 2) +  # Points at each time
    labs(
      title = paste("Ordinal Ranking of CH4 Fluxes for", sp, "in", loc),
      x = "Date",
      y = "Rank of CH4 Flux",
      color = "Tree ID"
    ) +
    scale_color_viridis_d() +  # Unique colors for each tree
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Store plot in list
  plots[[paste(sp, loc, sep = "_")]] <- p
}

# Display all plots
#plots
patchwork::wrap_plots(plots)




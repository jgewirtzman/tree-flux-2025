# Load necessary libraries
library(ggplot2)
library(dplyr)
library(patchwork)  # For combining plots

set.seed(123)

# Simulated data for Wetland
wetland_data <- data.frame(
  soil_ORP = runif(50, -200, 100)  # Wide range for ORP
)
wetland_data$CH4_flux <- 15 - 0.05 * wetland_data$soil_ORP + rnorm(50, 0, 2)  # Strong negative correlation
wetland_data$decay <- runif(50, 0, 5)  # Random decay (low correlation)
wetland_data$CH4_flux <- wetland_data$CH4_flux + rnorm(50, 0, 1)  # Adding noise

# Simulated data for Upland
upland_data <- data.frame(
  decay = runif(50, 1, 5)  # Strong decay variability
)
upland_data$CH4_flux <- 5 + 2 * upland_data$decay + rnorm(50, 0, 1)  # Strong positive correlation
upland_data$soil_ORP <- runif(50, -50, 50)  # Small ORP range (low variability)
upland_data$CH4_flux <- upland_data$CH4_flux + rnorm(50, 0, 1)  # Adding noise

# Wetland Hypothesis Plots
wetland_ORP_plot <- ggplot(wetland_data, aes(x = soil_ORP, y = CH4_flux)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  #labs(title = "Wetland: CH4 Flux vs. Soil ORP", x = "Soil ORP", y = "CH₄ Flux") +
  theme_minimal()

wetland_decay_plot <- ggplot(wetland_data, aes(x = decay, y = CH4_flux)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  #labs(title = "Wetland: CH4 Flux vs. Internal Decay", x = "Decay", y = "CH₄ Flux") +
  theme_minimal()

# Upland Hypothesis Plots
upland_ORP_plot <- ggplot(upland_data, aes(x = soil_ORP, y = CH4_flux)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  #labs(title = "Upland: CH4 Flux vs. Soil ORP", x = "Soil ORP", y = "CH4 Flux") +
  theme_minimal()

upland_decay_plot <- ggplot(upland_data, aes(x = decay, y = CH4_flux)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  #labs(title = "Upland: CH4 Flux vs. Internal Decay", x = "Decay", y = "CH4 Flux") +
  theme_minimal()

# Arrange hypothesis plots
hypothesis_plots <- (upland_ORP_plot + upland_decay_plot) / (wetland_ORP_plot + wetland_decay_plot)

# Combine with tree mean boxplots
final_plot <- tree_means_plot | hypothesis_plots  # Assuming 'tree_means_plot' is already defined

# Print final combined plot
print(final_plot)

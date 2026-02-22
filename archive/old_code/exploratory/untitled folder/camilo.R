# Load required libraries
library(tidyverse)
library(readxl)
library(lubridate)

# Read the data
flux_data <- read_excel('/Users/jongewirtzman/Downloads/Fluxes_mg (1).xlsx', sheet = "press")

# Clean and prepare the data
flux_clean <- flux_data %>%
  # Convert date to proper format
  mutate(Date = as.Date(Date)) %>%
  # Rename columns for easier handling
  rename(
    date = Date,
    hour = Hour,
    individual = Individual,
    measurement_type = type,
    ch4_flux_nmol = `CH4 Flux (nmol/m2/s)`,
    ch4_flux_mg = `CH4 flux (mg C-CH4 m-2 h-1)`
  ) %>%
  # Filter out rows without proper data
  filter(!is.na(individual), !is.na(measurement_type), !is.na(ch4_flux_mg)) %>%
  # Extract species from individual names
  mutate(
    species = case_when(
      str_detect(individual, "BlackGum") ~ "bg",
      str_detect(individual, "SwampBay") ~ "sb",
      TRUE ~ "unknown"
    ),
    species_full = case_when(
      species == "bg" ~ "Nyssa sylvatica",
      species == "sb" ~ "Persea palustris",
      TRUE ~ individual
    ),
    # Clean up measurement type names
    location = case_when(
      measurement_type == "dbh" ~ "DBH (~1.3m)",
      measurement_type == "bottom" ~ "Bottom (~0.4m)",
      measurement_type == "collar" ~ "Soil collar",
      TRUE ~ measurement_type
    )
  ) %>%
  # Filter out extreme outliers (optional - adjust thresholds as needed)
  filter(ch4_flux_mg > 0 & ch4_flux_mg < 50)

# Check the data structure
print("Species found:")
print(unique(flux_clean$species_full))
print("Measurement locations:")
print(unique(flux_clean$location))
print("Data summary:")
print(summary(flux_clean$ch4_flux_mg))

# Calculate means and standard errors by species and location
flux_summary <- flux_clean %>%
  group_by(species_full, location) %>%
  summarise(
    n = n(),
    mean_flux = mean(ch4_flux_mg, na.rm = TRUE),
    se_flux = sd(ch4_flux_mg, na.rm = TRUE) / sqrt(n()),
    median_flux = median(ch4_flux_mg, na.rm = TRUE),
    sd_flux = sd(ch4_flux_mg, na.rm = TRUE),
    min_flux = min(ch4_flux_mg, na.rm = TRUE),
    max_flux = max(ch4_flux_mg, na.rm = TRUE),
    .groups = 'drop'
  )

# Filter summary to only include species-location combinations present in data
flux_summary <- flux_summary %>%
  semi_join(flux_clean, by = c("species_full", "location"))

# Create faceted boxplot with greyscale colors and improved readability
flux_boxplot <- ggplot() +
  # Add boxplots of raw data
  geom_boxplot(data = flux_clean, 
               aes(x = species_full, y = ch4_flux_mg, fill = species_full),
               alpha = 0.8, outlier.shape = NA) +
  # Add individual data points
  geom_jitter(data = flux_clean,
              aes(x = species_full, y = ch4_flux_mg, color = species_full),
              alpha = 0.6, size = 2, width = 0.2) +
  # Add means as larger points
  geom_point(data = flux_summary,
             aes(x = species_full, y = mean_flux),
             size = 4, shape = 23, fill = "white", color = "black", stroke = 1.5) +
  # Add error bars for standard error
  geom_errorbar(data = flux_summary,
                aes(x = species_full, ymin = mean_flux - se_flux, ymax = mean_flux + se_flux),
                width = 0.3, size = 1, color = "black") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.01),
    breaks = c(0, 0.01, 0.1, 0.5, 1, 2, 5, 10, 20),
    labels = c("0", "0.01", "0.1", "0.5", "1", "2", "5", "10", "20")
  ) +
  scale_fill_manual(
    values = c(
      "Nyssa sylvatica" = "#808080",     # Medium grey
      "Persea palustris" = "#D3D3D3"     # Light grey
    )
  ) +
  scale_color_manual(
    values = c(
      "Nyssa sylvatica" = "#505050",     # Darker medium grey
      "Persea palustris" = "#A9A9A9"     # Darker light grey
    )
  ) +
  labs(
    y = expression("CH"[4] ~ "Flux (mg C-CH"[4] ~ m^-2 ~ h^-1 ~ ")"),
    x = NULL,
    fill = "Species",
    color = "Species"
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +  # Split species names to two lines
  theme_classic() +
  theme(
    axis.text.x = element_text(face = "italic", size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",  # Remove legend
    strip.text = element_text(size = 16, face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )

# Display the plot
print(flux_boxplot)

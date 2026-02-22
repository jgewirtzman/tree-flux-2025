# Load required libraries
library(tidyverse)

# Read the Howland Forest stem flux data
stem_data <- read.csv("/Users/jongewirtzman/Downloads/generic_stemflux_data_ho1 (1).csv")

# Clean and prepare the data
stem_clean <- stem_data %>%
  # Clean species names (remove newline characters)
  mutate(
    spp = str_trim(spp),
    species_full = case_when(
      spp == "RS" ~ "Picea rubens",      # Red spruce
      spp == "RM" ~ "Acer rubrum",       # Red maple
      spp == "EH" ~ "Tsuga canadensis",  # Eastern hemlock
      TRUE ~ spp
    ),
    # Convert units from nmol m-2 s-1 to mg C-CH4 m-2 h-1
    ch4_flux_mg = ch4.flux.nmol.m2.s * 43.236 / 1000,
    # Use actual height values as factor for plotting
    height_cm = as.factor(height.cm)
  ) %>%
  # Remove any missing values
  filter(!is.na(ch4_flux_mg), !is.na(species_full))

# Check the data structure
print("Species found:")
print(unique(stem_clean$species_full))
print("Height values:")
print(unique(stem_clean$height.cm))
print("Data summary:")
print(summary(stem_clean$ch4_flux_mg))

# Calculate summary statistics
stem_summary <- stem_clean %>%
  group_by(species_full) %>%
  summarise(
    n = n(),
    mean_flux = mean(ch4_flux_mg, na.rm = TRUE),
    median_flux = median(ch4_flux_mg, na.rm = TRUE),
    sd_flux = sd(ch4_flux_mg, na.rm = TRUE),
    se_flux = sd_flux / sqrt(n),
    min_flux = min(ch4_flux_mg, na.rm = TRUE),
    max_flux = max(ch4_flux_mg, na.rm = TRUE),
    .groups = 'drop'
  )

# Create the main boxplot with greyscale colors
howland_boxplot <- ggplot() +
  # Add boxplots of raw data
  geom_boxplot(data = stem_clean, 
               aes(x = species_full, y = ch4_flux_mg, fill = species_full),
               alpha = 0.8, outlier.shape = NA) +
  # Add individual data points
  geom_jitter(data = stem_clean,
              aes(x = species_full, y = ch4_flux_mg, color = species_full),
              alpha = 0.4, size = 0.8, width = 0.2) +
  # Add means as larger points
  geom_point(data = stem_summary,
             aes(x = species_full, y = mean_flux),
             size = 4, shape = 23, fill = "white", color = "black", stroke = 1.5) +
  # Add error bars for standard error
  geom_errorbar(data = stem_summary,
                aes(x = species_full, 
                    ymin = mean_flux - se_flux, 
                    ymax = mean_flux + se_flux),
                width = 0.3, size = 1, color = "black") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.01),
    breaks = c(-0.1, -0.01, 0, 0.01, 0.1, 0.5, 1, 2),
    labels = c("-0.1", "-0.01", "0", "0.01", "0.1", "0.5", "1", "2")
  ) +
  scale_fill_manual(
    values = c(
      "Acer rubrum" = "#D3D3D3",         # Light grey
      "Picea rubens" = "#808080",        # Medium grey
      "Tsuga canadensis" = "#404040"     # Dark grey
    )
  ) +
  scale_color_manual(
    values = c(
      "Acer rubrum" = "#A9A9A9",         # Darker light grey
      "Picea rubens" = "#505050",        # Darker medium grey
      "Tsuga canadensis" = "#202020"     # Darker dark grey
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
    plot.margin = margin(20, 20, 20, 20)
  )

# Display the plot
print(howland_boxplot)
# Load required libraries
library(tidyverse)
library(readxl)
library(lubridate)

# Read the data
tonzi_data <- read_excel("/Users/jongewirtzman/Downloads/tonzi_tree_ch4.xlsx")

# Clean and prepare the data
tonzi_clean <- tonzi_data %>%
  # Convert date to proper format
  mutate(date = as.Date(date)) %>%
  # Extract tree ID and height information from point column
  mutate(
    tree_id = case_when(
      str_detect(point, "^T[1-6][ABC]") ~ str_extract(point, "T[1-6]"),
      TRUE ~ point
    ),
    height_position = case_when(
      str_detect(point, "^T[1-6]A") ~ "A",
      str_detect(point, "^T[1-6]B") ~ "B", 
      str_detect(point, "^T[1-6]C") ~ "C",
      TRUE ~ "Other"
    ),
    measurement_type = case_when(
      str_detect(point, "^T[1-6][ABC]") ~ "Tree_stem",
      str_detect(point, "^P[1-5]") ~ "Soil",
      TRUE ~ "Other"
    )
  ) %>%
  # Filter to only tree stem measurements (T1-T6, A/B/C heights)
  filter(measurement_type == "Tree_stem", height_position %in% c("A", "B", "C")) %>%
  # Remove any missing values
  filter(!is.na(ch4_c_ug_m2_h)) %>%
  # Keep original units (μg C m-2 h-1) instead of converting to mg
  mutate(ch4_flux_ug = ch4_c_ug_m2_h) %>%
  # Filter out extreme outliers (adjusted for microgram scale)
  filter(ch4_flux_ug > -200 & ch4_flux_ug < 300) %>%
  # Create flux direction categories and colors
  mutate(
    flux_direction = case_when(
      ch4_flux_ug > 0 ~ "Positive",
      ch4_flux_ug < 0 ~ "Negative", 
      TRUE ~ "Zero"
    ),
    # Create height labels and numeric values for y-axis positioning
    height_numeric = case_when(
      height_position == "A" ~ 1,
      height_position == "B" ~ 2,
      height_position == "C" ~ 3,
      TRUE ~ NA_real_
    ),
    height_label = case_when(
      height_position == "A" ~ "Height A (lowest)",
      height_position == "B" ~ "Height B (~1.3m)",
      height_position == "C" ~ "Height C (highest)",
      TRUE ~ height_position
    )
  )

# Calculate mean positive and negative fluxes for arrows
mean_fluxes <- tonzi_clean %>%
  filter(flux_direction != "Zero") %>%
  group_by(height_numeric, height_position, flux_direction) %>%
  summarise(
    mean_flux = mean(ch4_flux_ug, na.rm = TRUE),
    .groups = 'drop'
  )

# Create the horizontal flux distribution plot
horizontal_flux_plot <- ggplot(tonzi_clean, aes(x = ch4_flux_ug, y = height_numeric)) +
  # Add points for individual measurements, colored by flux direction
  geom_point(aes(color = flux_direction), 
             alpha = 0.6, size = 2.5, 
             position = position_jitter(height = 0.15, width = 0)) +
  # Add vertical line at x = 0
  geom_vline(xintercept = 0, linetype = "solid", color = "black", size = 1.2) +
  # Add white outline arrows first (background)
  geom_segment(data = mean_fluxes %>% filter(flux_direction == "Positive"),
               aes(x = 0, xend = mean_flux, 
                   y = height_numeric + 0.2, yend = height_numeric + 0.2),
               arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
               color = "white", size = 4, linewidth = 4) +
  geom_segment(data = mean_fluxes %>% filter(flux_direction == "Negative"),
               aes(x = 0, xend = mean_flux, 
                   y = height_numeric - 0.2, yend = height_numeric - 0.2),
               arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
               color = "white", size = 4, linewidth = 4) +
  # Add colored arrows on top (foreground)
  geom_segment(data = mean_fluxes %>% filter(flux_direction == "Positive"),
               aes(x = 0, xend = mean_flux, 
                   y = height_numeric + 0.2, yend = height_numeric + 0.2),
               arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
               color = "#E63946", size = 3, linewidth = 3) +
  geom_segment(data = mean_fluxes %>% filter(flux_direction == "Negative"),
               aes(x = 0, xend = mean_flux, 
                   y = height_numeric - 0.2, yend = height_numeric - 0.2),
               arrow = arrow(length = unit(0.5, "cm"), type = "closed"),
               color = "#2E86AB", size = 3, linewidth = 3) +
  # Add boxplots for each height on top (foreground layer)
  geom_boxplot(aes(group = height_position), 
               alpha = 0.3, fill = "gray90", 
               outlier.shape = NA, width = 0.2) +
  # Color scheme: blue for negative (uptake), red for positive (emission)
  scale_color_manual(
    values = c(
      "Negative" = "#2E86AB",  # Blue for uptake/consumption
      "Positive" = "#E63946",  # Red for emission/production
      "Zero" = "#6C757D"       # Gray for zero
    ),
    labels = c(
      "Negative" = "Uptake (Consumption)",
      "Positive" = "Emission (Production)", 
      "Zero" = "Zero Flux"
    )
  ) +
  # Log scale for x-axis with symmetric transformation around zero
  scale_x_continuous(
    trans = scales::pseudo_log_trans(sigma = 1),
    breaks = c(-100, -10, -1, 0, 1, 10, 100),
    labels = c("-100", "-10", "-1", "0", "1", "10", "100")
  ) +
  # Set y-axis to show height labels
  scale_y_continuous(
    breaks = c(1, 2, 3),
    labels = c("Height A\n(lowest)", "Height B\n(~1.3m)", "Height C\n(highest)"),
    limits = c(0.5, 3.5)
  ) +
  # Labels and title (updated units)
  labs(
    x = expression("CH"[4] ~ "Flux (" * mu * "g C-CH"[4] ~ m^-2 ~ h^-1 ~ ")"),
    y = NULL,
    color = "Flux Direction",
    title = expression(italic("Quercus douglasii"))
  ) +
  # Add annotations for clarity
  annotate("text", x = min(tonzi_clean$ch4_flux_ug) * 0.5, y = 3.3, 
           label = "UPTAKE", color = "#2E86AB", fontface = "bold", size = 4) +
  annotate("text", x = max(tonzi_clean$ch4_flux_ug) * 0.5, y = 3.3, 
           label = "EMISSION", color = "#E63946", fontface = "bold", size = 4) +
  # Theme
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major.x = element_line(color = "gray90", size = 0.3),
    panel.grid.minor.x = element_line(color = "gray95", size = 0.2)
  )

# Display the plot
print(horizontal_flux_plot)

# Summary statistics for the horizontal plot (updated for microgram units)
summary_by_height <- tonzi_clean %>%
  group_by(height_position, height_label) %>%
  summarise(
    n_total = n(),
    n_positive = sum(ch4_flux_ug > 0, na.rm = TRUE),
    n_negative = sum(ch4_flux_ug < 0, na.rm = TRUE),
    n_zero = sum(ch4_flux_ug == 0, na.rm = TRUE),
    percent_positive = round(n_positive / n_total * 100, 1),
    percent_negative = round(n_negative / n_total * 100, 1),
    mean_flux = round(mean(ch4_flux_ug, na.rm = TRUE), 2),
    median_flux = round(median(ch4_flux_ug, na.rm = TRUE), 2),
    mean_positive = round(mean(ch4_flux_ug[ch4_flux_ug > 0], na.rm = TRUE), 2),
    mean_negative = round(mean(ch4_flux_ug[ch4_flux_ug < 0], na.rm = TRUE), 2),
    .groups = 'drop'
  )

print("\nSummary by height position:")
print(summary_by_height)
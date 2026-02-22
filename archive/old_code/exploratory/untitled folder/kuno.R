# Load required libraries
library(tidyverse)
library(readxl)
library(lme4)
library(emmeans)
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
  # Filter to only tree stem measurements (T1-T6, A/B heights only - remove C)
  filter(measurement_type == "Tree_stem", height_position %in% c("A", "B")) %>%
  # Remove any missing values
  filter(!is.na(ch4_c_ug_m2_h)) %>%
  # Convert units from μg C m-2 h-1 to mg C m-2 h-1 for consistency
  mutate(ch4_flux_mg = ch4_c_ug_m2_h / 1000) %>%
  # Filter out extreme outliers (optional - adjust thresholds as needed)
  filter(ch4_flux_mg > -200 & ch4_flux_mg < 300) %>%
  # Create height labels without species name
  mutate(
    height_label = case_when(
      height_position == "A" ~ "Height A\n(lowest)",
      height_position == "B" ~ "Height B*\n(~1.3m)",
      TRUE ~ height_position
    )
  )

# Check the data structure
print("Tree IDs:")
print(unique(tonzi_clean$tree_id))
print("Height positions:")
print(unique(tonzi_clean$height_position))
print("Data summary:")
print(summary(tonzi_clean$ch4_flux_mg))

# Fit mixed-effects model using linear mixed model
tonzi_model <- lmer(ch4_flux_mg ~ height_position + (1 | tree_id) + (1 | date),
                    data = tonzi_clean)

# Get estimated marginal means
marginal_means_tonzi <- emmeans(tonzi_model, ~ height_position)
marginal_means_tonzi_df <- as.data.frame(marginal_means_tonzi) %>%
  mutate(
    height_label = case_when(
      height_position == "A" ~ "Height A\n(lowest)",
      height_position == "B" ~ "Height B*\n(~1.3m)",
      TRUE ~ height_position
    )
  )

# Create the plot with greyscale colors and different patterns
tonzi_plot <- ggplot() +
  # Add boxplots with same fill color but different alpha for distinction
  geom_boxplot(data = tonzi_clean,
               aes(x = height_label, y = ch4_flux_mg, 
                   fill = height_position, alpha = height_position),
               outlier.shape = NA, color = "black", size = 0.8) +
  # Add individual points
  geom_jitter(data = tonzi_clean,
              aes(x = height_label, y = ch4_flux_mg, shape = height_position),
              alpha = 0.6, size = 2, width = 0.2, color = "#505050") +
  # Add marginal means
  geom_point(data = marginal_means_tonzi_df,
             aes(x = height_label, y = emmean),
             size = 4, shape = 23, fill = "white", color = "black", stroke = 1.5) +
  geom_errorbar(data = marginal_means_tonzi_df,
                aes(x = height_label, ymin = lower.CL, ymax = upper.CL),
                width = 0.3, size = 1, color = "black") +
  # Create custom transformation that expands around zero
  scale_y_continuous(
    trans = scales::trans_new(
      name = "expand_zero",
      transform = function(x) sign(x) * sqrt(abs(x) * 1000),
      inverse = function(x) sign(x) * (x^2) / 1000
    ),
    breaks = c(-0.1, -0.01, 0, 0.01, 0.1, 0.2),
    labels = c("-0.1", "-0.01", "0", "0.01", "0.1", "0.2")
  ) +
  # Same greyscale color for both heights, different alpha levels
  scale_fill_manual(
    values = c("A" = "#808080", "B" = "#808080"),
    guide = "none"
  ) +
  scale_alpha_manual(
    values = c("A" = 0.8, "B" = 0.4),
    guide = "none"
  ) +
  # Different shapes for distinction
  scale_shape_manual(
    values = c("A" = 16, "B" = 17),  # Circle vs triangle
    guide = "none"
  ) +
  # Add species name annotation below x-axis
  annotate("text", x = 1.5, y = -0.08, 
           label = expression(italic("Quercus douglasii")), 
           size = 4, hjust = 0.5) +
  labs(
    y = expression("CH"[4] ~ "Flux (mg C-CH"[4] ~ m^-2 ~ h^-1 ~ ")"),
    x = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, hjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.grid.major.y = element_line(color = "grey90", size = 0.3),
    plot.margin = margin(20, 20, 40, 20)
  )

# Display the plot
tonzi_plot

# Summary statistics by height
height_summary <- tonzi_clean %>%
  group_by(height_position) %>%
  summarise(
    n = n(),
    mean_flux = mean(ch4_flux_mg, na.rm = TRUE),
    median_flux = median(ch4_flux_mg, na.rm = TRUE),
    sd_flux = sd(ch4_flux_mg, na.rm = TRUE),
    min_flux = min(ch4_flux_mg, na.rm = TRUE),
    max_flux = max(ch4_flux_mg, na.rm = TRUE),
    .groups = 'drop'
  )

print("Summary by height position (A and B only):")
print(height_summary)
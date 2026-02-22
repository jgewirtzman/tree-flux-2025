# Load required libraries
library(tidyverse)
library(lubridate)

# Read the extracted cypress knee data
knee_data <- read.csv("/Users/jongewirtzman/Downloads/cypress_knee_data.csv")

# Clean and prepare the data
knee_clean <- knee_data %>%
  # Filter out rows without knee CH4 data
  filter(!is.na(kneeCH4_mean)) %>%
  # Create simplified site names
  mutate(
    site_short = case_when(
      str_detect(site, "Main Channel") ~ "Main Channel",
      str_detect(site, "Side Channel") ~ "Side Channel", 
      str_detect(site, "Kentucky Lake") ~ "Kentucky Lake",
      TRUE ~ site
    ),
    # Create season categories
    season_simple = case_when(
      str_detect(season, "Fall") ~ "Fall",
      str_detect(season, "Winter") ~ "Winter",
      str_detect(season, "Spring") ~ "Spring", 
      str_detect(season, "Summer") ~ "Summer",
      TRUE ~ season
    ),
    # Convert units from nmol m-2 s-1 to mg C-CH4 m-2 h-1
    ch4_flux_mg = kneeCH4_mean * 43.236 / 1000,
    ch4_flux_mg_sd = kneeCH4_sd * 43.236 / 1000
  ) %>%
  # Filter to positive values for pseudo-log scale and remove Kentucky Lake
  filter(ch4_flux_mg > 0, site_short != "Kentucky Lake")

# Create the main boxplot with greyscale colors
cypress_boxplot <- ggplot() +
  # Add boxplots of the mean data (simulating distribution)
  geom_boxplot(data = knee_clean, 
               aes(x = site_short, y = ch4_flux_mg, fill = site_short),
               alpha = 0.8, outlier.shape = NA) +
  # Add individual data points
  geom_jitter(data = knee_clean,
              aes(x = site_short, y = ch4_flux_mg, color = site_short),
              alpha = 0.6, size = 2, width = 0.2) +
  # Add means as larger points
  geom_point(data = knee_clean,
             aes(x = site_short, y = ch4_flux_mg),
             size = 4, shape = 23, fill = "white", color = "black", stroke = 1.5) +
  # Add error bars for standard deviation
  geom_errorbar(data = knee_clean %>% filter(!is.na(ch4_flux_mg_sd)),
                aes(x = site_short, 
                    ymin = pmax(0.001, ch4_flux_mg - ch4_flux_mg_sd), 
                    ymax = ch4_flux_mg + ch4_flux_mg_sd),
                width = 0.3, size = 1, color = "black") +
  facet_wrap(~ season_simple, scales = "free_x") +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.01),
    breaks = c(0, 0.01, 0.1, 0.5, 1, 5, 10, 20),
    labels = c("0", "0.01", "0.1", "0.5", "1", "5", "10", "20")
  ) +
  scale_fill_manual(
    values = c(
      "Main Channel" = "#808080",        # Medium grey
      "Side Channel" = "#D3D3D3"         # Light grey
    )
  ) +
  scale_color_manual(
    values = c(
      "Main Channel" = "#505050",        # Darker medium grey
      "Side Channel" = "#A9A9A9"         # Darker light grey
    )
  ) +
  labs(
    y = expression("CH"[4] ~ "Flux (mg C-CH"[4] ~ m^-2 ~ h^-1 ~ ")"),
    x = NULL,
    fill = "Site Location",
    color = "Site Location",
    title = expression(italic("Taxodium distichum") ~ "knee flux")
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +  # Split site names to two lines
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",  # Remove legend
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, hjust = 0.5),
    plot.margin = margin(20, 20, 20, 20)
  )

# Display the plot
cypress_boxplot

# Clean and prepare the data for overall plot
knee_clean_overall <- knee_data %>%
  # Filter out rows without knee CH4 data
  filter(!is.na(kneeCH4_mean)) %>%
  # Create simplified site names
  mutate(
    site_short = case_when(
      str_detect(site, "Main Channel") ~ "Main Channel",
      str_detect(site, "Side Channel") ~ "Side Channel", 
      str_detect(site, "Kentucky Lake") ~ "Kentucky Lake",
      TRUE ~ site
    ),
    # Convert units from nmol m-2 s-1 to mg C-CH4 m-2 h-1
    ch4_flux_mg = kneeCH4_mean * 43.236 / 1000,
    ch4_flux_mg_sd = kneeCH4_sd * 43.236 / 1000
  ) %>%
  # Filter to positive values for pseudo-log scale and remove Kentucky Lake
  filter(ch4_flux_mg > 0, site_short != "Kentucky Lake")

# Calculate overall means across all seasons for each site
knee_overall <- knee_clean_overall %>%
  group_by(site_short) %>%
  summarise(
    mean_flux = mean(ch4_flux_mg, na.rm = TRUE),
    sd_flux = sd(ch4_flux_mg, na.rm = TRUE),
    n_seasons = n(),
    .groups = 'drop'
  )

# Create the overall mean plot with greyscale colors
overall_mean_plot <- ggplot() +
  # Add boxplots of the data
  geom_boxplot(data = knee_clean_overall, 
               aes(x = site_short, y = ch4_flux_mg, fill = site_short),
               alpha = 0.8, outlier.shape = NA) +
  # Add individual data points
  geom_jitter(data = knee_clean_overall,
              aes(x = site_short, y = ch4_flux_mg, color = site_short),
              alpha = 0.6, size = 2, width = 0.2) +
  # Add overall means as larger points
  geom_point(data = knee_overall,
             aes(x = site_short, y = mean_flux),
             size = 4, shape = 23, fill = "white", color = "black", stroke = 1.5) +
  # Add error bars for standard deviation
  geom_errorbar(data = knee_overall,
                aes(x = site_short, 
                    ymin = pmax(0.001, mean_flux - sd_flux), 
                    ymax = mean_flux + sd_flux),
                width = 0.3, size = 1, color = "black") +
  
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.01),
    breaks = c(0, 0.01, 0.1, 0.5, 1, 5, 10, 20),
    labels = c("0", "0.01", "0.1", "0.5", "1", "5", "10", "20")
  ) +
  scale_fill_manual(
    values = c(
      "Main Channel" = "#808080",        # Medium grey
      "Side Channel" = "#D3D3D3"         # Light grey
    )
  ) +
  scale_color_manual(
    values = c(
      "Main Channel" = "#505050",        # Darker medium grey
      "Side Channel" = "#A9A9A9"         # Darker light grey
    )
  ) +
  labs(
    y = expression("CH"[4] ~ "Flux (mg C-CH"[4] ~ m^-2 ~ h^-1 ~ ")"),
    x = NULL,
    fill = "Site Location",
    color = "Site Location"
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +  # Clean site names only
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",  # Remove legend
    plot.margin = margin(20, 20, 40, 20)
  ) +
  # Add species name annotation below x-axis labels using annotation_custom
  annotate("text", x = 1.25, y = 0.01, 
           label = expression(italic("Taxodium distichum")), 
           size = 4, hjust = 0.5)

# Display the plot
overall_mean_plot

# Print the summary data
print("Overall means across all seasons (excluding Kentucky Lake):")
print(knee_overall)

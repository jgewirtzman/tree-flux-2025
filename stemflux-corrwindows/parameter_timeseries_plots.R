# Load necessary libraries if not already loaded
library(tidyverse)
library(lubridate)
library(patchwork) # Alternative to cowplot for combining plots

# Convert to long format for easier plotting
wtd_met_long <- wtd_met %>%
  select(datetime, bgs_wtd_cm, bvs_wtd_cm, tair_C, p_kPa, P_mm, 
         RH, PAR, rnet, slrr, s10t, VPD_kPa) %>%
  pivot_longer(cols = -datetime, 
               names_to = "variable", 
               values_to = "value")

# Define variable groups and their units for better organization
variable_info <- tribble(
  ~variable,    ~group,       ~unit,            ~color,
  "bgs_wtd_cm", "Water",      "Depth (cm)",     "darkblue",
  "bvs_wtd_cm", "Water",      "Depth (cm)",     "royalblue",
  "tair_C",     "Temperature", "°C",            "red",
  "s10t",       "Temperature", "°C",            "brown",
  "p_kPa",      "Atmospheric", "kPa",           "purple",
  "VPD_kPa",    "Atmospheric", "kPa",           "darkred",
  "P_mm",       "Precipitation", "mm",          "skyblue",
  "RH",         "Atmospheric", "%",             "darkgreen",
  "PAR",        "Radiation",   "μmol/m²/s",     "orange",
  "rnet",       "Radiation",   "W/m²",          "forestgreen",
  "slrr",       "Radiation",   "W/m²",          "gold"
)

# Join variable information to data
wtd_met_long <- wtd_met_long %>%
  left_join(variable_info, by = "variable")

# Method 1: Create faceted plot by variable group
grouped_facet_plot <- wtd_met_long %>%
  ggplot(aes(x = datetime, y = value, color = variable)) +
  geom_line() +
  facet_wrap(~ group, scales = "free_y") +
  scale_color_manual(values = setNames(variable_info$color, variable_info$variable)) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(title = "Time Series of All Variables by Group",
       x = "Date") 

# Display the faceted plot
grouped_facet_plot

# Save the plot
ggsave("grouped_facet_plot.pdf", grouped_facet_plot, width = 12, height = 10)

# Method 2: Create individual plots programmatically and combine
# Get unique groups
groups <- unique(variable_info$group)

# Create a list to store plots
plot_list <- list()

# Loop through each group and create a plot
for(group_name in groups) {
  # Filter variables for this group
  group_vars <- variable_info %>% 
    filter(group == group_name) %>% 
    pull(variable)
  
  # Get unit for this group (assume same unit per group)
  group_unit <- variable_info %>% 
    filter(group == group_name) %>% 
    pull(unit) %>% 
    first()
  
  # Create plot for this group
  p <- wtd_met_long %>%
    filter(variable %in% group_vars) %>%
    ggplot(aes(x = datetime, y = value, color = variable)) +
    geom_line(linewidth = 0.7) +
    scale_color_manual(values = setNames(
      variable_info %>% filter(variable %in% group_vars) %>% pull(color),
      variable_info %>% filter(variable %in% group_vars) %>% pull(variable)
    )) +
    labs(
      title = paste(group_name, "Variables"),
      y = group_unit,
      x = NULL
    ) +
    theme_bw() +
    theme(
      legend.title = element_blank(),
      legend.position = "right",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Add to plot list
  plot_list[[group_name]] <- p
}

# Combine all plots using patchwork
combined_plots <- wrap_plots(plot_list, ncol = 1) +
  plot_annotation(
    title = 'Time Series of Environmental Variables',
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
  )

# Display the combined plot
combined_plots

# Save the combined plot
#ggsave("combined_time_series.pdf", combined_plots, width = 10, height = 15)

# Method 3: For a more interactive approach with automatic looping
if(require(plotly)) {
  # Create interactive plot
  p <- plot_ly()
  
  # Loop through each variable to add traces
  for(i in 1:nrow(variable_info)) {
    var_name <- variable_info$variable[i]
    var_color <- variable_info$color[i]
    
    p <- p %>% add_trace(
      data = wtd_met,
      x = ~datetime,
      y = as.formula(paste0("~", var_name)),
      name = var_name,
      type = 'scatter',
      mode = 'lines',
      line = list(color = var_color)
    )
  }
  
  # Add layout
  p <- p %>% layout(
    title = "Interactive Time Series of All Variables",
    xaxis = list(title = "Date"),
    yaxis = list(title = "Values (various units)"),
    legend = list(orientation = "h")
  )
  
  # Display the interactive plot
  p
}
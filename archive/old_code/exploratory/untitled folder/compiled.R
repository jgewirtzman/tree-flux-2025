overall_mean_plot
howland_boxplot
flux_boxplot
tonzi_plot

library(patchwork)

# Arrange plots in 4x1 layout (4 rows, 1 column)
combined_plot <- howland_boxplot | overall_mean_plot | flux_boxplot | tonzi_plot

# Display the combined plot
combined_plot


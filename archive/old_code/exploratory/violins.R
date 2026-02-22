# 2. Create a half-violin plot with jittered points on log scale
half_violin_plot <- ggplot() +
  # Half violin on the left
  geom_half_violin(data = ymf_cleaned, 
                   aes(x = Species.Label, y = CH4_flux * 1000, fill = Species.Code), 
                   side = "l", alpha = 0.7) +
  # Add jittered points on the right
  geom_jitter(data = ymf_cleaned,
              aes(x = Species.Label, y = CH4_flux * 1000),
              width = 0.1, height = 0, 
              alpha = 0.7, size = 1.5, color = "darkblue") +
  # Add mean and SE
  geom_errorbar(data = summary_stats,
                aes(x = Species.Label, 
                    ymin = (mean_flux - se_flux) * 1000, 
                    ymax = (mean_flux + se_flux) * 1000),
                width = 0.1) +
  geom_point(data = summary_stats,
             aes(x = Species.Label, y = mean_flux * 1000),
             size = 2.5, shape = 18, color = "black") +
  # Horizontal line at y=0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Use viridis color palette for discrete values
  scale_fill_viridis_d(option = "D") +
  # Apply log transformation with sign preservation
  scale_y_continuous(
    trans = "symlog",  # Symmetric log transformation - handles negative values
    breaks = scales::trans_breaks("symlog", n = 6),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  theme_minimal() +
  labs(y = "CH4 Flux (nmol CH4/m²/s) [log scale]", 
       x = "Species", 
       title = "CH4 Flux by Tree Species (Half-Violin with Points) - Log Scale") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
print(half_violin_plot)

# Alternative approach if 'symlog' isn't available in your ggplot2 version
# You can use the "pseudo_log" transformation from scales package

# Load required package
library(scales)

# Create custom transformation function for pseudo-log (handles negative values)
half_violin_plot_alt <- ggplot() +
  # Half violin on the left
  geom_half_violin(data = ymf_cleaned, 
                   aes(x = Species.Label, y = CH4_flux * 1000, fill = Species.Code), 
                   side = "l", alpha = 0.7) +
  # Add jittered points on the right
  geom_jitter(data = ymf_cleaned,
              aes(x = Species.Label, y = CH4_flux * 1000),
              width = 0.1, height = 0, 
              alpha = 0.7, size = 1.5, color = "darkblue") +
  # Add mean and SE
  geom_errorbar(data = summary_stats,
                aes(x = Species.Label, 
                    ymin = (mean_flux - se_flux) * 1000, 
                    ymax = (mean_flux + se_flux) * 1000),
                width = 0.1) +
  geom_point(data = summary_stats,
             aes(x = Species.Label, y = mean_flux * 1000),
             size = 2.5, shape = 18, color = "black") +
  # Horizontal line at y=0
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Use viridis color palette for discrete values
  scale_fill_viridis_d(option = "D") +
  # Apply pseudo-log transformation
  scale_y_continuous(
    trans = pseudo_log_trans(base = 10),
    breaks = scales::pretty_breaks(n = 8),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  theme_minimal() +
  labs(y = "CH4 Flux (nmol CH4/m²/s) [log scale]", 
       x = "Species", 
       title = "CH4 Flux by Tree Species (Half-Violin with Points) - Log Scale") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
 print(half_violin_plot_alt) # Uncomment to use this alternative approach
 
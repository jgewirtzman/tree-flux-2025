# Filter out missing site or species values
cor_all_filtered <- cor_all %>%
  filter(!is.na(site))


# Reshape data for plotting multiple correlations separately
cor_all_long <- cor_all_filtered %>%
  pivot_longer(cols = starts_with("cor_"), names_to = "predictor", values_to = "correlation")

# Create separate plots for each predictor
ggplot(cor_all_long, aes(x = interval, y = correlation, color=site)) + 
  geom_point() +
  geom_line() +
  facet_wrap(~site) +
  labs(x = "Time Interval (Days)", y = "Correlation") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~predictor, scales = "free_y")

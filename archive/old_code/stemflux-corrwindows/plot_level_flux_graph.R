# Calculate mean and standard error (SE) by location and date
location_flux_summary <- fluxes %>%
  group_by(location, date) %>%
  summarise(
    mean_CH4 = mean(CH4_flux, na.rm = TRUE),
    se_CH4 = sd(CH4_flux, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

# Ensure both locations have a full date range
location_flux_summary <- location_flux_summary %>%
  complete(location, date = seq(min(date), max(date), by = "day")) %>%
  drop_na(mean_CH4)  # Remove rows where mean_CH4 is NA

# Create the plot with scaled flux values (x1000) and error bars
ggplot(location_flux_summary, aes(x = date, y = mean_CH4*1000, color = location)) +
  geom_point(size = 2) +  # Points for mean values
  geom_line() +  # Connect points with lines
  geom_ribbon(aes(ymin = (mean_CH4 - se_CH4)*1000, 
                  ymax = (mean_CH4 + se_CH4)*1000, 
                  fill = location), 
              alpha = 0.2, color = NA) +  # Ribbon for SE
  scale_color_manual(values = c("wetland" = "blue", "upland" = "darkgreen")) +
  scale_fill_manual(values = c("wetland" = "blue", "upland" = "darkgreen")) +
  labs(
    #title = "Mean CH4 Flux by Location Over Time",
    x = "Date",
    y = "CH4 Flux (nmol/m2/s)",
    color = "Location",
    fill = "Location"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14, face = "bold")
  )+
  xlim(as.Date("2023-01-01"), as.Date("2025-01-01"))

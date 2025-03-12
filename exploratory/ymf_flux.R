ymf<-read.csv("data/raw/ymf_dataset.csv")

library(ggplot2)

ggplot(ymf, aes(x=CH4_flux*1000))+
  geom_histogram()+facet_wrap(~Species.Code, scales="free")



library(ggplot2)
library(dplyr)
library(scales)

# Define pseudo-log transformation function
pseudo_log <- function(x) {
  sign(x) * log1p(abs(x))
}

# Compute mean CH4_flux for each species
mean_plot <- ymf %>%
  group_by(Species.Code) %>%
  summarise(mean_flux = mean(CH4_flux, na.rm = TRUE)) %>%
  ggplot(aes(x = Species.Code, y = mean_flux * 1000, fill = Species.Code)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(y = "Mean CH4 Flux (nmol CH4/m²/s)", x = "Species", title = "Mean CH4 Flux per Species") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Density plot with pseudo-log transformation
density_plot <- ggplot(ymf, aes(x = pseudo_log(CH4_flux * 1000), color = Species.Code)) +
  geom_density() +
  theme_minimal() +
  labs(x = "Pseudo-Log Transformed CH4 Flux", y = "Density", title = "Density of CH4 Flux (Transformed)") +
  theme(legend.position = "bottom")

# Histogram with colored points for values above/below zero
histogram_plot <- ggplot(ymf, aes(x = CH4_flux * 1000)) +
  geom_histogram(binwidth = 10, fill = "grey70", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(data = ymf %>% filter(CH4_flux < 0), aes(y = 0, x = CH4_flux * 1000), color = "red", size = 2) +
  geom_point(data = ymf %>% filter(CH4_flux > 0), aes(y = 0, x = CH4_flux * 1000), color = "blue", size = 2) +
  facet_wrap(~Species.Code, scales = "free") +
  theme_minimal() +
  labs(x = "CH4 Flux (nmol CH4/m²/s)", y = "Count", title = "Histogram of CH4 Flux by Species")

# Display the plots
print(mean_plot)
print(density_plot)
print(histogram_plot)



library(ggplot2)
library(dplyr)
library(viridis)

# Modify the dataset: merge BEPO into BEPA, remove CALA
ymf_cleaned <- ymf %>%
  mutate(Species.Code = ifelse(Species.Code %in% c("BEPA", "BEPO"), "BEPA", Species.Code)) %>%
  filter(Species.Code != "CALA")

# Compute mean, SE, and sample size for each species
summary_stats <- ymf_cleaned %>%
  group_by(Species.Code) %>%
  summarise(
    mean_flux = mean(CH4_flux, na.rm = TRUE),
    se_flux = sd(CH4_flux, na.rm = TRUE) / sqrt(n()),  # Standard Error
    n = n()
  ) %>%
  mutate(
    Species.Label = paste0(Species.Code, " (n=", n, ")"),  # Append n to species name
    Species.Label = factor(Species.Label, levels = Species.Label[order(-mean_flux)])  # Order by mean_flux descending
  )

# Create the plot
mean_plot <- ggplot(summary_stats, aes(x = Species.Label, y = mean_flux * 1000, fill = mean_flux)) +
  geom_bar(stat = "identity", color = "black") +
  geom_errorbar(aes(ymin = (mean_flux - se_flux) * 1000, ymax = (mean_flux + se_flux) * 1000),
                width = 0.2) +
  scale_fill_viridis(option = "D", direction = -1) +  # Use viridis scale, dark to light
  theme_minimal() +
  labs(y = "Mean CH4 Flux (nmol CH4/m²/s)", x = "Species", title = "Mean CH4 Flux per Species (with SE)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+theme(legend.position = "none")

# Display the plot
print(mean_plot)

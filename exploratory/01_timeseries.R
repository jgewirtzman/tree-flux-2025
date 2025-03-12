library(lubridate)
library(tidyverse)

#load data
fluxes<-read.csv("data/raw/flux_dataset.csv")

#add date with year
fluxes$date <- as.Date(fluxes$jday, origin = paste0(fluxes$year, "-01-01"))

# Add the location column based on the plot values
fluxes <- fluxes %>%  filter(!is.na(Plot)) %>%
  mutate(location = ifelse(Plot == "BGS", "wetland", "upland"))

#linear scale, all fluxes
ggplot(fluxes, aes(x=date, y=(CH4_flux*1000), color=species))+
  geom_point()+
  geom_smooth(se=FALSE)+
  facet_wrap(~location)

#log scale, all fluxes
ggplot(fluxes, aes(x=date, y=log(CH4_flux*1000), color=species))+
  geom_point()+
  geom_smooth(se=FALSE)+
  facet_wrap(~location)

# Calculate mean and standard error (SE) by species and date
flux_summary <- fluxes %>%
  group_by(species, date, location) %>%
  summarise(
    mean_CH4 = mean(CH4_flux, na.rm = TRUE),  # Replace 'value' with your column of interest
    se_CH4 = sd(CH4_flux, na.rm = TRUE) / sqrt(n()) # Standard Error
  ) %>%
  ungroup()

# Ensure all species-plot combinations have a full date range
flux_summary <- flux_summary %>%
  complete(species, location, date = seq(min(date), max(date), by = "day")) %>%
  drop_na(mean_CH4)  # Remove rows where mean_value is NA

# Create the plot
ggplot(flux_summary, aes(x = date, y = (mean_CH4), color = species)) +
  geom_point() +
  #geom_line() +
  geom_smooth(se=F)+
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 0.5) +
  facet_wrap(~ location) +  # Facet by plot
  theme_classic()

# Create the plot
ggplot(flux_summary, aes(x = date, y = mean_CH4*1000, color = location)) +
  geom_point() +
  geom_line() +
  #geom_smooth(se=F)+
  geom_errorbar(aes(ymin = mean_CH4 - se_CH4, ymax = mean_CH4 + se_CH4), width = 0.5) +
  facet_wrap(~ species) +  # Facet by plot
  theme_classic()

#linear model
summary(lm(CH4_flux~species*jday*location, data=fluxes))



library(lme4)
library(emmeans)
library(readr)

# Fit a mixed-effects model: Flux as a function of species and location, with date as a random effect
model <- lmer(CH4_flux*1000 ~ species * location + (1 | date), data = fluxes, REML = TRUE)

# Get estimated marginal means (marginal mean flux per species per location)
marginal_means <- emmeans(model, ~ species | location)

# Print results
print(marginal_means)

# Convert to a dataframe if needed
marginal_means_df <- as.data.frame(marginal_means)




plot(fluxes$CH4_flux~fluxes$CO2_flux)

ggplot(fluxes, aes(x=CH4_flux))+
  geom_histogram()+facet_wrap(~location*species)


ggplot(fluxes, aes(x=CH4_flux))+
  geom_histogram()+facet_wrap(~location*species, scales="free_x")


ggplot(fluxes, aes(x=CO2_flux))+
  geom_histogram()+facet_wrap(~location*species)


# Load necessary libraries
library(dplyr)
library(glmmTMB)
library(emmeans)
library(ggplot2)

# Ensure correct transformation of data
fluxes <- fluxes %>%
  filter(!is.na(location)) %>%
  mutate(
    CH4_flux_scaled = CH4_flux * 1000  # Scaling flux as in your original model
  )

# Filter out non-positive flux values
fluxes <- fluxes %>%
  filter(CH4_flux * 1000 > 0)  # Keep only positive values

# Fit the Gamma model
gamma_model <- glmmTMB(CH4_flux * 1000 ~ species * location + (1 | date),
                       family = Gamma(link = "log"), data = fluxes)



# Fit a Gamma model with a log link
gamma_model <- glmmTMB(CH4_flux_scaled ~ species * location + (1 | date),
                       family = Gamma(link = "log"), data = fluxes)

# Get estimated marginal means (on response scale)
marginal_means <- emmeans(gamma_model, ~ species | location, type = "response")

# Convert to dataframe for visualization
marginal_means_df <- as.data.frame(marginal_means)

# Plot estimated marginal means with confidence intervals
ggplot(marginal_means_df, aes(x = species, y = response, fill = species)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  facet_wrap(~ location) +  # Separate by location
  labs(
    title = "Marginal Mean CH4 Flux by Species & Location",
    x = "Species",
    y = "Estimated CH4 Flux (scaled)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")  # Hide redundant legend






# Load necessary libraries
library(ggplot2)
library(dplyr)

# Calculate mean CH₄ flux per tree (assuming TreeID uniquely identifies each tree)
tree_means <- fluxes %>%
  group_by(species, location, TreeID) %>%
  summarise(mean_CH4_flux = mean(CH4_flux, na.rm = TRUE)) %>%
  ungroup()

# Create the plot
tree_means_plot<-ggplot(tree_means, aes(y = mean_CH4_flux * 1000)) +
  geom_boxplot(aes(x = 1), fill = "lightblue", alpha = 0.5, outlier.shape = NA) +  # Boxplot of means
  geom_jitter(aes(x = 1), color = "black", size = 2, width = 0.2) +  # Individual tree means as points
  facet_wrap(~ location + species, scales = "free_y") +  # Separate panels for each species-plot combo
  labs(
    title = "CH4 Flux Distribution of Tree Means by Species & Plot",
    x = NULL,  # No x-axis label since species names are in facet titles
    y = "Mean CH4 Flux (scaled)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 10, face = "bold"),  # Format facet titles
    axis.text.x = element_blank(),  # Remove x-axis text (species names are in facets)
    axis.ticks.x = element_blank()   # Remove x-axis ticks
  )
tree_means_plot



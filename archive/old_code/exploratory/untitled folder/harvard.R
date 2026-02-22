# Convert methane flux: μmol CH4 m-2 s-1 to mg C-CH4/m2/h
convert_ch4_flux <- function(flux_umol_m2_s) {
  # Conversion factor: 10^-6 mol/μmol * 12.01 g C/mol * 1000 mg/g * 3600 s/h
  flux_umol_m2_s * 43.236
}

# Prepare data with conversions and full species names
fluxes_converted <- fluxes %>%
  filter(!is.na(location)) %>%
  filter(CH4_flux > 0) %>%  # Keep only positive values for Gamma model
  mutate(
    CH4_flux_mg_c = convert_ch4_flux(CH4_flux),
    species_full = case_when(
      species == "bg" ~ "Nyssa sylvatica",
      species == "hem" ~ "Tsuga canadensis", 
      species == "rm" ~ "Acer rubrum",
      species == "ro" ~ "Quercus rubra",
      TRUE ~ species
    )
  )

# Fit the Gamma model with converted units
gamma_model <- glmmTMB(CH4_flux_mg_c ~ species_full * location + (1 | date),
                       family = Gamma(link = "log"), data = fluxes_converted)

# Get estimated marginal means (on response scale)
marginal_means <- emmeans(gamma_model, ~ species_full | location, type = "response")

# Convert to dataframe for visualization
marginal_means_df <- as.data.frame(marginal_means)

# Filter marginal means to only include species-location combinations present in data
marginal_means_df <- marginal_means_df %>%
  semi_join(fluxes_converted, by = c("species_full", "location"))

# Create faceted boxplot
marginal_means_boxplot <- ggplot() +
  # Add boxplots of raw data
  geom_boxplot(data = fluxes_converted, 
               aes(x = species_full, y = CH4_flux_mg_c, fill = species_full),
               alpha = 0.8, outlier.shape = NA) +
  # Add individual data points
  geom_jitter(data = fluxes_converted,
              aes(x = species_full, y = CH4_flux_mg_c, color = species_full),
              alpha = 0.5, size = 1, width = 0.2) +
  # Add marginal means as larger points
  geom_point(data = marginal_means_df,
             aes(x = species_full, y = response),
             size = 4, shape = 23, fill = "white", color = "black", stroke = 1.5) +
  # Add error bars for marginal means confidence intervals
  geom_errorbar(data = marginal_means_df,
                aes(x = species_full, ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.3, size = 1, color = "black") +
  facet_wrap(~ location, scales = "free_x") +  # Added scales = "free_x"
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 0.01),
    breaks = c(-0.1, -0.01, 0, 0.01, 0.1, 0.5, 1, 2, 5, 10),
    labels = c("-0.1", "-0.01", "0", "0.01", "0.1", "0.5", "1", "2", "5", "10")
  ) +
  scale_fill_manual(
    values = c(
      "Acer rubrum" = "#E69F00",        # Orange (maple - warm autumn)
      "Nyssa sylvatica" = "#56B4E9",    # Sky blue (black gum - wetland association)
      "Quercus rubra" = "#CC79A7",      # Pink/magenta (oak - distinctive)
      "Tsuga canadensis" = "#009E73"    # Teal green (hemlock - evergreen)
    )
  ) +
  scale_color_manual(
    values = c(
      "Acer rubrum" = "#B8780A",        # Darker orange
      "Nyssa sylvatica" = "#3A8BC2",    # Darker sky blue
      "Quercus rubra" = "#A85D7A",      # Darker pink/magenta
      "Tsuga canadensis" = "#006B4F"    # Darker teal green
    )
  ) +
  labs(
    y = expression("CH"[4] ~ "Flux (mg C-CH"[4] ~ m^-2 ~ h^-1 ~ ")"),
    x = NULL,
    fill = "Species",
    color = "Species"
  ) +
  scale_x_discrete(labels = function(x) gsub(" ", "\n", x)) +  # Split species names to two lines
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "italic", size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",  # Remove legend
    strip.text = element_text(size = 16, face = "bold"),
    plot.margin = margin(20, 20, 20, 20)
  )

marginal_means_boxplot
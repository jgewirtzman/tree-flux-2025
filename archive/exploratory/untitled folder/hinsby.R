# Convert methane flux: μmol CH4 m-2 s-1 to mg C-CH4/m2/h
convert_ch4_flux <- function(flux_umol_m2_s) {
  # Conversion factor: 10^-6 mol/μmol * 12.01 g C/mol * 1000 mg/g * 3600 s/h
  flux_umol_m2_s * 43.236
}

# Convert flux data and add full species names
flux_summary_converted <- flux_summary %>%
  mutate(
    mean_CH4_mg_c = convert_ch4_flux(mean_CH4),
    se_CH4_mg_c = convert_ch4_flux(se_CH4),
    species_full = case_when(
      species == "bg" ~ "Nyssa sylvatica",
      species == "hem" ~ "Tsuga canadensis", 
      species == "rm" ~ "Acer rubrum",
      species == "ro" ~ "Quercus rubra",
      TRUE ~ species
    )
  )

# Create the plot
ggplot(flux_summary_converted, aes(x = date, y = mean_CH4_mg_c, color = location)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = mean_CH4_mg_c - se_CH4_mg_c, ymax = mean_CH4_mg_c + se_CH4_mg_c), width = 0.5) +
  facet_wrap(~ species_full, scales = "free_y") +
  labs(y = expression("Mean CH"[4] ~ "Flux (mg C-CH"[4] ~ m^-2 ~ h^-1 ~ ")"),
       x = "Date",
       color = "Location") +
  theme_classic()
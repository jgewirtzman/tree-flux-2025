# Find correlations between tree flux data and met data over windows
library(tidyverse)
library(patchwork)  # For combining plots

# Load flux data 
stem_flux = read_csv("data/raw/flux_dataset.csv")

# Format columns
stem_flux$datetime = lubridate::round_date(lubridate::force_tz(as.POSIXct(stem_flux$real_start-2190),
                                                               tzone="EST"),"hour")
stem_flux$date = lubridate::date(stem_flux$datetime)
stem_flux$ID = stem_flux$Tree
stem_flux$site = ifelse(stem_flux$Plot=="BGS","BGS","EMS")

# Fix missing species
stem_flux$species[stem_flux$ID==288] = "hem"
stem_flux$species[stem_flux$ID==153] = "rm"
stem_flux$species[stem_flux$ID==414] = "bg"
stem_flux$species[stem_flux$ID==452] = "bg"

# Table with number of measurements x site x species
n_bySiteSpecies = stem_flux %>%
  group_by(site,species) %>%
  summarize(CH4_tree = median(CH4_flux,na.rm=T),
            count = n())

# Options to loop through for days time period
days = c(2/24, 6/24, 1:30)  # Reduced to 30 days for faster processing, expand as needed

# Define all meteorological variables to analyze
met_vars <- c(
  "tair_C", "P_mm", "RH", "VPD_kPa", "bgs_wtd_cm", "bvs_wtd_cm", 
  "PAR", "rnet", "slrr", "s10t", "p_kPa"
)

# Define variable attributes for plots
var_info <- tribble(
  ~variable,    ~label,                    ~color,
  "tair_C",     "Air temperature",         "red",
  "P_mm",       "Precipitation",           "skyblue",
  "RH",         "Relative humidity",       "darkgreen",
  "VPD_kPa",    "Vapor pressure deficit",  "darkred",
  "bgs_wtd_cm", "BGS water table depth",   "darkblue",
  "bvs_wtd_cm", "BVS water table depth",   "royalblue",
  "PAR",        "PAR",                     "orange",
  "rnet",       "Net radiation",           "forestgreen",
  "slrr",       "Solar radiation",         "gold",
  "s10t",       "Soil temperature 10cm",   "brown",
  "p_kPa",      "Atmospheric pressure",    "purple"
)

# Create empty data frame to store results
cor_all <- data.frame()

# Loop through each time window
for(d in 1:length(days)) {
  
  # hours to calculate rolling sum/mean
  window = 24*days[d] 
  
  # Calculate rolling window stats for all variables
  met_cul <- wtd_met %>%
    mutate(date = lubridate::date(datetime),
           month = lubridate::month(date),
           jday = lubridate::yday(date),
           jday_d = ifelse(month>=11, -1*(365-jday), jday),
           year = lubridate::year(date), 
           wyear = ifelse(month>11,year+1,year)) %>%
    arrange(wyear,jday_d)
  
  # Apply rolling calculations to each met variable
  for (var in met_vars) {
    # For precipitation, use sum; for others, use mean
    if (var == "P_mm") {
      met_cul[[paste0(var, "_sum")]] <- RcppRoll::roll_sum(
        met_cul[[var]], window, align = "right", na.rm = TRUE, fill = NA
      )
    } else {
      met_cul[[paste0(var, "_mn")]] <- RcppRoll::roll_mean(
        met_cul[[var]], window, align = "right", na.rm = TRUE, fill = NA
      )
      # Also calculate SD for water table depth variables
      if (grepl("wtd", var)) {
        met_cul[[paste0(var, "_sd")]] <- RcppRoll::roll_sd(
          met_cul[[var]], window, align = "right", na.rm = TRUE, fill = NA
        )
      }
    }
  }
  
  # Join flux and meteorological data
  stem_met <- left_join(stem_flux, met_cul)
  
  # Initialize data frame for correlation results for this window
  cor_results <- data.frame(
    site = c("BGS", "EMS"),
    interval = days[d]
  )
  
  # Calculate correlations for all variables
  for (site_name in c("BGS", "EMS")) {
    site_data <- stem_met %>% filter(site == site_name)
    
    for (var in met_vars) {
      # Define variable suffix based on variable type
      suffix <- if(var == "P_mm") "_sum" else "_mn"
      var_name <- paste0(var, suffix)
      
      # Calculate correlation with CH4 flux
      cor_test <- tryCatch({
        test <- cor.test(site_data$CH4_flux, site_data[[var_name]], method = "pearson")
        c(estimate = test$estimate, p_value = test$p.value)
      }, error = function(e) {
        c(estimate = NA, p_value = NA)
      })
      
      # Add to results
      col_name_cor <- paste0("cor_", var)
      col_name_p <- paste0("p_", var)
      
      cor_results[cor_results$site == site_name, col_name_cor] <- cor_test["estimate"]
      cor_results[cor_results$site == site_name, col_name_p] <- cor_test["p_value"]
      
      # Add SD correlations for water table depth variables
      if (grepl("wtd", var)) {
        var_sd <- paste0(var, "_sd")
        cor_test_sd <- tryCatch({
          test <- cor.test(site_data$CH4_flux, site_data[[var_sd]], method = "pearson")
          c(estimate = test$estimate, p_value = test$p.value)
        }, error = function(e) {
          c(estimate = NA, p_value = NA)
        })
        
        col_name_cor_sd <- paste0("cor_", var, "var")
        col_name_p_sd <- paste0("p_", var, "var")
        
        cor_results[cor_results$site == site_name, col_name_cor_sd] <- cor_test_sd["estimate"]
        cor_results[cor_results$site == site_name, col_name_p_sd] <- cor_test_sd["p_value"]
      }
    }
  }
  
  # Append to master results dataframe
  cor_all <- bind_rows(cor_all, cor_results)
}

# Convert to long format for easier plotting
cor_all_long <- cor_all %>%
  # Select only columns that start with "cor_"
  select(site, interval, starts_with("cor_")) %>%
  # Pivot to long format
  pivot_longer(
    cols = starts_with("cor_"),
    names_to = "variable",
    values_to = "correlation"
  ) %>%
  # Clean up variable names
  mutate(variable = gsub("cor_", "", variable))

# Match with p-values
p_values <- cor_all %>%
  select(site, interval, starts_with("p_")) %>%
  pivot_longer(
    cols = starts_with("p_"),
    names_to = "variable",
    values_to = "p_value"
  ) %>%
  mutate(variable = gsub("p_", "", variable))

# Join correlation and p-values
cor_all_long <- left_join(cor_all_long, p_values, 
                          by = c("site", "interval", "variable"))

# Add significance flag
cor_all_long <- cor_all_long %>%
  mutate(significant = p_value < 0.05)

# Create individual plots for each variable
plot_list <- list()

# Get unique variables
unique_vars <- unique(cor_all_long$variable)

# Generate plots for each variable
for (var in unique_vars) {
  # Get display name from var_info if available
  var_label <- var_info %>% 
    filter(variable == gsub("var$", "", var)) %>% 
    pull(label)
  
  if (length(var_label) == 0) var_label <- var
  
  # Create plot
  p <- cor_all_long %>%
    filter(variable == var) %>%
    ggplot(aes(x = interval, y = correlation, color = site)) +
    geom_line(size = 1) +
    geom_point(aes(shape = significant), size = 3) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1)) +
    labs(
      title = paste("Correlation with CH4 flux:", var_label),
      x = "Time window (days)",
      y = "Pearson correlation coefficient"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(size = 12)
    ) +
    scale_x_log10() +  # Log scale for x axis to better show short time windows
    coord_cartesian(ylim = c(-1, 1))  # Fixed y-axis range
  
  plot_list[[var]] <- p
}

# Combine all plots
combined_plots <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(
    title = "Correlation between CH4 flux and meteorological variables across time windows",
    subtitle = "Solid points indicate significant correlations (p < 0.05)",
    theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
  )

# Display plot
combined_plots

# Save plot
#ggsave("CH4_correlations_all_variables.pdf", combined_plots, width = 15, height = 15)

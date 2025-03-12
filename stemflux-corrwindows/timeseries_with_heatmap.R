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

# Create empty dataframe to collect all results
all_results <- data.frame()

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
  
  # Loop through each site
  for (site_name in c("BGS", "EMS")) {
    site_data <- stem_met %>% filter(site == site_name)
    
    # Loop through each variable to calculate correlations
    for (var in met_vars) {
      # Define variable suffix based on variable type
      suffix <- if(var == "P_mm") "_sum" else "_mn"
      var_name <- paste0(var, suffix)
      
      # Calculate correlation with CH4 flux - with proper error handling
      tryCatch({
        # Make sure we have valid data
        valid_data <- site_data %>% 
          select(CH4_flux, !!sym(var_name)) %>% 
          filter(!is.na(CH4_flux), !is.na(!!sym(var_name)))
        
        # Only proceed if we have enough data points
        if(nrow(valid_data) > 3) {
          test <- cor.test(valid_data$CH4_flux, valid_data[[var_name]], method = "pearson")
          
          # Store result in a row of data
          result <- data.frame(
            site = site_name,
            interval = days[d],
            variable = var,
            correlation = as.numeric(test$estimate),
            p_value = as.numeric(test$p.value),
            n_samples = nrow(valid_data)
          )
          
          # Add to results
          all_results <- bind_rows(all_results, result)
        }
      }, error = function(e) {
        message(paste("Error with", site_name, var, "at interval", days[d], ":", e$message))
      })
      
      # Handle SD for water table depth variables
      if (grepl("wtd", var)) {
        var_sd <- paste0(var, "_sd")
        
        tryCatch({
          # Make sure we have valid data
          valid_data <- site_data %>% 
            select(CH4_flux, !!sym(var_sd)) %>% 
            filter(!is.na(CH4_flux), !is.na(!!sym(var_sd)))
          
          # Only proceed if we have enough data points
          if(nrow(valid_data) > 3) {
            test <- cor.test(valid_data$CH4_flux, valid_data[[var_sd]], method = "pearson")
            
            # Store result in a row of data
            result <- data.frame(
              site = site_name,
              interval = days[d],
              variable = paste0(var, "var"),
              correlation = as.numeric(test$estimate),
              p_value = as.numeric(test$p.value),
              n_samples = nrow(valid_data)
            )
            
            # Add to results
            all_results <- bind_rows(all_results, result)
          }
        }, error = function(e) {
          message(paste("Error with", site_name, var_sd, "at interval", days[d], ":", e$message))
        })
      }
    }
  }
  
  # Print progress
  cat("Completed interval", days[d], "\n")
}

# Check data structure after collection
print(paste("Total rows in results:", nrow(all_results)))
print(paste("Number of NAs in correlation:", sum(is.na(all_results$correlation))))

# Add significance flag
all_results <- all_results %>%
  mutate(significant = p_value < 0.05)

# Create individual plots for each variable
plot_list <- list()

# Get unique variables
unique_vars <- unique(all_results$variable)

# Generate plots for each variable
for (var in unique_vars) {
  # Get display name from var_info if available
  var_label <- var_info %>% 
    filter(variable == gsub("var$", "", var)) %>% 
    pull(label)
  
  if (length(var_label) == 0) var_label <- var
  
  var_data <- all_results %>% filter(variable == var)
  
  # Skip if no data
  if (nrow(var_data) == 0) next
  
  # Create plot
  p <- ggplot(var_data, aes(x = interval, y = correlation, color = site)) +
    geom_line(size = 1) +
    geom_point(aes(shape = significant), size = 3) +
    scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                       name = "p < 0.05") +
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

# Combine all plots if we have any
if (length(plot_list) > 0) {
  combined_plots <- wrap_plots(plot_list, ncol = 3) +
    plot_annotation(
      title = "Correlation between CH4 flux and meteorological variables across time windows",
      subtitle = "Solid points indicate significant correlations (p < 0.05)",
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5))
    )
  
  # Display plot
  print(combined_plots)
  
  # Save plot
  ggsave("CH4_correlations_all_variables.pdf", combined_plots, width = 15, height = 15)
} else {
  print("No plots were created. Check your data.")
}

# Also create a heatmap visualization for an overview
# Filter out variables with "var" suffix for cleaner heatmap
heatmap_data <- all_results %>%
  filter(!grepl("var$", variable)) %>%
  mutate(
    # Create readable variable names
    variable_clean = gsub("_mn$|_sum$", "", variable)
  )

# Join with variable info for better labels
heatmap_data <- heatmap_data %>%
  left_join(var_info, by = c("variable_clean" = "variable")) %>%
  mutate(variable_label = ifelse(is.na(label), variable_clean, label))

# Create heatmap
if (nrow(heatmap_data) > 0) {
  heatmap_plot <- ggplot(heatmap_data, aes(x = factor(interval), y = variable_label)) +
    facet_wrap(~ site) +
    geom_tile(aes(fill = correlation)) +
    geom_point(data = heatmap_data %>% filter(significant), 
               size = 1, color = "black") +
    scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", 
      midpoint = 0, limits = c(-1, 1),
      name = "Correlation"
    ) +
    labs(
      title = "Correlation between CH4 flux and meteorological variables",
      x = "Time window (days)",
      y = "Variable"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = element_blank()
    )
  
  # Display heatmap
  print(heatmap_plot)
  
  # Save heatmap
  ggsave("CH4_correlations_heatmap.pdf", heatmap_plot, width = 12, height = 8)
} else {
  print("No data available for heatmap. Check your data.")
}

# Print summary of correlation significance by variable
all_results %>%
  group_by(site, variable) %>%
  summarize(
    total_intervals = n(),
    significant_intervals = sum(significant),
    percent_significant = round(100 * significant_intervals / total_intervals, 1)
  ) %>%
  arrange(site, desc(percent_significant)) %>%
  print(n = 100)
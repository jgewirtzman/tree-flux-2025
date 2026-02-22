# Fixed Comprehensive AmeriFlux Variables Rolling Window Correlation Analysis
library(tidyverse)
library(patchwork)
library(RcppRoll)
library(lubridate)
library(googledrive)

# Download the AmeriFlux file using the direct file ID
drive_download("https://drive.google.com/uc?id=1NxqW8wRYJll8daI1CHIMA036EhazG9cb", 
               path = "data/raw/AMF_US-Ha1_BASE_HR_25-5_2324.csv",
               overwrite = TRUE)

# Load flux data 
stem_flux = read_csv("data/raw/flux_dataset.csv")

# Load AmeriFlux data
ameriflux_data = read_csv("data/raw/AMF_US-Ha1_BASE_HR_25-5_2324.csv")

cat("AmeriFlux data loaded successfully!\n")
cat("Dimensions:", nrow(ameriflux_data), "rows x", ncol(ameriflux_data), "columns\n")

# Format stem flux columns (same as before)
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

# Convert AmeriFlux timestamps - handle scientific notation
ameriflux_data <- ameriflux_data %>%
  mutate(
    # Convert scientific notation to character, then parse
    timestamp_char = sprintf("%.0f", TIMESTAMP_START),
    datetime = as.POSIXct(timestamp_char, format = "%Y%m%d%H%M", tz = "EST")
  ) %>%
  select(-timestamp_char)  # Remove temporary column

cat("AmeriFlux datetime range:", as.character(range(ameriflux_data$datetime, na.rm = TRUE)), "\n")
cat("Stem flux datetime range:", as.character(range(stem_flux$datetime, na.rm = TRUE)), "\n")

# Test join to see overlap
test_join <- left_join(stem_flux, ameriflux_data, by = "datetime")
matches <- sum(!is.na(test_join$TIMESTAMP_START))
cat("Number of successful datetime matches:", matches, "\n")
cat("Percentage of stem flux data matched:", round(100 * matches / nrow(stem_flux), 1), "%\n")

if(matches == 0) {
  stop("No datetime matches found! Please check the datetime alignment.")
}

# Define all AmeriFlux variables to analyze (excluding timestamps and datetime)
ameriflux_vars <- names(ameriflux_data)[!names(ameriflux_data) %in% 
                                          c("TIMESTAMP_START", "TIMESTAMP_END", "datetime")]

# Remove variables that are entirely missing (-9999) to speed up analysis
vars_with_data <- sapply(ameriflux_vars, function(var) {
  sum(ameriflux_data[[var]] != -9999, na.rm = TRUE)
})

# Keep only variables with at least some real data
useful_vars <- names(vars_with_data)[vars_with_data > 100]  # At least 100 non-missing values
cat("Variables with sufficient data:", length(useful_vars), "out of", length(ameriflux_vars), "\n")

print("Variables to analyze:")
print(useful_vars)

# Time windows (same as before)
days = c(2/24, 6/24, 12/24, 1:30)

# Create empty dataframe to collect all results
all_results <- data.frame()

# Progress tracking
total_combinations <- length(days) * 2 * length(useful_vars)
current_combination <- 0

# Loop through each time window
for(d in 1:length(days)) {
  window = 24*days[d] 
  cat("\nProcessing time window:", days[d], "days (", window, "hours)\n")
  
  # Calculate rolling window stats for AmeriFlux variables
  ameriflux_rolling <- ameriflux_data %>%
    mutate(date = lubridate::date(datetime),
           month = lubridate::month(date),
           jday = lubridate::yday(date),
           jday_d = ifelse(month>=11, -1*(365-jday), jday),
           year = lubridate::year(date), 
           wyear = ifelse(month>11,year+1,year)) %>%
    arrange(wyear,jday_d)
  
  # Replace -9999 with NA for all variables
  for(var in useful_vars) {
    ameriflux_rolling[[var]][ameriflux_rolling[[var]] == -9999] <- NA
  }
  
  # Apply rolling calculations to each AmeriFlux variable
  for (var in useful_vars) {
    # Skip if variable is all NA after cleaning
    if(all(is.na(ameriflux_rolling[[var]]))) {
      next
    }
    
    # For precipitation (P), use sum; for others, use mean
    if (var == "P") {
      ameriflux_rolling[[paste0(var, "_sum")]] <- RcppRoll::roll_sum(
        ameriflux_rolling[[var]], window, align = "right", na.rm = TRUE, fill = NA
      )
    } else {
      ameriflux_rolling[[paste0(var, "_mn")]] <- RcppRoll::roll_mean(
        ameriflux_rolling[[var]], window, align = "right", na.rm = TRUE, fill = NA
      )
      
      # Calculate SD for soil variables (SWC, TS) that might show important variability
      if (grepl("^(SWC|TS)_", var)) {
        ameriflux_rolling[[paste0(var, "_sd")]] <- RcppRoll::roll_sd(
          ameriflux_rolling[[var]], window, align = "right", na.rm = TRUE, fill = NA
        )
      }
    }
  }
  
  # Join flux and AmeriFlux data
  stem_ameriflux <- left_join(stem_flux, ameriflux_rolling, by = "datetime")
  
  # Loop through each site
  for (site_name in c("BGS", "EMS")) {
    site_data <- stem_ameriflux %>% filter(site == site_name)
    
    # Loop through each variable to calculate correlations
    for (var in useful_vars) {
      current_combination <- current_combination + 1
      
      # Progress update every 50 combinations
      if (current_combination %% 50 == 0) {
        cat("Progress:", round(100 * current_combination / total_combinations, 1), "%\n")
      }
      
      # Skip if variable is all NA
      if(all(is.na(site_data[[var]]))) {
        next
      }
      
      # Define variable suffix based on variable type
      suffix <- if(var == "P") "_sum" else "_mn"
      var_name <- paste0(var, suffix)
      
      # Calculate correlation with CH4 flux
      tryCatch({
        valid_data <- site_data %>% 
          select(CH4_flux, all_of(var_name)) %>% 
          filter(!is.na(CH4_flux), !is.na(.data[[var_name]]))
        
        if(nrow(valid_data) > 3) {
          test <- cor.test(valid_data$CH4_flux, valid_data[[var_name]], method = "pearson")
          
          result <- data.frame(
            site = site_name,
            interval = days[d],
            variable = var,
            correlation = as.numeric(test$estimate),
            p_value = as.numeric(test$p.value),
            n_samples = nrow(valid_data),
            variable_type = "mean"
          )
          
          all_results <- bind_rows(all_results, result)
        }
      }, error = function(e) {
        # Silent error handling
      })
      
      # Handle SD for soil variables
      if (grepl("^(SWC|TS)_", var)) {
        var_sd <- paste0(var, "_sd")
        
        if(var_sd %in% names(site_data)) {
          tryCatch({
            valid_data <- site_data %>% 
              select(CH4_flux, all_of(var_sd)) %>% 
              filter(!is.na(CH4_flux), !is.na(.data[[var_sd]]))
            
            if(nrow(valid_data) > 3) {
              test <- cor.test(valid_data$CH4_flux, valid_data[[var_sd]], method = "pearson")
              
              result <- data.frame(
                site = site_name,
                interval = days[d],
                variable = paste0(var, "_var"),
                correlation = as.numeric(test$estimate),
                p_value = as.numeric(test$p.value),
                n_samples = nrow(valid_data),
                variable_type = "variability"
              )
              
              all_results <- bind_rows(all_results, result)
            }
          }, error = function(e) {
            # Silent error handling
          })
        }
      }
    }
  }
  
  cat("Completed interval", days[d], "days\n")
}

# Check results
cat("\nAnalysis complete!\n")
cat("Total rows in results:", nrow(all_results), "\n")
cat("Number of unique variables:", length(unique(all_results$variable)), "\n")

if(nrow(all_results) == 0) {
  stop("No correlations calculated! Check data alignment and variable availability.")
}

# Add significance flag
all_results <- all_results %>%
  mutate(significant = p_value < 0.05)

# Create variable categories for better visualization
all_results <- all_results %>%
  mutate(
    category = case_when(
      grepl("^TS_", variable) | grepl("^T_SONIC", variable) ~ "Temperature",
      grepl("^SWC_", variable) ~ "Soil Water Content", 
      grepl("^TA_", variable) ~ "Air Temperature",
      grepl("^RH_", variable) ~ "Relative Humidity",
      grepl("^PPFD_", variable) ~ "Radiation",
      grepl("^(FC|CO2)_", variable) ~ "Carbon Fluxes",
      grepl("^(LE|FH2O|H2O)_", variable) ~ "Water/Energy Fluxes",
      grepl("^(WS|WD|USTAR|TAU)", variable) ~ "Wind/Turbulence",
      grepl("^(PA|P)$", variable) ~ "Atmospheric/Precip",
      grepl("^(H|NETRAD)", variable) ~ "Energy Balance",
      grepl("_SIGMA_", variable) ~ "Turbulence Statistics",
      TRUE ~ "Other"
    )
  )

# Print summary by category
cat("\nVariables by category:\n")
category_summary <- all_results %>%
  group_by(category) %>%
  summarise(n_variables = length(unique(variable))) %>%
  arrange(desc(n_variables))
print(category_summary)

# Create comprehensive heatmap
heatmap_data <- all_results %>%
  filter(interval <= 20) %>%  # Focus on shorter intervals for visibility
  mutate(
    interval_factor = factor(round(interval, 2)),
    correlation_binned = pmax(-1, pmin(1, correlation))
  )

# Create main heatmap
comprehensive_heatmap <- ggplot(heatmap_data, 
                                aes(x = interval_factor, y = variable)) +
  facet_grid(category ~ site, scales = "free_y", space = "free_y") +
  geom_tile(aes(fill = correlation_binned)) +
  geom_point(data = heatmap_data %>% filter(significant), 
             size = 0.5, color = "white", alpha = 0.8) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red", 
    midpoint = 0, limits = c(-1, 1),
    name = "Correlation"
  ) +
  labs(
    title = "Comprehensive CH4 Flux Correlations with AmeriFlux Variables",
    subtitle = "White dots indicate significant correlations (p < 0.05)",
    x = "Time window (days)",
    y = "Variable"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 7),
    panel.grid = element_blank(),
    strip.text.y = element_text(angle = 0, size = 9),
    strip.text.x = element_text(size = 10),
    legend.position = "bottom"
  )

# Display the heatmap
print(comprehensive_heatmap)

# Save the heatmap
ggsave("comprehensive_CH4_correlations_heatmap.pdf", comprehensive_heatmap, 
       width = 18, height = 20, limitsize = FALSE)

# Identify top correlations overall
top_correlations <- all_results %>%
  filter(!is.na(correlation)) %>%
  mutate(abs_correlation = abs(correlation)) %>%
  arrange(desc(abs_correlation)) %>%
  head(20)

cat("\nTop 20 strongest correlations:\n")
print(top_correlations %>% 
        select(site, variable, interval, correlation, p_value, significant), n = 20)

# Summary statistics by category
summary_stats <- all_results %>%
  group_by(category, site) %>%
  summarize(
    n_variables = length(unique(variable)),
    mean_abs_correlation = round(mean(abs(correlation), na.rm = TRUE), 3),
    max_abs_correlation = round(max(abs(correlation), na.rm = TRUE), 3),
    pct_significant = round(100 * mean(significant, na.rm = TRUE), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(max_abs_correlation))

cat("\nSummary by category and site:\n")
print(summary_stats, n = 50)

# Save all results
write_csv(all_results, "comprehensive_ameriflux_correlations_results.csv")

cat("\nAnalysis complete! Results saved.\n")
cat("Heatmap saved as: comprehensive_CH4_correlations_heatmap.pdf\n")
cat("Data saved as: comprehensive_ameriflux_correlations_results.csv\n")
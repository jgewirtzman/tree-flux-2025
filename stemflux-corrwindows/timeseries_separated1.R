# Modified code to create grouped plots with linear x-axis
library(tidyverse)
library(patchwork)  # For combining plots

# Create individual plots for each variable, but organize them in groups
# Assuming all_results data frame already contains your correlation results

# Define variable groups for better organization
var_groups <- list(
  "Temperature" = c("tair_C", "s10t"),
  "Water" = c("bgs_wtd_cm", "bvs_wtd_cm", "P_mm"),
  "Humidity" = c("RH", "VPD_kPa"),
  "Radiation" = c("PAR", "rnet", "slrr"),
  "Pressure" = c("p_kPa")
)

# Create reversed lookup from variable to group
var_to_group <- list()
for (group_name in names(var_groups)) {
  for (var in var_groups[[group_name]]) {
    var_to_group[[var]] <- group_name
  }
}

# Add group info to the results
all_results <- all_results %>%
  mutate(
    # Extract the base variable name by removing suffixes
    base_var = gsub("var$", "", variable),
    base_var = gsub("_mn$|_sum$", "", base_var),
    # Assign group
    group = sapply(base_var, function(v) {
      if (v %in% names(var_to_group)) {
        return(var_to_group[[v]])
      } else {
        return("Other")
      }
    })
  )

# Function to create plots for a group of variables
create_group_plot <- function(group_name, group_data) {
  # Get unique variables in this group
  group_vars <- unique(group_data$variable)
  
  # Generate plots for each variable in the group
  plot_list <- list()
  
  for (var in group_vars) {
    # Get display name from var_info if available
    var_label <- var_info %>% 
      filter(variable == gsub("var$", "", gsub("_mn$|_sum$", "", var))) %>% 
      pull(label)
    
    if (length(var_label) == 0) var_label <- var
    
    var_data <- group_data %>% filter(variable == var)
    
    # Skip if no data
    if (nrow(var_data) == 0) next
    
    # Create plot with linear scale
    p <- ggplot(var_data, aes(x = interval, y = correlation, color = site)) +
      geom_line(size = 1) +
      geom_point(aes(shape = significant), size = 3) +
      scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1),
                         name = "p < 0.05") +
      labs(
        title = paste(var_label),
        x = "Time window (days)",
        y = "Correlation"
      ) +
      theme_bw() +
      theme(
        legend.position = "bottom",
        plot.title = element_text(size = 12)
      ) +
      scale_x_continuous() +  # Linear scale instead of log
      coord_cartesian(ylim = c(-1, 1))  # Fixed y-axis range
    
    plot_list[[var]] <- p
  }
  
  # Combine all plots in this group
  if (length(plot_list) > 0) {
    combined_plots <- wrap_plots(plot_list, ncol = min(length(plot_list), 2)) +
      plot_annotation(
        title = paste("Group:", group_name),
        theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
      )
    
    return(combined_plots)
  } else {
    return(NULL)
  }
}

# Create plots for each group
for (group_name in unique(c(names(var_groups), "Other"))) {
  group_data <- all_results %>% filter(group == group_name)
  
  if (nrow(group_data) > 0) {
    group_plot <- create_group_plot(group_name, group_data)
    
    if (!is.null(group_plot)) {
      # Display plot
      print(group_plot)
      
      # Save plot
      ggsave(paste0("CH4_correlations_", gsub(" ", "_", tolower(group_name)), ".pdf"), 
             group_plot, width = 10, height = 8)
    }
  }
}

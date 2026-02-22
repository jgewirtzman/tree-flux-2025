# ============================================================
# 01_prepare_data.R
#
# Prepare flux and environmental data for analysis.
# - Apply unit corrections to flux data
# - Create asinh-transformed CH4 response
# - Join flux with environmental data
#
# Inputs:
#   - data/raw/HF_2023-2025_tree_flux.csv (or similar)
#   - data/processed/aligned_hourly_dataset.csv
#
# Outputs:
#   - data/processed/flux_cleaned.csv
#   - data/processed/flux_with_env.csv
#   - figures/01_data_overview.png
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})

# ============================================================
# CONFIGURATION
# ============================================================

PATHS <- list(
  flux_raw = "data/HF_2023-2025_tree_flux.csv",
  aligned = "data/processed/aligned_hourly_dataset.csv"
)

OUTPUT_DIR <- "figures/01_prepare_data"
DATA_DIR <- "data/processed"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

# Date range
DATE_MIN <- as.Date("2023-01-01")
DATE_MAX <- as.Date("2026-01-01")

# ============================================================
# LOAD AND CLEAN FLUX DATA
# ============================================================

message("Loading flux data...")

flux_raw <- read_csv(PATHS$flux_raw, show_col_types = FALSE)

message("  Raw observations: ", nrow(flux_raw))

# Apply corrections and create cleaned dataset
flux_cleaned <- flux_raw %>%
  # Fill missing species within tree
  group_by(Tree) %>%
  mutate(SPECIES = ifelse(is.na(SPECIES), first(na.omit(SPECIES)), SPECIES)) %>%
  ungroup() %>%
  # Parse datetime
  mutate(
    datetime = as.POSIXct(datetime_posx, tz = "UTC"),
    datetime = round_date(datetime, "hour"),
    date = as.Date(datetime),
    year = year(datetime),
    month = month(datetime),
    doy = yday(datetime),
    hour = hour(datetime)
  ) %>%
  # Filter date range and valid plots
  filter(
    date >= DATE_MIN, 
    date <= DATE_MAX,
    !is.na(PLOT)
  ) %>%
  # Create site variable
  mutate(
    site = factor(
      ifelse(PLOT == "BGS", "Wetland", "Upland"),
      levels = c("Wetland", "Upland")
    ),
    Tree = as.factor(Tree),
    species = as.factor(SPECIES)
  ) %>%
  # Apply unit correction to pre-2025 CH4 SE
  mutate(
    CH4_SE_corrected = ifelse(year < 2025, CH4_SE * 1000, CH4_SE)
  ) %>%
  # Create transformed response
  mutate(
    CH4_flux_asinh = asinh(CH4_flux_nmolpm2ps)
  ) %>%
  # Select key columns
  dplyr::select(
    datetime, date, year, month, doy, hour,
    site, Tree, species,
    CH4_flux_nmolpm2ps, CH4_flux_asinh,
    CH4_r2, CH4_SE = CH4_SE_corrected,
    CO2_flux_umolpm2ps, CO2_r2, CO2_SE
  )

message("  Cleaned observations: ", nrow(flux_cleaned))
message("  Date range: ", min(flux_cleaned$date), " to ", max(flux_cleaned$date))

# Summary by site
site_summary <- flux_cleaned %>%
  group_by(site) %>%
  summarize(
    n_obs = n(),
    n_trees = n_distinct(Tree),
    n_species = n_distinct(species),
    species = paste(unique(species), collapse = ", "),
    .groups = "drop"
  )

cat("\nSite summary:\n")
print(site_summary)

# Save cleaned flux data
write_csv(flux_cleaned, file.path(DATA_DIR, "flux_cleaned.csv"))
message("\nSaved: ", file.path(DATA_DIR, "flux_cleaned.csv"))

# ============================================================
# LOAD ENVIRONMENTAL DATA
# ============================================================

message("\nLoading environmental data...")

env_data <- read_csv(PATHS$aligned, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC")) %>%
  filter(as.Date(datetime) >= DATE_MIN, as.Date(datetime) <= DATE_MAX) %>%
  arrange(datetime)

message("  Environmental records: ", nrow(env_data))
message("  Variables: ", ncol(env_data) - 1)

# ============================================================
# JOIN FLUX WITH ENVIRONMENTAL DATA
# ============================================================

message("\nJoining flux with environmental data...")

flux_with_env <- flux_cleaned %>%
  left_join(env_data, by = "datetime")

# Check join success
n_with_env <- sum(!is.na(flux_with_env$TS_Ha2))
message("  Observations with environmental data: ", n_with_env, 
        " (", round(100 * n_with_env / nrow(flux_with_env), 1), "%)")

write_csv(flux_with_env, file.path(DATA_DIR, "flux_with_env.csv"))
message("Saved: ", file.path(DATA_DIR, "flux_with_env.csv"))

# ============================================================
# DATA OVERVIEW FIGURE
# ============================================================

message("\nCreating overview figure...")

# Temporal distribution
p1 <- flux_cleaned %>%
  count(date, site) %>%
  ggplot(aes(x = date, y = n, fill = site)) +
  geom_col(width = 1) +
  scale_fill_manual(values = c("Wetland" = "#1B9E77", "Upland" = "#D95F02")) +
  labs(title = "Observations per Day", x = NULL, y = "Count", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

# CH4 distribution
p2 <- flux_cleaned %>%
  ggplot(aes(x = CH4_flux_asinh, fill = site)) +
  geom_histogram(bins = 50, alpha = 0.7, position = "identity") +
  scale_fill_manual(values = c("Wetland" = "#1B9E77", "Upland" = "#D95F02")) +
  labs(title = "CH4 Flux Distribution", x = "asinh(CH4 flux)", y = "Count", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom")

# By species
p3 <- flux_cleaned %>%
  ggplot(aes(x = species, y = CH4_flux_asinh, fill = site)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = c("Wetland" = "#1B9E77", "Upland" = "#D95F02")) +
  labs(title = "CH4 by Species", x = NULL, y = "asinh(CH4 flux)", fill = NULL) +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

# Combine
p_overview <- p1 / (p2 | p3) +
  plot_annotation(title = "Data Overview")

# Need patchwork
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_overview <- p1 / (p2 | p3) +
    plot_annotation(title = "Data Overview")
  
  ggsave(file.path(OUTPUT_DIR, "data_overview.png"), p_overview,
         width = 12, height = 10, dpi = 300)
  message("Saved: ", file.path(OUTPUT_DIR, "data_overview.png"))
} else {
  ggsave(file.path(OUTPUT_DIR, "temporal_distribution.png"), p1, width = 10, height = 4, dpi = 300)
  ggsave(file.path(OUTPUT_DIR, "flux_distribution.png"), p2, width = 6, height = 4, dpi = 300)
  message("Saved individual plots (install patchwork for combined figure)")
}

# ============================================================
# SUMMARY
# ============================================================

message("\n", paste(rep("=", 50), collapse = ""))
message("DATA PREPARATION COMPLETE")
message(paste(rep("=", 50), collapse = ""))

cat("\nFlux data:\n")
cat("  Total observations:", nrow(flux_cleaned), "\n")
cat("  Wetland:", sum(flux_cleaned$site == "Wetland"), "\n")
cat("  Upland:", sum(flux_cleaned$site == "Upland"), "\n")
cat("  Trees:", n_distinct(flux_cleaned$Tree), "\n")
cat("  Species:", n_distinct(flux_cleaned$species), "\n")

cat("\nOutput files:\n")
cat("  data/processed/flux_cleaned.csv\n")
cat("  data/processed/flux_with_env.csv\n")

message("\nDone!")
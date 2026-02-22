# ============================================================
# 07c_correlation_matrix.R
# Generate correlation matrix for aligned dataset variables
# ============================================================

library(tidyverse)
library(lubridate)
library(corrplot)

# ============================================================
# CONFIGURATION
# ============================================================

DATA_PATH <- "data/processed/aligned_hourly_dataset.csv"
OUTPUT_DIR <- "figures"

# Date range (same as timeseries plots)
DATE_MIN <- as.Date("2023-01-01")
DATE_MAX <- as.Date("2026-01-01")

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# LOAD DATA
# ============================================================

message("Loading aligned dataset...")
aligned_data <- read_csv(DATA_PATH, show_col_types = FALSE) %>%
  mutate(datetime = as.POSIXct(datetime, tz = "UTC")) %>%
  filter(as.Date(datetime) >= DATE_MIN, as.Date(datetime) <= DATE_MAX)

message("  ", nrow(aligned_data), " rows after date filter")

# ============================================================
# SELECT VARIABLES FOR CORRELATION
# ============================================================

# Exclude time columns
time_cols <- c("datetime", "date", "year", "month", "doy", "hour", "wyear")
all_vars <- names(aligned_data)[!names(aligned_data) %in% time_cols]
all_vars <- all_vars[sapply(aligned_data[all_vars], is.numeric)]

message("  ", length(all_vars), " numeric variables")

# Get data matrix
cor_data <- aligned_data %>%
  select(all_of(all_vars))

# ============================================================
# COMPUTE CORRELATION MATRIX
# ============================================================

message("\nComputing correlation matrix...")

# Use pairwise complete observations
cor_matrix <- cor(cor_data, use = "pairwise.complete.obs", method = "pearson")

# Count valid pairs for each correlation
n_obs <- matrix(NA, nrow = ncol(cor_data), ncol = ncol(cor_data))
rownames(n_obs) <- names(cor_data)
colnames(n_obs) <- names(cor_data)

for (i in seq_along(names(cor_data))) {
  for (j in seq_along(names(cor_data))) {
    n_obs[i, j] <- sum(complete.cases(cor_data[, c(i, j)]))
  }
}

# ============================================================
# PLOT 1: FULL CORRELATION MATRIX
# ============================================================

message("\nGenerating correlation plots...")

# Full matrix heatmap
png(file.path(OUTPUT_DIR, "correlation_matrix_full.png"), 
    width = 14, height = 12, units = "in", res = 150)

corrplot(cor_matrix, 
         method = "color",
         type = "full",
         order = "hclust",  # Hierarchical clustering
         hclust.method = "ward.D2",
         addrect = 5,       # Add rectangles around clusters
         tl.col = "black",
         tl.cex = 0.6,
         tl.srt = 45,
         col = colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0", 
                                  "#FFFFFF", 
                                  "#FDDBC7", "#EF8A62", "#B2182B"))(100),
         title = "Correlation Matrix (Hierarchically Clustered)",
         mar = c(0, 0, 2, 0))

dev.off()
message("  Saved: correlation_matrix_full.png")

# ============================================================
# PLOT 2: LOWER TRIANGLE WITH COEFFICIENTS
# ============================================================

png(file.path(OUTPUT_DIR, "correlation_matrix_lower.png"), 
    width = 14, height = 12, units = "in", res = 150)

corrplot(cor_matrix, 
         method = "color",
         type = "lower",
         order = "hclust",
         hclust.method = "ward.D2",
         addCoef.col = "black",
         number.cex = 0.4,
         tl.col = "black",
         tl.cex = 0.6,
         tl.srt = 45,
         col = colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0", 
                                  "#FFFFFF", 
                                  "#FDDBC7", "#EF8A62", "#B2182B"))(100),
         title = "Correlation Matrix with Coefficients",
         mar = c(0, 0, 2, 0))

dev.off()
message("  Saved: correlation_matrix_lower.png")

# ============================================================
# PLOT 3: GROUPED BY VARIABLE TYPE
# ============================================================

# Define variable groups for ordering
var_order <- c(
  # Fluxes
  all_vars[grepl("^FC_", all_vars)],
  all_vars[grepl("^SC_", all_vars)],
  all_vars[grepl("^LE_", all_vars)],
  all_vars[grepl("^H_", all_vars)],
  all_vars[grepl("^G_", all_vars)],
  # Gas mixing ratios
  all_vars[grepl("^CH4_MR", all_vars)],
  all_vars[grepl("^CO2_MR", all_vars)],
  # Radiation
  all_vars[grepl("^PAR$|^slrr$|^rnet$", all_vars)],
  # Temperature
  all_vars[grepl("^tair|^T_CANOPY|^s10t|^TS_", all_vars)],
  # Atmospheric moisture
  all_vars[grepl("^RH$|^VPD", all_vars)],
  # Soil moisture
  all_vars[grepl("^NEON_SWC|^SWC_", all_vars)],
  all_vars[grepl("wtd", all_vars, ignore.case = TRUE)],
  # Precipitation
  all_vars[grepl("^P_mm$|^THROUGHFALL", all_vars)],
  # Pressure
  all_vars[grepl("^p_kPa$", all_vars)],
  # Wind
  all_vars[grepl("^WS_|^WD_|^USTAR", all_vars)],
  # Phenology
  all_vars[grepl("^gcc|^ndvi", all_vars)]
)

# Remove duplicates and keep only variables that exist
var_order <- unique(var_order)
var_order <- var_order[var_order %in% all_vars]

# Add any missing variables at the end
missing_vars <- setdiff(all_vars, var_order)
var_order <- c(var_order, missing_vars)

# Reorder correlation matrix
cor_matrix_ordered <- cor_matrix[var_order, var_order]

png(file.path(OUTPUT_DIR, "correlation_matrix_grouped.png"), 
    width = 14, height = 12, units = "in", res = 150)

corrplot(cor_matrix_ordered, 
         method = "color",
         type = "full",
         order = "original",  # Keep our custom order
         tl.col = "black",
         tl.cex = 0.6,
         tl.srt = 45,
         col = colorRampPalette(c("#2166AC", "#67A9CF", "#D1E5F0", 
                                  "#FFFFFF", 
                                  "#FDDBC7", "#EF8A62", "#B2182B"))(100),
         title = "Correlation Matrix (Grouped by Variable Type)",
         mar = c(0, 0, 2, 0))

dev.off()
message("  Saved: correlation_matrix_grouped.png")

# ============================================================
# SAVE CORRELATION VALUES
# ============================================================

# Convert to long format for easier inspection
cor_long <- as.data.frame(as.table(cor_matrix)) %>%
  rename(var1 = Var1, var2 = Var2, correlation = Freq) %>%
  filter(var1 != var2) %>%
  mutate(abs_cor = abs(correlation)) %>%
  arrange(desc(abs_cor))

# Add sample size
cor_long$n_obs <- mapply(function(v1, v2) n_obs[v1, v2], 
                         cor_long$var1, cor_long$var2)

# Remove duplicate pairs (keep only upper triangle)
cor_long <- cor_long %>%
  mutate(pair = paste(pmin(as.character(var1), as.character(var2)),
                      pmax(as.character(var1), as.character(var2)), sep = "_")) %>%
  distinct(pair, .keep_all = TRUE) %>%
  select(-pair)

write_csv(cor_long, file.path(OUTPUT_DIR, "correlation_values.csv"))
message("  Saved: correlation_values.csv")

# ============================================================
# PRINT TOP CORRELATIONS
# ============================================================

message("\n============================================================")
message("TOP 20 STRONGEST CORRELATIONS")
message("============================================================")

top_cors <- head(cor_long, 20)
for (i in seq_len(nrow(top_cors))) {
  message(sprintf("  %s <-> %s: r = %.3f (n = %d)",
                  top_cors$var1[i], top_cors$var2[i], 
                  top_cors$correlation[i], top_cors$n_obs[i]))
}

message("\n============================================================")
message("TOP 10 NEGATIVE CORRELATIONS")
message("============================================================")

neg_cors <- cor_long %>% 
  filter(correlation < 0) %>%
  arrange(correlation) %>%
  head(10)

for (i in seq_len(nrow(neg_cors))) {
  message(sprintf("  %s <-> %s: r = %.3f (n = %d)",
                  neg_cors$var1[i], neg_cors$var2[i], 
                  neg_cors$correlation[i], neg_cors$n_obs[i]))
}

message("\n============================================================")
message("DONE")
message("============================================================")
message("Outputs saved to: ", OUTPUT_DIR)
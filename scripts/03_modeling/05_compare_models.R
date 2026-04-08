# ============================================================
# 106_compare_ems_models.R
# 
# Compare EMS model variants:
#   A: Instantaneous (TS_Ha1 × NEON_SWC_shallow, 1h)
#   B: BGS drivers (TS_Ha2 × bvs_wtd_cm, 282h/132h)
# ============================================================

library(tidyverse)
library(lme4)

# ============================================================
# LOAD MODELS
# ============================================================

m_A <- readRDS("outputs/models/ems_instantaneous/m_final.rds")
m_B <- readRDS("outputs/models/ems_bgs_drivers/m_final.rds")

# ============================================================
# MODEL COMPARISON
# ============================================================

cat("══════════════════════════════════════════════════════════════\n")
cat("           EMS MODEL COMPARISON: A vs B\n")
cat("══════════════════════════════════════════════════════════════\n\n")

# Extract key metrics
get_model_stats <- function(m, name) {
  r2 <- var(predict(m, re.form = NA)) / var(m@frame[[1]])
  tibble(
    model = name,
    R2 = round(r2 * 100, 1),
    AIC = round(AIC(m), 1),
    BIC = round(BIC(m), 1),
    N = nrow(m@frame),
    n_trees = n_distinct(m@frame$Tree),
    residual_SD = round(sigma(m), 4)
  )
}

comparison <- bind_rows(
  get_model_stats(m_A, "A: Instantaneous (TS_Ha1 × SWC)"),
  get_model_stats(m_B, "B: BGS drivers (TS_Ha2 × WTD)")
)

print(comparison, width = 100)

cat("\n──────────────────────────────────────────────────────────────\n")
cat("Model formulas:\n")
cat("──────────────────────────────────────────────────────────────\n")
cat("A:", deparse(formula(m_A)), "\n")
cat("B:", deparse(formula(m_B)), "\n")

# ============================================================
# SPECIES SLOPES COMPARISON
# ============================================================

cat("\n══════════════════════════════════════════════════════════════\n")
cat("           SPECIES-SPECIFIC SLOPES\n")
cat("══════════════════════════════════════════════════════════════\n\n")

slopes_A <- read_csv("outputs/models/ems_instantaneous/species_slopes.csv",
                     show_col_types = FALSE) %>%
  mutate(model = "A")

slopes_B <- read_csv("outputs/models/ems_bgs_drivers/species_slopes.csv",
                     show_col_types = FALSE) %>%
  rename(wtd_slope = any_of(c("wtd_slope", "swc_slope"))) %>%
  mutate(model = "B")

# Standardize column names
if (!"wtd_slope" %in% names(slopes_A)) {
  slopes_A <- slopes_A %>% rename(wtd_slope = swc_slope)
}

slopes_combined <- bind_rows(slopes_A, slopes_B) %>%
  mutate(across(where(is.numeric), ~round(., 3)))

cat("Model A (Instantaneous):\n")
print(slopes_A %>% dplyr::select(-model))

cat("\nModel B (BGS drivers):\n")
print(slopes_B %>% dplyr::select(-model))

# ============================================================
# COEFFICIENT COMPARISON
# ============================================================

cat("\n══════════════════════════════════════════════════════════════\n")
cat("           FIXED EFFECTS SUMMARY\n")
cat("══════════════════════════════════════════════════════════════\n\n")

coef_A <- fixef(m_A)
coef_B <- fixef(m_B)

cat("Model A - significant effects (|t| > 2):\n")
summ_A <- summary(m_A)$coefficients
sig_A <- summ_A[abs(summ_A[, "t value"]) > 2, , drop = FALSE]
if (nrow(sig_A) > 0) {
  print(round(sig_A, 4))
} else {
  cat("  (none)\n")
}

cat("\nModel B - significant effects (|t| > 2):\n")
summ_B <- summary(m_B)$coefficients
sig_B <- summ_B[abs(summ_B[, "t value"]) > 2, , drop = FALSE]
if (nrow(sig_B) > 0) {
  print(round(sig_B, 4))
} else {
  cat("  (none)\n")
}

# ============================================================
# SUMMARY
# ============================================================

cat("\n══════════════════════════════════════════════════════════════\n")
cat("           SUMMARY\n")
cat("══════════════════════════════════════════════════════════════\n\n")

better <- if (comparison$R2[1] > comparison$R2[2]) "A" else "B"
r2_diff <- abs(comparison$R2[1] - comparison$R2[2])

cat(sprintf("Better fit (by R²): Model %s\n", better))
cat(sprintf("R² difference: %.1f percentage points\n", r2_diff))
cat(sprintf("AIC difference: %.1f (lower is better)\n", comparison$AIC[1] - comparison$AIC[2]))

cat("\n──────────────────────────────────────────────────────────────\n")
cat("Interpretation:\n")
cat("──────────────────────────────────────────────────────────────\n")
if (r2_diff < 2) {
  cat("→ Models explain similar variance - neither driver set is clearly better\n")
} else if (better == "A") {
  cat("→ Instantaneous conditions (TS_Ha1 × SWC) explain more variance\n")
  cat("→ Upland CH4 may respond to current rather than integrated conditions\n")
} else {
  cat("→ BGS drivers (TS_Ha2 × WTD) explain more variance\n
")
  cat("→ Same temporal integration as wetland applies to upland\n")
}

cat("\n══════════════════════════════════════════════════════════════\n")
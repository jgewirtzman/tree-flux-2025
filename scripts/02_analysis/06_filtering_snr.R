# ============================================================
# filtering_snr.R
#
# QC visualization: flux vs R², SNR flagging
# Shows CO2 and CH4 flux distributions with quality flags
#
# Outputs:
#   - 2x2 figure: flag-colored and density-colored scatter plots
#
# Requires: corrected flux file loaded as 'fluxes'
# ============================================================

library(dplyr)
library(ggplot2)
library(patchwork)
library(viridis)
library(ggpointdensity)

# ============================================================
# CONFIGURATION
# ============================================================

# QC thresholds
SNR_THR_CO2 <- 2
SNR_THR_CH4 <- 2
CO2_R2_THR  <- 0.8

# ============================================================
# LOAD DATA (if not already in environment)
# ============================================================

if (!exists("fluxes")) {
  message("Loading flux data...")
  fluxes <- read.csv("data/input/HF_2023-2025_tree_flux_corrected.csv")
}

# ============================================================
# PREP + UNIT CORRECTION
# ============================================================

# NOTE: The corrected flux file already has proper units.
#       We only need to correct pre-2025 CH4_SE (not flux)

df_qc <- fluxes %>%
  filter(
    !is.na(year),
    !is.na(CO2_flux_umolpm2ps), !is.na(CO2_r2), !is.na(CO2_SE),
    !is.na(CH4_flux_nmolpm2ps), !is.na(CH4_r2), !is.na(CH4_SE)
  ) %>%
  mutate(
    # Unit correction: pre-2025 CH4_SE only (flux already corrected in file)
    CH4_SE_corr = if_else(year < 2025, CH4_SE * 1000, CH4_SE),
    
    # SNR calculations
    CO2_snr = abs(CO2_flux_umolpm2ps) / CO2_SE,
    CH4_snr = abs(CH4_flux_nmolpm2ps) / CH4_SE_corr,
    
    # QC flags
    flag_CO2_snr = CO2_snr < SNR_THR_CO2,
    flag_CO2_r2  = CO2_r2 < CO2_R2_THR,
    flag_CO2_neg = CO2_flux_umolpm2ps < 0,
    flag_CH4_snr = CH4_snr < SNR_THR_CH4,
    
    # Y transforms for plotting
    CO2_flux_lin   = CO2_flux_umolpm2ps,
    CH4_flux_asinh = asinh(CH4_flux_nmolpm2ps)
  )

# ============================================================
# CLASSIFY INTO BUCKETS
# ============================================================

df_qc <- df_qc %>%
  mutate(
    nflags_co2 = (flag_CO2_snr + flag_CO2_r2 + flag_CO2_neg),
    nflags_ch4 = (flag_CH4_snr + flag_CO2_snr + flag_CO2_r2 + flag_CO2_neg),
    
    co2_class = case_when(
      nflags_co2 >= 2 ~ "Multiple flags",
      flag_CO2_snr    ~ "CO2 low SNR",
      flag_CO2_r2     ~ "CO2 low R2",
      flag_CO2_neg    ~ "CO2 < 0",
      TRUE            ~ "Neither"
    ),
    
    ch4_class = case_when(
      nflags_ch4 >= 2 ~ "Multiple flags",
      flag_CH4_snr    ~ "CH4 low SNR",
      flag_CO2_snr    ~ "CO2 low SNR",
      flag_CO2_r2     ~ "CO2 low R2",
      flag_CO2_neg    ~ "CO2 < 0",
      TRUE            ~ "Neither"
    ),
    
    co2_class = factor(co2_class, levels = c("Neither", "CO2 < 0", "CO2 low R2", 
                                             "CO2 low SNR", "Multiple flags")),
    ch4_class = factor(ch4_class, levels = c("Neither", "CO2 < 0", "CO2 low R2", 
                                             "CO2 low SNR", "CH4 low SNR", "Multiple flags"))
  )

# ============================================================
# COUNTS
# ============================================================

counts_co2 <- df_qc %>%
  summarise(
    n_total      = n(),
    n_co2_lowr2  = sum(flag_CO2_r2,  na.rm = TRUE),
    n_co2_lowsnr = sum(flag_CO2_snr, na.rm = TRUE),
    n_co2_neg    = sum(flag_CO2_neg, na.rm = TRUE),
    n_multiple   = sum(nflags_co2 >= 2, na.rm = TRUE),
    n_any        = sum(nflags_co2 >= 1, na.rm = TRUE)
  )

counts_ch4 <- df_qc %>%
  summarise(
    n_total      = n(),
    n_ch4_lowsnr = sum(flag_CH4_snr, na.rm = TRUE),
    n_co2_lowsnr = sum(flag_CO2_snr, na.rm = TRUE),
    n_co2_lowr2  = sum(flag_CO2_r2,  na.rm = TRUE),
    n_co2_neg    = sum(flag_CO2_neg, na.rm = TRUE),
    n_multiple   = sum(nflags_ch4 >= 2, na.rm = TRUE),
    n_any        = sum(nflags_ch4 >= 1, na.rm = TRUE)
  )

cat("=== CO2 QC Summary ===\n")
print(counts_co2)
cat("\n=== CH4 QC Summary ===\n")
print(counts_ch4)

# Annotation text
ann_co2 <- with(counts_co2, paste0(
  "Total: ", n_total,
  "\nCO2 low R² (<", CO2_R2_THR, "): ", n_co2_lowr2,
  "\nCO2 low SNR (<", SNR_THR_CO2, "): ", n_co2_lowsnr,
  "\nCO2 < 0: ", n_co2_neg,
  "\nMultiple: ", n_multiple,
  "\nAny flag: ", n_any
))

ann_ch4 <- with(counts_ch4, paste0(
  "Total: ", n_total,
  "\nCH4 low SNR (<", SNR_THR_CH4, "): ", n_ch4_lowsnr,
  "\nCO2 low SNR (<", SNR_THR_CO2, "): ", n_co2_lowsnr,
  "\nCO2 low R² (<", CO2_R2_THR, "): ", n_co2_lowr2,
  "\nCO2 < 0: ", n_co2_neg,
  "\nMultiple: ", n_multiple,
  "\nAny flag: ", n_any
))

# ============================================================
# COLOR PALETTES
# ============================================================

flag_cols_co2 <- c(
  "Neither"        = "black",
  "CO2 < 0"        = "#1B9E77",
  "CO2 low R2"     = "#2166AC",
  "CO2 low SNR"    = "#D55E00",
  "Multiple flags" = "#7A0177"
)

flag_cols_ch4 <- c(
  "Neither"        = "black",
  "CO2 < 0"        = "#1B9E77",
  "CO2 low R2"     = "#2166AC",
  "CO2 low SNR"    = "#D55E00",
  "CH4 low SNR"    = "#B2182B",
  "Multiple flags" = "#7A0177"
)

# Helper for density legend
dens_breaks <- function(lim) c(lim[1], (lim[1] + lim[2]) / 2, lim[2])

# ============================================================
# FLAG-COLORED PLOTS
# ============================================================

p_co2_flags <- ggplot(df_qc, aes(x = CO2_r2, y = CO2_flux_lin, color = co2_class)) +
  geom_point(alpha = 0.6, size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  coord_cartesian(xlim = c(-0.2, 1.02)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_manual(
    values = flag_cols_co2,
    limits = levels(df_qc$co2_class),
    breaks = levels(df_qc$co2_class),
    name = NULL,
    drop = FALSE
  ) +
  annotate("text", x = -0.18, y = Inf, label = ann_co2, hjust = 0, vjust = 1.1, size = 3.2) +
  labs(
    x = expression(CO[2]~R^2),
    y = expression(CO[2]~flux~(mu*mol~m^{-2}~s^{-1}))
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3), nrow = 2, byrow = TRUE))

p_ch4_flags <- ggplot(df_qc, aes(x = CH4_r2, y = CH4_flux_asinh, color = ch4_class)) +
  geom_point(alpha = 0.6, size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  coord_cartesian(xlim = c(-0.2, 1.02)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_manual(
    values = flag_cols_ch4,
    limits = levels(df_qc$ch4_class),
    breaks = levels(df_qc$ch4_class),
    name = NULL,
    drop = FALSE
  ) +
  annotate("text", x = -0.18, y = Inf, label = ann_ch4, hjust = 0, vjust = 1.1, size = 3.2) +
  labs(
    x = expression(CH[4]~R^2),
    y = expression(asinh(CH[4]~flux~(nmol~m^{-2}~s^{-1})))
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3), nrow = 2, byrow = TRUE))

# ============================================================
# DENSITY-COLORED PLOTS
# ============================================================

p_co2_dens <- ggplot(df_qc, aes(x = CO2_r2, y = CO2_flux_lin)) +
  ggpointdensity::geom_pointdensity(adjust = 1, size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  coord_cartesian(xlim = c(-0.2, 1.02)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_viridis(option = "viridis", name = "Local\ndensity", breaks = dens_breaks) +
  labs(
    x = expression(CO[2]~R^2),
    y = expression(CO[2]~flux~(mu*mol~m^{-2}~s^{-1}))
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

p_ch4_dens <- ggplot(df_qc, aes(x = CH4_r2, y = CH4_flux_asinh)) +
  ggpointdensity::geom_pointdensity(adjust = 1, size = 1.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  coord_cartesian(xlim = c(-0.2, 1.02)) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_viridis(option = "viridis", name = "Local\ndensity", breaks = dens_breaks) +
  labs(
    x = expression(CH[4]~R^2),
    y = expression(asinh(CH[4]~flux~(nmol~m^{-2}~s^{-1})))
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

# ============================================================
# COMBINED 2x2 FIGURE
# ============================================================

fig_qc <- (p_co2_flags | p_ch4_flags) /
  (p_co2_dens | p_ch4_dens)

print(fig_qc)

# Save
ggsave("outputs/figures/qc_flux_vs_r2.png", fig_qc, width = 12, height = 10, dpi = 300)
ggsave("outputs/figures/qc_flux_vs_r2.pdf", fig_qc, width = 12, height = 10)

message("Saved: figures/qc_flux_vs_r2.png/pdf")
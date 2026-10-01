# ============================================================
# 12_driver_timeseries_si.R
# SI figure: daily environmental drivers over the study period (June 2023 - October 2025),
# with flux sampling dates marked. Replaces a static image in the earlier draft.
# Input: data/processed/aligned_hourly_dataset.csv (EST clock), data/processed/flux_with_quality_flags.csv
# Output: outputs/figures/si_driver_timeseries.png / .pdf
# ============================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(lubridate) })
a <- read.csv("data/processed/aligned_hourly_dataset.csv") %>% mutate(date = as.Date(substr(datetime, 1, 10)))
fl <- read.csv("data/processed/flux_with_quality_flags.csv") %>% transmute(date = as.Date(date), site = location) %>% distinct()
START <- as.Date("2023-06-01"); END <- as.Date("2025-10-31")
vars <- c(tair_C = "Air temperature (°C)", TS_Ha2 = "Soil temperature, Ha2 (°C)",
          bvs_wtd_cm = "Water table, BVS well (cm)", bgs_wtd_cm = "Water table, BGS well (cm)",
          NEON_SWC_shallow = "Soil water content, NEON shallow (m3 m-3)", P_mm = "Precipitation (mm per day)",
          VPD_kPa = "Vapour pressure deficit (kPa)", LE_Ha1 = "Latent heat flux, Ha1 (W m-2)")
d <- a %>% filter(date >= START, date <= END) %>% select(date, all_of(names(vars))) %>%
  pivot_longer(-date, names_to = "var", values_to = "v") %>% group_by(var, date) %>%
  summarise(v = if (first(var) == "P_mm") sum(v, na.rm = TRUE) * (sum(!is.na(v)) >= 20) else mean(v, na.rm = TRUE),
            .groups = "drop") %>% mutate(v = ifelse(is.nan(v), NA, v), var = factor(vars[var], vars))
samp <- fl %>% filter(date >= START, date <= END)
g <- ggplot(d, aes(date, v)) +
  geom_vline(data = samp %>% filter(site == "Wetland"), aes(xintercept = date), colour = "#2A7F7A", alpha = 0.25, linewidth = 0.3) +
  geom_vline(data = samp %>% filter(site == "Upland"), aes(xintercept = date), colour = "#6E8B3D", alpha = 0.25, linewidth = 0.3, linetype = "22") +
  geom_line(linewidth = 0.3, colour = "grey15") +
  facet_wrap(~ var, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_x_date(date_breaks = "3 months", date_labels = "%b\n%Y", expand = c(0.005, 0)) +
  labs(x = NULL, y = NULL) + theme_bw(base_size = 9) +
  theme(strip.placement = "outside", strip.background = element_blank(), strip.text.y.left = element_text(angle = 0, hjust = 1, size = 7.5),
        panel.grid.minor = element_blank())
ggsave("outputs/figures/si_driver_timeseries.png", g, width = 8, height = 10, dpi = 300, bg = "white")
ggsave("outputs/figures/si_driver_timeseries.pdf", g, width = 8, height = 10, bg = "white")
message("Saved si_driver_timeseries.png/pdf")

library(tidyverse)

files <- tribble(
  ~camera,                 ~path,                                                                                         ~type,
  "harvardems2_DB",         "data/raw/phenocam/harvardems2_DB_1000_ndvi_1day.txt",       "txt",
)

read_camera <- function(camera, path, type) {
  
  df <- if (type == "txt") {
    readr::read_csv(path, comment = "#", show_col_types = FALSE)
  } else {
    readr::read_csv(path, show_col_types = FALSE)
  }
  
  # Harmonize GCC safely
  if ("gcc_mean" %in% names(df)) {
    df$gcc <- df$gcc_mean
  } else if ("gcc_90" %in% names(df)) {
    df$gcc <- df$gcc_90
  } else {
    df$gcc <- NA_real_
  }
  
  df %>%
    mutate(
      camera = camera,
      date = as.Date(date)
    ) %>%
    select(any_of(c("camera", "date", "gcc", "ndvi_90")))
}

dat <- purrr::pmap_dfr(files, read_camera)

dat_long <- dat %>%
  pivot_longer(cols = c(gcc, ndvi_90), names_to = "metric", values_to = "value") %>%
  filter(!is.na(value))

ggplot(dat_long, aes(x = date, y = value)) +
  geom_line(na.rm = TRUE) +
  facet_grid(metric ~ camera, scales = "free_y") +
  labs(x = "Date", y = NULL,
       title = "GCC (smooth_gcc_mean / gcc_90) and NDVI (ndvi_90) by camera") +
  theme_minimal()


# ---- Restrict to 2022–2026 ----
dat_2022_2026 <- dat %>%
  filter(date >= as.Date("2022-01-01"),
         date <= as.Date("2026-12-31"))

# ---- Long format ----
dat_long <- dat_2022_2026 %>%
  pivot_longer(
    cols = c(gcc, ndvi_90),
    names_to = "metric",
    values_to = "value"
  ) %>%
  filter(!is.na(value))

# ---- Plot ----
ggplot(dat_long, aes(x = date, y = value)) +
  geom_line(na.rm = TRUE) +
  facet_grid(metric ~ camera, scales = "free_y") +
  labs(
    x = "Date",
    y = NULL,
    title = "GCC and NDVI by camera (2022–2026)"
  ) +
  theme_minimal()

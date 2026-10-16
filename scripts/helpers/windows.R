# ============================================================
# windows.R
# Helper: right-aligned rolling windows and anomaly filtering of hourly drivers,
# shared by 4_analysis/05_window_screening.R and the 5_models/ scripts.
#
# Usage:
#   source("scripts/helpers/windows.R")
#   x_66h <- roll_mean(met$TS_Ha2, 66)
#   met   <- make_anomaly(met, "LE_Ha2")   # adds LE_Ha2_anom
#
# Requires: RcppRoll, dplyr, lubridate (yday, hour)
# ============================================================

# Rolling windows need at least half of their hours observed; otherwise NA
# (a mean over a mostly empty window is not the window it claims to be)
MIN_WINDOW_COMPLETE <- 0.5
.win_ok <- function(x, n) RcppRoll::roll_sum(as.numeric(!is.na(x)), n = n, align = "right", fill = NA) >= MIN_WINDOW_COMPLETE * n
roll_mean <- function(x, n) {
  m <- RcppRoll::roll_mean(x, n = n, align = "right", fill = NA, na.rm = TRUE); m[!.win_ok(x, n) %in% TRUE] <- NA; m
}

roll_sum <- function(x, n) {
  m <- RcppRoll::roll_sum(x, n = n, align = "right", fill = NA, na.rm = TRUE); m[!.win_ok(x, n) %in% TRUE] <- NA; m
}

# Anomaly filtering: remove day-of-year and hour-of-day main effects
# (adds column <var>_anom = x - (mean by DOY + mean by hour - grand mean))
make_anomaly <- function(df, var, time_col = "datetime") {
  if (!var %in% names(df)) return(df)
  x <- df[[var]]
  if (!is.numeric(x)) return(df)

  df <- df %>%
    mutate(
      .doy  = yday(.data[[time_col]]),
      .hour = hour(.data[[time_col]])
    )

  grand <- mean(x, na.rm = TRUE)

  doy_means <- df %>%
    group_by(.doy) %>%
    summarize(.m_doy = mean(.data[[var]], na.rm = TRUE), .groups = "drop")

  hr_means <- df %>%
    group_by(.hour) %>%
    summarize(.m_hr = mean(.data[[var]], na.rm = TRUE), .groups = "drop")

  df %>%
    left_join(doy_means, by = ".doy") %>%
    left_join(hr_means, by = ".hour") %>%
    mutate(
      .expected = .m_doy + .m_hr - grand,
      "{var}_anom" := .data[[var]] - .expected
    ) %>%
    dplyr::select(-.doy, -.hour, -.m_doy, -.m_hr, -.expected)
}

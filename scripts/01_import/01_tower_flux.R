Ha1<-read.csv(
  'data/raw/ameriflux/AMF_US-Ha1_BASE-BADM_26-5/AMF_US-Ha1_BASE_HR_26-5.csv',
  header=T, skip=2)
head(Ha1)

Ha2<-read.csv(
  'data/raw/ameriflux/AMF_US-Ha2_BASE-BADM_15-5/AMF_US-Ha2_BASE_HH_15-5.csv',
  header=T, skip=2)
head(Ha2)

xHA<-read.csv(
  'data/raw/ameriflux/AMF_US-xHA_BASE-BADM_11-5/AMF_US-xHA_BASE_HH_11-5.csv',
  header=T, skip=2)
head(xHA)

library(dplyr)
library(lubridate)
library(stringr)
library(tidyr)

to_posix_af <- function(x) {
  # AmeriFlux half-hourly/hourly TIMESTAMP_START is typically YYYYMMDDHHMM (12 digits)
  x_chr <- format(as.numeric(x), scientific = FALSE, trim = TRUE)
  x_chr <- str_pad(x_chr, width = 12, side = "left", pad = "0")
  ymd_hm(x_chr, tz = "UTC")
}

summarize_fc_coverage <- function(df, tower_name, start_year = 2022, end_year = 2026) {
  # Identify all FC columns (e.g., FC_1_1_1, FC_2_1_1, FC_1_1_2)
  fc_cols <- names(df)[grepl("^FC(_|$)", names(df))]
  if (length(fc_cols) == 0) {
    return(tibble(
      tower = tower_name,
      fc_col = character(),
      position = character(),
      n_total = integer(),
      n_nonmissing = integer(),
      pct_nonmissing = double(),
      has_any_data = logical()
    ))
  }
  
  df2 <- df %>%
    mutate(datetime = to_posix_af(TIMESTAMP_START)) %>%
    filter(year(datetime) >= start_year, year(datetime) <= end_year)
  
  n_total <- nrow(df2)
  
  # Summarize each FC column
  out <- lapply(fc_cols, function(col) {
    x <- df2[[col]]
    x <- ifelse(x == -9999, NA_real_, as.numeric(x))
    n_nonmissing <- sum(!is.na(x))
    tibble(
      tower = tower_name,
      fc_col = col,
      position = sub("^FC_", "", col), # "1_1_1", "2_1_1", etc. (for plain "FC" this will stay "FC")
      n_total = n_total,
      n_nonmissing = n_nonmissing,
      pct_nonmissing = if (n_total > 0) 100 * n_nonmissing / n_total else NA_real_,
      has_any_data = n_nonmissing > 0
    )
  }) %>% bind_rows()
  
  out %>% arrange(desc(has_any_data), desc(pct_nonmissing), fc_col)
}

# Run for your towers
ha1_cov <- summarize_fc_coverage(Ha1, "Ha1")
ha2_cov <- summarize_fc_coverage(Ha2, "Ha2")
xha_cov <- summarize_fc_coverage(xHA, "xHA")

coverage_all <- bind_rows(ha1_cov, ha2_cov, xha_cov)

coverage_all






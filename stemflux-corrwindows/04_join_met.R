ems_met<-read.csv("data/raw/wtd_ems-met.csv")
str(ems_met)

ggplot(ems_met, aes(x=datetime))+
  geom_point(aes(y=vwc_60cm, col="red"))+
  geom_point(aes(y=vwc_top, col="blue"))

ggplot(ems_met, aes(y=vwc_60cm, x=vwc_top))+
  geom_point()
  

# Convert both datetime and date in ems_met to the correct data types
ems_met <- ems_met %>%
  mutate(
    datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
    date = as.Date(Date)  # Note: using 'Date' column to create 'date'
  ) %>%
  select(-Date)  # Remove the original 'Date' column to avoid confusion

ems_met <- ems_met %>%
  select(datetime, vwc_top, vwc_60cm, LE, tair_C) %>%
  rename(t_airC_ems = tair_C)

# Now join the dataframes
wtd_met <- left_join(wtd_met, ems_met, by = "datetime")

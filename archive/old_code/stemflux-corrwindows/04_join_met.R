dt10 <- read.csv("https://harvardforest.fas.harvard.edu/data/p00/hf001/hf001-10-15min-m.csv")
dt4 <- read.csv("https://harvardforest.fas.harvard.edu/data/p07/hf070/hf070-04-15min.csv")


# Aggregate HF stream temp & discharge
# first run download-.R files

library(tidyverse)

## 15-min hydro data
hydro15 = dt4 %>%
  mutate(datetime = lubridate::ymd_hm(datetime),
         hour = lubridate::hour(datetime),
         date = lubridate::date(datetime),
         year = lubridate::year(datetime)) %>%
  #filter(hour>=8 & hour<=16, year>2022) %>% # stream filter
  filter(year>2022) 

# clean up
#rm(dt4)

# after download-15min-FisherMet.R
fishermet = dt10 %>%
  mutate(datetime = lubridate::ymd_hm(datetime),
         hour = lubridate::hour(datetime),
         date = lubridate::date(datetime),
         year = lubridate::year(datetime),
         month = lubridate::month(datetime))

hydromet = left_join(hydro15,fishermet)

# clean up
#rm(dt10)

# get BVS & BGS wtd at hourly timescale
wtd_ems = hydromet %>%
  select(datetime,jd,bgs.stg,bvs.stg,bgs.wt,bvs.wt) %>%
  rename(bgs_wtd_cm = bgs.stg, bvs_wtd_cm = bvs.stg,
         bgs_wtemp_C = bgs.wt, bvs_wtemp_C = bvs.wt) %>%
  mutate(hour = str_pad(lubridate::hour(datetime),pad=0,
                        width=2,side="left"),
         date = lubridate::date(datetime),
         datetime_old = datetime,
         datetime = lubridate::ymd_hms(paste0(date," ",hour,":00:00",tz="EST"))) %>%
  group_by(datetime) %>%
  summarize(bgs_wtd_cm = mean(bgs_wtd_cm,na.rm=T),
            bvs_wtd_cm = mean(bvs_wtd_cm,na.rm=T)) #%>%
#left_join(ems_vwc)

# Met from 15min to hourly
met = fishermet %>%
  group_by(date,hour) %>%
  summarize(tair_C = mean(airt,na.rm=T),
            p_kPa = mean(bar/10,na.rm=T),
            P_mm = sum(prec,na.rm=T),
            RH = mean(rh,na.rm=T),
            PAR = mean(parr,na.rm=T),
            rnet = mean(netr,na.rm=T),
            slrr = mean(slrr,na.rm=T),   # Added solar radiation
            s10t = mean(s10t,na.rm=T)) %>%   # Added soil temperature at 10cm 
  mutate(datetime = lubridate::ymd_hms(paste0(date," ",
                                              stringr::str_pad(hour,pad="0",side="left",width=2),
                                              ":00:00")),
         VPD_kPa = plantecophys::RHtoVPD(RH=RH,TdegC=tair_C,Pa = p_kPa*1000))

wtd_met = wtd_ems %>%
  left_join(met,by=c("datetime"))




#ems_met<-read.csv("data/raw/wtd_ems-met.csv")
#str(ems_met)

#ggplot(ems_met, aes(x=datetime))+
#  geom_point(aes(y=vwc_60cm, col="red"))+
#  geom_point(aes(y=vwc_top, col="blue"))

#ggplot(ems_met, aes(y=vwc_60cm, x=vwc_top))+
#  geom_point()
#   
# 
# # Convert both datetime and date in ems_met to the correct data types
# ems_met <- ems_met %>%
#   mutate(
#     datetime = as.POSIXct(datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"),
#     date = as.Date(Date)  # Note: using 'Date' column to create 'date'
#   ) %>%
#   select(-Date)  # Remove the original 'Date' column to avoid confusion
# 
# ems_met <- ems_met %>%
#   select(datetime, vwc_top, vwc_60cm, LE, tair_C) %>%
#   rename(t_airC_ems = tair_C)
# 
# # Now join the dataframes
# wtd_met <- left_join(wtd_met, ems_met, by = "datetime")

write_csv(wtd_met, "data/processed/wtd_met.csv")


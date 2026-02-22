# Harmonize soil moisture datasets to characterize 
# associations between soil moisture and eddy flux

library(tidyverse)
library(lubridate)

# Harmonize soil moisture data 1993-early 2000s from point measurements

# Warming experiment control plots, 1993-2002
inUrl1  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-hfr/5/37/f75235fbbc3fc2396171f4b1fa0b18c2" 
infile1 <- tempfile()
try(download.file(inUrl1,infile1,method="curl"))
if (is.na(file.size(infile1))) download.file(inUrl1,infile1,method="auto")
dt1 <-read.csv(infile1,header=F,skip=1,sep=",", 
               col.names=c("year","rep","date","treatment","block",     
                           "co2_flux","ch4_flux","n2o_flux","temp_2cm",     
                           "temp_4cm","grav_h2o_org","grav_h2o_min","vol_h2o"), 
               check.names=TRUE)

unlink(infile1)

warm_sm = dt1 %>%
  filter(treatment == 1) %>%
  filter(lubridate::month(date)>=5 & lubridate::month(date)<=10)%>%
  group_by(date) %>%
  summarize(#grav_h2o_mn = mean(grav_h2o_org),
    vwc = mean(vol_h2o, na.rm=T)) %>%
  filter(!is.na(vwc))

inUrl7  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-hfr/5/37/db1539e580ff99659be277667d9e34c1" 
infile7 <- tempfile()
try(download.file(inUrl7,infile7,method="curl"))
if (is.na(file.size(infile7))) download.file(inUrl7,infile7,method="auto")

dt7 <-read.csv(infile7,header=F,skip=1,sep=",", 
               col.names=c("date","p1_sm","p2_sm","p3_sm","p4_sm",     
                           "p5_sm","p6_sm","p7_sm","p8_sm","p9_sm","p10_sm",     
                           "p11_sm","p12_sm","p13_sm","p14_sm","p15_sm",     
                           "p16_sm","p17_sm", "p18_sm"), check.names=TRUE)

unlink(infile7)

# select only control plots
warm_ctlplots = c("date","p2_sm","p4_sm","p7_sm","p11_sm","p14_sm","p18_sm")
warm_sm2 = dt7[,warm_ctlplots] 
warm_sm2$vwc = apply(warm_sm2[,2:7],1,FUN = "mean",
                     na.action="na.rm")
warm_sm2 = select(warm_sm2, date, vwc)

# combine all warming soil moisture
warm_vwc = bind_rows(warm_sm, warm_sm2) %>%
  mutate(year = lubridate::year(date)) 
warm_vwc$vwc_anom = (warm_vwc$vwc-mean(warm_vwc$vwc,na.rm=T))/sd(warm_vwc$vwc,na.rm=T)
warm_vwc = warm_vwc %>% group_by(year) %>%
  summarize(vwc_anom = mean(vwc_anom, na.rm=T))

# EMS continuous soil moisture sensors 
# in biometry plots: 2010-present
inUrl5  <- "https://pasta.lternet.edu/package/data/eml/knb-lter-hfr/206/28/ff8708b51173fc2ba08abfef514e39c2" 
infile5 <- tempfile()
try(download.file(inUrl5,infile5,method="curl"))
if (is.na(file.size(infile5))) download.file(inUrl5,infile5,method="auto")

dt5 <-read.csv(infile5,header=F,skip=1,sep=",", col.names=c("year","datetime",     
                                                            "doy","dec_hour","dec_date","par_c3_ave","par_b2_ave","par_c2_ave",     
                                                            "par_e3_ave","par_f3_ave","par_f2_ave","tsoil1_n_ave","tsoil1_s_ave",     
                                                            "tsoil1_e_ave","tsoil1_w_ave","tsoil2_n_ave","tsoil2_s_ave",     
                                                            "tsoil2_e_ave","tsoil2_w_ave","tair_us_e_1m_ave",     
                                                            "tair_us_e_20cm_ave","tair_us_w_1m_ave","tair_us_w_20cm_ave",     
                                                            "tair_us_n_1m_ave","tair_us_n_20cm_ave","tair_us_s_1m_ave",     
                                                            "tair_us_s_20cm_ave","swc1_ave","swc2_ave","swc3_ave",     
                                                            "swc4_ave","par_c3_max","par_b2_max","par_c2_max",     
                                                            "par_e3_max","par_f3_max","par_f2_max","tsoil1_n_max",     
                                                            "tsoil1_s_max","tsoil1_e_max","tsoil1_w_max","tsoil2_n_max",     
                                                            "tsoil2_s_max","tsoil2_e_max","tsoil2_w_max","tair_us_e_1m_max",     
                                                            "tair_us_e_20cm_max","tair_us_w_1m_max","tair_us_w_20cm_max",     
                                                            "tair_us_n_1m_max","tair_us_n_20cm_max","tair_us_s_1m_max",     
                                                            "tair_us_s_20cm_max","swc1_max","swc2_max","swc3_max",     
                                                            "swc4_max","par_c3_min","par_b2_min","par_c2_min",     
                                                            "par_e3_min","par_f3_min","par_f2_min","tsoil1_n_min",     
                                                            "tsoil1_s_min","tsoil1_e_min","tsoil1_w_min","tsoil2_n_min",     
                                                            "tsoil2_s_min","tsoil2_e_min","tsoil2_w_min","tair_us_e_1m_min",     
                                                            "tair_us_e_20cm_min","tair_us_w_1m_min","tair_us_w_20cm_min",     
                                                            "tair_us_n_1m_min","tair_us_n_20cm_min","tair_us_s_1m_min",     
                                                            "tair_us_s_20cm_min","swc1_min","swc2_min","swc3_min",     
                                                            "swc4_min"), check.names=TRUE)

unlink(infile5)

# EMS continuous soil moisture sensors 
# in biometry plots: 2010-present
# EMS1: C3, C2, B2; EMS2: F3, F2, E3
ems_biom = dt5 %>%
  mutate(date = lubridate::date(datetime),
         month = lubridate::month(datetime),
         year = lubridate::year(datetime)) %>%
  filter(month>=5 & month<=10)%>%
  select(year, swc1_ave:swc2_ave) 

ems_biom$vwc = rowMeans(cbind(ems_biom$swc1_ave, ems_biom$swc2_ave),na.rm=T)

emsbiom_vwc = ems_biom %>%
  mutate(vwc_anom = (vwc-mean(ems_biom$vwc,na.rm=T))/sd(ems_biom$vwc,na.rm=T)) %>%
  filter(!is.na(vwc)) %>%
  group_by(year) %>%
  summarize(vwc_anom = mean(vwc_anom, na.rm=T))


### EMS TOWER MOISTURE
ems_swc_0717 = read_csv("~/Documents/field-data/hf069-18-swc-2007-2017.csv") %>%
  mutate(date = as.Date(doy, origin=paste0(year-1,"-12-31")),
         month = lubridate::month(date)) %>%
  filter(month >= 5 & month <= 9) 

ems_swc_0717$vwc = rowMeans(cbind(ems_swc_0717$vwc_15_a,
                                  ems_swc_0717$vwc_15_b))

ems_swc_1722 = read_csv("~/Documents/field-data/hf069-19-swc-since-2017.csv") %>%
  mutate(date = as.Date(doy, origin=paste0(year-1,"-12-31")),
         month = lubridate::month(date)) %>%
  filter(month >= 5 & month <= 9) 

ems_swc_1722$vwc = rowMeans(cbind(ems_swc_1722$vwc_15cm_avg_m_3_m_3_a,
                                  ems_swc_1722$vwc_15cm_avg_m_3_m_3_b))

ems_swc = bind_rows(ems_swc_0717, ems_swc_1722) %>%
  select(date, year, vwc)

ems_tower_anom = ems_swc %>% 
  mutate(vwc_anom = (vwc-mean(ems_swc$vwc, na.rm=T))/sd(ems_swc$vwc, na.rm=T)) %>%
  group_by(year) %>%
  summarize(vwc_anom = mean(vwc_anom, na.rm=T))

# HEM tower understory
hem_sm = read_csv("~/Documents/flux_data/met/hf206-03-HK-understory.csv") %>%
  mutate(date = lubridate::date(datetime)) %>%
  filter(lubridate::month(datetime)>=5 & lubridate::month(datetime)<=10)%>%
  select(date, swc1_ave, swc2_ave, swc3_ave, swc4_ave) 

hem_sm$vwc = rowMeans(cbind(hem_sm$swc1_ave, hem_sm$swc2_ave, hem_sm$swc3_ave, 
                            hem_sm$swc4_ave), na.rm=T)
hem_sm$vwc_anom = (hem_sm$vwc - mean(hem_sm$vwc,na.rm=T))/sd(hem_sm$vwc,na.rm=T)
hem_sm$year = lubridate::year(hem_sm$date)  

hem_sm = hem_sm %>% group_by(year) %>%
  summarize(vwc_anom = mean(vwc_anom, na.rm=T))

### SOIL MOISTURE ANOMALY BY SITE
colors = c("HEM Tower" = "darkgreen", 
           "EMS Bio Plots" = "black",
           "Soil Warm. Control" = "tomato4",
           "EMS Tower" = "skyblue3")
shapes = c("HEM Tower" = 15, 
           "EMS Bio Plots" = 16,
           "Soil Warm. Control" = 17,
           "EMS Tower" = 18)

ggplot(hem_sm) +
  geom_point(aes(year, vwc_anom, color = "HEM Tower",shape="HEM Tower")) +
  geom_smooth(aes(year, vwc_anom, color = "HEM Tower"),method="lm",se=F,lwd=0.75) +
  geom_point(data = emsbiom_vwc, aes(year, vwc_anom, color="EMS Bio Plots",shape="EMS Bio Plots"))+
  geom_smooth(data = emsbiom_vwc,aes(year, vwc_anom, color = "EMS Bio Plots"),method="lm",se=F,lwd=0.75) +
  geom_point(data = warm_vwc, aes(year, vwc_anom, color="Soil Warm. Control",shape="Soil Warm. Control")) +
  geom_smooth(data = warm_vwc,aes(year, vwc_anom, color = "Soil Warm. Control"),method="lm",se=F,lwd=0.75) +
  geom_point(data = ems_tower_anom, aes(year, vwc_anom,color="EMS Tower",shape="EMS Tower")) +
  geom_smooth(data = ems_tower_anom,aes(year, vwc_anom, color = "EMS Tower"),method="lm",se=F,lwd=0.75) +
  scale_color_manual(name = "Site",
                     values = colors) +
  scale_shape_manual(name = "Site", values = shapes) +
  theme_linedraw() +
  labs(y="Soil moisture anomaly (z-score)", x="Year",
       title="Long-term soil moisture") 

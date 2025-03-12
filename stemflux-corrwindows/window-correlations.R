# Find correlations between tree flux data and met data over windows
library(tidyverse)

# Load flux data 
stem_flux = read_csv("~/Documents/SPLATS-trees/stemflux-2324/HFtree-stemflux-2023-24.csv")

# Format columns
stem_flux$datetime = lubridate::round_date(lubridate::force_tz(as.POSIXct(stem_flux$real_start-2190),
                                                               tzone="EST"),"hour")
stem_flux$date = lubridate::date(stem_flux$datetime)
stem_flux$ID = stem_flux$Tree
stem_flux$site = ifelse(stem_flux$Plot=="BGS","BGS","EMS")

# Fix missing species
stem_flux$species[stem_flux$ID==288] = "hem"
stem_flux$species[stem_flux$ID==153] = "rm"
stem_flux$species[stem_flux$ID==414] = "bg"
stem_flux$species[stem_flux$ID==452] = "bg"

# Table with number of measurements x site x species
n_bySiteSpecies = stem_flux %>%
  group_by(site,species) %>%
  summarize(CH4_tree = median(CH4_flux,na.rm=T),
            count = n())

# Options to loop through for days time period
days = c(2/24, 6/24, 1:60)

# find running stats for windows
for(d in 1:length(days)){
  
  # hours to calculate rolling sum
  window = 24*days[d] 
  
  # Calculate rolling window from daily precip & temp
  met_cul = wtd_met %>%
    mutate(date = lubridate::date(datetime),
           month = lubridate::month(date),
           jday = lubridate::yday(date),
           jday_d = ifelse(month>=11, -1*(365-jday), jday),
           year = lubridate::year(date), wyear = ifelse(month>11,year+1,year)) %>%
    arrange(wyear,jday_d) %>%
    mutate(Ta_mn = RcppRoll::roll_mean(tair_C, window, 
                                       align = "right",na.rm=T, fill = NA),
           P_mm = RcppRoll::roll_sum(P_mm, window, 
                                     align = "right",na.rm=T, fill = NA),
           RH_mn = RcppRoll::roll_mean(RH, window, 
                                       align = "right",na.rm=T, fill = NA),
           VPD_mn = RcppRoll::roll_mean(VPD_kPa, window, 
                                        align = "right",na.rm=T, fill = NA),
           # LE_mn = RcppRoll::roll_mean(LE, window, 
           #                             align = "right",na.rm=T, fill = NA),
           wtd_mn = RcppRoll::roll_mean(bgs_wtd_cm, window, 
                                        align = "right",na.rm=T, fill = NA),
           wtd_sd = RcppRoll::roll_sd(bgs_wtd_cm, window, 
                                      align = "right",na.rm=T, fill = NA),
           PAR_mn = RcppRoll::roll_mean(PAR, window, 
                                        align = "right",na.rm=T, fill = NA))
  
  stem_met = left_join(stem_flux,met_cul)
  
  # do correlation tests for met windows
  cor_test = stem_met %>%
    group_by(site) %>%
    summarize(cor_ta = cor.test(CH4_flux,Ta_mn)$estimate,
              p_ta = cor.test(CH4_flux,Ta_mn)$p.value,
              cor_vpd = cor.test(CH4_flux,VPD_mn)$estimate,
              p_vpd = cor.test(CH4_flux,VPD_mn)$p.value,
              cor_wtd = cor.test(CH4_flux,wtd_mn)$estimate,
              p_wtd = cor.test(CH4_flux,wtd_mn)$p.value,
              cor_wtdvar = cor.test(CH4_flux,wtd_sd)$estimate,
              p_wtdvar = cor.test(CH4_flux,wtd_sd)$p.value,
              cor_Pmm = cor.test(CH4_flux,P_mm)$estimate,
              p_Pmm = cor.test(CH4_flux,P_mm)$p.value,
              cor_PAR = cor.test(CH4_flux,PAR_mn)$estimate,
              p_PAR = cor.test(CH4_flux,PAR_mn)$p.value,
              cor_RH = cor.test(CH4_flux,RH)$estimate,
              p_RH = cor.test(CH4_flux,RH)$p.value) %>%
    mutate(interval = days[d])
  if(d == 1){
    cor_all = cor_test
  } else {
    cor_all = bind_rows(cor_all,cor_test)
  }
}
# Plot interval x cor variable
ggplot(cor_all) + 
  geom_point(aes(interval,cor_wtd,color=site))


# corflux_bgs = filter(stem_met,site=="BGS") %>% 
#   select(CO2_flux,CH4_flux,Ta_mn,VPD_mn,LE_mn,
#          wtd_mn,wtd_sd,P_mm,PAR_mn,RH)
# corflux_bgs = corflux_bgs[complete.cases(corflux_bgs),]
# bgs_cor = cor(as.matrix(corflux_bgs))
# testRes = corrplot::cor.mtest(corflux_bgs, conf.level = 0.95)
# corrplot::corrplot(bgs_cor,method = 'number',
#                    tl.col = 'black',type = 'lower',
#                    p.mat = testRes$p, sig.level = 0.05,
#                    diag=F,pinsig = "n")
# 
# 
# corflux_bg = filter(stem_met,site=="BGS",species=="rm") %>% 
#   select(CO2_flux,CH4_flux,Ta_mn,VPD_mn,LE_mn,
#          wtd_mn,wtd_sd,P_mm,PAR_mn,RH)
# corflux_bgs = corflux_bg[complete.cases(corflux_bg),]
# bgs_cor = cor(as.matrix(corflux_bgs))
# testRes = corrplot::cor.mtest(corflux_bgs, conf.level = 0.95)
# corrplot::corrplot(bgs_cor,method = 'number',
#                    tl.col = 'black',type = 'lower',
#                    diag=F,p.mat = testRes$p, sig.level = 0.05)
# 
# corflux_ems = filter(stem_met,site=="EMS",species=="hem") %>% 
#   select(CO2_flux,CH4_flux,Ta_mn,VPD_mn,LE_mn,
#          wtd_mn,wtd_sd,P_mm,PAR_mn,RH)
# corflux_ems = corflux_ems[complete.cases(corflux_ems),]
# ems_cor = cor(as.matrix(corflux_ems))
# corrplot::corrplot(ems_cor,method = 'number',
#                    tl.col = 'black',type = 'lower',
#                    diag=F)
# 


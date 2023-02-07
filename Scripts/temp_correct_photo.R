
#Working on temperature correction for Photosynthesis package

library(tidyverse)
library(rpmodel)

pars_photo_dat <- read_csv("~/Documents/GitHub/DAT_Proj/Results/dat_fits_photo_pars_filt.csv")
pars_photo_trad <- read_csv("~/Documents/GitHub/DAT_Proj/Results/trad_fits_photo_pars.csv")

cmplt.rm_out <- read_csv("~/Documents/GitHub/DAT_Proj/Inputs/Aci_no_out.csv")

grp_curv <- cmplt.rm_out %>% 
  group_by(Data_point, unique) %>% 
  summarize(meanTleaf = mean(Tleaf))
grp_curv2 <- grp_curv %>% 
  rename(ID = unique)
grp_dat <- filter(grp_curv2, Data_point == "Before_DAT")
grp_trad <- filter(grp_curv2, Data_point == "Traditional")

pars_photo_dat$method <- "DAT"
pars_photo_trad$method <- "Traditional"

curv_dat_temp <- left_join(by = "ID", pars_photo_dat, grp_dat)
curv_trad_temp <- left_join(by = "ID", pars_photo_dat, grp_dat)

for (i in 1:length(curv_dat_temp)){
  tcleaf[i] <- curv_dat_temp$meanTleaf[i]
  tcref <- 25
  f[i] <- ftemp_inst_vcmax(tcleaf[i],tcref) #this function is in the rpmodel package
}
curv_dat_temp$fref <- f
curv_dat_temp_adj <- curv_dat_temp %>% 
  mutate(vcmax_25 = V_cmax * fref)

for (i in 1:length(curv_trad_temp)){
  tcleaf[i] <- curv_trad_temp$meanTleaf[i]
  tcref <- 25
  f[i] <- ftemp_inst_vcmax(tcleaf[i],tcref) #this function is in the rpmodel package
}
curv_trad_temp$fref <- f
curv_trad_temp_adj <- curv_trad_temp %>% 
  mutate(vcmax_25 = V_cmax * fref)


#Working on temperature correction for Photosynthesis package

library(tidyverse)
library(rpmodel)

#CHANGE THIS
pars_photo_dat <- read_csv("~/Documents/GitHub/DAT_Proj/Results/dat_fits_photo_pars_filt_no_TPU.csv")
pars_photo_trad <- read_csv("~/Documents/GitHub/DAT_Proj/Results/trad_fits_photo_pars_no_TPU.csv")

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

f_fact_dat <- data.frame()
for (i in 1:length(curv_dat_temp)){
  tcleaf <- curv_dat_temp$meanTleaf[i]
  tcref <- 25
  f <- ftemp_inst_vcmax(tcleaf,tcref) #this function is in the rpmodel package
  f_fact_dat <- rbind(f_fact_dat, f)
}
colnames(f_fact_dat) <- "f_fact"
curv_dat_temp_adj <- curv_dat_temp %>% 
  mutate(vcmax_25 = V_cmax * f_fact_dat$f_fact)

f_fact_trad <- data.frame()
for (i in 1:length(curv_trad_temp)){
  tcleaf <- curv_trad_temp$meanTleaf[i]
  tcref <- 25
  f <- ftemp_inst_vcmax(tcleaf,tcref) #this function is in the rpmodel package
  f_fact_trad <- rbind(f_fact_trad, f)
}
colnames(f_fact_trad) <- "f_fact"
curv_trad_temp_adj <- curv_trad_temp %>% 
  mutate(vcmax_25 = V_cmax * f_fact_trad$f_fact)


#Working on temperature correction for Photosynthesis package

library(tidyverse)
# library(rpmodel)
# 
# #CHANGE THIS
# pars_photo_dat <- read_csv("~/Documents/GitHub/DAT_Proj/Results/dat_fits_photo_pars_filt_no_TPU.csv")
# pars_photo_trad <- read_csv("~/Documents/GitHub/DAT_Proj/Results/trad_fits_photo_pars_no_TPU.csv")
# 
# cmplt.rm_out <- read_csv("~/Documents/GitHub/DAT_Proj/Inputs/Aci_no_out.csv")
# 
# grp_curv <- cmplt.rm_out %>% 
#   group_by(Data_point, unique) %>% 
#   summarize(meanTleaf = mean(Tleaf))
# grp_curv2 <- grp_curv %>% 
#   rename(ID = unique)
# grp_dat <- filter(grp_curv2, Data_point == "Before_DAT")
# grp_trad <- filter(grp_curv2, Data_point == "Traditional")
# 
# pars_photo_dat$method <- "DAT"
# pars_photo_trad$method <- "Traditional"
# 
# curv_dat_temp <- left_join(by = "ID", pars_photo_dat, grp_dat)
# curv_trad_temp <- left_join(by = "ID", pars_photo_dat, grp_dat)
# 
# f_fact_dat <- data.frame()
# for (i in 1:length(curv_dat_temp)){
#   tcleaf <- curv_dat_temp$meanTleaf[i]
#   tcref <- 25
#   f <- ftemp_inst_vcmax(tcleaf,tcref) #this function is in the rpmodel package
#   f_fact_dat <- rbind(f_fact_dat, f)
# }
# colnames(f_fact_dat) <- "f_fact"
# curv_dat_temp_adj <- curv_dat_temp %>% 
#   mutate(vcmax_25 = V_cmax * f_fact_dat$f_fact)
# 
# f_fact_trad <- data.frame()
# for (i in 1:length(curv_trad_temp)){
#   tcleaf <- curv_trad_temp$meanTleaf[i]
#   tcref <- 25
#   f <- ftemp_inst_vcmax(tcleaf,tcref) #this function is in the rpmodel package
#   f_fact_trad <- rbind(f_fact_trad, f)
# }
# colnames(f_fact_trad) <- "f_fact"
# curv_trad_temp_adj <- curv_trad_temp %>% 
#   mutate(vcmax_25 = V_cmax * f_fact_trad$f_fact)
# 



## Maquelle approach
cmplt.rm_out <- read_csv("~/Documents/GitHub/DAT_Proj/Inputs/Aci_no_out.csv")
pars_photo_dat <- read_csv("~/Documents/GitHub/DAT_Proj/Results/dat_fits_photo_pars_filt_no_TPU.csv")
pars_photo_trad <- read_csv("~/Documents/GitHub/DAT_Proj/Results/trad_fits_photo_pars_no_TPU.csv")

grp_curv <- cmplt.rm_out %>% 
    group_by(Data_point, unique) %>% 
    summarize(meanTleaf = mean(Tleaf)) %>% 
    mutate(meanTleafK = meanTleaf + 273.15)
grp_curv2 <- grp_curv %>% 
    rename(ID = unique)
grp_dat <- filter(grp_curv2, Data_point == "Before_DAT")
grp_trad <- filter(grp_curv2, Data_point == "Traditional")

pars_photo_dat$method <- "DAT"
pars_photo_trad$method <- "Traditional"

curv_dat_temp <- left_join(by = "ID", pars_photo_dat, grp_dat)
curv_trad_temp <- left_join(by = "ID", pars_photo_dat, grp_dat)

meanTleafK    <- mean(Tleafk)

###Constants used in the Farquhar´s model ***NOT*** considering mesophyll conductance
R             <- 0.008314    # Gas constant
R2            <- R * 1000
Kc   		      <- 40.49		   # Michaelis-Menten constant for CO2 (Pa) (Bernacchi et al 2001,2002) or 404.9 microbar von Caemmerer et al. (1994)
delta_Kc 		  <- 79430 	     # (J mol-1) from Medlyn et al 2002
Ko   		      <- 27.84		   # Michaelis-Menten constant for CO2 (kPa)(Bernacchi et al 2001,2002) or 278.4 mbar von Caemmerer et al. (1994)
delta_Ko 		  <- 36380	     # (J mol-1) from Medlyn et al 2002
gstar  		    <- 4.275		   # CO2 compensation point (Pa) (Bernacchi et al 2001,2002) or 36.9 microbar von Caemmerer et al. (1994)
delta_gstar 	<- 37830	     # (J mol-1)

#Other Constants
curv          <- 0.85        # Curvature for calculating J from Evans (1989)
Ambient.CO2   <- 400
O2  	      	<- 21 		     # Estimated Oxygen concentration at chloroplast - (kPa)
K0            <- 0.388       # for accounting for diffusion of CO2 through the gasket (Licor 6400 only) maybe we need to remove this part the new machine is alredy accounting for it

### Peaked function Medlin et al. 2002PCE
### Ea = Activation Energy, DeltaS = entropy factor = 200, Hd = Deactivation Energy)
### updated to Kumarathunge et al 2019, New Phytologist.(222: 768-784) doi: 10.1111/nph.15668
Ea_V     <- 82620.87 #58550 or 82620.87 or ---- these values are slight different from Kamarathunge
Delta_V  <- 645.1013   #629.26 or 645.1013 or
Ed       <- 200000  
b_V <- 1+exp((298.15 * Delta_V - Ed) / (298.15  *R2))

Ea_J     <- 39676.89 # or 39676.89 
Delta_J  <- 641.3615 # or 641.3615
#EdVJ    <- 200000
b_J <- 1 + exp((298.15 * Delta_J - Ed) / (298.15 * R2))

#Temp_Correction <- function (p1,p2,p3) {
    a <- p1/(exp(26.355 - (65.33 /(R*meanTleafK))))#Vcmax as in Bernacchi
    b <- p2/(exp( 17.71 - (43.9  /(R*meanTleafK))))#this is different from Bernacchi
    #c <- p3/((exp(21.46 - 53.1   /(R*meanTleafK)))/(1+exp((0.65*meanTleafK-201.8)/(R*meanTleafK))))
    d <- p3/(exp(18.715 - (46.39 /(R*meanTleafK))))#Rd as in bernacchi 
    
    return(c(a,b,d))
}
Vc_peaked_25C <- function (Vc_LeafTemp,Ea_V,Delta_V,Ed) {
    a <- exp(Ea_V * (meanTleafK - 298.15) / (298.15 * R2 * meanTleafK))
    c <- 1+exp((meanTleafK * Delta_V - Ed) / (meanTleafK * R2))
    Vc_25 <- Vc_LeafTemp / (a * (b_V / c))
}

J_peaked_25C <- function (J_LeafTemp,Ea_J,Delta_J,Ed) {
    a <- exp(Ea_J * (meanTleafK - 298.15)/(298.15 * R2 * meanTleafK))
    c <- 1 + exp((meanTleafK * Delta_J - Ed) / (meanTleafK * R2))
    Jmax_25 <- J_LeafTemp / (a * (b_J / c))
}  

meanTleafK <- curv_dat_temp$meanTleafK
    
best_Vcmax_25C        <- Vc_peaked_25C(       curv_dat_temp[[2]],Ea_V, Delta_V, Ed)
best_Jmax_25C         <- J_peaked_25C(       curv_dat_temp[[4]], Ea_J, Delta_J, Ed)

dat_corrected <- curv_dat_temp %>% 
    mutate(best_Vcmax_25C = best_Vcmax_25C,
           best_Jmax_25C = best_Jmax_25C)

write.csv(dat_corrected, "/Users/charlessouthwick/Documents/GitHub/DAT_Proj/Results/attempt_correct_photo_20220213.csv")

## This script is to graph, fit, and save the ACi data for our Tapajos 2022 campaign

## Load Packages
library(tidyverse)
library(ggpubr)


# Set working directory to DAT_proj
wd <- "C://Users/emmel/Desktop/DAT_proj/"
setwd(wd)



# Identify and Filter Outliers ---------------------------------------------------
## Load Data
complete_sp <- read.csv("Inputs/clean_aci_with_uniquecode.csv", sep = ",", header = TRUE,
                        fileEncoding="latin1") 
complete_sp <- filter(complete_sp, Data_QC == "OK") # All A/Ci curve data


## Identify Outliers and Write new data frame
cmplt.rm_out <- filter(complete_sp, Ci > -5 & A < 40 & A > -1)


# Data frame without outliers
write.csv(x = cmplt.rm_out, file = paste0(getwd(), "/Inputs/Aci_no_out.csv"), 
          row.names = FALSE)



# Load Filtered Data and Split into DAT and Trad  --------------------------
library(greekLetters)

## Load Filtered data
cmplt.rm_out <- read.csv(file = paste0(wd, "Inputs/Aci_no_out.csv"), header = TRUE, sep = ",")


### May not need
## Separate into DF of DAT and Trad
cmplt_DAT <- filter(cmplt.rm_out, Data_point == "Before_DAT") %>% 
  select(-contains(greeks("Delta"))) # removes the columns with deltas
head(cmplt_DAT)

write.csv(x = cmplt_DAT, file = paste0(wd, "/Inputs/all_DAT_data_no_out.csv"), row.names = FALSE)


cmplt_trad <- filter(cmplt.rm_out, Data_point == "Traditional") %>% 
  select(-contains(greeks("Delta")))
cmplt_trad <- as.data.frame(cmplt_trad)
head(cmplt_trad)

write.csv(x = cmplt_trad, file = paste0(wd, "/Inputs/all_trad_data_no_out.csv"), row.names = FALSE)


# Make a function to remove backwards points ------------------------------

## Make function to find min Ci and exclude points with Anet below the Anet for that min Ci
exclude_backwardsCi <- function(data, givedf){
  min_Ci_ind <- which(data$Ci == min(data$Ci))
  data_new <- slice(data, -which(data$A < data$A[min_Ci_ind]))
  if(givedf =="TRUE"){
    data_new <- as.data.frame(data_new)
  }
  return(data_new)
}


# # Apply function to tibble. if .keep = TRUE this throws an error
# # Additional info: https://stackoverflow.com/questions/63412850/managing-dplyr-group-by-function-to-keep-the-grouping-variable-when-used-in-comb
# DAT_filt_ex <- DAT_filt %>%
#   group_by(unique) %>%
#   group_modify(~exclude_backwardsCi(data = .x, givedf = TRUE), .keep = FALSE)
# DAT_filt_ex <- as.data.frame(DAT_filt_ex)





# Plotting ACi Curves -----------------------------------------------------

cmplt.grp <- group_by(cmplt.rm_out, k67.id)


## Plot all ACi curves on one graph by tree
ggplot(cmplt.grp, mapping = aes(x = Ci, y = A, color = unique)) +
  geom_point(mapping = aes(pch = Data_point)) +
  theme_classic()


# Make and save plots for each leaf
for (id in unique(cmplt.grp$unique)) {
    df1 <- cmplt.grp %>% filter(unique == id)
    gg1 <- ggplot(data = df1, mapping = aes(x = Ci, y = A, shape = Data_point, color = Data_point)) +
        geom_point(size = 3) +
        theme_classic() +
        labs(y = expression("A"[net]*" ??mol "*m^{-2}*" "*s^{-1}), 
             x = expression("C"[i]*" ??mol "* mol^{-1})) +
        scale_shape_manual(name = "Method", labels = c("DAT", "Steady-State"), values = c(19, 17)) +
        scale_color_viridis_d(begin = 0.3, name = "Method", labels = c("DAT", "Steady-State")) +
        ggtitle(id)
    plot(gg1)
    filename1 <- paste0("plot_", id, ".png")
    ggsave(filename1, gg1, path = paste0(getwd(), "/Figures/"))
}



# Fitting ACi Curves with Plantecophys ------------------------------------------------------


library(plantecophys)


cmplt.rm_out <- read.csv(file = paste0(wd, "Inputs/Aci_no_out.csv"), header = TRUE, sep = ",")
cmplt_DAT <- filter(cmplt.rm_out, Data_point == "Before_DAT")
cmplt_trad <- filter(cmplt.rm_out, Data_point == "Traditional")


DAT_filt <- cmplt_DAT %>%
  group_by(unique) %>%
  group_modify(~exclude_backwardsCi(data = .x, givedf = TRUE), .keep = FALSE)
DAT_filt <- as.data.frame(DAT_filt)




## Remove K6714L1, since it acts weird with fitacis
DAT_filt_ex <- filter(DAT_filt_ex, unique != "K6714L1") %>% 
    as.data.frame()


## Fit the ACi curves for each species for DAT using fitacis
DAT_fits_ecophys <- fitacis(DAT_filt_ex, group = "unique", fitmethod = "bilinear",
                            varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                            PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(DAT_fits_ecophys[[23]], main = coef(DAT_fits_ecophys)$unique[[23]])
coef(DAT_fits_ecophys)

### Run K6706L1 separately, since it gives a weird curve
k6714l1 <- filter(DAT_filt, unique == "K6714L1")
k6714l1_fit <- fitaci(k6714l1, fitmethod = "bilinear", 
                      varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                      PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE,
                      citransition = 200) # The Ci transition is specified as 200, as per Sharkey's
                                          # recommendations and to avoid an unreasonable Jmax value
plot(k6714l1_fit)
coef(k6714l1_fit)



# PDF file of all plots
pdf(file = paste0(wd,"Figures/dataci_ecophys_no_TPU.pdf"), height=10, width=20)
plot.new()
for (curve in 1:33){
  title <- coef(DAT_fits_ecophys)$unique[[curve]]
  plot(DAT_fits_ecophys[[curve]], main = title)
}
plot(k6714l1_fit, main = "K6714L1")
dev.off()


## Make a dataframe out of coefficients
par_dat <- as.data.frame(coef(DAT_fits_ecophys), row.names = NULL)
par_dat[23,] <- c("K6714L1", 20.8393492, 13.5285784, 0.1658439, 0.47818072, NA, 0.02406953)
#par_dat <- par_dat[-1,]
par_dat <- par_dat %>%
  add_column(method = "dat")
par_dat$Vcmax <- as.double(par_dat$Vcmax)
par_dat$Jmax <- as.double(par_dat$Jmax)
par_dat$Rd <- as.double(par_dat$Rd)
par_dat$Vcmax_SE <- as.double(par_dat$Vcmax_SE)
par_dat$Jmax_SE <- as.double(par_dat$Jmax_SE)
par_dat$Rd_SE <- as.double(par_dat$Rd_SE)
table(par_dat$method) ## number of initial DAT curves






# Fit the ACi curves for Traditional using fitacis
trad_fits_ecophys <- fitacis(cmplt_trad, group = "unique", fitmethod = "bilinear",
                             varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                             PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(trad_fits_ecophys[[14]], main = coef(trad_fits_ecophys)$unique[14])
coef(trad_fits_ecophys)


# Save all plots to a pdf
pdf(file = paste0(wd,"Figures/tradaci_ecophys.pdf"), height=10, width=20)
plot.new()
for (curve in 1:28){
  title <- coef(trad_fits_ecophys)$unique[[curve]]
  plot(trad_fits_ecophys[[curve]], main = title)
}
dev.off()





## Make a dataframe of coefficients
par_trad <- as.data.frame(coef(trad_fits_ecophys), row.names = NULL)
#par_trad <- par_trad[-1,]
par_trad <- par_trad %>% 
  add_column(method = "trad")





# Merge DAT and Trad coef dfs
par_join <- bind_rows(par_dat, par_trad)
head(par_join)


write.csv(x = par_join, file = paste0(wd, "/Results/params_ecophys_no_TPU.csv"), row.names = FALSE)





# Fitting A/Ci curves with photosynthesis ----------------------------------

library(photosynthesis)

cmplt.rm_out <- read.csv(file = paste0(wd, "Inputs/Aci_no_out.csv"), header = TRUE, sep = ",")
cmplt_DAT <- filter(cmplt.rm_out, Data_point == "Before_DAT")
cmplt_trad <- filter(cmplt.rm_out, Data_point == "Traditional")


DAT_filt_ex <- cmplt_DAT %>%
  group_by(unique) %>%
  group_modify(~exclude_backwardsCi(data = .x, givedf = TRUE), .keep = FALSE)
DAT_filt_ex <- as.data.frame(DAT_filt_ex)


# Convert leaf temperature to Kelvin
cmplt_DAT$Tleaf <- cmplt_DAT$Tleaf + 273
cmplt_trad$Tleaf <- cmplt_trad$Tleaf + 273
DAT_filt_ex$Tleaf <- DAT_filt_ex$Tleaf + 273
cmplt_trad <- as.data.frame(cmplt_trad)



# fitting many curves at a time: NO TPU --------------------------------------------------
trad_fits_photo_noTPU <- fit_many(data = cmplt_trad, fitTPU = FALSE,
                            varnames = list(A_net = "A", T_leaf = "Tleaf", C_i = "Ci", PPFD = "Qin"), 
                            funct = fit_aci_response,
                            group = "unique")
trad_fits_photo_graphs_noTPU <- compile_data(trad_fits_photo_noTPU,
                                       list_element = 2)
trad_fits_photo_pars_noTPU <- compile_data(trad_fits_photo_noTPU,
                                     output_type = "dataframe",
                                     list_element = 1)


dat_fits_photo_noTPU <- fit_many(data = DAT_filt_ex, fitTPU = FALSE,
                           varnames = list(A_net = "A", T_leaf = "Tleaf", C_i = "Ci", PPFD = "Qin"), 
                           funct = fit_aci_response,
                           group = "unique")
dat_fits_photo_graphs_noTPU <- compile_data(dat_fits_photo_noTPU,
                                      list_element = 2)
dat_fits_photo_pars_noTPU <- compile_data(dat_fits_photo_noTPU,
                                    output_type = "dataframe",
                                    list_element = 1)

write.csv(x = trad_fits_photo_pars_noTPU, file = paste0(wd, "Results/trad_fits_photo_pars_no_TPU.csv"),
          row.names = FALSE)

write.csv(x = dat_fits_photo_pars_noTPU, file = paste0(wd, "Results/dat_fits_photo_pars_filt_no_TPU.csv"),
          row.names = FALSE)



#Write PDF of outputs

pdf(file = paste0(wd,"Figures/trad_fits_photo_figs_no_TPU.pdf"))
plot.new()
for (curve in 1:nrow(trad_fits_photo_pars_noTPU)){
    plot(trad_fits_photo_graphs_noTPU[[curve]])
}
dev.off()


pdf(file = paste0(wd,"Figures/dat_fits_photo_figs_filt_no_TPU.pdf"))
plot.new()
for (curve in 1:nrow(dat_fits_photo_pars_noTPU)){
    plot(dat_fits_photo_graphs_noTPU[[curve]])
}
dev.off()


# fitting many curves at a time: WITH TPU --------------------------------------------------
trad_fits_photo <- fit_many(data = cmplt_trad, fitTPU = TRUE,
                            varnames = list(A_net = "A", T_leaf = "Tleaf", C_i = "Ci", PPFD = "Qin"), 
                            funct = fit_aci_response,
                            group = "unique")
trad_fits_photo_graphs <- compile_data(trad_fits_photo,
                                       list_element = 2)
trad_fits_photo_pars <- compile_data(trad_fits_photo,
                                     output_type = "dataframe",
                                     list_element = 1)


dat_fits_photo <- fit_many(data = DAT_filt_ex, fitTPU = TRUE,## uses the back-filtered data
                           varnames = list(A_net = "A", T_leaf = "Tleaf", C_i = "Ci", PPFD = "Qin"), 
                           funct = fit_aci_response,
                           group = "unique")
dat_fits_photo_graphs <- compile_data(dat_fits_photo,
                                      list_element = 2)
dat_fits_photo_pars <- compile_data(dat_fits_photo,
                                    output_type = "dataframe",
                                    list_element = 1)



write.csv(x = trad_fits_photo_pars, file = paste0(wd, "Results/trad_fits_photo_pars.csv"),
          row.names = FALSE)

write.csv(x = dat_fits_photo_pars, file = paste0(wd, "Results/dat_fit_ex_photo_pars.csv"),
          row.names = FALSE)



#Write PDF of outputs

pdf(file = paste0(wd,"Figures/trad_fits_photo_figs_with_TPU.pdf"))
plot.new()
for (curve in 1:nrow(trad_fits_photo_pars)){
    plot(trad_fits_photo_graphs[[curve]])
}
dev.off()

pdf(file = paste0(wd,"Figures/dat_fits_photo_figs_filt_with_TPU.pdf"))
plot.new()
for (curve in 1:nrow(dat_fits_photo_pars)){ 
    plot(dat_fits_photo_graphs[[curve]])
}
dev.off()


#Temperature correction: For NO TPU data -------------------------------
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
curv_trad_temp <- left_join(by = "ID", pars_photo_trad, grp_trad)

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

best_Vcmax_25C_dat        <- Vc_peaked_25C(curv_dat_temp[[2]],Ea_V, Delta_V, Ed)
best_Jmax_25C_dat         <- J_peaked_25C(curv_dat_temp[[4]], Ea_J, Delta_J, Ed)

meanTleafK <- curv_trad_temp$meanTleafK

best_Vcmax_25C_trad        <- Vc_peaked_25C(curv_trad_temp[[2]],Ea_V, Delta_V, Ed)
best_Jmax_25C_trad         <- J_peaked_25C(curv_trad_temp[[4]], Ea_J, Delta_J, Ed)

dat_corrected <- curv_dat_temp %>% 
    mutate(Best_Vcmax_25C = best_Vcmax_25C_dat,
           Best_Jmax_25C = best_Jmax_25C_dat)

trad_corrected <- curv_trad_temp %>% 
    mutate(Best_Vcmax_25C = best_Vcmax_25C_trad,
           Best_Jmax_25C = best_Jmax_25C_trad)

#Write csvs

write.csv(x = dat_corrected, file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_no_TPU.csv"),
          row.names = FALSE)

write.csv(x = trad_corrected, file = paste0(wd, "Results/trad_fits_photo_pars_correct_no_TPU.csv"),
          row.names = FALSE)


#Temperature correction: For WITH TPU data -------------------------------
cmplt.rm_out <- read_csv("~/Documents/GitHub/DAT_Proj/Inputs/Aci_no_out.csv")
pars_photo_dat <- read_csv("~/Documents/GitHub/DAT_Proj/Results/dat_fit_ex_photo_pars.csv")
pars_photo_trad <- read_csv("~/Documents/GitHub/DAT_Proj/Results/trad_fits_photo_pars.csv")

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
curv_trad_temp <- left_join(by = "ID", pars_photo_trad, grp_trad)

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

best_Vcmax_25C_dat        <- Vc_peaked_25C(curv_dat_temp[[2]],Ea_V, Delta_V, Ed)
best_Jmax_25C_dat         <- J_peaked_25C(curv_dat_temp[[4]], Ea_J, Delta_J, Ed)

meanTleafK <- curv_trad_temp$meanTleafK

best_Vcmax_25C_trad        <- Vc_peaked_25C(curv_trad_temp[[2]],Ea_V, Delta_V, Ed)
best_Jmax_25C_trad         <- J_peaked_25C(curv_trad_temp[[4]], Ea_J, Delta_J, Ed)

dat_corrected <- curv_dat_temp %>% 
    mutate(Best_Vcmax_25C = best_Vcmax_25C_dat,
           Best_Jmax_25C = best_Jmax_25C_dat)

trad_corrected <- curv_trad_temp %>% 
    mutate(Best_Vcmax_25C = best_Vcmax_25C_trad,
           Best_Jmax_25C = best_Jmax_25C_trad)

#Write csvs

write.csv(x = trad_corrected, file = paste0(wd, "Results/trad_fits_photo_pars_correct_with_TPU.csv"),
          row.names = FALSE)

write.csv(x = dat_corrected, file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_with_TPU.csv"),
          row.names = FALSE)


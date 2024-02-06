## This script is to graph, fit, and save the ACi data for our Tapajos 2022 campaign

## Load Packages
library(tidyverse)
library(ggpubr)
library(here)


# # Set working directory to DAT_proj -- use here package instead
# wd <- "C://Users/emmel/Desktop/DAT_proj/"
# setwd(wd)




# Load Data and Filter DAT and SS ------------------------------------------------

cmplt.rm_out <- read.csv(file = here("3_Clean_data/clean_aci_noOutliers.csv"),
                         header = TRUE, sep = ",", fileEncoding = "latin1")
cmplt_DAT <- filter(cmplt.rm_out, curv_meth == "DAT")
cmplt_SS <- filter(cmplt.rm_out, curv_meth == "SS")



# Make and use a function to remove backwards points ------------------------------

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
#   group_by(unique_id) %>%
#   group_modify(~exclude_backwardsCi(data = .x, givedf = TRUE), .keep = FALSE)
# DAT_filt_ex <- as.data.frame(DAT_filt_ex)



### Now filter out the 'backwards' points in the DAT data

DAT_filt <- cmplt_DAT %>%
    group_by(unique_id) %>%
    group_modify(~exclude_backwardsCi(data = .x, givedf = TRUE), .keep = FALSE)
DAT_filt <- as.data.frame(DAT_filt)



# Fitting ACi Curves with Plantecophys ----------------------------------------------------

library(plantecophys)


## Remove K6714L1, since it acts weird with fitacis
DAT_filt_ex <- filter(DAT_filt, unique_id != "K6714L1") %>% 
    as.data.frame()


## Fit the ACi curves for each species for DAT using fitacis
DAT_fits_ecophys <- fitacis(DAT_filt_ex, group = "unique_id", fitmethod = "bilinear",
                            varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                            PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)


### Run K6706L1 separately, since it gives a weird curve
k6714l1 <- filter(DAT_filt, unique_id == "K6714L1")
k6714l1_fit <- fitaci(k6714l1, fitmethod = "bilinear", 
                      varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                      PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE,
                      # The Ci transition is specified as 200, as per Sharkey's
                      # recommendations and to avoid an unreasonable Jmax value
                      citransition = 200) 



# PDF file of all plots
pdf(file = here("6_Figures/DAT_ecophys_figs_noTPU.pdf"))
plot.new()
for (curve in 1:length(DAT_fits_ecophys)){
  title <- coef(DAT_fits_ecophys)$unique_id[[curve]]
  plot(DAT_fits_ecophys[[curve]], main = title)
}
plot(k6714l1_fit, main = "K6714L1")
dev.off()


## Make a dataframe out of coefficients
par_dat <- as.data.frame(coef(DAT_fits_ecophys), row.names = NULL)
k6714l1_pars <- matrix(data = c(k6714l1_fit$pars[,1], k6714l1_fit$pars[,2]), 
                       nrow = 1, ncol = 6)
par_dat[33,] <- c("K6714L1", k6714l1_pars)
par_dat <- par_dat %>%
  add_column(curv_meth = "DAT")
par_dat$Vcmax <- as.double(par_dat$Vcmax)
par_dat$Jmax <- as.double(par_dat$Jmax)
par_dat$Rd <- as.double(par_dat$Rd)
par_dat$Vcmax_SE <- as.double(par_dat$Vcmax_SE)
par_dat$Jmax_SE <- as.double(par_dat$Jmax_SE)
par_dat$Rd_SE <- as.double(par_dat$Rd_SE)
table(par_dat$curv_meth) ## number of initial DAT curves






# Fit the ACi curves for SS using fitacis

SS_fits_ecophys <- fitacis(cmplt_SS, group = "unique_id", fitmethod = "bilinear",
                             varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                             PPFD = "Qin"),
                           fitTPU = FALSE, Tcorrect = TRUE)

# Save all plots to a pdf
pdf(file = here("6_Figures/SS_ecophys_figs_noTPU.pdf"), height=10, width=20)
plot.new()
for (curve in 1:28){
  title <- coef(SS_fits_ecophys)$unique_id[[curve]]
  plot(SS_fits_ecophys[[curve]], main = title)
}
dev.off()





## Make a dataframe of coefficients
par_SS <- as.data.frame(coef(SS_fits_ecophys), row.names = NULL)
par_SS <- par_SS %>% 
  add_column(curv_meth = "SS")




# Merge DAT and SS coef dfs
par_join <- bind_rows(par_dat, par_SS)
head(par_join)


write.csv(x = par_join, file = here("5_Results/pars_ecophys_noTPU.csv"),
          row.names = FALSE)


### Since plantecophys does not work for DAT curves, we did not run with TPU



# Fitting A/Ci curves with photosynthesis ----------------------------------

library(photosynthesis)

### Warning: fitting all DAT curves with photosynthesis takes several hours; including TPU fitting takes longer

# Convert leaf temperature to Kelvin
cmplt_SS$Tleaf <- cmplt_SS$Tleaf + 273.15
cmplt_SS <- as.data.frame(cmplt_SS)
DAT_filt$Tleaf <- DAT_filt$Tleaf + 273.15
DAT_filt <- as.data.frame(DAT_filt)



# Fitting all curves WITHOUT TPU --------------------------------------------------

### Fit SS
SS_fits_photo_noTPU <- fit_many(data = cmplt_SS, 
                                fitTPU = FALSE,
                            varnames = list(A_net = "A", T_leaf = "Tleaf",
                                            C_i = "Ci", PPFD = "Qin"), 
                            funct = fit_aci_response,
                            group = "unique_id")

### Save the output in case we need to do some editing
saveRDS(SS_fits_photo_noTPU, file = here("5_Results/SS_photo_fits_noTPU.RData"))
#### readRDS() to load in


### Pull out the parameters and figures
SS_fits_photo_graphs_noTPU <- compile_data(SS_fits_photo_noTPU,
                                       list_element = 2)

SS_fits_photo_pars_noTPU <- compile_data(SS_fits_photo_noTPU,
                                     output_type = "dataframe",
                                     list_element = 1)




### Fit DAT
dat_fits_photo_noTPU <- fit_many(data = DAT_filt,
                                 fitTPU = FALSE,
                           varnames = list(A_net = "A", T_leaf = "Tleaf",
                                           C_i = "Ci", PPFD = "Qin"), 
                           funct = fit_aci_response,
                           group = "unique_id")

### Save the output in case we need to do some editing
saveRDS(dat_fits_photo_noTPU, file = here("5_Results/DAT_photo_fits_noTPU.RData"))




### Pull out the parameters and figures
dat_fits_photo_graphs_noTPU <- compile_data(dat_fits_photo_noTPU,
                                      list_element = 2)

dat_fits_photo_pars_noTPU <- compile_data(dat_fits_photo_noTPU,
                                    output_type = "dataframe",
                                    list_element = 1)




# Make .pdf files with all fitted figures

### Add the species names back in
ids <- read.csv(here("3_Clean_data/id_codebook.csv")) %>% 
    rename(treeid = ï..treeid)

dat_fits_photo_pars_noTPU <- dat_fits_photo_pars_noTPU %>% 
    mutate(treeid = substring(ID, 1, 5)) %>% 
    left_join(., ids, by = "treeid")

SS_fits_photo_pars_noTPU <- SS_fits_photo_pars_noTPU %>% 
    mutate(treeid = substring(ID, 1, 5)) %>% 
    left_join(., ids, by = "treeid")


### Make the PDF files with all fitted figures
pdf(file = here("6_Figures/SS_photo_fits_noTPU.pdf"))
plot.new()
for (curve in 1:nrow(SS_fits_photo_pars_noTPU)){
    plot(SS_fits_photo_graphs_noTPU[[curve]])
    title(paste0(SS_fits_photo_pars_noTPU$gen_spec_id[curve], 
                 ", Leaf ", 
                 substring(SS_fits_photo_pars_noTPU$ID[curve], 6)),
          adj = 0,
          line = -27,
          cex.main = 1.5)
}
dev.off()


pdf(file = here("6_Figures/DAT_photo_figs_noTPU.pdf"))
plot.new()
for (curve in 1:nrow(dat_fits_photo_pars_noTPU)){
    plot(dat_fits_photo_graphs_noTPU[[curve]])
    title(paste0(dat_fits_photo_pars_noTPU$gen_spec_id[curve], 
                 ", Leaf ", 
                 substring(dat_fits_photo_pars_noTPU$ID[curve], 6)),
          adj = 0,
          line = -27,
          cex.main = 1.5)
}
dev.off()


# Fitting all curves WITH TPU --------------------------------------------------

### SS curves
SS_fits_photo_TPU <- fit_many(data = cmplt_SS, fitTPU = TRUE,
                            varnames = list(A_net = "A", T_leaf = "Tleaf",
                                            C_i = "Ci", PPFD = "Qin"), 
                            funct = fit_aci_response,
                            group = "unique_id")
SS_fits_photo_graphs_TPU <- compile_data(SS_fits_photo_TPU,
                                       list_element = 2)
SS_fits_photo_pars_TPU <- compile_data(SS_fits_photo_TPU,
                                     output_type = "dataframe",
                                     list_element = 1)

### Save the output in case we need to do some editing
saveRDS(SS_fits_photo_TPU, file = here("5_Results/SS_photo_fits_TPU.RData"))




### DAT curves
dat_fits_photo_TPU <- fit_many(data = DAT_filt_ex, fitTPU = TRUE,
                           varnames = list(A_net = "A", T_leaf = "Tleaf",
                                           C_i = "Ci", PPFD = "Qin"), 
                           funct = fit_aci_response,
                           group = "unique_id")
dat_fits_photo_graphs_TPU <- compile_data(dat_fits_photo_TPU,
                                      list_element = 2)
dat_fits_photo_pars_TPU <- compile_data(dat_fits_photo_TPU,
                                    output_type = "dataframe",
                                    list_element = 1)


### Save the output in case we need to do some editing
saveRDS(dat_fits_photo_TPU, file = here("5_Results/DAT_photo_fits_TPU.RData"))




# Make .pdf files with all fitted figures

### Add the species names back in
ids <- read.csv(here("3_Clean_data/id_codebook.csv")) %>% 
    rename(treeid = ï..treeid)

dat_fits_photo_pars_TPU <- dat_fits_photo_pars_TPU %>% 
    mutate(treeid = substring(ID, 1, 5)) %>% 
    left_join(., ids, by = "treeid")

SS_fits_photo_pars_TPU <- SS_fits_photo_pars_TPU %>% 
    mutate(treeid = substring(ID, 1, 5)) %>% 
    left_join(., ids, by = "treeid")


### Make the PDF files with all fitted figures
pdf(file = here("6_Figures/SS_photo_fits_TPU.pdf"))
plot.new()
for (curve in 1:nrow(SS_fits_photo_pars_TPU)){
    plot(SS_fits_photo_graphs_TPU[[curve]])
    title(paste0(SS_fits_photo_pars_TPU$gen_spec_id[curve], 
                 ", Leaf ", 
                 substring(SS_fits_photo_pars_TPU$ID[curve], 6)),
          adj = 0,
          line = -27,
          cex.main = 1.5)
}
dev.off()


pdf(file = here("6_Figures/DAT_photo_figs_TPU.pdf"))
plot.new()
for (curve in 1:nrow(dat_fits_photo_pars_TPU)){
    plot(dat_fits_photo_graphs_TPU[[curve]])
    title(paste0(dat_fits_photo_pars_TPU$gen_spec_id[curve], 
                 ", Leaf ", 
                 substring(dat_fits_photo_pars_TPU$ID[curve], 6)),
          adj = 0,
          line = -27,
          cex.main = 1.5)
}
dev.off()


# Photosynthesis Temperature Correction; Constants and Functions ----------

### Photosynthesis does not perform it's own temperature correction, so we apply our own

### Constants used in the Farquhar´s model ***NOT*** considering mesophyll conductance
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


# Temperature correction: For NO TPU data -------------------------------

### Calculate the mean Tleaf for each curve

grp_dat <- DAT_filt %>% 
    group_by(unique_id) %>% 
    summarize(meanTleafK = mean(Tleaf)) %>% 
    rename(ID = unique_id)

grp_SS <- cmplt_SS %>% 
    group_by(unique_id) %>% 
    summarize(meanTleafK = mean(Tleaf)) %>% 
    rename(ID = unique_id)

# grp_curv <- cmplt.rm_out %>% 
#     ### Group by curve type and unique_id
#     group_by(curv_meth, unique_id) %>% 
#     ### Calculate the mean Tleaf per curve
#     summarize(meanTleaf = mean(Tleaf)) %>% 
#     ### Convert to Kelvin
#     mutate(meanTleafK = meanTleaf + 273.15) %>% 
#     rename(ID = unique_id)
# ### Filter by curve type
# grp_dat <- filter(grp_curv, curv_meth == "DAT")
# grp_SS <- filter(grp_curv, curv_meth == "SS")

### Joining the mean Tleaf in Kelvin to fitted parameters



curv_dat_temp <- left_join(by = "ID",
                           dat_fits_photo_pars_noTPU,
                           grp_dat)
curv_SS_temp <- left_join(by = "ID", 
                          SS_fits_photo_pars_noTPU, 
                          grp_SS)



 
### Run the temperature correction on the DAT data
meanTleafK <- curv_dat_temp$meanTleafK

best_Vcmax_25C_dat <- Vc_peaked_25C(curv_dat_temp[[2]], Ea_V, Delta_V, Ed)
best_Jmax_25C_dat <- J_peaked_25C(curv_dat_temp[[4]], Ea_J, Delta_J, Ed)

dat_corrected <- curv_dat_temp %>% 
    mutate(Best_Vcmax_25C = best_Vcmax_25C_dat,
           Best_Jmax_25C = best_Jmax_25C_dat,
           curv_meth = "DAT") ####### Check whether you use curv_meth or method in following scripts


### Run the temperature correction on the SS data
meanTleafK <- curv_SS_temp$meanTleafK

best_Vcmax_25C_SS <- Vc_peaked_25C(curv_SS_temp[[2]],Ea_V, Delta_V, Ed)
best_Jmax_25C_SS <- J_peaked_25C(curv_SS_temp[[4]], Ea_J, Delta_J, Ed)

SS_corrected <- curv_SS_temp %>% 
    mutate(Best_Vcmax_25C = best_Vcmax_25C_SS,
           Best_Jmax_25C = best_Jmax_25C_SS,
           curv_meth = "SS")




# Write .csv files for corrected data
write.csv(x = dat_corrected, 
          file = here("5_Results/DAT_photo_pars_crct_noTPU.csv"),
          row.names = FALSE)

write.csv(x = SS_corrected, 
          file = here("5_Results/SS_photo_pars_crct_noTPU.csv"),
          row.names = FALSE)




# Temperature correction: For WITH TPU data -------------------------------

### Calculate the mean Tleaf for each curve

grp_dat <- DAT_filt %>% 
    group_by(unique_id) %>% 
    summarize(meanTleafK = mean(Tleaf)) %>% 
    rename(ID = unique_id)

grp_SS <- cmplt_SS %>% 
    group_by(unique_id) %>% 
    summarize(meanTleafK = mean(Tleaf)) %>% 
    rename(ID = unique_id)



### Joining the mean Tleaf in Kelvin to fitted parameters

curv_dat_temp <- left_join(by = "ID",
                           DAT_fits_photo_pars_TPU, 
                           grp_dat)
curv_SS_temp <- left_join(by = "ID",
                          SS_fits_photo_pars_TPU,
                          grp_SS)



### Run the temperature correction on the SS data

meanTleafK <- curv_SS_temp$meanTleafK

best_Vcmax_25C_SS <- Vc_peaked_25C(curv_SS_temp[[2]],Ea_V, Delta_V, Ed)
best_Jmax_25C_SS <- J_peaked_25C(curv_SS_temp[[4]], Ea_J, Delta_J, Ed)

SS_corrected <- curv_SS_temp %>% 
    mutate(Best_Vcmax_25C = best_Vcmax_25C_SS,
           Best_Jmax_25C = best_Jmax_25C_SS,
           curv_meth = "SS")



### Run the temperature correction on the DAT data

meanTleafK <- curv_dat_temp$meanTleafK

best_Vcmax_25C_dat <- Vc_peaked_25C(curv_dat_temp[[2]], Ea_V, Delta_V, Ed)
best_Jmax_25C_dat <- J_peaked_25C(curv_dat_temp[[4]], Ea_J, Delta_J, Ed)

dat_corrected <- curv_dat_temp %>% 
    mutate(Best_Vcmax_25C = best_Vcmax_25C_dat,
           Best_Jmax_25C = best_Jmax_25C_dat,
           curv_meth = "DAT") ####### Check whether you use curv_meth or method in following scripts




# Write .csv files for corrected data
write.csv(x = dat_corrected, 
          file = here("5_Results/DAT_photo_pars_crct_TPU.csv"),
          row.names = FALSE)

write.csv(x = SS_corrected, 
          file = here("5_Results/SS_photo_pars_crct_TPU.csv"),
          row.names = FALSE)


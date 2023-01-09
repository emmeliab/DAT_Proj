## This script is to graph, fit, and save the ACi data for our Tapajos 2022 campaign

## Load Packages
library(tidyverse)
library(greekLetters)
library(ggpubr)


# Set working directory to DAT_proj
wd <- "C://Users/emmel/Desktop/DAT_proj/"
setwd(wd)



# Identify and Filter Outliers ---------------------------------------------------
## Load Data
complete_sp <- read.csv("Inputs/clean_aci_with_uniquecode.csv", sep = ",", header = TRUE,
                        fileEncoding="latin1") %>% 
  filter(complete_sp, Data_QC == "OK") # All A/Ci curve data
unique_ids <- read.csv("Inputs/unique_ids.csv") # ID data table


## Identify Outliers and Write new data frame
which(complete_sp$Ci < -50) #4 Ci values are super negative: 394, 395, 396, 398
which(complete_sp$Ci < -5) # Add in 16 more values, all for MACA1
complete_sp[c(394,395,396,398),19] # What are those values? Should they be outliers?
which(complete_sp$A > 40)
which(complete_sp$A < -1)
cmplt.rm_out1 <- filter(complete_sp, Ci > -5)
cmplt.rm_out2 <- filter(cmplt.rm_out1, A < 40) ## A < 40
cmplt.rm_out <- filter(cmplt.rm_out2, A > -1)

# Data frame without outliers
write.csv(x = cmplt.rm_out, file = paste0(getwd(), "/Inputs/Aci_no_out.csv"), 
          row.names = FALSE)



# Load Filtered Data and Split into DAT and Trad --------------------------
## Load Filtered data
cmplt.rm_out <- read.csv(file = paste0(wd, "Inputs/Aci_no_out.csv"), header = TRUE, sep = ",")


## Separate into df of DAT and Trad
cmplt_DAT <- filter(cmplt.rm_out, Data_point == "Before_DAT") %>% 
  select(-contains(greeks("Delta"))) # removes the columns with deltas
head(cmplt_DAT)

cmplt_trad <- filter(cmplt.rm_out, Data_point == "Traditional") %>% 
  select(-contains(greeks("Delta")))
cmplt_trad <- as.data.frame(cmplt_trad)
head(cmplt_trad)


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

cmplt.grp <- groupby(cmplt.rm_out, k67.id)


## Plot all ACi curves on one graph by Species
ggplot(cmplt.grp, mapping = aes(x = Ci, y = A, color = k67.id)) +
  geom_point(mapping = aes(pch = Data_point)) +
  theme_classic()

## Make and save plots for each individual species
for (code in unique(cmplt.grp$k67.id)) {
  df1 <- cmplt.grp %>% filter(k67.id == code)
  gg1 <- ggplot(data = df1, mapping = aes(x = Ci, y = A, color = Data_point)) +
    geom_point() +
    theme_classic() +
    scale_color_viridis_d() +
    ggtitle(code)
  plot(gg1)
  filename1 <- paste("plot_", code, ".png")
  ggsave(filename1, gg1, path = paste0(getwd(), "/Figures/"))
}

#The same, but for each leaf
for (id in unique(cmplt.grp$unique)) {
  df1 <- cmplt.grp %>% filter(unique == id)
  gg1 <- ggplot(data = df1, mapping = aes(x = Ci, y = A, color = Data_point)) +
    geom_point() +
    theme_classic() +
    scale_color_viridis_d() +
    ggtitle(id)
  plot(gg1)
  filename1 <- paste0("plot_", id, ".png")
  ggsave(filename1, gg1, path = paste0(getwd(), "/Figures/"))
}



# Fitting ACi Curves with Plantecophys ------------------------------------------------------


library(plantecophys)


cmplt.rm_out <- read.csv(file = paste0(wd, "Inputs/Aci_no_out.csv"), header = TRUE, sep = ",")
cmply_DAT <- filter(cmplt.rm_out, Data_point == "Before_DAT")
cmplt_trad <- filter(cmplt.rm_out, Data_point == "Traditional")



## Remove K6709L2-2 cause it is giving us trouble; will fit separately
DAT_filt <- filter(cmplt_DAT, unique != "K6709L2-2" & unique != "K6714L2" & unique != "K6718L2") %>% 
  as.data.frame()
write.csv(x = DAT_filt, file = paste0(getwd(), "/Inputs/filtered_DAT_data.csv"), 
          row.names = FALSE)
DAT09l22 <- filter(cmplt_DAT, unique == "K6709L2-2") %>% 
  as.data.frame()
DAT14L2 <- filter(cmplt_DAT, unique == "K6714L2") %>% 
  as.data.frame()
DAT18L2 <- filter(cmplt_DAT, unique == "K6718L2") %>% 
  as.data.frame()


DAT_filt_ex <- DAT_filt %>%
  group_by(unique) %>%
  group_modify(~exclude_backwardsCi(data = .x, givedf = TRUE), .keep = FALSE)
DAT_filt_ex <- as.data.frame(DAT_filt_ex)




## Fit the ACi curves for each species for DAT using fitacis
DAT_fits_ecophys <- fitacis(DAT_filt_ex, group = "unique", id = "unique",
                            varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                            PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(DAT_fits_ecophys[[23]], main = coef(DAT_fits_ecophys)$unique[[23]]) ##keep an eye on #7 as the example
coef(DAT_fits_ecophys)

#For loop to save all the plots
for (curve in 1:33){
  title <- coef(DAT_fits_ecophys)$unique[[curve]]
  png(filename = paste0(getwd(), "/Figures/", title,"_dataci_curve.png"))
  plot(DAT_fits_ecophys[[curve]], main = title)
  dev.off()
}


#strange curves: #7, 9, 10, 11, 12, 15, 17, 22, 23, 24? 

## Make a dataframe out of coefficients
par_dat <- as.data.frame(coef(DAT_fits_ecophys), row.names = NULL)
par_dat <- par_dat[-1,]
par_dat <- par_dat %>%
  add_column(method = "dat")
table(par_dat$method) ## number of initial DAT curves


plot(DAT_fits_ecophys[[15]], main = par_dat$unique[14])






# Fit the ACi curves for Traditional using fitacis
trad_fits_ecophys <- fitacis(cmplt_trad, group = "unique", #id = "unique",
                             varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                             PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(trad_fits_ecophys[[14]], main = coef(trad_fits_ecophys)$unique[14]) # 20 is pretty ugly
coef(trad_fits_ecophys)

# For loop to plot all the curves
for (curve in 1:28){
  title <- coef(trad_fits_ecophys)$unique[[curve]]
  png(filename = paste0(getwd(), "/Figures/", title,"_tradaci_curve.png"))
  plot(trad_fits_ecophys[[curve]], main = title)
  dev.off()
}


pdf(file = paste0(wd,"Figures/tradaci_ecophys.pdf"), height=10, width=20)
plot.new()
for (curve in 1:28){
  title <- coef(trad_fits_ecophys)$unique[[curve]]
  plot(trad_fits_ecophys[[curve]], main = title)
  #text(30, 5, labels = as.character(title))
}
dev.off()

write.csv(x = trad_fits_ecophys, file = paste0(wd, "Results/trad_fits_photo_pars.csv"),
          row.names = FALSE)



## Make a dataframe of coefficients
par_trad <- as.data.frame(coef(trad_fits_ecophys), row.names = NULL)
par_trad <- par_trad[-1,]
par_trad <- par_trad %>% 
  add_column(method = "trad")



# Merge DAT and Trad coeff dfs
par_join <- bind_rows(par_dat, par_trad)
#unique_ids <- rename(unique_ids, "unique" = "?..unique") # I was having some weird problems reading it in
#par_species <- left_join(par_join, unique_ids, by = "unique")
head(par_species)


write.csv(x = par_species, file = paste0(wd, "/Results/params_ecophys.csv"), row.names = FALSE)





# Fitting A/Ci curves with photosynthesis ----------------------------------

library(photosynthesis)

cmplt.rm_out <- read.csv(file = paste0(wd, "Inputs/Aci_no_out.csv"), header = TRUE, sep = ",")
DAT_filt <- filter(cmplt.rm_out, Data_point == "Before_DAT")
cmplt_trad <- filter(cmplt.rm_out, Data_point == "Traditional")



# Try filtering out K6709L2-2, K6712L2, and K6718L2
DAT_filt <- filter(DAT_filt, unique != "K6709L2-2" & unique != "K6714L2" & unique != "K6718L2") %>% 
  as.data.frame()



# Convert leaf temperature to Kelvin
DAT_filt$Tleaf <- DAT_filt$Tleaf + 273
cmplt_trad$Tleaf <- cmplt_trad$Tleaf + 273




# Fitting one curve at a time
k6708l1_dat_fit_photo <- fit_aci_response(data = DAT_filt[DAT_filt$unique == "K6708L1", ],
                                          varnames = list(A_net = "A", T_leaf = "Tleaf", 
                                                          C_i = "Ci", PPFD = "Qin"), 
                                          fitTPU = TRUE)
k6708l1_dat_fit_photo[[1]]
k6708l1_dat_fit_photo[[2]]



k6708l1_trad_fit_photo <- fit_aci_response(data = cmplt_trad[cmplt_trad$unique == "K6708L1", ],
                                           varnames = list(A_net = "A", T_leaf = "Tleaf", 
                                                           C_i = "Ci", PPFD = "Qin"), 
                                           fitTPU = TRUE)

k6708l1_trad_fit_photo[[1]]
k6708l1_trad_fit_photo[[2]]




# fitting many curves at a time
trad_fits_photo <- fit_many(data = cmplt_trad, 
                            varnames = list(A_net = "A", T_leaf = "Tleaf", C_i = "Ci", PPFD = "Qin"), 
                            funct = fit_aci_response,
                            group = "unique")
trad_fits_photo_graphs <- compile_data(trad_fits_photo,
                                       list_element = 2)
trad_fits_photo_pars <- compile_data(trad_fits_photo,
                                     output_type = "dataframe",
                                     list_element = 1)


pdf(file = paste0(wd,"Figures/trad_fits_photo_figs.pdf"), height=10, width=20)
plot.new()
for (curve in 1:28){
  title <- trad_fits_photo_pars$ID[[curve]]
  plot(trad_fits_photo_graphs[[curve]], main = title)
  text(30, 5, labels = as.character(title))
}
dev.off()

write.csv(x = trad_fits_photo_pars, file = paste0(wd, "Results/trad_fits_photo_pars.csv"),
          row.names = FALSE)






dat_fits_photo <- fit_many(data = DAT_filt, 
                           varnames = list(A_net = "A", T_leaf = "Tleaf", C_i = "Ci", PPFD = "Qin"), 
                           funct = fit_aci_response,
                           group = "unique")
dat_fits_photo_graphs <- compile_data(dat_fits_photo,
                                      list_element = 2)
dat_fits_photo_pars <- compile_data(dat_fits_photo,
                                    output_type = "dataframe",
                                    list_element = 1)


pdf(file = paste0(wd,"Figures/dat_fit_photo_figs.pdf"), height=10, width=20)
plot.new()
for (curve in 1:33){
  title <- dat_fits_photo_pars$ID[[curve]]
  plot(dat_fits_photo_graphs[[curve]], main = title)
  text(30, 5, labels = as.character(title))
}
dev.off()

write.csv(x = dat_fits_photo_pars, file = paste0(wd, "Results/dat_fits_photo_pars.csv"),
          row.names = FALSE)






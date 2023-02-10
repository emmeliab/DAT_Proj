## This script is to graph, fit, and save the ACi data for our Tapajos 2022 campaign

## Load Packages
library(tidyverse)
library(greekLetters)
library(ggpubr)


# Set working directory to DAT_proj
wd <- "~/Documents/GitHub/DAT_Proj/"
setwd(wd)



# Identify and Filter Outliers ---------------------------------------------------
## Load Data
complete_sp <- read.csv("Inputs/clean_aci_with_uniquecode.csv", sep = ",", header = TRUE,
                        fileEncoding="latin1") 
complete_sp <- filter(complete_sp, Data_QC == "OK") # All A/Ci curve data
unique_ids <- read.csv("Inputs/unique_ids.csv") # ID data table


## Identify Outliers and Write new data frame -------
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



# Load Filtered Data and Split into DAT and Trad  --------------------------
## Load Filtered data
cmplt.rm_out <- read.csv(file = paste0(wd, "Inputs/Aci_no_out.csv"), header = TRUE, sep = ",")

## Separate into DF of DAT and Trad
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
ggplot(cmplt.grp, mapping = aes(x = Ci, y = A, color = unique)) +
  geom_point(mapping = aes(pch = Data_point)) +
  theme_classic()

## Make and save plots for each individual species
for (code in unique(cmplt.grp$unique)) {
  df1 <- cmplt.grp %>% filter(unique == code)
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
cmplt_DAT <- filter(cmplt.rm_out, Data_point == "Before_DAT")
cmplt_trad <- filter(cmplt.rm_out, Data_point == "Traditional")


DAT_filt_ex <- cmplt_DAT %>%
  group_by(unique) %>%
  group_modify(~exclude_backwardsCi(data = .x, givedf = TRUE), .keep = FALSE)
DAT_filt_ex <- as.data.frame(DAT_filt_ex)




## Fit the ACi curves for each species for DAT using fitacis
DAT_fits_ecophys <- fitacis(DAT_filt_ex, group = "unique", fitmethod = "bilinear",
                            varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                            PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(DAT_fits_ecophys[[23]], main = coef(DAT_fits_ecophys)$unique[[23]])
coef(DAT_fits_ecophys)

### Run K6706L1 separately, since it gives a weird curve
k6714l1 <- filter(DAT_filt_ex, unique == "K6714L1")
k6714l1_fit <- fitaci(k6714l1, fitmethod = "bilinear", varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                               PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE, 
                      citransition = 200) # The Ci transition is specified as 200, as per Sharkey's
                                          # recommendations and to avoid an unreasonable Jmax value
plot(k6714l1_fit)
coef(k6714l1_fit)

# #For loop to save all the plots
# for (curve in 1:33){
#   title <- coef(DAT_fits_ecophys)$unique[[curve]]
#   png(filename = paste0(getwd(), "/Figures/", title,"_dataci_curve.png"))
#   plot(DAT_fits_ecophys[[curve]], main = title)
#   dev.off()
# }

# PDF file of all plots
pdf(file = paste0(wd,"Figures/dataci_ecophys_no_TPU.pdf"), height=10, width=20)
plot.new()
for (curve in 1:33){
  title <- coef(DAT_fits_ecophys)$unique[[curve]]
  plot(DAT_fits_ecophys[[curve]], main = title)
}
plot(k6714l1_fit, main = "K6714L1 Fixed")
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


plot(DAT_fits_ecophys[[15]], main = par_dat$unique[14])






# Fit the ACi curves for Traditional using fitacis
trad_fits_ecophys <- fitacis(cmplt_trad, group = "unique", fitmethod = "bilinear",
                             varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                             PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(trad_fits_ecophys[[14]], main = coef(trad_fits_ecophys)$unique[14]) # 20 is pretty ugly
coef(trad_fits_ecophys)

# # For loop to plot all the curves
# for (curve in 1:28){
#   title <- coef(trad_fits_ecophys)$unique[[curve]]
#   png(filename = paste0(getwd(), "/Figures/", title,"_tradaci_curve.png"))
#   plot(trad_fits_ecophys[[curve]], main = title)
#   dev.off()
# }


# Save all plots to a pdf
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
#par_trad <- par_trad[-1,]
par_trad <- par_trad %>% 
  add_column(method = "trad")





# Merge DAT and Trad coef dfs
par_join <- bind_rows(par_dat, par_trad)
head(par_join)


write.csv(x = par_join, file = paste0(wd, "/Results/params_ecophys_no_TPU.csv"), row.names = FALSE)
## Note, this does not contain the fixed K6706L1 DAT curve. This is fixed in Tapajos_stat_analysis.R




# Fitting A/Ci curves with photosynthesis ----------------------------------

library(photosynthesis)
library(tidyverse)
library(rpmodel)

cmplt.rm_out <- read.csv(file = paste0(wd, "Inputs/Aci_no_out.csv"), header = TRUE, sep = ",")

DAT_filt <- filter(cmplt.rm_out, Data_point == "Before_DAT")
cmplt_trad <- filter(cmplt.rm_out, Data_point == "Traditional")


DAT_filt_ex <- DAT_filt %>%
  group_by(unique) %>%
  group_modify(~exclude_backwardsCi(data = .x, givedf = TRUE), .keep = FALSE)
DAT_filt_ex <- as.data.frame(DAT_filt_ex)


# Convert leaf temperature to Kelvin
DAT_filt$Tleaf <- DAT_filt$Tleaf + 273
cmplt_trad$Tleaf <- cmplt_trad$Tleaf + 273
DAT_filt_ex$Tleaf <- DAT_filt_ex$Tleaf + 273
cmplt_trad <- as.data.frame(cmplt_trad)




# Fitting one curve at a time
k6708l1_dat_fit_photo <- fit_aci_response(data = DAT_filt[DAT_filt$unique == "K6708L1", ],
                                          varnames = list(A_net = "A", T_leaf = "Tleaf", 
                                                          C_i = "Ci", PPFD = "Qin"), 
                                          fitTPU = TRUE)
k6718l2_dat_fit_photo[[1]]
k6718l2_dat_fit_photo[[2]]



k6708l1_trad_fit_photo <- fit_aci_response(data = cmplt_trad[cmplt_trad$unique == "K6708L1", ],
                                           varnames = list(A_net = "A", T_leaf = "Tleaf", 
                                                           C_i = "Ci", PPFD = "Qin"), 
                                           fitTPU = TRUE)

k6708l1_trad_fit_photo[[1]]




# fitting many curves at a time
trad_fits_photo <- fit_many(data = cmplt_trad, fitTPU = FALSE,
                            varnames = list(A_net = "A", T_leaf = "Tleaf", C_i = "Ci", PPFD = "Qin"), 
                            funct = fit_aci_response,
                            group = "unique")
trad_fits_photo_graphs <- compile_data(trad_fits_photo,
                                       list_element = 2)
trad_fits_photo_pars <- compile_data(trad_fits_photo,
                                     output_type = "dataframe",
                                     list_element = 1)


dat_fits_photo <- fit_many(data = DAT_filt_ex, fitTPU = FALSE,## uses the back-filtered data
                           varnames = list(A_net = "A", T_leaf = "Tleaf", C_i = "Ci", PPFD = "Qin"), 
                           funct = fit_aci_response,
                           group = "unique")
dat_fits_photo_graphs <- compile_data(dat_fits_photo,
                                      list_element = 2)
dat_fits_photo_pars <- compile_data(dat_fits_photo,
                                    output_type = "dataframe",
                                    list_element = 1)

#Temperature correction. Important: watch celsius vs fahrenheit

grp_curv <- cmplt.rm_out %>% 
  group_by(Data_point, unique) %>% 
  summarize(meanTleaf = mean(Tleaf))
grp_curv2 <- grp_curv %>% 
  rename(ID = unique)
grp_dat <- filter(grp_curv2, Data_point == "Before_DAT")
grp_trad <- filter(grp_curv2, Data_point == "Traditional")

dat_fits_photo_pars$method <- "DAT"
trad_fits_photo_pars$method <- "Traditional"

curv_dat_temp <- left_join(by = "ID", dat_fits_photo_pars, grp_dat)
curv_trad_temp <- left_join(by = "ID", trad_fits_photo_pars, grp_trad)

f_fact_dat <- data.frame()
for (i in 1:nrow(curv_dat_temp)){
  tcleaf <- curv_dat_temp$meanTleaf[i]
  tcref <- 25
  f <- ftemp_inst_vcmax(tcleaf,tcref) #this function is in the rpmodel package
  f_fact_dat <- rbind(f_fact_dat, f)
}
colnames(f_fact_dat) <- "f_fact"
curv_dat_temp_adj <- curv_dat_temp %>% 
  mutate(vcmax_25 = V_cmax * f_fact_dat$f_fact)

f_fact_trad <- data.frame()
for (i in 1:nrow(curv_trad_temp)){
  tcleaf <- curv_trad_temp$meanTleaf[i]
  tcref <- 25
  f <- ftemp_inst_vcmax(tcleaf,tcref) #this function is in the rpmodel package
  f_fact_trad <- rbind(f_fact_trad, f)
}
colnames(f_fact_trad) <- "f_fact"
curv_trad_temp_adj <- curv_trad_temp %>% 
  mutate(vcmax_25 = V_cmax * f_fact_trad$f_fact)


#Write PDF of outputs

pdf(file = paste0(wd,"Figures/trad_fits_photo_figs_no_TPU.pdf"), height=10, width=20)
plot.new()
for (curve in 1:28){
  title <- trad_fits_photo_pars$ID[[curve]]
  plot(trad_fits_photo_graphs[[curve]], main = title)
  text(30, 5, labels = as.character(title))
}
dev.off()

pdf(file = paste0(wd,"Figures/dat_fit_photo_figs_filt_no_TPU.pdf"), height=10, width=20)
plot.new()
for (curve in 1:33){ ### Change this depending on the number of curves!
  title <- dat_fits_photo_pars$ID[[curve]]
  plot(dat_fits_photo_graphs[[curve]], main = title)
  text(30, 5, labels = as.character(title))
}
dev.off()

#Write csvs

write.csv(x = curv_trad_temp_adj, file = paste0(wd, "Results/trad_fits_photo_pars_correct_no_TPU.csv"),
          row.names = FALSE)

write.csv(x = curv_dat_temp_adj, file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_no_TPU.csv"),
          row.names = FALSE)






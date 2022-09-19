
# Purpose -----------------------------------------------------------------
## This script is to graph, fit, and save the ACi data for our Tapajos 2022 campaign

####A bit more processing and some exploring


##### Loren, please take a look at the fitacis function, and the traditional fitaci functions


# Load Packages/Data ------------------------------------------------------
## Load Packages
library(tidyverse)
library(plantecophys)


# Load Data
## Make sure working directory is DAT_proj
getwd()
#setwd()
complete_sp <- read.csv("Inputs/clean_aci_data_one_file.csv", 
                        sep = ",", 
                        fileEncoding ="latin1")


# Identify Outliers -------------------------------------------------------

which(complete_sp$Ci < -50) #4 Ci values are super negative: 394, 395, 396, 398
complete_sp[c(394,395,396,398),19] # What are those values? Should they be outliers?
which(complete_sp$A > 31) # 2 A values greater than 31: 399 and 400
complete_sp[c(397,398,399,400,401),17]
cmplt.rm_out <- filter(complete_sp, Ci > -50 & A < 31)



# Plotting ACi Curves -----------------------------------------------------

## Group by Four Letter Code
cmplt.grp <- group_by(cmplt.rm_out, fourlettercode) %>% 
  group_by(Leaf_number)



## Plot all ACi curves on one graph by Species
ggplot(cmplt.grp, mapping = aes(x = Ci, y = A, color = fourlettercode)) +
  geom_point(mapping = aes(pch = Data_point)) +
  theme_classic()




## Make and save plots for each individual species
for (code in unique(cmplt.grp$fourlettercode)) {
  df1 <- cmplt.grp %>% filter(fourlettercode == code)
  gg1 <- ggplot(data = df1, mapping = aes(x = Ci, y = A, color = Data_point)) +
    geom_point() +
    theme_classic() +
    scale_color_viridis_d() +
    ggtitle(code)
  plot(gg1)
  filename1 <- paste("plot_", code, ".png")
  ggsave(filename1, gg1, path = paste0(getwd(), "/Outputs/"))
}


## The same, but for each individual tree
for (id in unique(cmplt.grp$k67.id)) {
  df1 <- cmplt.grp %>% filter(k67.id == id)
  gg1 <- ggplot(data = df1, mapping = aes(x = Ci, y = A, color = Data_point)) +
    geom_point() +
    theme_classic() +
    scale_color_viridis_d() +
    ggtitle(id)
  plot(gg1)
  filename1 <- paste("plot_", id, ".png")
  ggsave(filename1, gg1, path = paste0(getwd(), "/Outputs/"))
}





# Fitting ACi Curves ------------------------------------------------------
library(greekLetters)

## Separate by DAT and trad, and convert to dataframe
cmplt_DAT <- filter(cmplt.grp, Data_point == "Before_DAT") %>% 
  select(-contains(greeks("Delta"))) #removes the columns with deltas
cmplt_DAT <- as.data.frame(cmplt_DAT)
head(cmplt_DAT)

cmplt_trad <- filter(cmplt.grp, Data_point == "Traditional") %>% 
  select(-contains(greeks("Delta")))
cmplt_trad <- as.data.frame(cmplt_trad)
head(cmplt_trad)



# # Fit the ACi curves for each species for DAT and Traditional using fitacis
# DAT_fits <- fitacis(cmplt_DAT, group = "k67.id",
#                     varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
#                                     PPFD = "Qin"), fitTPU = TRUE, Tcorrect = TRUE,
#                     id = Leaf_number)
# plot(DAT_fits[[1]])
# 
# 
# trad_fits <- fitacis(cmplt_trad, group = "k67.id",
#                     varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
#                                     PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE,
#                     id = Leaf_number)
#plot(trad_fits[[1]])




## For loop to fit aci curves for each leaf individually using fitaci

# DAT
for (code in length(unique(cmplt_DAT$k67.id))) {
  ind <- cmplt_DAT$k67.id == cmplt_DAT$k67.id[code]
  dat <- cmplt_DAT[ind, ]
  #print(head(dat))
  for (lf in length(unique(dat$Leaf_number))){
    indlf <- dat$Leaf_number == dat$Leaf_number[lf]
    datlf <- dat[indlf, ]
    fit <- fitaci(data = datlf, varnames = list(ALEAF = "Adyn", Tleaf = "Tleaf", 
                                                Ci = "Ci", PPFD = "Qin"),
           fitTPU = FALSE, Tcorrect = TRUE)
    assign(paste0("DAT_fit_", cmplt_DAT$k67.id[code], dat$Leaf_number[lf]), fit)
  }
  plot(fit[[lf]])
  }

# Traditional
remove(trad_pars) #remove the data frame if you've ran it before
trad_pars <- data.frame()
for (code in 1:length(unique(cmplt_trad$k67.id))) {
  # First, separate for each tree
  ind <- cmplt_trad$k67.id == cmplt_trad$k67.id[code]
  trad <- cmplt_trad[ind, ]
  print(trad$k67.id)
  for (lf in 1:length(unique(trad$Leaf_number))) {
    # Next, separate for each leaf
    indlf <- trad$Leaf_number == trad$Leaf_number[lf]
    tradlf <- trad[indlf, ]
    print(tradlf$Leaf_number)
    fit <- fitaci(data = tradlf, varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
                                                 Ci = "Ci", PPFD = "Qin"),
                  fitTPU = FALSE, Tcorrect = TRUE)
    print(fit)
    plot(fit, main = paste(cmplt_trad$k67.id[code], "Leaf", trad$Leaf_number[lf]))
    assign(paste0("trad_fit_", cmplt_trad$k67.id[code], "_lf", trad$Leaf_number[lf]), fit)
    summary(fit)
    #dev.print(png, file =  paste0("trad_fit_", cmplt_trad$k67.id[code], "_lf",
    #                       trad$Leaf_number[lf], ".png"))
  }
  # Third, pull the parameters for each and append to a single dataframe
  par <- fit$pars
  par.df <- as.data.frame(par, row.names = NULL)
  par.df2 <- mutate(par.df, Result = c("Vcmax", "Jmax", "Rd"))
  par.df3 <- mutate(par.df2, Curve = as.character(paste0(cmplt_trad$k67.id[code], "_lf", 
                                          trad$Leaf_number[lf])))
  trad_pars <- bind_rows(trad_pars, par.df3)
}
print(trad_pars)
trad_fit_list <- ls(pattern = "trad_fit_K67")










# Examine Fit ACi Results -------------------------------------------------
library(ggpubr)

## DAT




## Traditional
trad_pars_grp <- group_by(trad_pars, Curve)

ggboxplot(data = trad_pars_grp, x = "Result", y = "Estimate", bxp.errorbar = TRUE,
          bxp.errorbar.width = "Std. Error", facet.by = "Curve")


overview(`trad_fit_K67-WT-08_lf1`$nlsfit)
overview(`trad_fit_K67-WT-09_lf2`$nlsfit)











#### Below is Charlie's code to run Acis on each leaf----------------------------


# ### plot Vcmax for one leaf from a given species, in this case ABMA.
# ### having a hard time getting the plural "fitacis" function to work, so iteratively used "fitaci"
# ### Leaf 1 DAT
# abma <- filter(cmplt.grp, fourlettercode == "ABMA")
# abma_dat <- filter(abma, Data_point == "Before_DAT")
# abma_dat_l1 <- filter(abma_dat, Leaf_number == 1)
# abma_dat_aci_l1 <- fitaci(data = abma_dat_l1, 
#                          varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", 
#                                          PPFD = "Qin"),
#                          fitTPU = TRUE,
#                          Tcorrect = TRUE)
# plot(abma_dat_aci_l1)
abma_dat_aci_l1$pars # the pars gives the parameter estimates as well as the standard error




### Leaf 2 DAT
abma_dat_l2 <- filter(abma_dat, Leaf_number == 2)
abma_dat_aci_l2 <- fitaci(data = abma_dat_l2, 
                          varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                          fitTPU = TRUE,
                          Tcorrect = TRUE)
plot(abma_dat_aci_l2)
abma_dat_aci_l2$pars

### Leaf 1 Trad
abma <- filter(complete_sp, fourlettercode == "ABMA")
abma_trad <- filter(abma, Data_point == "Traditional")
abma_trad_l1 <- filter(abma_trad, Leaf_number == 1)
abma_trad_aci_l1 <- fitaci(data = abma_trad_l1, 
                          varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                          fitTPU = TRUE,
                          Tcorrect = TRUE)
plot(abma_trad_aci_l1)
abma_trad_aci_l1$pars

### Leaf 2 Trad
abma_trad_l2 <- filter(abma_trad, Leaf_number == 2)
abma_trad_aci_l2 <- fitaci(data = abma_trad_l2, 
                          varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                          fitTPU = TRUE,
                          Tcorrect = TRUE)
plot(abma_trad_aci_l2)
abma_trad_aci_l2$pars

### trying to get multiple summary data tables into one data frame for further comparison


### the plantecophys creators recommend working with the nlstools packaging to analyze data
### but I haven't figured it out yet
install.packages('nlstools')
library(nlstools)

nls_abma_tracl2 <- abma_trad_aci_l2$nlsfit
overview(nls_abma_tracl2)

### plot Vcmax for a given species from multiple leaves ### this hasn't worked yet for some reason
## abma <- filter(complete_species, fourlettercode == "ABMA")
##abma_dat <- filter(abma, Data_point == "Before_DAT")
## abma_dat_multi <- fitacis(data = abma_dat, 
                       # group = "Leaf_number",
                       # varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                       # fitTPU = TRUE,
                       # Tcorrect = TRUE)


# Purpose -----------------------------------------------------------------
## This script is to graph, fit, and save the ACi data for our Tapajos 2022 campaign

####A bit more processing and some exploring


##### Loren, please take a look at the fitacis function, and the traditional fitaci functions


# Load Packages/Data ------------------------------------------------------
## Load Packages
library(tidyverse)
library(plantecophys)


# Load Data
## I'm gonna work on a way to make this universal
complete_sp <- read.csv("~/Documents/PhD/DAT_Tapajos/Inputs/clean_aci_data_one_file.csv")




# Identify Outliers -------------------------------------------------------

which(complete_sp$Ci < -50) #4 Ci values are super negative: 394, 395, 396, 398
which(complete_sp$A > 31) # 2 A values greater than 31: 395 and 396
cmplt.rm_out <- filter(complete_sp, Ci > -50 & A < 31)



# Plotting ACi Curves -----------------------------------------------------

## Group by Four Letter Code
cmplt.grp <- group_by(cmplt.rm_out, fourlettercode)



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
    scale_color_viridis_d()
  filename1 <- paste("plot_", code, ".png")
  ggsave(filename1, gg1)
}






# Fitting ACi Curves ------------------------------------------------------

# Fit the ACi curves for each species for DAT and Traditional, respectively
######## This is the part that is broken :(
cmplt_DAT <- filter(cmplt.grp, Data_point == "Before_DAT")
DAT_fits <- fitacis(cmplt_DAT, group = "fourlettercode", 
                    varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                    PPFD = "Qin"), fitTPU = TRUE, Tcorrect = TRUE,
                    id = Leaf_number)
plot(DAT_fits[[1]])


### Attempts to make a for loop to fit curves, ignore
#for (code in unique(cmplt.grp$fourlettercode)) {
  each <- cmplt.grp %>% filter(fourlettercode == code)
  if (each$Data_point == "Before_DAT") {
    fitacis(each, group = Leaf_number, varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                                       PPFD = "Qin"), id = fourlettercode,
            fitTPU = TRUE, Tcorrect = TRUE, how = "oneplot")
  }
}

#for (code in unique(cmplt.grp$fourlettercode)) {
  each <- cmplt.grp %>% filter(Data_point == "Before_DAT")
  fitacis(each, group = "fourlettercode", varnames = list(ALEAF = "A", Tleaf = "Tleaf",
                                                          Ci = "Ci",
                                                          PPFD = "Qin"), id = code, fitTPU = TRUE,
          Tcorrect = TRUE)
}



#### Below is Charlie's code to run Acis on each leaf----------------------------
#### This part was also giving me trouble, possibly because of the columns with deltas

### plot Vcmax for one leaf from a given species, in this case ABMA.
### having a hard time getting the plural "fitacis" function to work, so iteratively used "fitaci"
### Leaf 1 DAT
abma <- filter(cmplt.grp, fourlettercode == "ABMA")
abma_dat <- filter(abma, Data_point == "Before_DAT")
abma_dat_l1 <- filter(abma_dat, Leaf_number == 1)
abma_dat_aci_l1 <- fitaci(data = abma_dat_l1, 
                         varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", 
                                         PPFD = "Qin"),
                         fitTPU = TRUE,
                         Tcorrect = TRUE)
plot(abma_dat_aci_l1)
abma_dat_aci_l1$pars # the pars gives the parameter estimates as well as the standard error




### Attempts to remove columns with deltas, ignore
#rmv_delta <- function(file_name) {
  file_name <- file_name[,!grepl(expression(delta), colnames(file_name))]
}
#cmplt_DAT_no_delta <- rmv_delta(cmplt_DAT)




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

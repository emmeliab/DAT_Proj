####A bit more processing and some exploring

library(tidyverse)
library(plantecophys)

### adding a column for a four-letter species code and a column for species name
### need to add to this once we get all of the data into Emmy's original csv file
complete <- read.csv("~/Documents/PhD/DAT_Tapajos/Inputs/clean_aci_data_one_file.csv")
complete <- as_tibble(complete)
complete_species <- complete %>% mutate(fourlettercode = Tree_Identifier, SciName = Tree_Identifier)

complete_species$fourlettercode <- recode(complete_species$fourlettercode,
                                          'Maca1' = 'MAEL',
                                          'Tree3' = 'CHTU',
                                          'tree8' = 'COSP',
                                          'tree9' = 'APCO',
                                          'Tree10' = 'VICA',
                                          'tree11' = 'COST',
                                          'tree12' = 'UNKN',
                                          'tree22' = 'ABMA')
complete_species$SciName <- recode(complete_species$SciName,
                                   'Maca1' = 'Manilkara elata',
                                   'Tree3' = 'Chimaris turbinata',
                                   'tree8' = 'Coussarea sp',
                                   'tree9' = 'Aparisthmium cordatum',
                                   'Tree10' = 'Vismia cayennensis',
                                   'tree11' = 'Couratari stellata',
                                   'tree12' = 'Unknown sp',
                                   'tree22' = 'Abarema mataybifolia')

### outliers?
Ci_outlier_index <- which(complete_species$Ci < -50) #  4 values are super negative, 394, 395, 396, 398
A_outlier_index <- which(complete_species$A > 31) #  2 values greater than 31, 395 and 396
complete_species <- filter(complete_species, Ci > -50 & A < 31)

### Group data by species four letter code to see them all in one place
p1 <- ggplot(complete_species, aes(x = Ci, y = A, colour = fourlettercode))
p1 + geom_point() + theme_classic()

### plot a set of curves for a given species (2 DATs and 2 Trads)
abma <- filter(complete_species, fourlettercode == "ABMA")
p2 <- ggplot(abma, aes(x = Ci, y = A, colour = Data_point))
p2 + geom_point() + theme_classic()


### plot Vcmax for one leaf from a given species, in this case ABMA.
### having a hard time getting the plural "fitacis" function to work, so iteratively used "fitaci"
### Leaf 1 DAT
abma <- filter(complete_species, fourlettercode == "ABMA")
abma_dat <- filter(abma, Data_point == "Before_DAT")
abma_dat_l1 <- filter(abma_dat, Leaf_number == 1)
abma_dat_aci_l1 <- fitaci(data = abma_dat_l1, 
                         varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin"),
                         fitTPU = TRUE,
                         Tcorrect = TRUE)
plot(abma_dat_aci_l1)
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
abma <- filter(complete_species, fourlettercode == "ABMA")
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

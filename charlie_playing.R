# Purpose -----------------------------------------------------------------
## This script is to graph, fit, and save the ACi data for our Tapajos 2022 campaign

####A bit more processing and some exploring


##### Loren, please take a look at the fitacis function, and the traditional fitaci functions


# Load Packages/Data ------------------------------------------------------
## Load Packages
library(tidyverse)
library(plantecophys)
library(greekLetters)
library(ggpubr)

# Load Data
## Make sure working directory is DAT_proj
#getwd()
#setwd()

complete_sp <- read.csv("/Users/charlessouthwick/Documents/GitHub/DAT_Proj/Inputs/clean_aci_with_uniquecode.csv",
                        sep = ",", 
                        fileEncoding="latin1")
#Create an id data table. This includes the species code so we can merge it later
unique_ids <- read.csv("/Users/charlessouthwick/Documents/GitHub/DAT_Proj/Inputs/unique_ids.csv")

# Identify Outliers and filtering ---------------------------------------------------

which(complete_sp$Ci < -50) #4 Ci values are super negative: 394, 395, 396, 398
which(complete_sp$Ci < 0) # Add in 16 more values, all for MACA1
complete_sp[c(394,395,396,398),19] # What are those values? Should they be outliers?
which(complete_sp$A > 40) # 2 A values greater than 31: 399 and 400
which(complete_sp$A < 1)
complete_sp[c(397,398,399,400,401),17]
cmplt.rm_out <- filter(complete_sp, Ci > 0 & A < 31)  ## Ci > 0 (we had been doing Ci > -50)

#filter again for A above 1 to try and fix the DAT curves that look weird
which(complete_sp$A < 1)
cmplt.rm_out <- filter(complete_sp, A > 1)

# Plotting ACi Curves -----------------------------------------------------

## Group by unique instead, Emmy did by leaf number.
cmplt.grp <- group_by(cmplt.rm_out, fourlettercode) %>% 
  group_by(unique)

## Plot all ACi curves on one graph by Species
#ggplot(cmplt.grp, mapping = aes(x = Ci, y = A, color = fourlettercode)) +
 # geom_point(mapping = aes(pch = Data_point)) +
 # theme_classic()

## Make and save plots for each individual species
#for (code in unique(cmplt.grp$fourlettercode)) {
  #df1 <- cmplt.grp %>% filter(fourlettercode == code)
 # gg1 <- ggplot(data = df1, mapping = aes(x = Ci, y = A, color = Data_point)) +
  #  geom_point() +
  #  theme_classic() +
  #  scale_color_viridis_d() +
  #  ggtitle(code)
 # plot(gg1)
 # filename1 <- paste("plot_", code, ".png")
 # ggsave(filename1, gg1, path = paste0(getwd(), "/Outputs/"))
#}

## The same, but for each individual tree
#for (id in unique(cmplt.grp$k67.id)) {
 # df1 <- cmplt.grp %>% filter(k67.id == id)
  #gg1 <- ggplot(data = df1, mapping = aes(x = Ci, y = A, color = Data_point)) +
   # geom_point() +
    #theme_classic() +
    #scale_color_viridis_d() +
    #ggtitle(id)
  #plot(gg1)
  #filename1 <- paste("plot_", id, ".png")
  #ggsave(filename1, gg1, path = paste0(getwd(), "/Outputs/"))
#}


# Fitting ACi Curves ------------------------------------------------------
# install.packages("greekLetters") # have to reinstall if you clean the env to get 
# the vectors you need
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
DAT_fits <- fitacis(cmplt_DAT, group = "unique",
                    varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                    PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(DAT_fits[[7]]) # I was just exploring here
coef(DAT_fits)
#strange curves: #7, 9, 10, 11, 12, 15, 17, 22, 23, 24?, 32? --> at least 11 are a weird fit!
# I wonder, should we try Loren's old approach from her grad school days?

## make a dataframe out of coefficients
par_dat <- as.data.frame(coef(DAT_fits), row.names = NULL)
par_dat <- par_dat[-1,]
par_dat <- par_dat %>%
  add_column(method = "dat")
table(par_dat$method) ## number of initial DAT curves

### note that some of these estimates are way incorrect. I think it is a problem of ecophys
### which doesn't seem to  be able to figure out the initial 'back correction'

## now let's run the traditional
trad_fits <- fitacis(cmplt_trad, group = "unique",
                     varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                     PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(trad_fits[[20]]) # 20 is pretty ugly
coef(trad_fits)

par_trad <- as.data.frame(coef(trad_fits), row.names = NULL)
par_trad <- par_trad[-1,]
par_trad <- par_trad %>% 
  add_column(method = "trad")

#merge trad and dat outputs into one
par_join <- bind_rows(par_dat, par_trad)
par_species <- left_join(par_join, unique_ids, by = "unique")
par_species

#filter some outliers (at least until we can figure out our plantecophys issues) 
which(par_species$Vcmax > 100)
which(par_species$Vcmax < 0)
filt_par_species <- par_species %>% 
  filter(Vcmax < 100 & Vcmax > 0)
filt_par_species

## some exploration
table(filt_par_species$method) ## number of each type of curve

group_by(filt_par_species, method) %>%
  summarise(
    count = n(),
    mean = mean(Vcmax, na.rm = TRUE),
    sd = sd(Vcmax, na.rm = TRUE)
  )

#qqplot for visual estimation of normality
ggqqplot(filt_par_species$Vcmax)
ggqqplot(filt_par_species$Jmax)

#shapiro-wilk test for normality
shapiro.test(filt_par_species$Vcmax)
shapiro.test(filt_par_species$Jmax)

#these results indicate a super low p-value, which means the data ARE significantly
# different from normal.
#since many statistical tests assume a normal distribution, we have to think about this.
#But the Central Limit Theorem says that if a sample size n > 30 we can ignore normality
#But how does this apply to paired samples? Do we need 30 for both trad and dat?
# Or just 30 total?

#I think a wilcox signed-rank test uses non-parametric, paired data.
#Unfortunately we still need the same number of before and afters...
wilcox.test(filt_par_species$Vcmax ~ filt_par_species$method, paired = TRUE)

#here's a histogram to prove it
hist_pars <- filt_par_species %>% 
  ggplot() +
  geom_histogram(aes(x = Vcmax))
hist_pars

#t-test doesn't work as paired data because we don't have the same number of curves for each method
# ttest_vcmax <- t.test(Vcmax ~ method, data = filt_par_species, paired = TRUE)
#ttest_vcmax


box_both_vcmax <- filt_par_species %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = Vcmax))
box_both_vcmax

box_both_jmax <- filt_par_species %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = Jmax))
box_both_jmax

scat_vcmax <- filt_par_species %>%
  ggplot() +
  geom_point(aes(x = method, y = Vcmax, color = unique))
scat_vcmax

box_both_vcmaxse <- par_species %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = Vcmax_SE))
box_both_vcmaxse

box_both_jmax <- par_species %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = Jmax))
box_both_jmax

box_both_jmaxse <- par_species %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = Jmax_SE))
box_both_jmaxse


## Emmy's For Loop below
## For loop to fit aci curves for each leaf individually using fitaci
# 
# # DAT
# for (code in length(unique(cmplt_DAT$k67.id))) {
#   ind <- cmplt_DAT$k67.id == cmplt_DAT$k67.id[code]
#   dat <- cmplt_DAT[ind, ]
#   #print(head(dat))
#   for (lf in length(unique(dat$Leaf_number))){
#     indlf <- dat$Leaf_number == dat$Leaf_number[lf]
#     datlf <- dat[indlf, ]
#     fit <- fitaci(data = datlf, varnames = list(ALEAF = "Adyn", Tleaf = "Tleaf", 
#                                                 Ci = "Ci", PPFD = "Qin"),
#                   fitTPU = FALSE, Tcorrect = TRUE)
#     assign(paste0("DAT_fit_", cmplt_DAT$k67.id[code], dat$Leaf_number[lf]), fit)
#   }
#   plot(fit[[lf]])
# }

# # Traditional
# remove(trad_pars) #remove the data frame if you've run it before
# trad_pars <- data.frame()
# for (code in 1:length(unique(cmplt_trad$k67.id))) {
#   # First, separate for each tree
#   ind <- cmplt_trad$k67.id == cmplt_trad$k67.id[code]
#   trad <- cmplt_trad[ind, ]
#   print(trad$k67.id)
#   for (lf in 1:length(unique(trad$Leaf_number))) {
#     # Next, separate for each leaf
#     indlf <- trad$Leaf_number == trad$Leaf_number[lf]
#     tradlf <- trad[indlf, ]
#     print(tradlf$Leaf_number)
#     fit <- fitaci(data = tradlf, varnames = list(ALEAF = "Asty", Tleaf = "Tleaf", 
#                                                  Ci = "Ci", PPFD = "Qin"),
#                   fitTPU = FALSE, Tcorrect = TRUE)
#     print(fit)
#     plot(fit, main = paste(cmplt_trad$k67.id[code], "Leaf", trad$Leaf_number[lf]))
#     assign(paste0("trad_fit_", cmplt_trad$k67.id[code], "_lf", trad$Leaf_number[lf]), fit)
#     summary(fit)
#     #dev.print(png, file =  paste0("trad_fit_", cmplt_trad$k67.id[code], "_lf",
#     #                       trad$Leaf_number[lf], ".png"))
#   }
#   # Third, pull the parameters for each and append to a single dataframe
#   par <- fit$pars
#   par.df <- as.data.frame(par, row.names = NULL)
#   par.df2 <- mutate(par.df, Result = c("Vcmax", "Jmax", "Rd"))
#   par.df3 <- mutate(par.df2, Curve = as.character(paste0(cmplt_trad$k67.id[code], "_lf", 
#                                                          trad$Leaf_number[lf])))
#   trad_pars <- bind_rows(trad_pars, par.df3)
# }
# # print(trad_pars)
# # trad_fit_list <- ls(pattern = "trad_fit_K67")
# # print(trad_fit_list)
# 
# 
# # Examine Fit ACi Results -------------------------------------------------
# library(ggpubr)
# 
# ## DAT
# 
# 
# ## Traditional
# trad_pars_grp <- group_by(trad_pars, Curve)
# 
# ggboxplot(data = trad_pars_grp, x = "Result", y = "Estimate", bxp.errorbar = TRUE,
#           bxp.errorbar.width = "Std. Error", facet.by = "Curve")
# 
# 
# overview(`trad_fit_K67-WT-08_lf1`$nlsfit)
# overview(`trad_fit_K67-WT-09_lf2`$nlsfit)
# 
# 


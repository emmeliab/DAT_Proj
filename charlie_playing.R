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

# Set working directory to DAT_proj
#getwd()
#setwd()

# Load Data
complete_sp <- read.csv("Inputs/clean_aci_with_uniquecode.csv", sep = ",", 
                        fileEncoding="latin1") #changed this to work on anyone's computer
#Create an id data table. This includes the species code so we can merge it later
unique_ids <- read.csv("Inputs/unique_ids.csv") # same here




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

## Group by unique
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

# The same, but for each unique tree ## Charlie changed this on 11/1/22
# for (id in unique(cmplt.grp$unique)) {
# df1 <- cmplt.grp %>% filter(unique == id)
# gg1 <- ggplot(data = df1, mapping = aes(x = Ci, y = A, color = Data_point)) +
# geom_point() +
# theme_classic() +
# scale_color_viridis_d() +
# ggtitle(id)
# plot(gg1)
# filename1 <- paste("plot_", id, ".png")
# ggsave(filename1, gg1, path = paste0(getwd(), "/Outputs/"))
# }



# Fitting ACi Curves ------------------------------------------------------

## Separate by DAT and trad, and convert to dataframe
cmplt_DAT <- filter(cmplt.grp, Data_point == "Before_DAT") %>% 
  select(-contains(greeks("Delta"))) #removes the columns with deltas
cmplt_DAT <- as.data.frame(cmplt_DAT)
head(cmplt_DAT)

cmplt_trad <- filter(cmplt.grp, Data_point == "Traditional") %>% 
  select(-contains(greeks("Delta")))
cmplt_trad <- as.data.frame(cmplt_trad)
head(cmplt_trad)

# Fit the ACi curves for each species for DAT using fitacis
DAT_fits <- fitacis(cmplt_DAT, group = "unique", id = "unique",
                    varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                    PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(DAT_fits[[15]], main = coef(DAT_fits)$unique[15])
coef(DAT_fits)


#jpeg("K6707L2-2.jpg")
#plot(trad_fits[[10]], main = coef(trad_fits)$unique[10])
#dev.off()
#strange curves: #7, 9, 10, 11, 12, 15, 17, 22, 23, 24? --> at least 11 are a weird fit!
# I wonder, should we try Loren's old approach from her grad school days?

## Make a dataframe out of coefficients
par_dat <- as.data.frame(coef(DAT_fits), row.names = NULL)
par_dat <- par_dat[-1,]
par_dat <- par_dat %>%
  add_column(method = "dat")
table(par_dat$method) ## number of initial DAT curves

### note that some of these estimates are way incorrect. I think it is a problem of ecophys
### which doesn't seem to  be able to figure out the initial 'back correction'



# Fit the ACi curves for Traditional using fitacis
trad_fits <- fitacis(cmplt_trad, group = "unique", #id = "unique",
                     varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                     PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(trad_fits[[10]], main = coef(trad_fits)$unique[10]) # 20 is pretty ugly
coef(trad_fits)

## Make a dataframe of coefficients
par_trad <- as.data.frame(coef(trad_fits), row.names = NULL)
par_trad <- par_trad[-1,]
par_trad <- par_trad %>% 
  add_column(method = "trad")

# Merge DAT and Trad dfs
par_join <- bind_rows(par_dat, par_trad)
unique_ids <- rename(unique_ids, "unique" = "?..unique") # I was having some weird problems reading it in
par_species <- left_join(par_join, unique_ids, by = "unique")
head(par_species)


# Filter some outliers (at least until we can figure out our plantecophys issues) 
which(par_species$Vcmax > 100)
which(par_species$Vcmax < 0)
filt_par_species <- par_species %>% 
  filter(Vcmax < 100 & Vcmax > 0)
head(filt_par_species)



## some exploration
table(filt_par_species$method) ## number of each type of curve

group_by(filt_par_species, method) %>%
  summarise(
    count = n(),
    mean = mean(Vcmax, na.rm = TRUE),
    sd = sd(Vcmax, na.rm = TRUE))




# Testing Normality -------------------------------------------------------

# Q-Q plots
ggqqplot(filt_par_species$Vcmax)
ggqqplot(filt_par_species$Jmax)

# Shapiro-Wilk test
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


# This fun little function plots a curve on top of the graph, just for visualization
library(rcompanion)
plotNormalHistogram(filt_par_species$Vcmax, breaks = 30)
plotNormalHistogram(filt_par_species$Jmax, breaks = 30)

#t-test doesn't work as paired data because we don't have the same number of curves for each method
# ttest_vcmax <- t.test(Vcmax ~ method, data = filt_par_species, paired = TRUE)
#ttest_vcmax


box_both_vcmax <- filt_par_species %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = Vcmax)) +
  theme_light()
box_both_vcmax

box_both_jmax <- filt_par_species %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = Jmax)) +
  theme_light()
box_both_jmax


filt_par_dummy <- mutate(.data = filt_par_species,# makes a dummy variable to plot
               dummy = if_else(filt_par_species$method == "dat", 0,1))

scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = Vcmax,
                                                          color = unique)) + 
  geom_line() + 
  geom_point() +
  theme_light() +
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Trad")) +
  xlab("Method")
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

scat_jmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = Jmax,
                                                          color = unique)) + 
  geom_line() + 
  geom_point() +
  theme_light() +
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Trad")) +
  xlab("Method")
scat_jmax




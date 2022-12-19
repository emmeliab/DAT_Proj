# Purpose -----------------------------------------------------------------
## This script is to graph, fit, and save the ACi data for our Tapajos 2022 campaign


# Load Packages/Data ------------------------------------------------------
## Load Packages
library(tidyverse)
library(plantecophys)
library(greekLetters)
library(ggpubr)


# Set working directory to DAT_proj
wd <- "C://Users/emmel/Desktop/DAT_proj/"
setwd(wd)

# Load Data
complete_sp <- read.csv("Stomata Test/clean_data_sr.csv", sep = ",", header = TRUE,
                        fileEncoding="latin1")
complete_sp <- filter(complete_sp, Data_QC == "OK")
#Create an id data table. This includes the species code so we can merge it later
#unique_ids <- read.csv("Inputs/unique_ids.csv") # same here

complete_sp <- mutate(complete_sp, unique = paste0("K67", substring(complete_sp$k67.id, 8, 9), 
                                                   "L", round(complete_sp$Leaf_number, 1)))
unique(complete_sp$unique)
#test <- read.csv("Inputs/clean_aci_with_uniquecode.csv", sep = ",", fileEncoding = "latin1")


# Identify Outliers and filtering ---------------------------------------------------

which(complete_sp$Ci < -50) #4 Ci values are super negative: 394, 395, 396, 398
which(complete_sp$Ci < -5) # Add in 16 more values, all for MACA1
complete_sp[c(394,395,396,398),19] # What are those values? Should they be outliers?
which(complete_sp$A > 40)
which(complete_sp$A < -1)
cmplt.rm_out1 <- filter(complete_sp, Ci > -5)
cmplt.rm_out2 <- filter(cmplt.rm_out1, A < 40) ## A < 40
cmplt.rm_out <- filter(cmplt.rm_out2, A > -1)
# Data frame without outliers
write.csv(x = cmplt.rm_out, file = paste0(getwd(), "/Stomata Test/Aci_no_out_sr.csv"), 
          row.names = FALSE)


## Group by unique
cmplt.grp <- group_by(cmplt.rm_out, fourlettercode) %>% 
  group_by(unique)

## Separate by DAT and trad, and convert to dataframe
cmplt_DAT <- filter(cmplt.grp, Data_point == "Before_DAT") %>% 
  select(-contains(greeks("Delta"))) #removes the columns with deltas
head(cmplt_DAT)
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


cmplt_trad <- filter(cmplt.grp, Data_point == "Traditional") %>% 
  select(-contains(greeks("Delta")))
cmplt_trad <- as.data.frame(cmplt_trad)
head(cmplt_trad)


# Make function to find min Ci and exclude points with Anet below the Anet for that min Ci
exclude_backwardsCi <- function(data, givedf){
  min_Ci_ind <- which(data$Ci == min(data$Ci))
  data_new <- slice(data, -which(data$A < data$A[min_Ci_ind]))
  if(givedf =="TRUE"){
    data_new <- as.data.frame(data_new)
  }
  return(data_new)
}

# # Can delete: Apply function and output a list (difficult to work with)
# test <- group_map(.data = cmplt_DAT_grp, .f = exclude_backwardsCi, .keep = TRUE)
# test2 <- purrr::map(test, tibble::as_tibble)
# list2env(test2, envir = .GlobalEnv)

# Apply function to tibble. if .keep = TRUE this throws an error
# Additional info: https://stackoverflow.com/questions/63412850/managing-dplyr-group-by-function-to-keep-the-grouping-variable-when-used-in-comb
DAT_filt_ex <- DAT_filt %>%
  group_by(unique) %>%
  group_modify(~exclude_backwardsCi(data = .x, givedf = TRUE), .keep = FALSE)
DAT_filt_ex <- as.data.frame(DAT_filt_ex)


#important to slice in descending order!!
#complete_sp <- slice(complete_sp, -(7992:8025)) # This one didn't work
#complete_sp <- slice(complete_sp, -(7589:7606)) # This one (somewhat) worked
#complete_sp <- slice(complete_sp, -(7196:7216)) # This one didn't work
#complete_sp <- slice(complete_sp, -(1564:1595)) # This one worked
#complete_sp <- slice(complete_sp, -(799:825)) # Check the points on this one... may not have worked.
#complete_sp <- slice(complete_sp, -(397:402)) # This one worked

##What I'm finding is that removing the back-correction sometimes helps
##Other times it breaks the plantecophys code and doesn't run



# Plotting ACi Curves -----------------------------------------------------

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
 # ggsave(filename1, gg1, path = paste0(getwd(), "/Figures/"))
#}

#The same, but for each unique tree ## Charlie changed this on 11/1/22
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



# Fitting ACi Curves ------------------------------------------------------

# Fit the ACi curves for each species for DAT using fitacis
DAT_fits <- fitacis(DAT_filt_ex, group = "unique", id = "unique",
                    varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                    PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(DAT_fits[[23]], main = coef(DAT_fits)$unique[[23]]) ##keep an eye on #7 as the example
coef(DAT_fits)
#For loop to save all the plots
for (curve in 1:30){
  title <- coef(DAT_fits)$unique[[curve]]
  png(filename = paste0(getwd(), "/Outputs/", title,"_dataci_curve.png"))
  plot(DAT_fits[[curve]], main = title)
  dev.off()
}


#strange curves: #7, 9, 10, 11, 12, 15, 17, 22, 23, 24? 
## Make a dataframe out of coefficients
par_dat <- as.data.frame(coef(DAT_fits), row.names = NULL)
par_dat <- par_dat[-1,]
par_dat <- par_dat %>%
  add_column(method = "dat")
table(par_dat$method) ## number of initial DAT curves

plot(DAT_fits[[15]], main = par_dat$unique[14])
### note that some of these estimates are way incorrect. I think it is a problem of ecophys
### which doesn't seem to  be able to figure out the initial 'back correction'



# Fit the ACi curves for Traditional using fitacis
trad_fits <- fitacis(cmplt_trad, group = "unique", #id = "unique",
                     varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci",
                                     PPFD = "Qin"), fitTPU = FALSE, Tcorrect = TRUE)
plot(trad_fits[[14]], main = coef(trad_fits)$unique[14]) # 20 is pretty ugly
coef(trad_fits)
# For loop to plot all the curves
for (curve in 1:28){
  title <- coef(trad_fits)$unique[[curve]]
  png(filename = paste0(getwd(), "/Outputs/", title,"_tradaci_curve.png"))
  plot(trad_fits[[curve]], main = title)
  dev.off()
}

## Make a dataframe of coefficients
par_trad <- as.data.frame(coef(trad_fits), row.names = NULL)
par_trad <- par_trad[-1,]
par_trad <- par_trad %>% 
  add_column(method = "trad")

# Merge DAT and Trad dfs
par_join <- bind_rows(par_dat, par_trad)
#unique_ids <- rename(unique_ids, "unique" = "?..unique") # I was having some weird problems reading it in
#par_species <- left_join(par_join, unique_ids, by = "unique")
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
  geom_boxplot(aes(x = method, y = Vcmax))
box_both_vcmax

box_both_jmax <- filt_par_species %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = Jmax))
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










# Testing -----------------------------------------------------------------

#Made this separate dataframe and plot to make it easy to inspect any given leaf
# k6706l1 <- complete_sp %>% subset(unique == 'K6706L1')
# ggplot() +
#    geom_point(data = k6706l1, aes(x = Ci, y = A, color = Data_point))

#complete_sp[c(397:402),]  ## K6709L2-2 outliers (back-correction); coef_DAT id: [[16]]
#complete_sp[c(799:825),] ## K6706L1 outliers (back-correction); coef_DAT id: [[7]]
#complete_sp[c(1564:1595),] ## K6709L6 outliers (back-correction); coef_DAT id: [[18]]
#complete_sp[c(7196:7216),] ## K6707L2 outliers (back correction); coef_DAT id: [[11]]
#complete_sp[c(7589:7606),] ## K6707L2-2 outliers (back correction); coef_DAT id: [[12]]


# # Compare sequential Cis
# first_obs <- k6706l1$obs[which(k6706l1$obs==min(k6706l1$obs))]
# k6706l1_new <- k6706l1
# for(i in 1:201){
#   if(k6706l1$obs[i] >= first_obs+1){
#     print(k6706l1$Ci[i])
#     print(k6706l1$Ci[i] - k6706l1$Ci[i+1])
#     if(k6706l1$Ci[i] - k6706l1$Ci[i+1] > 0){
#       k6706l1_new$Data_QC[i] <- "exclude"
#     } 
#   }
# }
# k6706l1_new <- slice(k6706l1_new, -which(k6706l1_new$Data_QC=="exclude"))
# 
# ggplot() +
#   geom_point(data = k6706l1_new, aes(x = Ci, y = A, color = Data_point))
# 
# # Find min Ci and exclude points with Anet below the Anet for that min Ci
# min_Ci_ind <- which(k6706l1$Ci == min(k6706l1$Ci))
# k6706l1_new <- slice(k6706l1, -which(k6706l1$A < k6706l1$A[min_Ci_ind]))
# 
# ggplot() +
#   geom_point(data = k6706l1_new, aes(x = Ci, y = A, color = Data_point))

# Quick test
# test<-fitacis(k6706l1_new, varnames = list(ALEAF = "A", Tleaf = "Tleaf", Ci = "Ci", PPFD = "Qin", Rd =
#                                              "Rd"), Tcorrect = TRUE, fitTPU = FALSE, group = "Data_point")
# plot(test[[1]])
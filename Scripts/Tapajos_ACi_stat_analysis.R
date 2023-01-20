# Statistical Analysis of A/Ci curve data

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(hydroGOF)

wd <- "C://Users/emmel/Desktop/DAT_proj/"
setwd(wd)


# Comparing fits from photosynthesis, plantecophys, and MG code -----------
## Note: for comparisons of plantecophys and MG results, I used Vcmax and Jmax corrected to 25 degrees

## Load in the datasets
params_ecophys <- read.csv(file = paste0(wd, "Results/params_ecophys.csv"), sep = ",", 
                           header = TRUE) %>% 
  filter(method == "dat")
params_ecophys <- add_row(.data = params_ecophys, unique = "K6702L1", Vcmax = 15.361861, 
                          Jmax = 27.750436, Rd = -1.010354,
                          TPU = 1.78906, Vcmax_SE = NA, Jmax_SE = NA, Rd_SE = NA, TPU_SE = NA, 
                          unique.1 = "K6702L1", method = "dat") %>% 
  arrange(unique)
params_ecophys[7,3] <- NA # since it says it is 800000
params_photo <- read.csv(file = paste0(wd, "Results/dat_fit_ex_photo_pars.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000) ## TPU values at 1000 are coded as NA
params_mg <- read.csv(file = paste0(wd, "Results/curve_fitting_MG_out.csv"), sep = ",",
                      header = TRUE) %>% 
  filter(DAT == "Before_DAT", back_filt == "back_filtered")





# Vcmax
rmse(params_ecophys$Vcmax, params_photo$V_cmax)
ecovphoto_vcmax <- ggplot(mapping = aes(x = params_ecophys$Vcmax, y = params_photo$V_cmax, 
                                        color = params_photo$ID)) +
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys Vcmax", y="Photosynthesis Vcmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 150)) + 
  scale_y_continuous(limits = c(1, 150)) +
  annotate(geom = "text", label = "RMSE = 33.32", x = 125, y = 50)
ecovphoto_vcmax


rmse(params_ecophys$Vcmax, params_mg$Best_Vcmax_25C)
ecovMG_vcmax <- ggplot(mapping = aes(x = params_ecophys$Vcmax, y = params_mg$Best_Vcmax_25C, 
                                     color = params_mg$Tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys Vcmax", y="MG Vcmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 170)) + 
  scale_y_continuous(limits = c(1, 170)) +
  annotate(geom = "text", label = "RMSE = 16.63", x = 125, y = 50)
ecovMG_vcmax


rmse(params_photo$V_cmax, params_mg$vcmax_Best_Model)
photovMG_vcmax <- ggplot(mapping = aes(x = params_photo$V_cmax, y = params_mg$vcmax_Best_Model,
                                       color = params_photo$ID))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x = "Photosynthesis Vcmax", y = "MG Vcmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 170)) + 
  scale_y_continuous(limits = c(1, 170)) +
  annotate(geom = "text", label = "RMSE = 22.06", x = 125, y = 50)
photovMG_vcmax


#Jmax
rmse(params_ecophys$Jmax, params_mg$Best.Jmax_25C) ## see if we can fix the 800000 and try again
ecovMG_jmax <- ggplot(mapping = aes(x = params_ecophys$Jmax, y = params_mg$Best.Jmax_25C,
                                    color = params_mg$Tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys Jmax", y="MG Jmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 200)) + 
  scale_y_continuous(limits = c(1, 200)) +
  annotate(geom = "text", label = "RMSE = 23.27", x = 125, y = 50)
ecovMG_jmax
# Note: K6706L1 jmax for plantecophys is 800000, and is not included in the graph


rmse(params_ecophys$Jmax, params_photo$J_max)
ecovphoto_jmax <- ggplot(mapping = aes(x = params_ecophys$Jmax, y = params_photo$J_max, 
                                       color = params_photo$ID))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys Jmax", y="Photosynthesis Jmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 170)) + 
  scale_y_continuous(limits = c(1, 170)) +
  annotate(geom = "text", label = "RMSE = 12.74", x = 125, y = 50)
ecovphoto_jmax
# Note: K6706L1 jmax for plantecophys is 800000, and is not included in the graph


rmse(params_photo$J_max, params_mg$Jmax_Best)
photovMG_jmax <- ggplot(mapping = aes(x = params_photo$J_max, y = params_mg$Jmax_Best,
                                      color = params_photo$ID))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Photosynthesis Jmax", y="MG Jmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 200)) + 
  scale_y_continuous(limits = c(1, 200)) +
  annotate(geom = "text", label = "RMSE = 37.64", x = 125, y = 50)
photovMG_jmax


#TPU
rmse(params_ecophys$TPU, params_mg$TPU_Best)
ecovMG_tpu <- ggplot(mapping = aes(x = params_ecophys$TPU, y = params_mg$TPU_Best,
                                   color = params_mg$Tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys TPU", y="MG TPU", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 11)) + 
  scale_y_continuous(limits = c(1, 11)) +
  annotate(geom = "text", label = "RMSE = 1.03", x = 7, y = 3)
ecovMG_tpu


rmse(params_ecophys$TPU, params_photo$V_TPU)
ecovphoto_tpu <- ggplot(mapping = aes(x = params_ecophys$TPU, y = params_photo$V_TPU,
                                      color = params_photo$ID))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys TPU", y="Photosynthesis TPU", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 10)) + 
  scale_y_continuous(limits = c(1, 10)) +
  annotate(geom = "text", label = "RMSE = 0.32", x = 7, y = 3)
ecovphoto_tpu
# Note: TPU for several of the curves via the photosynthesis package are at 1000 (likely meaning it
# wasn't fit). These are not included in the chart


rmse(params_photo$V_TPU, params_mg$TPU_Best)
photovMG_tpu <- ggplot(mapping = aes(x = params_photo$V_TPU, y = params_mg$TPU_Best, 
                                     color = params_photo$ID))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Photosynthesis TPU", y="MG TPU", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 11)) + 
  scale_y_continuous(limits = c(1, 11)) +
  annotate(geom = "text", label = "RMSE = 1.47", x = 7, y = 3)
photovMG_tpu





# Testing Normality of Plantecophys results -------------------------------------------------------

## read in the data
par_species <- read.csv(file = paste0(wd, "/Results/params_ecophys.csv"), header = TRUE, sep = ",")

## Filter some outliers
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



# Exploration of Fits from Maquelle's code --------------------------------

setwd(paste0(wd, "/Figures"))

curves_df <- read.csv("~/Documents/GitHub/DAT_Proj/Results/curve_fitting_MG_out.csv")


## Separate the concatenated tree ID column
curve_split <- unlist(str_split(curves_df$Tree_id, "_", n=2))
curve_sub <- subset(curve_split, curve_split != "Before_DAT" & curve_split != "Traditional")
curves_df$leaf_id <- curve_sub
curves_df_fixed <- subset(curves_df, select = -c(Tree_id))

leaf_split <- unlist(str_split(curves_df_fixed$leaf_id, "L", n=2))
leaf_sub <- subset(leaf_split, leaf_split != "1" & leaf_split != "2" & leaf_split != "1-1"
                   & leaf_split != 3 & leaf_split != "6" & leaf_split != "8" & leaf_split != "2-2")
curves_df_fixed$tree_id <- leaf_sub
curves_final2 <- curves_df_fixed

#Add in the relative canopy position and four letter code
codes <- read.csv("~/Documents/GitHub/DAT_Proj/Inputs/unique_ids.csv")
canopy_pos <- read.csv("~/Documents/GitHub/DAT_Proj/Inputs/rel_canopy_pos.csv")
names(codes)[1] ="leaf_id"
codes_and_can <- left_join(codes, canopy_pos, by = "k67.id")
names(codes_and_can)[3]="code4let"
codes_and_can <- subset(codes_and_can, select = -code.y)

curves_final <- left_join(curves_final2, codes_and_can, by = "leaf_id")






## Boxplots of filtered vs. nonfiltered data
dat_all <- subset(curves_final, DAT=="Before_DAT")

label_backfilt <- c('Back Filtered', 'Original')
b1 <- ggplot(dat_all, aes(x=back_filt, y=vcmax_Best_Model)) +
  geom_boxplot()+
  labs(x="Dataset", y = "Vcmax")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.position="none")+
  scale_x_discrete(labels=label_backfilt)
b1

b2 <- ggplot(dat_all, aes(x=back_filt, y=Jmax_Best)) +
  geom_boxplot()+
  labs(x="Dataset", y = "Jmax")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.position="none")+
  scale_x_discrete(labels=label_backfilt)
b2

b3 <- ggplot(dat_all, aes(x=back_filt, y=TPU_Best)) +
  geom_boxplot()+
  labs(x="Dataset", y = "TPU")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.position="none")+
  scale_x_discrete(labels=label_backfilt)
b3

# Figures from MG data ----------------------------------------------------


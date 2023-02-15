# Statistical Analysis of A/Ci curve data
rm(list = ls())

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(hydroGOF)
library(Publish) 
library(moments) 
library(vcd)
library(effsize)
library(car)
library(pwr)

wd <- "/Users/charlessouthwick/Documents/GitHub/DAT_Proj/"
setwd(wd)

## Read in the datasets -- These are files without TPU
params_ecophys <- read.csv(file = paste0(wd, "Results/params_ecophys_no_TPU.csv"), sep = ",", 
                           header = TRUE) %>% 
  filter(method == "dat")
params_photo <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_no_TPU.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000) ## TPU values at 1000 are coded as NA
init_mg <- read.csv(file = paste0(wd, "Results/MG_fixed_aci_fits_230213.csv"), sep = ",",
                      header = TRUE)
params_mg <- init_mg %>%
    select(-c(X)) %>% 
    filter(DAT == "Before_DAT") #filters out traditional



# Comparing fits across photosynthesis, plantecophys, and MG code -----------
## Note:

# photosynth vs ecophys ----------
rmse(params_ecophys$Vcmax, params_photo$V_cmax)
ecovphoto_vcmax <- ggplot(mapping = aes(x = params_ecophys$Vcmax, y = params_photo$vcmax_25, 
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
    scale_x_continuous(limits = c(1, 500)) + 
    scale_y_continuous(limits = c(1, 500)) +
    annotate(geom = "text", label = "RMSE", x = 125, y = 50)
ecovphoto_vcmax

#ecophys vs MG ---------------
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
  scale_x_continuous(limits = c(1, 500)) + 
  scale_y_continuous(limits = c(1, 500))
 # annotate(geom = "text", label = "RMSE = 16.63", x = 125, y = 50)
ecovMG_vcmax

#Photosynthesis vs MG temp corrected ----------
cor(params_photo$Best_Vcmax_25C, params_mg$Best_Vcmax_25C)
photovMG_vcmax <- ggplot(mapping = aes(x = params_photo$Best_Vcmax_25C, y = params_mg$Best_Vcmax_25C,
                                       color = params_photo$ID))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x = "Photosynthesis Vcmax", y = "MG Vcmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=14, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 90)) + 
  scale_y_continuous(limits = c(1, 90))+
  annotate(geom = "text", label = "r = 0.963", x = 65, y = 20)
photovMG_vcmax

#Jmax
# Ecophys vs MG
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

# Ecophys vs Photosynthesis
rmse(params_ecophys$Jmax, params_photo$Best_Jmax_25C)
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

#Photosynthesis vs MG
cor(params_photo$Best_Jmax_25C, params_mg$Best.Jmax_25C)
photovMG_jmax <- ggplot(mapping = aes(x = params_photo$Best_Jmax_25C, y = params_mg$Best.Jmax_25C,
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
  scale_x_continuous(limits = c(1, 150)) + 
  scale_y_continuous(limits = c(1, 150)) +
  annotate(geom = "text", label = "r = 0.996", x = 125, y = 50)
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




# Plantecophys result visualization (fix code to make sense with new variables)---------------------------------------
## Boxplots

### Vcmax
box_both_vcmax <- filt_par_species %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = Vcmax))
box_both_vcmax


### Jmax
box_both_jmax <- filt_par_species %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = Jmax))
box_both_jmax



# Stacked Scatters
filt_par_dummy <- mutate(.data = filt_par_species,# makes a dummy variable to plot
                         dummy = if_else(filt_par_species$method == "dat", 0,1))



## Vcmax
scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = Vcmax,
                                                          color = unique)) + 
  geom_line() + 
  geom_point() +
  theme_light() +
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Trad")) +
  xlab("Method")
scat_vcmax


## Jmax
scat_jmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = Jmax,
                                                         color = unique)) + 
  geom_line() + 
  geom_point() +
  theme_light() +
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Trad")) +
  xlab("Method")
scat_jmax


# Testing Assumptions of Plantecophys results -------------------------------------------------------


#### Delete this after fixing code v
# ## read in the data
# par_species <- read.csv(file = paste0(wd, "/Results/params_ecophys.csv"), header = TRUE, sep = ",")
# 
# ## Filter some outliers
# which(par_species$Vcmax > 100)
# which(par_species$Vcmax < 0)
# filt_par_species <- par_species %>% 
#   filter(Vcmax < 100 & Vcmax > 0)
# head(filt_par_species)
###### ^

## Summary Stats
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



#I think a wilcox signed-rank test uses non-parametric, paired data.
wilcox.test(filt_par_species$Vcmax ~ filt_par_species$method, paired = TRUE)






# Photosynthesis results visualization ------------------------------------

##### NEED INITIAL PROCESSING!
#CHANGE ALL VARIABLE NAMES
#Merge the params_photo DAT with traditional params

# Add in the trad: "Results/trad_fits_photo_pars_correct_no_TPU.csv"
photo_trad <- read.csv(file = paste0(wd, "Results/trad_fits_photo_pars_correct_no_TPU.csv"), sep = ",", header = TRUE, na.strings = 1000)
photo_both <- rbind(photo_trad, params_photo)

#change all mg_leaf to photo_leaf
photo_leaf <- photo_both %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           DAT = Data_point,
           leaf_id = ID)

#group data for further analysis -----------------
grp_pho_leaf <- photo_leaf %>% group_by(DAT, leaf_unique) %>%
    summarise(mean_vcmax=mean(Best_Vcmax_25C),
              mean_jmax= mean(Best_Jmax_25C),
              #mean_tpumax= mean(TPU_Best),
              leaf_unique=leaf_unique,
              leaf_id=leaf_id) %>%
    as.data.frame()
summary(grp_pho_leaf)

#Subset variables of interest
pho_stat <- select(photo_leaf, 'DAT', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'leaf_id', 'leaf_unique')
pho_stat$DAT[pho_stat$DAT == "Before_DAT"] <- "DAT"
pho_stat <- rename(pho_stat,
                    method = DAT,
                    vcmax = Best_Vcmax_25C,
                    jmax = Best_Jmax_25C,
                    #tpu = TPU_Best,
)
#Describe factor levels. 0 is traditional, 1 is DAT
pho_stat$method <- factor(pho_stat$method)

#photo stat analysis for vcmax ---------------
#displays grouped summary
pho_summary <- pho_stat %>%
    group_by(method) %>%
    summarise(n_vcmax = length(vcmax),
              mean_vcmax = round(mean(vcmax),3),
              sd_vcmax = round(sd(vcmax),3),
              se_vcmax = sd(vcmax)/sqrt(n())) %>% 
    mutate(low_ci_vcmax = mean_vcmax - qt(1 - (0.05 / 2), n_vcmax - 1) * se_vcmax,
           up_ci_vcmax = mean_vcmax + qt(1 - (0.05 / 2), n_vcmax - 1) * se_vcmax)
print(pho_summary)

#Visualize Vcmax by Method
ci.mean(vcmax ~ method, data = pho_stat)
ci1<-ci.mean(vcmax~method, data=pho_stat)
plot(ci1,title.labels="Method")

#Histogram to visualize
pho_hist<-ggplot(pho_stat, aes(x=vcmax)) + 
    geom_histogram(color="black", fill="white", bins = 8)+
    geom_vline(aes(xintercept=mean(vcmax)),
               color="red", linetype="dashed", linewidth=0.5)+
    theme_classic()
pho_hist
#data are positively skewed
skewness(pho_stat$vcmax)
#value close to 1, should be okay
kurtosis(pho_stat$vcmax)
#It's a bit high

#Test the assumption of equal variances for each group for t-test with Levene's
leveneTest(vcmax ~ method, data = pho_stat)
#Non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances

# Shapiro-Wilk normality test for Vcmax for the one-sample t-test
with(pho_stat, shapiro.test(vcmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(pho_stat, shapiro.test(vcmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(pho_stat, shapiro.test(vcmax[method == "Traditional"])) 
#All are significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

grp_pho_dat <- pho_stat %>%
    filter(method == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_vcmax = length(vcmax),
              mean_vcmax = mean(vcmax),
              sd_vcmax = sd(vcmax),
              se_vcmax = sd(vcmax)/sqrt(n())) %>% 
    mutate(method = "DAT")

grp_pho_trad <- pho_stat %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_vcmax = length(vcmax),
              mean_vcmax = mean(vcmax),
              sd_vcmax = sd(vcmax),
              se_vcmax = sd(vcmax)/sqrt(n())) %>% 
    mutate(method = "Traditional")

grp_pho_all <- rbind(grp_pho_dat, grp_pho_trad)

with(grp_pho_all, shapiro.test(mean_vcmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(grp_pho_all, shapiro.test(mean_vcmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(grp_pho_all, shapiro.test(mean_vcmax[method == "Traditional"])) 
#All are significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

wilcox.test(mean_vcmax ~ method, data = grp_pho_all, paired = TRUE)
#This result is significant

#Effect size for the independent sample t-test:
cohen.ES(test = "t", size = "small") # To remind oneself
cohen.d(mean_vcmax ~ method | Subject(leaf_unique), data=grp_pho_all, paired = TRUE)
#small effect size

#### Power Analysis

d <- cohen.d(mean_vcmax ~ method | Subject(leaf_unique), data=grp_pho_all, paired = TRUE)
d[["estimate"]]

#How many samples to achieve a certain power?
pwr1 <- pwr.t.test(n = , d = d[["estimate"]], power = 0.7, sig.level = 0.05, type = "paired", alternative = "two.sided")
plot(pwr1)

#What was the power of our study?
pwr.t.test(n = 28, d = d[["estimate"]], power = , sig.level = 0.05, type = "paired", alternative = "two.sided")
plot(pwr2)

#photo stat analysis for jmax -------------------------------
pho_jmax_summary <- pho_stat %>%
    group_by(method) %>%
    summarise(n_jmax = length(jmax),
              mean_jmax = round(mean(jmax),3),
              sd_jmax = round(sd(jmax),3),
              se_jmax = sd(jmax)/sqrt(n())) %>% 
    mutate(low_ci_jmax = mean_jmax - qt(1 - (0.05 / 2), n_jmax - 1) * se_jmax,
           up_ci_jmax = mean_jmax + qt(1 - (0.05 / 2), n_jmax - 1) * se_jmax)
print(pho_jmax_summary)

#Visualize Vcmax by Method
ci.mean(jmax ~ method, data = pho_stat)
ci1<-ci.mean(jmax~method, data=pho_stat)
plot(ci1,title.labels="Method")

#Histogram to visualize
pho_jmax_hist<-ggplot(pho_stat, aes(x=jmax)) + 
    geom_histogram(color="black", fill="white", bins = 8)+
    geom_vline(aes(xintercept=mean(jmax)),
               color="red", linetype="dashed", linewidth=0.5)+
    theme_classic()
pho_jmax_hist
#data are positively skewed
skewness(pho_stat$jmax)
#value close to 1, should be okay
kurtosis(pho_stat$jmax)
#It's a bit high

#Test the assumption of equal variances for each group for t-test with Levene's
leveneTest(jmax ~ method, data = pho_stat)
#Non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances

# Shapiro-Wilk normality test for Vcmax for the one-sample t-test
with(pho_stat, shapiro.test(jmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(pho_stat, shapiro.test(jmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(pho_stat, shapiro.test(jmax[method == "Traditional"])) 
#2/3 are significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

grp_pho_jmax_dat <- pho_stat %>%
    filter(method == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_jmax = length(jmax),
              mean_jmax = mean(jmax),
              sd_jmax = sd(jmax),
              se_jmax = sd(jmax)/sqrt(n())) %>% 
    mutate(method = "DAT")

grp_pho_jmax_trad <- pho_stat %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_jmax = length(jmax),
              mean_jmax = mean(jmax),
              sd_jmax = sd(jmax),
              se_jmax = sd(jmax)/sqrt(n())) %>% 
    mutate(method = "Traditional")

grp_pho_jmax_all <- rbind(grp_pho_jmax_dat, grp_pho_jmax_trad)

with(grp_pho_jmax_all, shapiro.test(mean_jmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(grp_pho_jmax_all, shapiro.test(mean_jmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(grp_pho_jmax_all, shapiro.test(mean_jmax[method == "Traditional"])) 
#2/3 are significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

wilcox.test(mean_jmax ~ method, data = grp_pho_jmax_all, paired = TRUE)
#This result is significant

#Effect size for the independent sample t-test:
cohen.ES(test = "t", size = "small") # To remind oneself
cohen.d(mean_jmax ~ method | Subject(leaf_unique), data=grp_pho_jmax_all, paired = TRUE)
#small effect size

#### Power Analysis

d <- cohen.d(mean_jmax ~ method | Subject(leaf_unique), data=grp_pho_jmax_all, paired = TRUE)
d[["estimate"]]

#How many samples to achieve a certain power?
pwr1 <- pwr.t.test(n = , d = d[["estimate"]], power = 0.7, sig.level = 0.05, type = "paired", alternative = "two.sided")
plot(pwr1)

#What was the power of our study?
pwr.t.test(n = 28, d = d[["estimate"]], power = , sig.level = 0.05, type = "paired", alternative = "two.sided")
plot(pwr2)


# Visualization of DAT vs Traditional in Photosynthesis package -- NEED TO CHANGE -----

lab_DATTrad <- c('DAT', 'Traditional')
b4 <- ggplot(photo_leaf, aes(x=DAT, y=Best_Vcmax_25C)) +
    geom_boxplot()+
    labs(x="Method", y = "Vcmax")+
    theme_classic()+
    theme(axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
b4

b5 <- ggplot(photo_leaf, aes(x=DAT, y=Best_Jmax_25C)) +
    geom_boxplot()+
    labs(x="Method", y = "Jmax")+
    theme_classic()+
    theme(axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
b5

# b6 <- ggplot(mg_leaf, aes(x=DAT, y=TPU_Best)) +
#     geom_boxplot()+
#     labs(x="Method", y = "TPU")+
#     theme_classic()+
#     theme(axis.title.x=element_text(size=18, family = "serif"),
#           axis.title.y=element_text(size=18, family = "serif"),
#           axis.text.x=element_text(size=15, family = "serif"),
#           axis.text.y=element_text(size=15, family = "serif"),
#           legend.position="none")+
#     scale_x_discrete(labels=lab_DATTrad)
# b6


# Stacked scatters for photosynthesis package
filt_par_dummy <- mutate(.data = grp_pho_leaf,# makes a dummy variable to plot
                         dummy = if_else(grp_pho_leaf$DAT == "Before_DAT", 0,1))

scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_vcmax,
                                                          color = leaf_unique,
                                                          label = leaf_unique)) +
    geom_line() + 
    geom_point() +
    theme_classic() +
    geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
                    size          = 2.4,
                    box.padding   = 0.25,
                    point.padding = 0.25,
                    segment.size  = 0.2,
                    segment.linetype = 3,
                    direction     = "y",
                    nudge_x = -0.2)+ ## Playing around with this to help visualize
    scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                       expand = expansion(mult=0.3)) +
    labs(x="Method", y="Vcmax")+
    theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"),
          legend.text=element_text(size=9, family = "serif"),
          legend.title=element_text(size=11, family = "serif"),
          legend.position="none")
scat_vcmax

scat_jmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_jmax,
                                                         color = leaf_unique,
                                                         label = leaf_unique)) +
    geom_line() + 
    geom_point() +
    theme_classic() +
    geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
                    size          = 2.8,
                    box.padding   = 0.25,
                    point.padding = 0.25,
                    segment.size  = 0.2,
                    direction     = "y",
                    nudge_x = -0.2)+ ## Playing around with this to help visualize
    scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                       expand = expansion(mult=0.3)) +
    labs(x="Method", y="Jmax")+
    theme(aspect.ratio = 1.5,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"),
          legend.text=element_text(size=9, family = "serif"),
          legend.title=element_text(size=11, family = "serif"),
          legend.position="none")+
    guides(color = guide_legend(title = "Leaf Identifier"))
scat_jmax

## 1:1 plots DAT vs Trad
# By leaf

#Vcmax
leaf_sub_vcmax <- select(grp_pho_leaf, mean_vcmax, DAT, leaf_unique)
leaf_wide_vcmax <- reshape(leaf_sub_vcmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_vcmax)[2:3]=c("vcmax_DAT", "vcmax_Trad")
photo_leaf_vcmax <- ggplot(data = leaf_wide_vcmax, mapping = aes(x = vcmax_Trad,
                                                               y = vcmax_DAT,
                                                               color = leaf_unique))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1)+
    theme_classic()+
    labs(x="Traditional Vcmax", y="DAT Vcmax", col = "Unique Leaf")+
    theme(axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"),
          legend.text=element_text(size=7, family = "serif"),
          legend.title=element_text(size=11, family = "serif"))+
    scale_x_continuous(limits = c(1, 100)) + 
    scale_y_continuous(limits = c(1, 100))
photo_leaf_vcmax

#Jmax
leaf_sub_jmax <- select(grp_pho_leaf, mean_jmax, DAT, leaf_unique)
leaf_wide_jmax <- reshape(leaf_sub_jmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_jmax)[2:3]=c("jmax_DAT", "jmax_Trad")
photo_leaf_jmax <- ggplot(data = leaf_wide_jmax, mapping = aes(x = jmax_Trad,
                                                             y = jmax_DAT,
                                                             color = leaf_unique))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1)+
    theme_classic()+
    labs(x="Traditional Jmax", y="DAT Jmax", col = "Unique Leaf")+
    theme(axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"),
          legend.text=element_text(size=7, family = "serif"),
          legend.title=element_text(size=11, family = "serif"))+
    scale_x_continuous(limits = c(1, 130)) + 
    scale_y_continuous(limits = c(1, 130))
photo_leaf_jmax

# MG data processing ------------------------------------------------

## Separate the concatenated tree ID column
mg_split <- unlist(str_split(init_mg$Tree_id, "_", n=2))
mg_sub <- subset(mg_split, mg_split != "Before_DAT" & mg_split != "Traditional")
init_mg$leaf_id <- mg_sub
mg_complete <- subset(init_mg, select = -c(Tree_id))

leaf_split <- unlist(str_split(mg_complete$leaf_id, "L", n=2))
leaf_sub <- subset(leaf_split, leaf_split != "1" & leaf_split != "2" & leaf_split != "1-1"
                   & leaf_split != 3 & leaf_split != "6" & leaf_split != "8" & leaf_split != "2-2")
mg_complete$tree_id <- leaf_sub

#Add in the relative canopy position
codes <- read.csv("~/Documents/GitHub/DAT_Proj/Inputs/unique_ids.csv") # need to keep this in for leaf_id
canopy_pos <- read.csv("~/Documents/GitHub/DAT_Proj/Inputs/rel_canopy_pos.csv")
names(codes)[1] ="leaf_id"
codes_and_can <- left_join(codes, canopy_pos, by = "k67.id")
names(codes_and_can)[3]="code4let"
codes_and_can <- subset(codes_and_can, select = -code.y)

mg_complete <- left_join(mg_complete, codes_and_can, by = "leaf_id") %>% 
  select(-c(k67.id,code4let))

# For back-filter vs no_back analysis
#mg_all_dat <- subset(mg_complete, DAT=="Before_DAT")

#This is the data with the back correction filtered out
#mg_no_back <- subset(mg_complete, back_filt == "back_filtered")
mg_leaf <- mg_complete %>%
  mutate(leaf_unique = substring(leaf_id, 1, 7))

#group data for further analysis
grp_leaf <- mg_leaf %>% group_by(DAT, leaf_unique) %>%
  summarise(mean_vcmax=mean(Best_Vcmax_25C),
            mean_jmax= mean(Best.Jmax_25C),
            #mean_tpumax= mean(TPU_Best),
            tree_id = tree_id,
            leaf_unique=leaf_unique,
            leaf_id=leaf_id,
            rel_can_pos=rel_can_pos) %>%
  as.data.frame()
summary(grp_leaf)

grp_tree <- mg_leaf %>% group_by(DAT,tree_id) %>% 
  summarise(mean_vcmax=mean(Best_Vcmax_25C),
            mean_jmax= mean(Best.Jmax_25C),
            #mean_tpumax= mean(TPU_Best),
            tree_id = tree_id,
            leaf_unique=leaf_unique,
            leaf_id=leaf_id,
            rel_can_pos=rel_can_pos,
            .groups = 'drop') %>%
  as.data.frame()


### MG Visualizations -------------------------------
## Boxplots of filtered vs. nonfiltered data
# label_backfilt <- c('Back Filtered', 'Original')
# b1 <- ggplot(mg_all_dat, aes(x=back_filt, y=Best_Vcmax_25C)) +
#   geom_boxplot()+
#   labs(x="Dataset", y = "Vcmax")+
#   scale_fill_manual(values=c("#E69F00","#7bccc4"))+
#   theme_classic()+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.position="none")+
#   scale_x_discrete(labels=label_backfilt)
# b1
# 
# b2 <- ggplot(mg_all_dat, aes(x=back_filt, y=Best.Jmax_25C)) +
#   geom_boxplot()+
#   labs(x="Dataset", y = "Jmax")+
#   scale_fill_manual(values=c("#E69F00","#7bccc4"))+
#   theme_classic()+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.position="none")+
#   scale_x_discrete(labels=label_backfilt)
# b2
# 
# b3 <- ggplot(mg_all_dat, aes(x=back_filt, y=TPU_Best)) +
#   geom_boxplot()+
#   labs(x="Dataset", y = "TPU")+
#   scale_fill_manual(values=c("#E69F00","#7bccc4"))+
#   theme_classic()+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.position="none")+
#   scale_x_discrete(labels=label_backfilt)
# b3


#### Visualization of DAT vs Traditional

## Boxplots
lab_DATTrad <- c('DAT', 'Traditional')
b4 <- ggplot(mg_leaf, aes(x=DAT, y=Best_Vcmax_25C)) +
  geom_boxplot()+
  labs(x="Method", y = "Vcmax")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.position="none")+
  scale_x_discrete(labels=lab_DATTrad)
b4

b5 <- ggplot(mg_leaf, aes(x=DAT, y=Best.Jmax_25C)) +
  geom_boxplot()+
  labs(x="Method", y = "Jmax")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.position="none")+
  scale_x_discrete(labels=lab_DATTrad)
b5

b6 <- ggplot(mg_leaf, aes(x=DAT, y=TPU_Best)) +
  geom_boxplot()+
  labs(x="Method", y = "TPU")+
  theme_classic()+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.position="none")+
  scale_x_discrete(labels=lab_DATTrad)
b6

#Just a scatter to understand the spread a bit
# g1 <- ggplot(mg_leaf, aes(x = tree_id, y = Best_Vcmax_25C)) +
#   geom_point() +
#   xlab("Tree")+
#   ylab("Vcmax") +
#   theme_classic() +
#   theme(axis.title.x=element_text(size=11, family = "serif"),
#         axis.title.y=element_text(size=11, family = "serif"),
#         axis.text.x=element_text(size=11, family = "serif"),
#         axis.text.y=element_text(size=11, family = "serif"))
# g1


## Stacked Scatters by leaf
filt_par_dummy <- mutate(.data = grp_leaf,# makes a dummy variable to plot
                         dummy = if_else(grp_leaf$DAT == "Before_DAT", 0,1))

scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_vcmax,
                                                          color = leaf_unique,
                                                          label = tree_id)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
                  size          = 2.4,
                  box.padding   = 0.25,
                  point.padding = 0.25,
                  segment.size  = 0.2,
                  segment.linetype = 3,
                  direction     = "y",
                  nudge_x = -0.2)+ ## Playing around with this to help visualize
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                     expand = expansion(mult=0.3)) +
  labs(x="Method", y="Vcmax")+
  theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
        axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")
scat_vcmax

scat_jmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_jmax,
                                                         color = leaf_unique,
                                                         label = tree_id)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
                  size          = 2.8,
                  box.padding   = 0.25,
                  point.padding = 0.25,
                  segment.size  = 0.2,
                  direction     = "y",
                  nudge_x = -0.2)+ ## Playing around with this to help visualize
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                     expand = expansion(mult=0.3)) +
  labs(x="Method", y="Jmax")+
  theme(aspect.ratio = 1.5,
        axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")+
  guides(color = guide_legend(title = "Leaf Identifier"))
scat_jmax

# scat_tpu <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_tpumax,
#                                                         color = leaf_unique,
#                                                         label = tree_id)) +
#   geom_line() + 
#   geom_point() +
#   theme_classic() +
#   geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
#                   size          = 2.8,
#                   box.padding   = 0.25,
#                   point.padding = 0.25,
#                   segment.size  = 0.2,
#                   direction     = "y",
#                   nudge_x = -0.2)+ ## Playing around with this to help visualize
#   scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
#                      expand = expansion(mult=0.3)) +
#   labs(x="Method", y="TPU")+
#   theme(aspect.ratio = 1.5,
#         axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=9, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"),
#         legend.position ="none")+
#   guides(color = guide_legend(title = "Leaf Identifier"))
# scat_tpu


## Stacked scatter by tree
# filt_par_dummy2 <- mutate(.data = grp_tree,# makes a dummy variable to plot
#                           dummy = if_else(grp_tree$DAT == "Before_DAT", 0,1))
# 
# scat_vcmax2 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_vcmax,
#                                                             color = tree_id,
#                                                             label = tree_id)) +
#   geom_line() + 
#   geom_point() +
#   theme_classic() +
#   geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
#                   size          = 2.8,
#                   box.padding   = 0.25,
#                   point.padding = 0.25,
#                   segment.size  = 0.2,
#                   direction     = "y",
#                   nudge_x = -0.2)+
#   scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
#                      expand = expansion(mult=0.3)) +
#   labs(x="Method", y="Vcmax")+
#   theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
#         axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=9, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"),
#         legend.position="none")
# scat_vcmax2
# 
# scat_jmax2 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_jmax,
#                                                            color = tree_id,
#                                                            label = tree_id)) +
#   geom_line() + 
#   geom_point() +
#   theme_classic() +
#   geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
#                   size          = 2.8,
#                   box.padding   = 0.25,
#                   point.padding = 0.25,
#                   segment.size  = 0.2,
#                   direction     = "y",
#                   nudge_x = -0.2)+
#   scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
#                      expand = expansion(mult=0.3)) +
#   labs(x="Method", y="Jmax")+
#   theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
#         axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=9, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"),
#         legend.position="none")
# scat_jmax2
# 
# scat_tpu2 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_tpumax,
#                                                           color = tree_id,
#                                                           label = tree_id)) +
#   geom_line() + 
#   geom_point() +
#   theme_classic() +
#   geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
#                   size          = 2.8,
#                   box.padding   = 0.25,
#                   point.padding = 0.25,
#                   segment.size  = 0.2,
#                   direction     = "y",
#                   nudge_x = -0.2)+
#   scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
#                      expand = expansion(mult=0.3)) +
#   labs(x="Method", y="TPU")+
#   theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
#         axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=9, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"),
#         legend.position="none")
# scat_tpu2

#grouped scatters labelled with relative canopy heights
# scat_vcmax3 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_vcmax,
#                                                             color = tree_id,
#                                                             label = rel_can_pos)) +
#   geom_line() + 
#   geom_point() +
#   theme_classic() +
#   geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
#                   size          = 2.8,
#                   box.padding   = 0.25,
#                   point.padding = 0.25,
#                   segment.size  = 0.2,
#                   direction     = "y",
#                   nudge_x = -0.2)+
#   scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
#                      expand = expansion(mult=0.3)) +
#   labs(x="Method", y="Vcmax")+
#   theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
#         axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=9, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"),
#         legend.position="none")
# scat_vcmax3
# 
# scat_jmax3 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_jmax,
#                                                            color = tree_id,
#                                                            label = rel_can_pos)) +
#   geom_line() + 
#   geom_point() +
#   theme_classic() +
#   geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
#                   size          = 2.8,
#                   box.padding   = 0.25,
#                   point.padding = 0.25,
#                   segment.size  = 0.2,
#                   direction     = "y",
#                   nudge_x = -0.2)+
#   scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
#                      expand = expansion(mult=0.3)) +
#   labs(x="Method", y="Jmax")+
#   theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
#         axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=9, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"),
#         legend.position="none")
# scat_jmax3
# 
# scat_tpu3 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_tpumax,
#                                                           color = tree_id,
#                                                           label = rel_can_pos)) +
#   geom_line() + 
#   geom_point() +
#   theme_classic() +
#   geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
#                   size          = 2.8,
#                   box.padding   = 0.25,
#                   point.padding = 0.25,
#                   segment.size  = 0.2,
#                   direction     = "y",
#                   nudge_x = -0.2)+
#   scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
#                      expand = expansion(mult=0.3)) +
#   labs(x="Method", y="TPU")+
#   theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
#         axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=9, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"),
#         legend.position="none")
# scat_tpu3


## 1:1 plots DAT vs Trad
# By leaf

#Vcmax
leaf_sub_vcmax <- select(grp_leaf, mean_vcmax, DAT, leaf_unique, tree_id)
leaf_wide_vcmax <- reshape(leaf_sub_vcmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_vcmax)[2:4]=c("vcmax_DAT", "tree_id", "vcmax_Trad")
leaf_wide_vcmax <- subset(leaf_wide_vcmax, select = -tree_id.Traditional)
mng_leaf_vcmax <- ggplot(data = leaf_wide_vcmax, mapping = aes(x = vcmax_Trad,
                                                               y = vcmax_DAT,
                                                               color = leaf_unique))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional Vcmax", y="DAT Vcmax", col = "Unique Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))+
    scale_x_continuous(limits = c(1, 80)) + 
    scale_y_continuous(limits = c(1, 80))
mng_leaf_vcmax

#Jmax
leaf_sub_jmax <- select(grp_leaf, mean_jmax, DAT, leaf_unique, tree_id)
leaf_wide_jmax <- reshape(leaf_sub_jmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_jmax)[2:4]=c("jmax_DAT", "tree_id", "jmax_Trad")
leaf_wide_jmax <- subset(leaf_wide_jmax, select = -tree_id.Traditional)
mng_leaf_jmax <- ggplot(data = leaf_wide_jmax, mapping = aes(x = jmax_Trad,
                                                             y = jmax_DAT,
                                                             color = leaf_unique))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional Jmax", y="DAT Jmax", col = "Unique Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))+
    scale_x_continuous(limits = c(1, 120)) + 
    scale_y_continuous(limits = c(1, 120))
mng_leaf_jmax

# #TPU
# leaf_sub_tpu <- select(grp_leaf, mean_tpumax, DAT, leaf_unique, tree_id)
# leaf_wide_tpu <- reshape(leaf_sub_tpu, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
# names(leaf_wide_tpu)[2:4]=c("tpu_DAT", "tree_id", "tpu_Trad")
# leaf_wide_tpu <- subset(leaf_wide_tpu, select = -tree_id.Traditional)
# mng_leaf_tpu <- ggplot(data = leaf_wide_tpu, mapping = aes(x = tpu_Trad,
#                                                            y = tpu_DAT,
#                                                            color = tree_id))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Traditional TPU", y="DAT TPU", col = "Tree ID")+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"))
# mng_leaf_tpu
# 
# 
# # by tree
# 
# #Vcmax
# tree_sub_vcmax <- select(grp_tree, mean_vcmax, DAT, tree_id)
# tree_wide_vcmax <- reshape(tree_sub_vcmax, idvar = "tree_id", timevar = "DAT", direction = "wide")
# names(tree_wide_vcmax)[2:4]=c("vcmax_DAT", "tree_id", "vcmax_Trad")
# tree_wide_vcmax <- subset(tree_wide_vcmax, select = -tree_id.Traditional)
# mng_tree_vcmax <- ggplot(data = tree_wide_vcmax, mapping = aes(x = vcmax_Trad,
#                                                                y = vcmax_DAT,
#                                                                color = tree_id))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Traditional Vcmax", y="DAT Vcmax", col = "Tree ID")+
#   theme(axis.title.x=element_text(size=11, family = "serif"),
#         axis.title.y=element_text(size=11, family = "serif"),
#         axis.text.x=element_text(size=11, family = "serif"),
#         axis.text.y=element_text(size=11, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"))
# mng_tree_vcmax
# 
# #Jmax
# tree_sub_jmax <- select(grp_tree, mean_jmax, DAT, tree_id)
# tree_wide_jmax <- reshape(tree_sub_jmax, idvar = "tree_id", timevar = "DAT", direction = "wide")
# names(tree_wide_jmax)[2:4]=c("jmax_DAT", "tree_id", "jmax_Trad")
# tree_wide_jmax <- subset(tree_wide_jmax, select = -tree_id.Traditional)
# mng_tree_jmax <- ggplot(data = tree_wide_jmax, mapping = aes(x = jmax_Trad,
#                                                              y = jmax_DAT,
#                                                              color = tree_id))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Traditional Jmax", y="DAT Jmax", col = "Tree ID")+
#   theme(axis.title.x=element_text(size=11, family = "serif"),
#         axis.title.y=element_text(size=11, family = "serif"),
#         axis.text.x=element_text(size=11, family = "serif"),
#         axis.text.y=element_text(size=11, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"))
# mng_tree_jmax
# 
# #TPU
# tree_sub_tpu <- select(grp_tree, mean_tpumax, DAT, tree_id)
# tree_wide_tpu <- reshape(tree_sub_tpu, idvar = "tree_id", timevar = "DAT", direction = "wide")
# names(tree_wide_tpu)[2:4]=c("tpu_DAT", "tree_id", "tpu_Trad")
# tree_wide_tpu <- subset(tree_wide_vcmax, select = -tree_id.Traditional)
# mng_tree_tpu <- ggplot(data = tree_wide_tpu, mapping = aes(x = tpu_Trad,
#                                                            y = tpu_DAT,
#                                                            color = tree_id))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Traditional TPU", y="DAT TPU", col = "Tree ID")+
#   theme(axis.title.x=element_text(size=11, family = "serif"),
#         axis.title.y=element_text(size=11, family = "serif"),
#         axis.text.x=element_text(size=11, family = "serif"),
#         axis.text.y=element_text(size=11, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"))
# mng_tree_tpu


# MG stat analysis for vcmax --------------------------------

#Subset variables of interest
leaf_stat <- select(mg_leaf, 'DAT', 'Best_Vcmax_25C', 'Best.Jmax_25C', 'leaf_id', 'tree_id', 'leaf_unique', 'rel_can_pos')
leaf_stat$DAT[leaf_stat$DAT == "Before_DAT"] <- "DAT"
leaf_stat <- rename(leaf_stat,
                    method = DAT,
                    vcmax = Best_Vcmax_25C,
                    jmax = Best.Jmax_25C,
                    #tpu = TPU_Best,
                    )
#Describe factor levels. 0 is traditional, 1 is DAT
leaf_stat$method <- factor(leaf_stat$method)

#displays grouped summary
leaf_summary <- leaf_stat %>%
  group_by(method) %>%
  summarise(n_vcmax = length(vcmax),
            mean_vcmax = round(mean(vcmax),3),
            sd_vcmax = round(sd(vcmax),3),
            se_vcmax = sd(vcmax)/sqrt(n())) %>% 
  mutate(low_ci_vcmax = mean_vcmax - qt(1 - (0.05 / 2), n_vcmax - 1) * se_vcmax,
         up_ci_vcmax = mean_vcmax + qt(1 - (0.05 / 2), n_vcmax - 1) * se_vcmax)
print(leaf_summary)

#Visualize Vcmax by Method
ci.mean(vcmax ~ method, data = leaf_stat)
ci1<-ci.mean(vcmax~method, data=leaf_stat)
plot(ci1,title.labels="Method")

#Histogram to visualize
leaf_hist<-ggplot(leaf_stat, aes(x=vcmax)) + 
  geom_histogram(color="black", fill="white", bins = 8)+
  geom_vline(aes(xintercept=mean(vcmax)),
             color="red", linetype="dashed", linewidth=0.5)+
  theme_classic()
leaf_hist
#data are positively skewed
skewness(leaf_stat$vcmax)
#value close to 1, should be okay
kurtosis(leaf_stat$vcmax)
#It's a bit high

#Test the assumption of equal variances for each group for t-test with Levene's
leveneTest(vcmax ~ method, data = leaf_stat)
#Non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances

# Shapiro-Wilk normality test for Vcmax for the one-sample t-test
with(leaf_stat, shapiro.test(vcmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(leaf_stat, shapiro.test(vcmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(leaf_stat, shapiro.test(vcmax[method == "Traditional"])) 
#All are significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

grp_dat <- leaf_stat %>%
  filter(method == "DAT") %>% 
  group_by(leaf_unique) %>% 
  summarise(n_vcmax = length(vcmax),
            mean_vcmax = mean(vcmax),
            sd_vcmax = sd(vcmax),
            se_vcmax = sd(vcmax)/sqrt(n())) %>% 
  mutate(method = "DAT")
  
grp_trad <- leaf_stat %>%
  filter(method == "Traditional") %>% 
  group_by(leaf_unique) %>% 
  summarise(n_vcmax = length(vcmax),
            mean_vcmax = mean(vcmax),
            sd_vcmax = sd(vcmax),
            se_vcmax = sd(vcmax)/sqrt(n())) %>% 
  mutate(method = "Traditional")

grp_all <- rbind(grp_dat, grp_trad)

with(grp_all, shapiro.test(mean_vcmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(grp_all, shapiro.test(mean_vcmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(grp_all, shapiro.test(mean_vcmax[method == "Traditional"])) 
#All are significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

wilcox.test(mean_vcmax ~ method, data = grp_all, paired = TRUE)
#This result is significant

#Effect size for the independent sample t-test:
cohen.ES(test = "t", size = "large") # To remind oneself
cohen.d(mean_vcmax ~ method | Subject(leaf_unique), data=grp_all, paired = TRUE)
#small effect size

#### Power Analysis

d <- cohen.d(mean_vcmax ~ method | Subject(leaf_unique), data=grp_all, paired = TRUE)
d[["estimate"]]

#How many samples to achieve a certain power?
pwr1 <- pwr.t.test(n = , d = d[["estimate"]], power = 0.7, sig.level = 0.05, type = "paired", alternative = "two.sided")
plot(pwr1)

#What was the power of our study?
pwr.t.test(n = 28, d = d[["estimate"]], power = , sig.level = 0.05, type = "paired", alternative = "two.sided")
plot(pwr2)


# MG stat analysis for jmax ------------------------------

leaf_jmax_summary <- leaf_stat %>%
    group_by(method) %>%
    summarise(n_jmax = length(jmax),
              mean_jmax = round(mean(jmax),3),
              sd_jmax = round(sd(jmax),3),
              se_jmax = sd(jmax)/sqrt(n())) %>% 
    mutate(low_ci_jmax = mean_jmax - qt(1 - (0.05 / 2), n_jmax - 1) * se_jmax,
           up_ci_jmax = mean_jmax + qt(1 - (0.05 / 2), n_jmax - 1) * se_jmax)
print(leaf_jmax_summary)

#Visualize Vcmax by Method
ci.mean(jmax ~ method, data = leaf_stat)
ci1<-ci.mean(jmax~method, data=leaf_stat)
plot(ci1,title.labels="Method")

#Histogram to visualize
leaf_jmax_hist<-ggplot(leaf_stat, aes(x=jmax)) + 
    geom_histogram(color="black", fill="white", bins = 8)+
    geom_vline(aes(xintercept=mean(vcmax)),
               color="red", linetype="dashed", linewidth=0.5)+
    theme_classic()
leaf_jmax_hist
#data are positively skewed
skewness(leaf_stat$jmax)
#value close to 1, should be okay
kurtosis(leaf_stat$jmax)
#It's a bit high

#Test the assumption of equal variances for each group for t-test with Levene's
leveneTest(jmax ~ method, data = leaf_stat)
#Non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances

# Shapiro-Wilk normality test for Vcmax for the one-sample t-test
with(leaf_stat, shapiro.test(jmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(leaf_stat, shapiro.test(jmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(leaf_stat, shapiro.test(jmax[method == "Traditional"])) 
#2/3 significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

grp_jmax_dat <- leaf_stat %>%
    filter(method == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_jmax = length(jmax),
              mean_jmax = mean(jmax),
              sd_jmax = sd(jmax),
              se_jmax = sd(jmax)/sqrt(n())) %>% 
    mutate(method = "DAT")

grp_jmax_trad <- leaf_stat %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_jmax = length(jmax),
              mean_jmax = mean(jmax),
              sd_jmax = sd(jmax),
              se_jmax = sd(jmax)/sqrt(n())) %>% 
    mutate(method = "Traditional")

grp_jmax_all <- rbind(grp_jmax_dat, grp_jmax_trad)

with(grp_jmax_all, shapiro.test(mean_jmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(grp_jmax_all, shapiro.test(mean_jmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(grp_jmax_all, shapiro.test(mean_jmax[method == "Traditional"])) 
#2/3 significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

wilcox.test(mean_jmax ~ method, data = grp_jmax_all, paired = TRUE)
#This result is significant

#Effect size for the independent sample t-test:
cohen.ES(test = "t", size = "large") # To remind oneself
cohen.d(mean_jmax ~ method | Subject(leaf_unique), data=grp_jmax_all, paired = TRUE)
#small effect size

#### Power Analysis

d <- cohen.d(mean_jmax ~ method | Subject(leaf_unique), data=grp_jmax_all, paired = TRUE)
d[["estimate"]]

#How many samples to achieve a certain power?
pwr1 <- pwr.t.test(n = , d = d[["estimate"]], power = 0.7, sig.level = 0.05, type = "paired", alternative = "two.sided")
plot(pwr1)

#What was the power of our study?
pwr.t.test(n = 28, d = d[["estimate"]], power = , sig.level = 0.05, type = "paired", alternative = "two.sided")
plot(pwr2)


### comparing photosynthesis and MG, barplots and ANOVA? -----------------------

mg_results <- leaf_stat %>%
    mutate(fit_type = "MG") %>% 
    select(-c(rel_can_pos, tree_id))
photo_results <- pho_stat %>% 
    mutate(fit_type = "photo")
all_results <- rbind(mg_results, photo_results)

#vcmax
all_res_sum_vcmax <- all_results %>%
    group_by(method, fit_type) %>%
    summarise(
        sd = sd(vcmax),
        vcmax_mean = mean(vcmax))
all_res_sum_vcmax

ggplot(all_res_sum_vcmax, aes(method, vcmax_mean)) +
    geom_errorbar(
        aes(ymin = vcmax_mean - sd, ymax = vcmax_mean + sd, color = fit_type),
        position = position_dodge(0.3), width = 0.2
    )+
    geom_point(aes(color = fit_type), position = position_dodge(0.3)) +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))+
    scale_y_continuous(limits = c(1, 100))+
    theme_classic()+
    labs(x="Method", y = "Vcmax")+
    theme(axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"))

#jmax
all_res_sum_jmax <- all_results %>%
    group_by(method, fit_type) %>%
    summarise(
        sd = sd(jmax),
        jmax_mean = mean(jmax))
all_res_sum_jmax

ggplot(all_res_sum_jmax, aes(method, jmax_mean)) +
    geom_errorbar(
        aes(ymin = jmax_mean - sd, ymax = jmax_mean + sd, color = fit_type),
        position = position_dodge(0.3), width = 0.2
    )+
    geom_point(aes(color = fit_type), position = position_dodge(0.3)) +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))+
    scale_y_continuous(limits = c(1, 100))+
    theme_classic()+
    labs(x="Method", y = "Jmax")+
    theme(axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"))

#Analysis of variance
library(AICcmodavg)

#vcmax
lm_method <- lm(vcmax ~ method, data = all_results)
lm_fit <- lm(vcmax ~ fit_type, data = all_results)
lm_both <- lm(vcmax ~ method + fit_type, data = all_results)
lm_null <- lm(vcmax ~ 1, data = all_results)

mod_names <- c("Method", "Fit Type", "Method + Fit Type", "Intercept Only")

mod_table <- aictab(list(lm_method, lm_fit, lm_both, lm_null), modnames = mod_names)
mod_table
summary(lm_null)


#jmax
lm2_method <- lm(jmax ~ method, data = all_results)
lm2_fit <- lm(jmax ~ fit_type, data = all_results)
lm2_both <- lm(jmax ~ method + fit_type, data = all_results)
lm2_null <- lm(jmax ~ 1, data = all_results)

mod2_names <- c("Method", "Fit Type", "Method + Fit Type", "Intercept Only")

mod2_table <- aictab(list(lm2_method, lm2_fit, lm2_both, lm2_null), modnames = mod_names)
mod2_table
summary(lm2_method)
kruskal.test(jmax ~ method, data = all_results) #non-parametric ANOVA

leveragePlots(lm2_method)
outlierTest(lm2_method)


cutoff <- 4/122
plot(lm2_method, which=4, cook.levels=cutoff)
abline(h = 4/122,  lty = 4, col = "red")
#We see there's 3 points of influence above our Cook's threshold

# Saving cook's distance values as a separate object
cd <- cooks.distance(lm2_method)
cd.idx <- cd[cd > cutoff] #creating an index with high CD values
view(cd.idx)

#Now let's filter out Case #70, which was present in all 3 outlier tests
all_res_filt <- all_results[-70,]

#create a new model with the new dataset
new_lm2_method <- lm(jmax ~ method, data = all_res_filt)

summary(lm2_method) #Now compare original model
summary(new_lm2_method) #Now compare original model to the new model




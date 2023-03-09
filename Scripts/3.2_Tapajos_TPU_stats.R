### Investigating TPU more fully

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

## Read in the datasets ---------------------

# These are 'photosynthesis' files WITHOUT TPU
pho_dat <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_no_TPU.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000) ## TPU values at 1000 are coded as NA
pho_trad <- read.csv(file = paste0(wd, "Results/trad_fits_photo_pars_correct_no_TPU.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000)
pho_both <- rbind(pho_dat, pho_trad)

pho_leaf <- pho_both %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           leaf_id = ID)

grp_pho <- pho_leaf %>% group_by(method, leaf_unique) %>%
    summarise(mean_vcmax=mean(Best_Vcmax_25C),
              mean_jmax= mean(Best_Jmax_25C),
              mean_tpu= mean(V_TPU),
              leaf_unique=leaf_unique,
              leaf_id=leaf_id) %>%
    as.data.frame()
summary(grp_pho)

#Subset variables of interest
pho_stat <- select(pho_leaf, 'method', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU', 'leaf_id', 'leaf_unique')
pho_stat <- rename(pho_stat,
                       vcmax = Best_Vcmax_25C,
                       jmax = Best_Jmax_25C,
                       tpu = V_TPU,
)
#Describe factor levels. 0 is traditional, 1 is DAT
pho_stat$method <- factor(pho_stat$method)





#These are the files WITH TPU
pho_dat_tpu <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_with_TPU.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000) ## TPU values at 1000 are coded as NA
pho_trad_tpu <- read.csv(file = paste0(wd, "Results/trad_fits_photo_pars_correct_with_TPU.csv"), sep = ",", 
                          header = TRUE, na.strings = 1000)

pho_both_tpu <- rbind(pho_dat_tpu, pho_trad_tpu)

pho_leaf_tpu <- pho_both_tpu %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           leaf_id = ID)

grp_pho_tpu <- pho_leaf_tpu %>% group_by(method, leaf_unique) %>%
    summarise(mean_vcmax=mean(Best_Vcmax_25C),
              mean_jmax= mean(Best_Jmax_25C),
              mean_tpu= mean(V_TPU),
              leaf_unique=leaf_unique,
              leaf_id=leaf_id) %>%
    as.data.frame()
summary(grp_pho_tpu)

#Subset variables of interest
pho_stat_tpu <- select(pho_leaf_tpu, 'method', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU', 'leaf_id', 'leaf_unique')
pho_stat_tpu <- rename(pho_stat_tpu,
                   vcmax = Best_Vcmax_25C,
                   jmax = Best_Jmax_25C,
                   tpu = V_TPU,
)
#Describe factor levels. 0 is traditional, 1 is DAT
pho_stat_tpu$method <- factor(pho_stat_tpu$method)


#photo stat analysis for vcmax ---------------

#Histogram to visualize
pho_hist_tpu<-ggplot(pho_stat_tpu, aes(x=vcmax)) + 
    geom_histogram(color="black", fill="white", bins = 8)+
    geom_vline(aes(xintercept=mean(vcmax)),
               color="red", linetype="dashed", linewidth=0.5)+
    theme_classic()
pho_hist_tpu

grp_pho_dat_tpu <- pho_stat_tpu %>%
    filter(method == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_vcmax = length(vcmax),
              mean_vcmax = mean(vcmax),
              sd_vcmax = sd(vcmax),
              se_vcmax = sd(vcmax)/sqrt(n())) %>% 
    mutate(method = "DAT")

grp_pho_trad_tpu <- pho_stat_tpu %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_vcmax = length(vcmax),
              mean_vcmax = mean(vcmax),
              sd_vcmax = sd(vcmax),
              se_vcmax = sd(vcmax)/sqrt(n())) %>% 
    mutate(method = "Traditional")

grp_pho_all_tpu <- rbind(grp_pho_dat_tpu, grp_pho_trad_tpu)

wilcox.test(mean_vcmax ~ method, data = grp_pho_all_tpu, conf.int = TRUE, paired = TRUE)

#Effect size for the independent sample t-test:
cohen.d(mean_vcmax ~ method | Subject(leaf_unique), data=grp_pho_all_tpu, paired = TRUE)
#small effect size

#### Power Analysis

d <- cohen.d(mean_vcmax ~ method | Subject(leaf_unique), data=grp_pho_all_tpu, paired = TRUE)
d[["estimate"]]

#What was the power?
pwr.t.test(n = 28, d = d[["estimate"]], power = , sig.level = 0.05, type = "paired", alternative = "two.sided")

#photo stat analysis for jmax -------------------------------
#Histogram to visualize
pho_jmax_hist_tpu<-ggplot(pho_stat_tpu, aes(x=jmax)) + 
    geom_histogram(color="black", fill="white", bins = 8)+
    geom_vline(aes(xintercept=mean(jmax)),
               color="red", linetype="dashed", linewidth=0.5)+
    theme_classic()
pho_jmax_hist_tpu

grp_pho_jmax_dat_tpu <- pho_stat_tpu %>%
    filter(method == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_jmax = length(jmax),
              mean_jmax = mean(jmax),
              sd_jmax = sd(jmax),
              se_jmax = sd(jmax)/sqrt(n())) %>% 
    mutate(method = "DAT")

grp_pho_jmax_trad_tpu <- pho_stat_tpu %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_jmax = length(jmax),
              mean_jmax = mean(jmax),
              sd_jmax = sd(jmax),
              se_jmax = sd(jmax)/sqrt(n())) %>% 
    mutate(method = "Traditional")

grp_pho_jmax_all_tpu <- rbind(grp_pho_jmax_dat_tpu, grp_pho_jmax_trad_tpu)

wilcox.test(mean_jmax ~ method, data = grp_pho_jmax_all_tpu, conf.int = TRUE, paired = TRUE)


#Effect size for the independent sample t-test:
cohen.d(mean_jmax ~ method | Subject(leaf_unique), data=grp_pho_jmax_all_tpu, paired = TRUE)
#small effect size

#### Power Analysis

d <- cohen.d(mean_jmax ~ method | Subject(leaf_unique), data=grp_pho_jmax_all_tpu, paired = TRUE)

#What was the power of our study?
pwr.t.test(n = 28, d = d[["estimate"]], power = , sig.level = 0.05, type = "paired", alternative = "two.sided")

#photo stat analysis for TPU ---------------------------------------
grp_pho_tpu_dat_tpu <- pho_stat_tpu %>%
    filter(method == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_tpu = length(tpu),
              mean_tpu = mean(tpu),
              sd_tpu = sd(tpu),
              se_tpu = sd(tpu)/sqrt(n())) %>% 
    mutate(method = "DAT")

grp_pho_tpu_trad_tpu <- pho_stat_tpu %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_tpu = length(tpu),
              mean_tpu = mean(tpu),
              sd_tpu = sd(tpu),
              se_tpu = sd(tpu)/sqrt(n())) %>% 
    mutate(method = "Traditional")

grp_pho_tpu_all_tpu <- rbind(grp_pho_tpu_dat_tpu, grp_pho_tpu_trad_tpu)

#different samples with TPU, therefore can't compare with wilcox test.

# Visualization of DAT vs Traditional in Photosynthesis package -----

lab_DATTrad <- c('DAT', 'Traditional')
b4 <- ggplot(pho_leaf_tpu, aes(x=method, y=Best_Vcmax_25C)) +
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

b5 <- ggplot(pho_leaf_tpu, aes(x=method, y=Best_Jmax_25C)) +
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

b6 <- ggplot(pho_leaf_tpu, aes(x=method, y=V_TPU)) +
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


## 1:1 plots DAT vs Trad
# By leaf

#Vcmax
leaf_sub_tpu <- select(grp_pho_tpu, mean_tpu, method, leaf_unique)
leaf_wide_tpu <- reshape(leaf_sub_tpu, idvar = "leaf_unique", timevar = "method", direction = "wide")
names(leaf_wide_tpu)[2:3]=c("tpu_DAT", "tpu_Trad")
photo_leaf_tpu <- ggplot(data = leaf_wide_tpu, mapping = aes(x = tpu_Trad,
                                                                y = tpu_DAT,
                                                                color = leaf_unique))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1)+
    theme_classic()+
    labs(x="Traditional TPU", y="DAT TPU", col = "Unique Leaf")+
    theme(axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"),
          legend.text=element_text(size=7, family = "serif"),
          legend.title=element_text(size=11, family = "serif"))+
    scale_x_continuous(limits = c(1, 15)) + 
    scale_y_continuous(limits = c(1, 15))
photo_leaf_tpu

### comparing photosynthesis with and without TPU using ANOVA? -----------------

tpu_results <- pho_stat_tpu %>%
    mutate(fit_type = "TPU")
notpu_results <- pho_stat %>% 
    mutate(fit_type = "No_TPU")
all_results <- rbind(tpu_results, notpu_results)

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
    scale_y_continuous(limits = c(1, 75))+
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


kruskal.test(vcmax ~ fit_type, data = all_results) #non-parametric ANOVA
#chi-squared = 0.029, df = 1, p-value = 0.8654. Not significant.


#REPORT THESE!!!!! THIS IS POOLED DATA
wilcox.test(vcmax ~ fit_type, data = all_results, conf.int = TRUE)
wilcox.test(vcmax ~ method, data = all_results, conf.int = TRUE)
#REPORT THESE!!!! THIS IS POOLED DATA


aov_both <- aov(vcmax ~ method + fit_type, data = all_results)
summary(aov_both)

lm_both <- lm(vcmax ~ method + fit_type, data = all_results)
summary(lm_both)
#So there's an increase in Vcmax with TPU relative to no_TPU
#and a decrease in Vcmax with Traditional relative to no_TPU
#This relationship is not significant

#jmax
lm2_method <- lm(jmax ~ method, data = all_results)
lm2_fit <- lm(jmax ~ fit_type, data = all_results)
lm2_both <- lm(jmax ~ method + fit_type, data = all_results)
lm2_null <- lm(jmax ~ 1, data = all_results)

mod2_names <- c("Method", "Fit Type", "Method + Fit Type", "Intercept Only")

mod2_table <- aictab(list(lm2_method, lm2_fit, lm2_both, lm2_null), modnames = mod_names)
mod2_table
summary(lm2_method)


kruskal.test(jmax ~ fit_type, data = all_results) #non-parametric ANOVA
#chi-squared = 0.029, df = 1, p-value = 0.8654. Not significant.
kruskal.test(jmax ~ method, data = all_results) #Not significant with pooled TPU and no-TPU data.

#REPORT THESE!!! THIS IS POOLED DATA
wilcox.test(jmax ~ fit_type, data = all_results, conf.int = TRUE)
wilcox.test(jmax ~ method, data = all_results, conf.int = TRUE)
#REPORT THESE!!! THIS IS POOLED DATA


aov_jmax_both <- aov(jmax ~ method + fit_type, data = all_results)
summary(aov_jmax_both)

lm_jmax_both <- lm(jmax ~ method + fit_type, data = all_results)
summary(lm_jmax_both)
#So there's an increase in Jmax with TPU relative to no_TPU
#and an in Jmax with Traditional relative to no_TPU
#This relationship is not significant


# 
# leveragePlots(lm2_method)
# outlierTest(lm2_method)
# 
# cutoff <- 4/122
# plot(lm2_method, which=4, cook.levels=cutoff)
# abline(h = 4/122,  lty = 4, col = "red")
# #We see there's 3 points of influence above our Cook's threshold
# 
# # Saving cook's distance values as a separate object
# cd <- cooks.distance(lm2_method)
# cd.idx <- cd[cd > cutoff] #creating an index with high CD values
# view(cd.idx)
# 
# #Now let's filter out Case #70, which was present in all 3 outlier tests
# all_res_filt <- all_results[-70,]
# 
# #create a new model with the new dataset
# new_lm2_method <- lm(jmax ~ method, data = all_res_filt)
# 
# summary(lm2_method) #Now compare original model
# summary(new_lm2_method) #Now compare original model to the new model
# 


all_res_flip_vcmax <- all_results %>%
    group_by(fit_type, method) %>%
    summarise(
        sd = sd(vcmax),
        vcmax_mean = mean(vcmax))
all_res_flip_vcmax

ggplot(all_res_flip_vcmax, aes(fit_type, vcmax_mean)) +
    geom_errorbar(
        aes(ymin = vcmax_mean - sd, ymax = vcmax_mean + sd, color = method),
        position = position_dodge(0.3), width = 0.2
    )+
    geom_point(aes(color = method), position = position_dodge(0.3)) +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))+
    scale_y_continuous(limits = c(1, 75))+
    theme_classic()+
    labs(x="Fitting Type", y = "Vcmax")+
    theme(axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"))


all_res_flip_jmax <- all_results %>%
    group_by(fit_type, method) %>%
    summarise(
        sd = sd(jmax),
        jmax_mean = mean(jmax))
all_res_flip_jmax

ggplot(all_res_flip_jmax, aes(fit_type, jmax_mean)) +
    geom_errorbar(
        aes(ymin = jmax_mean - sd, ymax = jmax_mean + sd, color = method),
        position = position_dodge(0.3), width = 0.2
    )+
    geom_point(aes(color = method), position = position_dodge(0.3)) +
    scale_color_manual(values = c("#00AFBB", "#E7B800"))+
    scale_y_continuous(limits = c(1, 100))+
    theme_classic()+
    labs(x="Fitting Type", y = "Jmax")+
    theme(axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"))

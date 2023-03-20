### Investigating TPU more fully

rm(list = ls())

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(hydroGOF)
library(Publish) 
library(moments) 
library(vcd)
library(car)
library(rcompanion) #for wilcoxonPairedRC function

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


#photo stat analysis for vcmax by method, with TPU fit = TRUE, all data ---------------

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
              mean_vcmax = mean(vcmax)) %>% 
    mutate(method = "DAT")

grp_pho_trad_tpu <- pho_stat_tpu %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_vcmax = length(vcmax),
              mean_vcmax = mean(vcmax)) %>% 
    mutate(method = "Traditional")


grp_pho_all_tpu <- rbind(grp_pho_dat_tpu, grp_pho_trad_tpu)

grp_pho_all_tpu %>% group_by(method) %>% summarize(mean = mean(mean_vcmax),
                              median = median(mean_vcmax),
                              sd = sd(mean_vcmax),
                              min = min(mean_vcmax),
                              max = max(mean_vcmax))

boxplot(mean_vcmax~ method, data = grp_pho_all_tpu)

#STATS
wt <- wilcox.test(mean_vcmax ~ method, data = grp_pho_all_tpu, conf.int = TRUE, paired = TRUE)
wt
zval <- qnorm(wt$p.value/2) #z-score applied to a normal distribution
zval

set.seed(67)
wilcoxonPairedRC(x = grp_pho_all_tpu$mean_vcmax,
                 g = grp_pho_all_tpu$method,
                 ci = TRUE,
                 R = 1000) # see King B. M., Rosopa P. J., Minium E. W. (2011) Statistical reasoning in the behavioral sciences (6th ed.). Hoboken, NJ: John Wiley. This is the matched-pairs rank biserial correlation coefficient (rc).


#photo stat analysis for jmax -------------------------------
#Histogram to visualize
pho_jmax_hist_tpu<-ggplot(pho_stat_tpu, aes(x=jmax)) + 
    geom_histogram(color="black", fill="white", bins = 8)+
    geom_vline(aes(xintercept=mean(jmax)),
               color="red", linetype="dashed", linewidth=0.5)+
    theme_classic()

grp_pho_jmax_dat_tpu <- pho_stat_tpu %>%
    filter(method == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_jmax = length(jmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(method = "DAT")

max(grp_pho_jmax_dat_tpu$mean_jmax)

grp_pho_jmax_trad_tpu <- pho_stat_tpu %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_jmax = length(jmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(method = "Traditional")

max(grp_pho_jmax_trad_tpu$mean_jmax)

grp_pho_jmax_all_tpu <- rbind(grp_pho_jmax_dat_tpu, grp_pho_jmax_trad_tpu)

grp_pho_jmax_all_tpu %>% group_by(method) %>% summarize(mean = mean(mean_jmax),
                                                   median = median(mean_jmax),
                                                   sd = sd(mean_jmax),
                                                   min = min(mean_jmax),
                                                   max = max(mean_jmax))

boxplot(mean_jmax~ method, data = grp_pho_jmax_all_tpu)
#STATS
wt1 <- wilcox.test(mean_jmax ~ method, data = grp_pho_jmax_all_tpu, conf.int = TRUE, paired = TRUE)
wt1
zval1 <- qnorm(wt1$p.value/2) #z-score applied to a normal distribution
zval1

set.seed(67)
wilcoxonPairedRC(x = grp_pho_jmax_all_tpu$mean_jmax,
                 g = grp_pho_jmax_all_tpu$method,
                 ci = TRUE,
                 R = 1000)

#photo stat analysis for TPU ---------------------------------------

pho_stat_tpu %>% filter(method == "DAT") %>% summary() #11 curves without TPU
pho_stat_tpu %>% filter(method == "Traditional") %>% summary() #22 curves without TPU

grp_narm_pho_tpu <- pho_stat_tpu %>%
    na.omit(.$tpu) #This leaves 6 trad and 22 DAT

grp_tpu_6trad <- grp_narm_pho_tpu %>% filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarize(vcmax = vcmax,
              jmax = jmax,
              tpu = tpu,
              method = "Traditional")

grp_tpu_6dat <- grp_narm_pho_tpu %>% filter(method == "DAT") %>%
    filter(leaf_unique %in% c("K6706L1", "K6707L1", "K6707L2", "K6709L6", "K6714L1", "K6714L2")) %>%
    group_by(leaf_unique) %>% 
    summarize(vcmax = mean(vcmax),
              jmax = mean(jmax),
              tpu = mean(tpu),
              method = "DAT")


tpu_just6_all <- rbind(grp_tpu_6trad, grp_tpu_6dat)

boxplot(tpu ~ method, data = tpu_just6_all)

#summary stats by group
tpu_just6_all %>% group_by(method) %>% summarize(median = median(tpu),
                                                 mean = mean(tpu),
                                                 min = min(tpu),
                                                 max = max(tpu),
                                                 sd = sd(tpu))


wt2 <- wilcox.test(tpu ~ method, data = tpu_just6_all, conf.int = TRUE, paired = TRUE)
wt2
zval2 <- qnorm(wt2$p.value/2) #z-score applied to a normal distribution
zval2
set.seed(67)
wilcoxonPairedRC(x = tpu_just6_all$tpu,
                 g = tpu_just6_all$method,
                 ci = TRUE,
                 R = 1000)



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

#TRYING THIS OUT ON 3/20 -- grouping by leaf!
tpu_results_grp <- pho_stat_tpu %>%
    group_by(method, leaf_unique) %>% 
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(fit_type = "TPU")

notpu_results_grp <- pho_stat %>%
    group_by(method, leaf_unique) %>% 
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>%
    mutate(fit_type = "No_TPU")
all_results2 <- rbind(tpu_results_grp, notpu_results_grp)

all_res_summ <- all_results2 %>%
    group_by(method, fit_type) %>%
    summarise(
        sd_vcmax = sd(mean_vcmax),
        vcmax_mean = mean(mean_vcmax),
        vcmax_median = median(mean_vcmax),
        sd_jmax = sd(mean_jmax),
        jmax_mean = mean(mean_jmax),
        jmax_median = median(mean_jmax))
all_res_summ


#vcmax
all_res_sum_vcmax <- all_results2 %>%
    group_by(method, fit_type) %>%
    summarise(
        sd = sd(mean_vcmax),
        vcmax_mean = mean(mean_vcmax),
        vcmax_median = median(mean_vcmax))
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
all_res_sum_jmax <- all_results2 %>%
    group_by(method, fit_type) %>%
    summarise(
        sd = sd(jmax),
        jmax_mean = mean(mean_jmax),
        jmax_median = median(mean_jmax))
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


#STATS ---------

#vcmax by fit_type, pooled data
wt3 <- wilcox.test(mean_vcmax ~ fit_type, data = all_results2, conf.int = TRUE, paired = TRUE, exact = FALSE) #Exact = False because can't compute exact confidence interval with zeroes
wt3
zval3 <- qnorm(wt3$p.value/2) #z-score applied to a normal distribution
zval3

set.seed(67)
wilcoxonPairedRC(x = all_results2$mean_vcmax,
                 g = all_results2$fit_type,
                 ci = TRUE,
                 R = 1000)


#jmax by fit_type, pooled data
wt4 <- wilcox.test(mean_jmax ~ fit_type, data = all_results2, conf.int = TRUE, exact = FALSE, paired = TRUE)
wt4
zval4 <- qnorm(wt4$p.value/2) #z-score applied to a normal distribution
zval4
set.seed(67)
wilcoxonPairedRC(x = all_results2$mean_jmax,
                 g = all_results2$fit_type,
                 ci = TRUE,
                 R = 1000)


#What about JUST the DAT data with TPU fitting
dat_all_results <- all_results2 %>% filter(method == "DAT")

wt7 <- wilcox.test(mean_vcmax ~ fit_type, data = dat_all_results, conf.int = TRUE, exact = FALSE, paired = TRUE)
wt7
zval7 <- qnorm(wt7$p.value/2) #z-score applied to a normal distribution
zval7
set.seed(67)
wilcoxonPairedRC(x = dat_all_results$mean_vcmax,
                 g = dat_all_results$fit_type,
                 ci = TRUE,
                 R = 1000)

#Not significant
wt8 <- wilcox.test(mean_jmax ~ fit_type, data = dat_all_results, conf.int = TRUE, exact = FALSE, paired = TRUE)
wt8
zval8 <- qnorm(wt8$p.value/2) #z-score applied to a normal distribution
zval8
abs(zval8)/sqrt(2*20)

set.seed(67)
wilcoxonPairedRC(x = dat_all_results$mean_jmax,
                 g = dat_all_results$fit_type,
                 ci = TRUE,
                 R = 1000)


#Analysis of variance
library(AICcmodavg)

#vcmax
lm_method <- lm(mean_vcmax ~ method, data = all_results2)
lm_fit <- lm(mean_vcmax ~ fit_type, data = all_results2)
lm_both <- lm(mean_vcmax ~ method + fit_type, data = all_results2)
lm_null <- lm(mean_vcmax ~ 1, data = all_results2)

mod_names <- c("Method", "Fit Type", "Method + Fit Type", "Intercept Only")

mod_table <- aictab(list(lm_method, lm_fit, lm_both, lm_null), modnames = mod_names)
mod_table
summary(lm_null)

lm_both <- lm(mean_vcmax ~ method + fit_type, data = all_results2)
summary(lm_both)
lm_fittype <- lm(mean_vcmax ~ fit_type, data = all_results2)
summary(lm_fittype)
#So there's an increase in Vcmax with TPU relative to no_TPU
#and a decrease in Vcmax with Traditional relative to no_TPU
#This relationship is not significant

#jmax
lm2_method <- lm(mean_jmax ~ method, data = all_results2)
lm2_fit <- lm(mean_jmax ~ fit_type, data = all_results2)
lm2_both <- lm(mean_jmax ~ method + fit_type, data = all_results2)
lm2_null <- lm(mean_jmax ~ 1, data = all_results2)

mod2_names <- c("Method", "Fit Type", "Method + Fit Type", "Intercept Only")

mod2_table <- aictab(list(lm2_method, lm2_fit, lm2_both, lm2_null), modnames = mod_names)
mod2_table
summary(lm2_method)


lm_jmax_both <- lm(mean_jmax ~ method + fit_type, data = all_results2)
summary(lm_jmax_both)
#So there's an increase in Jmax with TPU relative to no_TPU
#and an in Jmax with Traditional relative to no_TPU
#This relationship is not significant


# Quick look at TPU fit = TRUE with No overshoot ------
pho_dat_tpu <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_with_TPU.csv"), sep = ",", 
                        header = TRUE, na.strings = 1000) ## TPU values at 1000 are coded as NA
pho_trad_tpu <- read.csv(file = paste0(wd, "Results/trad_fits_photo_pars_correct_with_TPU.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000)

pho_dat_nd_tpu <- pho_dat_tpu %>% subset(ID != "K6707L2" & ID != "K6707L2-2" & ID != "K6707L1" & ID != "K6707L1-1" & ID != "K6709L6" & ID != "K6714L2" & ID != "K6714L1" & ID != "K6702L1" & ID != "K6706L2" & ID != "K6706L1" & ID != "K6709L2")

pho_nd_both <- rbind(pho_trad_tpu, pho_dat_nd_tpu)
pho_nd_leaf <- pho_nd_both %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           DAT = Data_point,
           leaf_id = ID)
pho_nd_stat <- select(pho_nd_leaf, 'DAT', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'leaf_id', 'leaf_unique')
pho_nd_stat$DAT[pho_nd_stat$DAT == "Before_DAT"] <- "DAT"
pho_nd_stat <- rename(pho_nd_stat,
                      method = DAT,
                      vcmax = Best_Vcmax_25C,
                      jmax = Best_Jmax_25C)
pho_nd_stat$method <- factor(pho_nd_stat$method)

grp_pho_nd_dat <- pho_nd_stat %>%
    filter(method == "DAT") %>%
    group_by(leaf_unique) %>%
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(method = "DAT")

sd(grp_pho_nd_dat$mean_vcmax)
sd(grp_pho_nd_dat$mean_jmax)

grp_pho_nd_trad <- pho_nd_stat %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    subset(leaf_unique != "K6702L1"
           & leaf_unique != "K6706L1"
           & leaf_unique != "K6706L2"
           & leaf_unique != "K6707L1"
           & leaf_unique != "K6707L2"
           & leaf_unique != "K6709L6"
           & leaf_unique != "K6714L1"
           & leaf_unique != "K6714L2") %>%
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(method = "Traditional")

sd(grp_pho_nd_trad$mean_vcmax)
sd(grp_pho_nd_trad$mean_jmax)

grp_pho_nd_all <- rbind(grp_pho_nd_dat, grp_pho_nd_trad)

grp_pho_nd_all %>% filter(method == "DAT") %>% summary()
grp_pho_nd_all %>% filter(method == "Traditional") %>% summary()

wt9 <- wilcox.test(mean_vcmax ~ method, data = grp_pho_nd_all, conf.int = TRUE, paired = TRUE)
wt9
zval9 <- qnorm(wt9$p.value/2) #z-score applied to a normal distribution
zval9

set.seed(67)
wilcoxonPairedRC(x = grp_pho_nd_all$mean_vcmax,
                 g = grp_pho_nd_all$method,
                 ci = TRUE,
                 R = 1000)

wt10 <- wilcox.test(mean_jmax ~ method, data = grp_pho_nd_all, conf.int = TRUE, paired = TRUE)
wt10
zval10 <- qnorm(wt4$p.value/2) #z-score applied to a normal distribution
zval10

set.seed(67)
wilcoxonPairedRC(x = grp_pho_nd_all$mean_jmax,
                 g = grp_pho_nd_all$method,
                 ci = TRUE,
                 R = 1000)


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


# Statistical Analysis of A/Ci curve data
rm(list = ls())

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(hydroGOF)
library(Publish)
library(moments) 
library(vcd)
library(car)
library(rstatix) #For wilcox_effsize function

wd <- "C://Users/emmel/Desktop/DAT_proj/"
setwd(wd)

##Read in Datasets -------------------------------------------

pho_dat <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_no_TPU.csv"), sep = ",", 
                    header = TRUE, na.strings = 1000) ## TPU values at 1000 are coded as NA
pho_trad <- read.csv(file = paste0(wd, "Results/trad_fits_photo_pars_correct_no_TPU.csv"), sep = ",", 
                     header = TRUE, na.strings = 1000)


pho_dat_tpu <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_with_TPU.csv"), sep = ",", 
                        header = TRUE, na.strings = 1000) ## TPU values at 1000 are coded as NA
pho_trad_tpu <- read.csv(file = paste0(wd, "Results/trad_fits_photo_pars_correct_with_TPU.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000)

pho_both <- rbind(pho_dat, pho_trad)
pho_both_tpu <- rbind(pho_dat_tpu, pho_trad_tpu)

#No overshoot data

pho_nd <- pho_dat %>% subset(ID != "K6707L2" & ID != "K6707L2-2" & ID != "K6707L1" & ID != "K6707L1-1" & ID != "K6709L6" & ID != "K6714L2" & ID != "K6714L1" & ID != "K6702L1" & ID != "K6706L2" & ID != "K6706L1" & ID != "K6709L2")
pho_nd_tpu <- pho_dat_tpu %>% subset(ID != "K6707L2" & ID != "K6707L2-2" & ID != "K6707L1" & ID != "K6707L1-1" & ID != "K6709L6" & ID != "K6714L2" & ID != "K6714L1" & ID != "K6702L1" & ID != "K6706L2" & ID != "K6706L1" & ID != "K6709L2")

pho_nd_both <- rbind(pho_nd, pho_trad)
pho_nd_both_tpu <- rbind(pho_trad_tpu, pho_nd_tpu)


#Initial Processing/Grouping ----------------------------------

#Full data, no TPU
pho_leaf <- pho_both %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           curv_meth = Data_point,
           leaf_id = ID)
#CHECK! MAY NOT NEED 'V_TPU'
pho_stat <- select(pho_leaf, 'curv_meth', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU', 'leaf_id', 'leaf_unique')
pho_stat$curv_meth[pho_stat$curv_meth == "Before_DAT"] <- "DAT"
pho_stat <- rename(pho_stat,
                   vcmax = Best_Vcmax_25C,
                   jmax = Best_Jmax_25C,
                   tpu = V_TPU) %>% mutate(fit_type = "no_tpu")
#Describe factor levels: 0 is traditional, 1 is DAT
pho_stat$curv_meth <- factor(pho_stat$curv_meth)


#Full data, with TPU
pho_leaf_tpu <- pho_both_tpu %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           curv_meth = Data_point,
           leaf_id = ID)
pho_stat_tpu <- select(pho_leaf_tpu, 'curv_meth', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU', 'leaf_id', 'leaf_unique')
pho_stat_tpu$curv_meth[pho_stat_tpu$curv_meth == "Before_DAT"] <- "DAT"
pho_stat_tpu <- rename(pho_stat_tpu,
                       vcmax = Best_Vcmax_25C,
                       jmax = Best_Jmax_25C,
                       tpu = V_TPU) %>% mutate(fit_type = "tpu")
#Describe factor levels: 0 is traditional, 1 is DAT
pho_stat_tpu$curv_meth <- factor(pho_stat_tpu$curv_meth)


all_results <- rbind(pho_stat, pho_stat_tpu)


#No overshoot, no TPU
pho_nd_leaf <- pho_nd_both %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           curv_meth = Data_point,
           leaf_id = ID)
pho_nd_stat <- select(pho_nd_leaf, 'curv_meth', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU', 'leaf_id', 'leaf_unique')
pho_nd_stat$curv_meth[pho_nd_stat$curv_meth == "Before_DAT"] <- "DAT"
pho_nd_stat <- rename(pho_nd_stat,
                      vcmax = Best_Vcmax_25C,
                      jmax = Best_Jmax_25C,
                      tpu = V_TPU) %>% mutate(fit_type = "no_tpu")
#Describe factor levels: 0 is traditional, 1 is DAT
pho_nd_stat$curv_meth <- factor(pho_nd_stat$curv_meth)


#No overshoot, with TPU
pho_nd_leaf_tpu <- pho_nd_both_tpu %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           curv_meth = Data_point,
           leaf_id = ID)
pho_nd_stat_tpu <- select(pho_nd_leaf_tpu, 'curv_meth', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU', 'leaf_id', 'leaf_unique')
pho_nd_stat_tpu$curv_meth[pho_nd_stat_tpu$curv_meth == "Before_DAT"] <- "DAT"
pho_nd_stat_tpu <- rename(pho_nd_stat_tpu,
                          vcmax = Best_Vcmax_25C,
                          jmax = Best_Jmax_25C,
                          tpu = V_TPU) %>% mutate(fit_type = "tpu")
#Describe factor levels: 0 is traditional, 1 is DAT
pho_nd_stat_tpu$curv_meth <- factor(pho_nd_stat_tpu$curv_meth)



#group data for further analysis -------------------------------------------

#Merge no-TPU with TPU
tpu_results_grp <- pho_stat_tpu %>%
    group_by(curv_meth, leaf_unique) %>% 
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>% 
    mutate(fit_type = "tpu")

notpu_results_grp <- pho_stat %>%
    group_by(curv_meth, leaf_unique) %>% 
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>%
    mutate(fit_type = "no_tpu")

all_results2 <- rbind(tpu_results_grp, notpu_results_grp)
all_results2$fit_type <- factor(all_results2$fit_type)

#No-overshoot data, no TPU
grp_pho_nd_dat <- pho_nd_stat %>%
    filter(curv_meth == "DAT") %>%
    group_by(leaf_unique) %>%
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>% 
    mutate(curv_meth = "DAT",
           fit_type = "no_tpu")

grp_pho_nd_trad <- pho_nd_stat %>%
    filter(curv_meth == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    subset(leaf_unique != "K6702L1"
           & leaf_unique != "K6706L1"
           & leaf_unique != "K6706L2"
           & leaf_unique != "K6707L1"
           & leaf_unique != "K6707L2"
           & leaf_unique != "K6709L6"
           & leaf_unique != "K6714L1"
           & leaf_unique != "K6714L2") %>%
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>% 
    mutate(curv_meth = "Traditional",
           fit_type = "no_tpu")

grp_pho_nd_all <- rbind(grp_pho_nd_dat, grp_pho_nd_trad)


#No-overshoot data, with TPU
grp_pho_nd_dat_tpu <- pho_nd_stat_tpu %>%
    filter(curv_meth == "DAT") %>%
    group_by(leaf_unique) %>%
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>% 
    mutate(curv_meth = "DAT",
           fit_type = "tpu")

grp_pho_nd_trad_tpu <- pho_nd_stat_tpu %>%
    filter(curv_meth == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    subset(leaf_unique != "K6702L1"
           & leaf_unique != "K6706L1"
           & leaf_unique != "K6706L2"
           & leaf_unique != "K6707L1"
           & leaf_unique != "K6707L2"
           & leaf_unique != "K6709L6"
           & leaf_unique != "K6714L1"
           & leaf_unique != "K6714L2") %>%
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>% 
    mutate(curv_meth = "Traditional",
           fit_type = "tpu")

grp_pho_nd_all_tpu <- rbind(grp_pho_nd_dat_tpu, grp_pho_nd_trad_tpu)


nd_complete <- rbind(grp_pho_nd_all, grp_pho_nd_all_tpu)
nd_complete$curv_meth <- factor(nd_complete$curv_meth)
nd_complete$fit_type <- factor(nd_complete$fit_type)

#Just for the TPU analysis

grp_narm_pho_tpu <- pho_stat_tpu %>%
    na.omit(.$tpu) #This leaves 6 trad and 22 DAT

grp_tpu_6trad <- grp_narm_pho_tpu%>%
    filter(curv_meth == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarize(vcmax = vcmax,
              jmax = jmax,
              tpu = tpu,
              curv_meth = "Traditional")

grp_tpu_6dat <- grp_narm_pho_tpu %>%
    filter(curv_meth == "DAT") %>%
    filter(leaf_unique %in% c("K6706L1", "K6707L1", "K6707L2", "K6709L6", "K6714L1", "K6714L2")) %>%
    group_by(leaf_unique) %>% 
    summarize(vcmax = mean(vcmax),
              jmax = mean(jmax),
              tpu = mean(tpu),
              curv_meth = "DAT")

tpu_just6_all <- rbind(grp_tpu_6trad, grp_tpu_6dat)
tpu_just6_all$curv_meth <- factor(tpu_just6_all$curv_meth)


#Summary statistics ----------------------

all_res_summ <- all_results2 %>%
    group_by(fit_type, curv_meth) %>%
    summarise(
        vcmax_mean = mean(vcmax),
        vcmax_median = median(vcmax),
        vcmax_sd = sd(vcmax),
        vcmax_min = min(vcmax),
        vcmax_max = max(vcmax),
        jmax_mean = mean(jmax),
        jmax_median = median(jmax),
        jmax_sd = sd(jmax),
        jmax_min = min(jmax),
        jmax_max = max(jmax))
all_res_summ

#No TPU

all_results2 %>% filter(fit_type == "no_tpu") %>%
ggplot(aes(x=vcmax)) + 
    geom_histogram()

all_results2 %>% filter(fit_type == "no_tpu") %>%
ggplot(aes(x=jmax)) + 
    geom_histogram()

library(reshape2)


#Levene's for homogeneity of variance
all_results2 %>% filter(fit_type == "no_tpu") %>% leveneTest(vcmax ~ curv_meth, data = .)
all_results2 %>% filter(fit_type == "no_tpu") %>% leveneTest(jmax ~ curv_meth, data = .)
# Both non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances

#Shapiro-Wilk test for normality
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax))
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax[curv_meth == "DAT"]))
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax[curv_meth == "Traditional"]))

all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax))
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "Traditional"]))
#All but one deviates from normal

# Yes TPU

all_results2 %>% filter(fit_type == "tpu") %>%
    ggplot(aes(x=vcmax)) + 
    geom_histogram()

all_results2 %>% filter(fit_type == "tpu") %>%
    ggplot(aes(x=jmax)) + 
    geom_histogram()

#Levene's for homogeneity of variance
all_results2 %>% filter(fit_type == "tpu") %>% leveneTest(vcmax ~ curv_meth, data = .)
all_results2 %>% filter(fit_type == "tpu") %>% leveneTest(jmax ~ curv_meth, data = .)
# Both non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances

#Shapiro-Wilk test for normality
all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax))
all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax[curv_meth == "DAT"]))
all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax[curv_meth == "Traditional"]))

all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax))
all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "Traditional"]))
#Most deviate from normal


#No overshoot, No TPU
nd_complete_summ <- nd_complete %>%
    group_by(fit_type, curv_meth) %>%
    summarise(
        vcmax_mean = mean(vcmax),
        vcmax_median = median(vcmax),
        vcmax_sd = sd(vcmax),
        vcmax_min = min(vcmax),
        vcmax_max = max(vcmax),
        jmax_mean = mean(jmax),
        jmax_median = median(jmax),
        jmax_sd = sd(jmax),
        jmax_min = min(jmax),
        jmax_max = max(jmax))
nd_complete_summ

#Levene's test
nd_complete %>% filter(fit_type == "tpu") %>% leveneTest(vcmax ~ curv_meth, data = .)
nd_complete %>% filter(fit_type == "tpu") %>% leveneTest(jmax ~ curv_meth, data = .)
nd_complete %>% filter(fit_type == "no_tpu") %>% leveneTest(vcmax ~ curv_meth, data = .)
nd_complete %>% filter(fit_type == "no_tpu") %>% leveneTest(jmax ~ curv_meth, data = .)
#Met assumption of homogeneity of variances

#Shapiro Wilk Test
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax[curv_meth == "Traditional"]))

nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "Traditional"]))

nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax[curv_meth == "Traditional"]))

nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "Traditional"]))
#They all deviate from normal.


#JUST THE TPU

tpu_just6_all %>% group_by(curv_meth)%>% summarize(median = median(tpu),
                                                mean = mean(tpu),
                                                min = min(tpu),
                                                max = max(tpu),
                                                sd = sd(tpu))


#Checking the symmetric distribution of Wilcoxon -----------------------------------------

#To use wilcoxon paired, we assume the differences between paired samples are distributed symmetrically about the median

#Vcmax
all_res_wide_vcmax <- all_results2 %>%
    dcast(., leaf_unique + fit_type ~ curv_meth, value.var="vcmax") %>%
    mutate(differences = Traditional - DAT)

all_res_wide_vcmax %>% filter(fit_type == "no_tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#That should be fine

all_res_wide_vcmax %>% filter(fit_type == "tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#okay

#Jmax
all_res_wide_jmax <- all_results2 %>%
    dcast(., leaf_unique + fit_type ~ curv_meth, value.var="jmax") %>%
    mutate(differences = Traditional - DAT)

all_res_wide_jmax %>% filter(fit_type == "no_tpu") %>%
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#This is a little rough, but the high value is real data... Perhaps we run both the sign test and the Wilcoxon?

all_res_wide_jmax %>% filter(fit_type == "tpu") %>%
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#This is also a little rough, but the high value is real data... Perhaps we run both the sign test and the Wilcoxon?


# No Overshoot data
#Vcmax
nd_comp_wide_vcmax <- nd_complete %>%
    dcast(., leaf_unique + fit_type ~ curv_meth, value.var="vcmax") %>%
    mutate(differences = Traditional - DAT)

nd_comp_wide_vcmax %>% filter(fit_type == "no_tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#That should be fine

nd_comp_wide_vcmax %>% filter(fit_type == "tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#okay

#Jmax
nd_comp_wide_jmax <- nd_complete %>%
    dcast(., leaf_unique + fit_type ~ curv_meth, value.var="jmax") %>%
    mutate(differences = Traditional - DAT)

nd_comp_wide_jmax %>% filter(fit_type == "no_tpu") %>%
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#This is fine

nd_comp_wide_jmax %>% filter(fit_type == "tpu") %>%
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#This is fine


#So it's just the Jmax where we don't totally meet the Wilcoxon assumptions. We will run both Sign and Wilcoxon for both of these, and make a note of this.



#Stat analysis, clean ------------------------------------------

#When we compared based on fit_type, is this paired data???

#all results, Wilcoxon signed rank test on paired samples (and a few sign tests)

#Vcmax
all_results2 %>%
    group_by(fit_type) %>%
    wilcox_test(data =., vcmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    add_significance()

all_results2 %>%
    group_by(fit_type) %>% wilcox_effsize(data = ., vcmax ~ curv_meth, paired = TRUE)


all_results2 %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>% 
    add_significance()

all_results2 %>%
    group_by(curv_meth) %>% wilcox_effsize(data = ., vcmax ~ fit_type, paired = TRUE)



#Jmax
all_results2 %>%
    group_by(fit_type) %>%
    wilcox_test(data =., jmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    add_significance()

#Note we're running ths sign test here in addition!
all_results2 %>%
    group_by(fit_type) %>%
    sign_test(data =., jmax ~ curv_meth, detailed = TRUE) %>%
    add_significance()

all_results2 %>%
    group_by(fit_type) %>% wilcox_effsize(data = ., jmax ~ curv_meth, paired = TRUE)


all_results2 %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., jmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    add_significance()

#Sign test here. Note we have a different result.
all_results2 %>%
    group_by(curv_meth) %>%
    sign_test(data =., jmax ~ fit_type, detailed = TRUE) %>%
    add_significance()

all_results2 %>%
    group_by(curv_meth) %>% wilcox_effsize(data = ., jmax ~ fit_type, paired = TRUE)

#no overshoot, Wilcoxon

#Vcmax
nd_complete %>%
    group_by(fit_type) %>%
    wilcox_test(data =., vcmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    add_significance()

nd_complete %>%
    group_by(fit_type) %>% wilcox_effsize(data = ., vcmax ~ curv_meth, paired = TRUE)

nd_complete %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    add_significance()

nd_complete %>%
    group_by(curv_meth) %>% wilcox_effsize(data = ., vcmax ~ fit_type, paired = TRUE)

#Jmax
nd_complete %>%
    group_by(fit_type) %>%
    wilcox_test(data =., jmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    add_significance()

nd_complete %>%
    group_by(fit_type) %>% wilcox_effsize(data = ., jmax ~ curv_meth, paired = TRUE)

nd_complete %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., jmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    add_significance()

nd_complete %>%
    group_by(curv_meth) %>% wilcox_effsize(data = ., jmax ~ fit_type, paired = TRUE)


#TPU comparison

tpu_just6_all %>% 
    wilcox_test(data = ., tpu ~ curv_meth, paired = TRUE, detailed = TRUE) %>% 
    add_significance()

tpu_just6_all %>% 
    wilcox_effsize(data = ., tpu ~ curv_meth, paired = TRUE)


#Boxplots ------------------------------------------------------------------


## TPU v No TPU
#Full data
vcmax_all_TPUvNoTPU <- ggplot(all_results2, aes(x = fit_type, y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), 
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("No TPU", "TPU")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(V[cmax]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=13, family = "serif", colour = "gray10"))
vcmax_all_TPUvNoTPU
ggsave(plot = vcmax_all_TPUvNoTPU, "Figures/box_vcmax_all_TPUvNoTPU.png")


jmax_all_TPUvNoTPU <- ggplot(all_results2, aes(x = fit_type, y = jmax, 
                                               fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("No TPU", "TPU"))+
    theme_classic()+
    labs(x="Fit Type",
         y = expression(J[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=13, family = "serif", color = "grey10"))
jmax_all_TPUvNoTPU
ggsave(plot = jmax_all_TPUvNoTPU, "Figures/box_jmax_all_TPUvNoTPU.png")


#No Overshoot
vcmax_nOS_TPUvNoTPU <- ggplot(nd_complete, aes(x = fit_type, y = vcmax,
                                               fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("No TPU", "TPU")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(V[cmax]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=13, family = "serif", color = "grey10"))
vcmax_nOS_TPUvNoTPU
ggsave(plot = vcmax_nOS_TPUvNoTPU, "Figures/box_vcmax_nOS_TPUvNoTPU.png")

jmax_nOS_TPUvNoTPU <- ggplot(nd_complete, aes(x = fit_type, y = jmax,
                                              fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("No TPU", "TPU")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(J[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=13, family = "serif", color = "grey10"))
jmax_nOS_TPUvNoTPU
ggsave(plot = jmax_nOS_TPUvNoTPU, "Figures/box_jmax_nOS_TPUvNoTPU.png")




### Test TPU v No TPU plot
test_vcmax_all_TPUvNoTPU <- ggplot(all_results2, aes(x = curv_meth, y = vcmax, fill = fit_type)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Fit Method", labels = c("Without TPU", "With TPU"), 
      values = c("skyblue", "red"))+
    scale_x_discrete(labels = c("DAT", "Steady-State")) +
    theme_classic()+
    labs(x="Curve Method",
         y = expression(V[cmax]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=13, family = "serif", colour = "gray10"))
test_vcmax_all_TPUvNoTPU
ggsave(plot = test_vcmax_all_TPUvNoTPU, "Figures/test_box_vcmax_all_TPUvNoTPU.png")







## Double Boxplots: All v No OS

pho_nd_stat <- mutate(pho_nd_stat, subset = "nOS")
pho_stat <- mutate(pho_stat, subset = "wOS")
all_and_nOS_noTPU <- rbind(pho_stat, pho_nd_stat)


## No TPU
vcmax_AllvNoOS_noTPU <- all_and_nOS_noTPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), 
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(V[cmax]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=13, family = "serif", colour = "gray10")) +
    scale_y_continuous(limits = c(0, 85))
vcmax_AllvNoOS_noTPU
ggsave(plot = vcmax_AllvNoOS_noTPU, "Figures/box_vcmax_AllvNoOS_noTPU.png")


jmax_AllvNoOS_noTPU <- all_and_nOS_noTPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = jmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("All", "No Overshoot"))+
    theme_classic()+
    labs(x="Fit Type",
         y = expression(J[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=13, family = "serif", color = "grey10")) +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_AllvNoOS_noTPU
ggsave(plot = jmax_AllvNoOS_noTPU, "Figures/box_jmax_AllvNoOS_noTPU.png")





pho_nd_stat_tpu <- mutate(pho_nd_stat_tpu, subset = "nOS")
pho_stat_tpu <- mutate(pho_stat_tpu, subset = "wOS")
all_and_nOS_TPU <- rbind(pho_stat_tpu, pho_nd_stat_tpu)


## With TPU
vcmax_AllvNoOS_TPU <- all_and_nOS_TPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), 
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(V[cmax]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=13, family = "serif", colour = "gray10")) +
    scale_y_continuous(limits = c(0, 85))
vcmax_AllvNoOS_TPU
ggsave(plot = vcmax_AllvNoOS_TPU, "Figures/box_vcmax_AllvNoOS_TPU.png")


jmax_AllvNoOS_TPU <- all_and_nOS_TPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = jmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("All", "No Overshoot"))+
    theme_classic()+
    labs(x="Fit Type",
         y = expression(J[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=13, family = "serif", color = "grey10")) + 
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_AllvNoOS_TPU
ggsave(plot = jmax_AllvNoOS_TPU, "Figures/box_jmax_AllvNoOS_TPU.png")







## Single Plots
lab_DATTrad <- c('DAT', 'Steady-State')

#### All data, No TPU
vcmax_all_noTPU <- ggplot(pho_leaf, aes(x=curv_meth, y=Best_Vcmax_25C)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
vcmax_all_noTPU
ggsave(plot = vcmax_all_noTPU, "Figures/pho_box_DvT_vcmax_noTPU.png")


jmax_all_noTPU <- ggplot(pho_leaf, aes(x=curv_meth, y=Best_Jmax_25C, fill = DAT)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(J[max]*" "*(mu*mol~m^{-2}~s^{-1} *"")))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad) +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_all_noTPU
ggsave(plot = jmax_all_noTPU, "Figures/photo_box_DvT_jmax_noTPU.png")




### All data, with TPU
vcmax_all_TPU <- ggplot(pho_leaf_tpu, aes(x=curv_meth, y=Best_Vcmax_25C)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
vcmax_all_TPU
ggsave(plot = vcmax_all_TPU, "Figures/pho_box_DvT_vcmax_TPU.png")


jmax_all_TPU <- ggplot(pho_leaf_tpu, aes(x=curv_meth, y=Best_Jmax_25C, fill = DAT)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(J[max]*" "*(mu*mol~m^{-2}~s^{-1} *"")))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad) +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_all_TPU
ggsave(plot = jmax_all_TPU, "Figures/photo_box_DvT_jmax_TPU.png")





## No Overshoot, without TPU
nd_vcmax_box_noTPU <- ggplot(grp_pho_nd_all, aes(x=curv_meth, y=vcmax, fill = method)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
nd_vcmax_box_noTPU
ggsave(plot = nd_vcmax_box_noTPU, "Figures/pho_box_noOS_DvT_vcmax_noTPU.png")

nd_jmax_box_noTPU <- ggplot(grp_pho_nd_all, aes(x=curv_meth, y=jmax, fill = method)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(J[max]*" "*(mu*mol~m^{-2}~s^{-1} *"")))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
nd_jmax_box_noTPU
ggsave(plot = nd_jmax_box_noTPU, "Figures/photo_box_nOS_DvT_jmax_noTPU.png")




## 1:1 Plots -----------------------------------------------------------------------------------

## 1:1 plots DAT vs Trad, NO TPU

#Vcmax -- work in progress!!
leaf_sub_vcmax <- select(pho_stat, vcmax, curv_meth, leaf_unique)
leaf_wide_vcmax <- reshape(leaf_sub_vcmax, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_vcmax)[2:3]=c("vcmax_DAT", "vcmax_Trad")
cor1 <- round(cor(leaf_wide_vcmax$vcmax_DAT, leaf_wide_vcmax$vcmax_Trad), 3)
pho_1to1_vcmax_NoTPU <- ggplot(data = leaf_wide_vcmax, mapping = aes(x = vcmax_Trad,
                                                                 y = vcmax_DAT,
                                                                 color = leaf_unique))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1, linetype = 5, cex = 0.6)+
    theme_classic()+
    labs(x = expression("Steady-State "*V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})),
         y = expression("DAT "*V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})),
         col = "Unique Leaf") +
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=12, family = "serif"),
          axis.title.y=element_text(size=12, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=8, family = "serif", color = "gray10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 100)) + 
    scale_y_continuous(limits = c(1, 100)) +
    annotate(geom = "text", label = paste0("r = ", cor1), x = 25, y = 75)
pho_1to1_vcmax_NoTPU
ggsave(plot = pho_1to1_vcmax_NoTPU, "Figures/pho_1to1_vcmax_NoTPU.png")

#Jmax
leaf_sub_jmax <- select(pho_stat, jmax, curv_meth, leaf_unique)
leaf_wide_jmax <- reshape(leaf_sub_jmax, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_jmax)[2:3]=c("jmax_DAT", "jmax_Trad")
cor2 <- round(cor(leaf_wide_jmax$jmax_DAT, leaf_wide_jmax$jmax_Trad), 3)
pho_1to1_jmax_noTPU <- ggplot(data = leaf_wide_jmax, mapping = aes(x = jmax_Trad,
                                                               y = jmax_DAT,
                                                               color = leaf_unique))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6)+
    theme_classic()+
    labs(x = expression("Steady-State "*J[max]*" "*(mu*mol~m^{-2}~s^{-1})),
         y = expression("DAT "*J[max]*" "*(mu*mol~m^{-2}~s^{-1})),
         col = "Unique Leaf")+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=12, family = "serif"),
          axis.title.y=element_text(size=12, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=8, family = "serif", color = "grey10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 130)) + 
    scale_y_continuous(limits = c(1, 130)) +
    annotate(geom = "text", label = paste0("r = ", cor2), x = 100, y = 30)
pho_1to1_jmax_noTPU
ggsave(plot = pho_1to1_jmax_noTPU, "Figures/pho_1to1_jmax_NoTPU.png")



## 1:1 plots DAT vs Trad WITH TPU

#Vcmax
leaf_sub_vcmax_tpu <- select(pho_stat_tpu, vcmax, curv_meth, leaf_unique)
leaf_wide_vcmax_tpu <- reshape(leaf_sub_vcmax_tpu, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_vcmax_tpu)[2:3]=c("vcmax_DAT", "vcmax_Trad")
cor3 <- round(cor(leaf_wide_vcmax_tpu$vcmax_DAT, leaf_wide_vcmax_tpu$vcmax_Trad), 3)
#leaf_wide_vcmax_tpu <- subset(leaf_wide_vcmax_tpu, select = -tree_id.Traditional)
pho_1to1_vcmax_tpu <- ggplot(data = leaf_wide_vcmax_tpu, mapping = aes(x = vcmax_Trad,
                                                               y = vcmax_DAT,
                                                               color = leaf_unique))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1, linetype = 5, cex = 0.6)+
    theme_classic()+
    labs(x = expression("Steady-State "*V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})),
         y = expression("DAT "*V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})),
         col = "Unique Leaf") +
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=12, family = "serif"),
          axis.title.y=element_text(size=12, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=8, family = "serif", color = "gray10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 100)) + 
    scale_y_continuous(limits = c(1, 100)) +
    annotate(geom = "text", label = paste0("r = ", cor3), x = 25, y = 75)
pho_1to1_vcmax_tpu
ggsave(plot = pho_1to1_vcmax_tpu, "Figures/pho_1to1_vcmax_tpu.png")


#Jmax
leaf_sub_jmax_tpu <- select(pho_stat_tpu, jmax, curv_meth, leaf_unique)
leaf_wide_jmax_tpu <- reshape(leaf_sub_jmax_tpu, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_jmax_tpu)[2:3]=c("jmax_DAT", "jmax_Trad")
cor4 <- round(cor(leaf_wide_jmax_tpu$jmax_DAT, leaf_wide_jmax_tpu$jmax_Trad), 3)
pho_1to1_jmax_tpu <- ggplot(data = leaf_wide_jmax_tpu, mapping = aes(x = jmax_Trad,
                                                                       y = jmax_DAT,
                                                                       color = leaf_unique))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6)+
    theme_classic()+
    labs(x = expression("Steady-State "*J[max]*" "*(mu*mol~m^{-2}~s^{-1})),
         y = expression("DAT "*J[max]*" "*(mu*mol~m^{-2}~s^{-1})),
         col = "Unique Leaf")+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=12, family = "serif"),
          axis.title.y=element_text(size=12, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=8, family = "serif", color = "grey10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 130)) + 
    scale_y_continuous(limits = c(1, 130)) +
    annotate(geom = "text", label = paste0("r = ", cor4), x = 100, y = 30)
pho_1to1_jmax_tpu
ggsave(plot = pho_1to1_jmax_tpu, "Figures/pho_1to1_datvtrad_jmax_tpu.png")


#TPU
leaf_sub_tpu <- select(pho_stat_tpu, tpu, curv_meth, leaf_unique)
leaf_wide_tpu <- reshape(leaf_sub_tpu, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_tpu)[2:3]=c("tpu_DAT", "tpu_Trad")
cor5 <- round(cor(leaf_wide_tpu$tpu_DAT, leaf_wide_tpu$tpu_Trad), 3)
pho_1to1_tpu_tpu <- ggplot(data = leaf_wide_tpu, mapping = aes(x = tpu_Trad,
                                                             y = tpu_DAT,
                                                             color = leaf_unique))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1, linetype = 5, cex = 0.6)+
    theme_classic()+
    labs(x="Steady-State TPU", y="DAT TPU", col = "Unique Leaf")+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=12, family = "serif"),
          axis.title.y=element_text(size=12, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=8, family = "serif", color = "grey10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 15)) + 
    scale_y_continuous(limits = c(1, 15)) +
    annotate(geom = "text", label = paste0("r = ", cor5), x = 4, y = 12)
pho_1to1_tpu_tpu
ggsave(plot = pho_1to1_tpu_tpu, "Figures/pho_1to1_datvtrad_tpu.png")



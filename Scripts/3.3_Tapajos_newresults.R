
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
library(rcompanion) #For the point biserial effect size

wd <- "/Users/charlessouthwick/Documents/GitHub/DAT_Proj/"
setwd(wd)

#Read in Datasets -------------------------------------------

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
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(fit_type = "tpu")

notpu_results_grp <- pho_stat %>%
    group_by(DAT, leaf_unique) %>% 
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>%
    mutate(fit_type = "no_tpu")

all_results2 <- rbind(tpu_results_grp, notpu_results_grp) %>%
    rename(curv_meth = method)


#No-overshoot data, no TPU
grp_pho_nd_dat <- pho_nd_stat %>%
    filter(curv_meth == "DAT") %>%
    group_by(leaf_unique) %>%
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(curv_meth = "DAT")

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
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(curv_meth = "Traditional")

grp_pho_nd_all <- rbind(grp_pho_nd_dat, grp_pho_nd_trad)


#No-overshoot data, with TPU
grp_pho_nd_dat_tpu <- pho_nd_stat_tpu %>%
    filter(curv_meth == "DAT") %>%
    group_by(leaf_unique) %>%
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
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
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(curv_meth = "Traditional",
           fit_type = "tpu")

grp_pho_nd_all_tpu <- rbind(grp_pho_nd_dat_tpu, grp_pho_nd_trad_tpu)


nd_complete <- rbind(grp_pho_nd_all, grp_pho_nd_all_tpu) %>%
    rename(curv_meth = method)


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


#photo stat analysis for vcmax -------------------------------------------------------------
#displays grouped summary
pho_summary <- pho_stat %>%
    group_by(curv_meth) %>%
    summarise(n_vcmax = length(vcmax),
              mean_vcmax = round(mean(vcmax),3),
              sd_vcmax = round(sd(vcmax)))
print(pho_summary)

wt <- wilcox.test(mean_vcmax ~ curv_meth, data = grp_pho_all, conf.int = TRUE, paired = TRUE)
wt
zval <- qnorm(wt$p.value/2) #z-score applied to a normal distribution
zval

set.seed(67)
wilcoxonPairedRC(x = grp_pho_all$mean_vcmax,
                 g = grp_pho_all$curv_meth,
                 ci = TRUE,
                 R = 1000)

#photo stat analysis for jmax -------------------------------


wt2 <- wilcox.test(mean_jmax ~ curv_meth, data = grp_pho_jmax_all, conf.int = TRUE, paired = TRUE)
wt2
zval2 <- qnorm(wt2$p.value/2) #z-score applied to a normal distribution
zval2

set.seed(67)
wilcoxonPairedRC(x = grp_pho_jmax_all$mean_jmax,
                 g = grp_pho_jmax_all$curv_meth,
                 ci = TRUE,
                 R = 1000)

# Quick analysis of methods WITHOUT THE OVERSHOOT CURVES --------------------------------------


#Summary statistics

all_res_summ <- all_results2 %>%
    group_by(curv_meth, fit_type) %>%
    summarise(
        sd_vcmax = sd(mean_vcmax),
        vcmax_mean = mean(mean_vcmax),
        vcmax_median = median(mean_vcmax),
        sd_jmax = sd(mean_jmax),
        jmax_mean = mean(mean_jmax),
        jmax_median = median(mean_jmax))
all_res_summ

#No TPU

ci.mean(vcmax ~ curv_meth, data = pho_stat)
ci.mean(jmax ~ curv_meth, data = pho_stat)

ggplot(pho_stat, aes(x=vcmax)) + 
    geom_histogram(color="black", fill="white")

ggplot(pho_stat, aes(x=jmax)) + 
    geom_histogram(color="black", fill="white")

grp_pho_dat <- pho_stat %>%
    filter(curv_meth == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(curv_meth = "DAT")

grp_pho_trad <- pho_stat %>%
    filter(curv_meth == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(curv_meth = "Traditional")

grp_pho_all <- rbind(grp_pho_dat, grp_pho_trad)

grp_pho_all %>% group_by(curv_meth) %>% summarize(mean_vcmax = mean(mean_vcmax),
                                               med_vcmax = median(mean_vcmax),
                                               sd_vcmax = sd(mean_vcmax),
                                               min_vcmax = min(mean_vcmax),
                                               max_vcmax = max(mean_vcmax),
                                               mean_jmax = mean(mean_jmax),
                                               med_jmax = median(mean_jmax),
                                               sd_jmax = sd(mean_jmax),
                                               min_jmax = min(mean_jmax),
                                               max_jmax = max(mean_jmax))

grp_pho_all %>% filter(curv_meth == "DAT") %>% summary()
grp_pho_all %>% filter(curv_meth == "Traditional") %>% summary()

with(grp_pho_all, shapiro.test(mean_vcmax))
with(grp_pho_all, shapiro.test(mean_vcmax[curv_meth == "DAT"]))
with(grp_pho_all, shapiro.test(mean_vcmax[curv_meth == "Traditional"])) 

with(grp_pho_all, shapiro.test(mean_jmax))
with(grp_pho_all, shapiro.test(mean_jmax[curv_meth == "DAT"]))
with(grp_pho_all, shapiro.test(mean_jmax[curv_meth == "Traditional"]))


# Yes TPU

ci.mean(vcmax ~ curv_meth, data = pho_stat_tpu)
ci.mean(jmax ~ curv_meth, data = pho_stat_tpu)

ggplot(pho_stat_tpu, aes(x=vcmax)) + 
    geom_histogram(color="black", fill="white")

ggplot(pho_stat_tpu, aes(x=jmax)) + 
    geom_histogram(color="black", fill="white")

grp_pho_dat_tpu <- pho_stat_tpu %>%
    filter(curv_meth == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(curv_meth = "DAT")

grp_pho_trad_tpu <- pho_stat %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax)) %>% 
    mutate(curv_meth = "Traditional")

grp_pho_all_tpu <- rbind(grp_pho_dat_tpu, grp_pho_trad_tpu)

grp_pho_all_tpu %>% group_by(curv_meth) %>% summarize(mean_vcmax = mean(mean_vcmax),
                                               med_vcmax = median(mean_vcmax),
                                               sd_vcmax = sd(mean_vcmax),
                                               min_vcmax = min(mean_vcmax),
                                               max_vcmax = max(mean_vcmax),
                                               mean_jmax = mean(mean_jmax),
                                               med_jmax = median(mean_jmax),
                                               sd_jmax = sd(mean_jmax),
                                               min_jmax = min(mean_jmax),
                                               max_jmax = max(mean_jmax))

grp_pho_all_tpu %>% filter(curv_meth == "DAT") %>% summary()
grp_pho_all_tpu %>% filter(curv_meth == "Traditional") %>% summary()

with(grp_pho_all_tpu, shapiro.test(mean_vcmax))
with(grp_pho_all_tpu, shapiro.test(mean_vcmax[curv_meth == "DAT"]))
with(grp_pho_all_tpu, shapiro.test(mean_vcmax[curv_meth == "Traditional"])) 

with(grp_pho_all_tpu, shapiro.test(mean_jmax))
with(grp_pho_all_tpu, shapiro.test(mean_jmax[curv_meth == "DAT"]))
with(grp_pho_all_tpu, shapiro.test(mean_jmax[curv_meth == "Traditional"]))



#No overshoot, No TPU
sd(grp_pho_nd_dat$mean_vcmax)
sd(grp_pho_nd_dat$mean_jmax)

sd(grp_pho_nd_trad$mean_vcmax)
sd(grp_pho_nd_trad$mean_jmax)

grp_pho_nd_all %>% filter(curv_meth == "DAT") %>% summary()
grp_pho_nd_all %>% filter(curv_meth == "Traditional") %>% summary()

ci.mean(grp_pho_nd_all$mean_vcmax, data = grp_pho_nd_all)
ci.mean(mean_vcmax ~ curv_meth, data = grp_pho_nd_all)

ci.mean(grp_pho_nd_all$mean_jmax, data = grp_pho_nd_all)
ci.mean(mean_jmax ~ curv_meth, data = grp_pho_nd_all)


#No overshoot, WITH TPU
sd(grp_pho_nd_dat_tpu$mean_vcmax)
sd(grp_pho_nd_dat_tpu$mean_jmax)

sd(grp_pho_nd_trad_tpu$mean_vcmax)
sd(grp_pho_nd_trad_tpu$mean_jmax)

grp_pho_nd_all_tpu %>% filter(curv_meth == "DAT") %>% summary()
grp_pho_nd_all_tpu %>% filter(curv_meth == "Traditional") %>% summary()

ci.mean(grp_pho_nd_all_tpu$mean_vcmax, data = grp_pho_nd_all_tpu)
ci.mean(mean_vcmax ~ curv_meth, data = grp_pho_nd_all_tpu)

ci.mean(grp_pho_nd_all_tpu$mean_jmax, data = grp_pho_nd_all)
ci.mean(mean_jmax ~ curv_meth, data = grp_pho_nd_all_tpu)


#JUST THE TPU

tpu_just6_all %>% group_by(curv_meth)%>% summarize(median = median(tpu),
                                                mean = mean(tpu),
                                                min = min(tpu),
                                                max = max(tpu),
                                                sd = sd(tpu))


wt2 <- wilcox.test(tpu ~ curv_meth, data = tpu_just6_all, conf.int = TRUE, paired = TRUE)
wt2
zval2 <- qnorm(wt2$p.value/2) #z-score applied to a normal distribution
zval2




#Stat analysis, clean ------------------------------------------

#all results, Wilcoxon

all_results2 %>%
    group_by(fit_type) %>%
    wilcox_test(data =., mean_vcmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")

all_results2 %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., mean_vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")

all_results2 %>%
    group_by(fit_type) %>%
    wilcox_test(data =., mean_jmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")

all_results2 %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., mean_jmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")

#no overshoot, Wilcoxon

nd_complete %>%
    group_by(fit_type) %>%
    wilcox_test(data =., mean_vcmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")

nd_complete %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., mean_vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")

nd_complete %>%
    group_by(fit_type) %>%
    wilcox_test(data =., mean_jmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")

nd_complete %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., mean_jmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    adjust_pvalue(method = "bonferroni") %>%
    add_significance("p.adj")

#Effect sizes

set.seed(67)
wilcoxonPairedRC(x = grp_pho_all$mean_vcmax,
                 g = grp_pho_all$curv_meth,
                 ci = TRUE,
                 R = 1000) # see King B. M., Rosopa P. J., Minium E. W. (2011) Statistical reasoning in the behavioral sciences (6th ed.). Hoboken, NJ: John Wiley. This is the matched-pairs rank biserial correlation coefficient (rc).

set.seed(67)
wilcoxonPairedRC(x = grp_pho_jmax_all$mean_jmax,
                 g = grp_pho_jmax_all$curv_meth,
                 ci = TRUE,
                 R = 1000)

set.seed(67)
wilcoxonPairedRC(x = grp_pho_all_tpu$mean_vcmax,
                 g = grp_pho_all_tpu$curv_meth,
                 ci = TRUE,
                 R = 1000)

set.seed(67)
wilcoxonPairedRC(x = grp_pho_jmax_all_tpu$mean_jmax,
                 g = grp_pho_jmax_all_tpu$curv_meth,
                 ci = TRUE,
                 R = 1000)

set.seed(67)
wilcoxonPairedRC(x = grp_pho_nd_all$mean_vcmax,
                 g = grp_pho_nd_all$curv_meth,
                 ci = TRUE,
                 R = 1000)

set.seed(67)
wilcoxonPairedRC(x = grp_pho_nd_all$mean_jmax,
                 g = grp_pho_nd_all$curv_meth,
                 ci = TRUE,
                 R = 1000)


set.seed(67)
wilcoxonPairedRC(x = grp_pho_nd_all_tpu$mean_vcmax,
                 g = grp_pho_nd_all_tpu$curv_meth,
                 ci = TRUE,
                 R = 1000)

set.seed(67)
wilcoxonPairedRC(x = grp_pho_nd_all_tpu$mean_jmax,
                 g = grp_pho_nd_all_tpu$curv_meth,
                 ci = TRUE,
                 R = 1000)

#TPU comparison
wilcox.test(tpu ~ method, data = tpu_just6_all, conf.int = TRUE, paired = TRUE)
set.seed(67)
wilcoxonPairedRC(x = tpu_just6_all$tpu,
                 g = tpu_just6_all$curv_meth,
                 ci = TRUE,
                 R = 1000)


#Boxplots ------------------------------------------------------------------

#Full data
ggplot(all_results2, aes(x = fit_type, y = mean_vcmax, fill = curv_meth)) +
    geom_boxplot()+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), values = c("#00AFBB", "#E7B800"))+
    scale_x_discrete(labels = c("No TPU", "TPU")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(Vc[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif"),
          axis.text.y=element_text(size=13, family = "serif"))
ggsave("Figures/box_vcmax_fulldata.png")


ggplot(all_results2, aes(x = fit_type, y = mean_jmax, fill = curv_meth)) +
    geom_boxplot()+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), values = c("#00AFBB", "#E7B800"))+
    scale_x_discrete(labels = c("No TPU", "TPU"))+
    theme_classic()+
    labs(x="Fit Type",
         y = expression(J[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif"),
          axis.text.y=element_text(size=13, family = "serif"))
ggsave("Figures/box_jmax_fulldata.png")


#No Overshoot
ggplot(nd_complete, aes(x = fit_type, y = mean_vcmax, fill = curv_meth)) +
    geom_boxplot()+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), values = c("#00AFBB", "#E7B800"))+
    scale_x_discrete(labels = c("No TPU", "TPU")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(Vc[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif"),
          axis.text.y=element_text(size=13, family = "serif"))
ggsave("Figures/box_vcmax_nd.png")

ggplot(nd_complete, aes(x = fit_type, y = mean_jmax, fill = curv_meth)) +
    geom_boxplot()+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), values = c("#00AFBB", "#E7B800"))+
    scale_x_discrete(labels = c("No TPU", "TPU")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(J[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif"),
          axis.text.y=element_text(size=13, family = "serif"))
ggsave("Figures/box_jmax_nd.png")

## 1:1 Plots -----------------------------------------------------------------------------------

## 1:1 plots DAT vs Trad, NO TPU

#Vcmax -- work in progress!!
leaf_sub_vcmax <- select(pho_stat, vcmax, curv_meth, leaf_unique)
leaf_wide_vcmax <- reshape(leaf_sub_vcmax, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_vcmax)[2:3]=c("vcmax_DAT", "vcmax_Trad")
cor1 <- round(cor(leaf_wide_vcmax$vcmax_DAT, leaf_wide_vcmax$vcmax_Trad), 3)
pho_1to1_vcmax <- ggplot(data = leaf_wide_vcmax, mapping = aes(x = vcmax_Trad,
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
    scale_y_continuous(limits = c(1, 100)) +
    annotate(geom = "text", label = paste0("r = ", cor1), x = 75, y = 25)
pho_1to1_vcmax
ggsave("Figures/pho_1to1_datvtrad_vcmax.png")

#Jmax
leaf_sub_jmax <- select(pho_stat, jmax, curv_meth, leaf_unique)
leaf_wide_jmax <- reshape(leaf_sub_jmax, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_jmax)[2:3]=c("jmax_DAT", "jmax_Trad")
cor2 <- round(cor(leaf_wide_jmax$jmax_DAT, leaf_wide_jmax$jmax_Trad), 3)
pho_1to1_jmax <- ggplot(data = leaf_wide_jmax, mapping = aes(x = jmax_Trad,
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
    scale_y_continuous(limits = c(1, 130)) +
    annotate(geom = "text", label = paste0("r = ", cor2), x = 100, y = 30)
photo_1to1_jmax
ggsave("Figures/pho_1to1_datvtrad_jmax.png")



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
    scale_y_continuous(limits = c(1, 80)) +
    annotate(geom = "text", label = paste0("r = ", cor3), x = 100, y = 30)
pho_1to1_vcmax_tpu
ggsave("Figures/pho_1to1_datvtrad_vcmax_tpu.png")

#Jmax
leaf_sub_jmax_tpu <- select(pho_stat_tpu, jmax, curv_meth, leaf_unique)
leaf_wide_jmax_tpu <- reshape(leaf_sub_jmax_tpu, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_jmax_tpu)[2:3]=c("jmax_DAT", "jmax_Trad")
cor4 <- round(cor(leaf_wide_jmax_tpu$jmax_DAT, leaf_wide_jmax_tpu$jmax_Trad), 3)
pho_1to1_jmax_tpu <- ggplot(data = leaf_wide_jmax_tpu, mapping = aes(x = jmax_Trad,
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
    scale_y_continuous(limits = c(1, 130)) +
    annotate(geom = "text", label = paste0("r = ", cor4), x = 100, y = 30)
pho_1to1_jmax_tpu
ggsave("Figures/pho_1to1_datvtrad_jmax_tpu.png")

#TPU
leaf_sub_tpu <- select(pho_stat_tpu, tpu, curv_meth, leaf_unique)
leaf_wide_tpu <- reshape(leaf_sub_tpu, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_tpu)[2:3]=c("tpu_DAT", "tpu_Trad")
cor5 <- round(cor(leaf_wide_tpu$tpu_DAT, leaf_wide_tpu$tpu_Trad), 3)
pho_1to1_tpu_tpu <- ggplot(data = leaf_wide_tpu, mapping = aes(x = tpu_Trad,
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
    scale_y_continuous(limits = c(1, 15)) +
    annotate(geom = "text", label = paste0("r = ", cor5), x = 100, y = 30)
pho_1to1_tpu_tpu

ggsave("Figures/pho_1to1_datvtrad_tpu_tpu.png")



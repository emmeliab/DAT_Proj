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

## Read in the datasets

#params_ecophys <- read.csv(file = paste0(wd, "Results/params_ecophys_no_TPU.csv"), sep = ",", 
#                           header = TRUE) %>% 
#  filter(method == "dat")

canopy_pos <- read.csv(file = paste0(wd, "Inputs/rel_canopy_pos.csv"))


params_photo <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_no_TPU.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000) %>% ## TPU values at 1000 are coded as NA
    mutate(tree_id = substr(ID, 1, 5)) %>% 
    left_join(., canopy_pos, by = "tree_id") %>% 
    select(-code)

params_photo_tpu <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_with_TPU.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000) %>%
    mutate(tree_id = substr(ID, 1, 5)) %>% 
    left_join(., canopy_pos, by = "tree_id") %>% 
    select(-code)

photo_trad <- read.csv(file = paste0(wd, "Results/trad_fits_photo_pars_correct_no_TPU.csv"), sep = ",", header = TRUE, na.strings = 1000) %>%
    mutate(tree_id = substr(ID, 1, 5)) %>% 
    left_join(., canopy_pos, by = "tree_id") %>% 
    select(-code)


# Photosynthesis results visualization ------------------------------------

#Merge the params_photo DAT with traditional params

photo_both <- rbind(photo_trad, params_photo)


photo_leaf <- photo_both %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           DAT = Data_point,
           leaf_id = ID)

#group data for further analysis -----------------
grp_pho_leaf <- photo_leaf %>% group_by(DAT, leaf_unique) %>%
    summarise(mean_vcmax=mean(Best_Vcmax_25C),
              mean_jmax= mean(Best_Jmax_25C),
              leaf_unique=leaf_unique,
              leaf_id=leaf_id,
              rel_can_pos=rel_can_pos) %>%
    as.data.frame()
summary(grp_pho_leaf)

#Subset variables of interest
pho_stat <- select(photo_leaf, 'DAT', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'leaf_id', 'leaf_unique', 'rel_can_pos')
pho_stat$DAT[pho_stat$DAT == "Before_DAT"] <- "DAT"
pho_stat <- rename(pho_stat,
                    method = DAT,
                    vcmax = Best_Vcmax_25C,
                    jmax = Best_Jmax_25C
)
#Describe factor levels. 0 is traditional, 1 is DAT
pho_stat$method <- factor(pho_stat$method)

#Create data frame with NONE of the overshoots ------------------------------------------------

pho_nd <- params_photo %>% subset(ID != "K6707L2" & ID != "K6707L2-2" & ID != "K6707L1" & ID != "K6707L1-1" & ID != "K6709L6" & ID != "K6714L2" & ID != "K6714L1" & ID != "K6702L1" & ID != "K6706L2" & ID != "K6706L1" & ID != "K6709L2")

pho_nd_both <- rbind(photo_trad, pho_nd)
pho_nd_leaf <- pho_nd_both %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           DAT = Data_point,
           leaf_id = ID)
pho_nd_stat <- select(pho_nd_leaf, 'DAT', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'leaf_id', 'leaf_unique', 'rel_can_pos')
pho_nd_stat$DAT[pho_nd_stat$DAT == "Before_DAT"] <- "DAT"
pho_nd_stat <- rename(pho_nd_stat,
                   method = DAT,
                   vcmax = Best_Vcmax_25C,
                   jmax = Best_Jmax_25C)
pho_nd_stat$method <- factor(pho_nd_stat$method)

#photo stat analysis for vcmax -------------------------------------------------------------
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

grp_pho_dat <- pho_stat %>%
    filter(method == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_vcmax = length(vcmax),
              mean_vcmax = mean(vcmax),
              rel_can_pos = mean(rel_can_pos)) %>% 
    mutate(method = "DAT")

grp_pho_trad <- pho_stat %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_vcmax = length(vcmax),
              mean_vcmax = mean(vcmax),
              rel_can_pos = mean(rel_can_pos)) %>% 
    mutate(method = "Traditional")

grp_pho_all <- rbind(grp_pho_dat, grp_pho_trad)

grp_pho_all %>% group_by(method) %>% summarize(mean = mean(mean_vcmax),
                                               median = median(mean_vcmax),
                                               sd = sd(mean_vcmax),
                                               min = min(mean_vcmax),
                                               max = max(mean_vcmax))

ci.mean(mean_vcmax ~ method, data = grp_pho_all)
ci1<-ci.mean(mean_vcmax~method, data=grp_pho_all)
plot(ci1,title.labels="Method")

grp_pho_all %>% filter(method == "DAT") %>% summary()
grp_pho_all %>% filter(method == "Traditional") %>% summary()

with(grp_pho_all, shapiro.test(mean_vcmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(grp_pho_all, shapiro.test(mean_vcmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(grp_pho_all, shapiro.test(mean_vcmax[method == "Traditional"])) 
#All are significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

wt <- wilcox.test(mean_vcmax ~ method, data = grp_pho_all, conf.int = TRUE, paired = TRUE)
wt
zval <- qnorm(wt$p.value/2) #z-score applied to a normal distribution
zval

set.seed(67)
wilcoxonPairedRC(x = grp_pho_all$mean_vcmax,
                 g = grp_pho_all$method,
                 ci = TRUE,
                 R = 1000)

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

grp_pho_jmax_dat <- pho_stat %>%
    filter(method == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_jmax = length(jmax),
              mean_jmax = mean(jmax),
              rel_can_pos = mean(rel_can_pos)) %>% 
    mutate(method = "DAT")

max(grp_pho_jmax_dat$mean_jmax)

grp_pho_jmax_trad <- pho_stat %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_jmax = length(jmax),
              mean_jmax = mean(jmax),
              rel_can_pos = mean(rel_can_pos)) %>% 
    mutate(method = "Traditional")

max(grp_pho_jmax_trad$mean_jmax)

grp_pho_jmax_all <- rbind(grp_pho_jmax_dat, grp_pho_jmax_trad)

grp_pho_jmax_all %>% group_by(method) %>% summarize(mean = mean(mean_jmax),
                                               median = median(mean_jmax),
                                               sd = sd(mean_jmax),
                                               min = min(mean_jmax),
                                               max = max(mean_jmax))

grp_pho_jmax_all %>% filter(method == "DAT") %>% summary()
grp_pho_jmax_all %>% filter(method == "Traditional") %>% summary()

with(grp_pho_jmax_all, shapiro.test(mean_jmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(grp_pho_jmax_all, shapiro.test(mean_jmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(grp_pho_jmax_all, shapiro.test(mean_jmax[method == "Traditional"])) 
#2/3 are significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

wt2 <- wilcox.test(mean_jmax ~ method, data = grp_pho_jmax_all, conf.int = TRUE, paired = TRUE)
wt2
zval2 <- qnorm(wt2$p.value/2) #z-score applied to a normal distribution
zval2

set.seed(67)
wilcoxonPairedRC(x = grp_pho_jmax_all$mean_jmax,
                 g = grp_pho_jmax_all$method,
                 ci = TRUE,
                 R = 1000)

# Quick analysis of methods WITHOUT THE OVERSHOOT CURVES --------------------------------------

grp_pho_nd_dat <- pho_nd_stat %>%
    filter(method == "DAT") %>%
    group_by(leaf_unique) %>%
    summarise(mean_vcmax = mean(vcmax),
              mean_jmax = mean(jmax),
              rel_can_pos = mean(rel_can_pos)) %>% 
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
              mean_jmax = mean(jmax),
              rel_can_pos = mean(rel_can_pos)) %>% 
    mutate(method = "Traditional")

sd(grp_pho_nd_trad$mean_vcmax)
sd(grp_pho_nd_trad$mean_jmax)

grp_pho_nd_all <- rbind(grp_pho_nd_dat, grp_pho_nd_trad)

grp_pho_nd_all %>% filter(method == "DAT") %>% summary()
grp_pho_nd_all %>% filter(method == "Traditional") %>% summary()

ci.mean(grp_pho_nd_all$mean_vcmax, data = grp_pho_nd_all)
ci <- ci.mean(grp_pho_nd_all$mean_vcmax, data = grp_pho_nd_all)
ci.mean(mean_jmax ~ method, data = grp_pho_nd_all)

wt3 <- wilcox.test(mean_vcmax ~ method, data = grp_pho_nd_all, conf.int = TRUE, paired = TRUE)
wt3
zval3 <- qnorm(wt3$p.value/2) #z-score applied to a normal distribution
zval3

set.seed(67)
wilcoxonPairedRC(x = grp_pho_nd_all$mean_vcmax,
                 g = grp_pho_nd_all$method,
                 ci = TRUE,
                 R = 1000)

wt4 <- wilcox.test(mean_jmax ~ method, data = grp_pho_nd_all, conf.int = TRUE, paired = TRUE)
wt4
zval4 <- qnorm(wt4$p.value/2) #z-score applied to a normal distribution
zval4

set.seed(67)
wilcoxonPairedRC(x = grp_pho_nd_all$mean_jmax,
                 g = grp_pho_nd_all$method,
                 ci = TRUE,
                 R = 1000)


# Visualization of DAT vs Traditional in Photosynthesis package -----

lab_DATTrad <- c('DAT', 'Steady-State')
b4 <- ggplot(photo_leaf, aes(x=DAT, y=Best_Vcmax_25C)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression("Vcmax "*(mu*mol~m^{-2}~s^{-1} *"")))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
b4
ggsave(plot = b4, "Figures/photo_box_datvtrad_vcmax.png")

b5 <- ggplot(photo_leaf, aes(x=DAT, y=Best_Jmax_25C, fill = DAT)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression("Jmax "*(mu*mol~m^{-2}~s^{-1} *"")))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
b5
ggsave(plot = b5, "Figures/photo_box_datvtrad_jmax.png")


nd_vcmax_box <- ggplot(grp_pho_nd_all, aes(x=method, y=mean_vcmax, fill = method)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression("Vcmax "*(mu*mol~m^{-2}~s^{-1} *"")))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
nd_vcmax_box
ggsave(plot = nd_vcmax_box, "Figures/photo_box_nodip_datvtrad_vcmax.png")

nd_jmax_box <- ggplot(grp_pho_nd_all, aes(x=method, y=mean_jmax, fill = method)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression("Jmax "*(mu*mol~m^{-2}~s^{-1} *"")))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
nd_jmax_box
ggsave(plot = nd_jmax_box, "Figures/photo_box_nodip_datvtrad_jmax.png")




## 1:1 plots DAT vs Trad
# By leaf

#Vcmax
leaf_sub_vcmax <- select(grp_pho_leaf, mean_vcmax, DAT, leaf_unique)
leaf_wide_vcmax <- reshape(leaf_sub_vcmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_vcmax)[2:3]=c("vcmax_DAT", "vcmax_Trad")
cor1 <- round(cor(leaf_wide_vcmax$vcmax_DAT, leaf_wide_vcmax$vcmax_Trad), 3)
photo_leaf_vcmax <- ggplot(data = leaf_wide_vcmax, mapping = aes(x = vcmax_Trad,
                                                               y = vcmax_DAT,
                                                               color = leaf_unique))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1)+
    theme_classic()+
    labs(x = expression("Steady-State Vcmax "*(mu*mol~m^{-2}~s^{-1} *"")),
         y = expression("DAT Vcmax "*(mu*mol~m^{-2}~s^{-1} *"")),
         col = "Unique Leaf")+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=14, family = "serif"),
          axis.title.y=element_text(size=14, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=8, family = "serif", color = "gray10"),
          #legend.text=element_text(size=7, family = "serif"),
          #legend.title=element_text(size=11, family = "serif"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 100)) + 
    scale_y_continuous(limits = c(1, 100)) +
    annotate(geom = "text", label = paste0("r = ", cor1), x = 25, y = 75)
photo_leaf_vcmax
ggsave(plot = photo_leaf_vcmax, "Figures/photo_1to1_datvtrad_vcmax_all.png")



#Jmax
leaf_sub_jmax <- select(grp_pho_leaf, mean_jmax, DAT, leaf_unique)
leaf_wide_jmax <- reshape(leaf_sub_jmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_jmax)[2:3]=c("jmax_DAT", "jmax_Trad")
cor2 <- round(cor(leaf_wide_jmax$jmax_DAT, leaf_wide_jmax$jmax_Trad), 3)
photo_leaf_jmax <- ggplot(data = leaf_wide_jmax, mapping = aes(x = jmax_Trad,
                                                             y = jmax_DAT,
                                                             color = leaf_unique))+
    geom_point()+
    geom_abline(intercept = 0, slope = 1)+
    theme_classic()+
    labs(x = expression("Steady-State Jmax "*(mu*mol~m^{-2}~s^{-1} *"")),
         y = expression("DAT Jmax "*(mu*mol~m^{-2}~s^{-1} *"")), col = "Unique Leaf") +
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=14, family = "serif"),
          axis.title.y=element_text(size=14, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif"),
          axis.text.y=element_text(size=8, family = "serif"),
          #legend.text=element_text(size=7, family = "serif"),
          #legend.title=element_text(size=11, family = "serif"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 130)) + 
    scale_y_continuous(limits = c(1, 130)) +
    annotate(geom = "text", label = paste0("r = ", cor2), x = 30, y = 100)
photo_leaf_jmax
ggsave(plot = photo_leaf_jmax, "Figures/photo_1to1_datvtrad_jmax_all.png")


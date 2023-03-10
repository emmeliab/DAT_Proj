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
params_ecophys <- read.csv(file = paste0(wd, "Results/params_ecophys_no_TPU.csv"), sep = ",", 
                           header = TRUE) %>% 
  filter(method == "dat")
params_photo <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_no_TPU.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000) ## TPU values at 1000 are coded as NA
params_photo_tpu <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_with_TPU.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000)

init_mg <- read.csv(file = paste0(wd, "Results/MG_fixed_aci_fits_230213.csv"), sep = ",",
                      header = TRUE)
params_mg <- init_mg %>%
    select(-c(X)) %>% 
    filter(DAT == "Before_DAT") #filters out traditional



# Comparing fits across photosynthesis, plantecophys, and MG code -----------

# photosynth vs ecophys
# rmse(params_ecophys$Vcmax, params_photo$V_cmax)
# ecovphoto_vcmax <- ggplot(mapping = aes(x = params_ecophys$Vcmax, y = params_photo$vcmax_25, 
#                                         color = params_photo$ID)) +
#     geom_point()+
#     geom_abline(intercept = 0, slope = 1)+
#     theme_classic()+
#     labs(x="Plantecophys Vcmax", y="Photosynthesis Vcmax", col = "Leaf")+
#     theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif")) +
#     scale_x_continuous(limits = c(1, 500)) + 
#     scale_y_continuous(limits = c(1, 500)) +
#     annotate(geom = "text", label = "RMSE", x = 125, y = 50)
# ecovphoto_vcmax
# 
# #ecophys vs MG 
# rmse(params_ecophys$Vcmax, params_mg$Best_Vcmax_25C)
# ecovMG_vcmax <- ggplot(mapping = aes(x = params_ecophys$Vcmax, y = params_mg$Best_Vcmax_25C, 
#                                      color = params_mg$Tree_id))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Plantecophys Vcmax", y="MG Vcmax", col = "Leaf")+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif")) +
#   scale_x_continuous(limits = c(1, 500)) + 
#   scale_y_continuous(limits = c(1, 500))
#  # annotate(geom = "text", label = "RMSE = 16.63", x = 125, y = 50)
# ecovMG_vcmax

#Photosynthesis temp corrected vs MG temp corrected
# cor(params_photo$Best_Vcmax_25C, params_mg$Best_Vcmax_25C)
# photovMG_vcmax <- ggplot(mapping = aes(x = params_photo$Best_Vcmax_25C, y = params_mg$Best_Vcmax_25C,
#                                        color = params_photo$ID))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x = "Photosynthesis Vcmax", y = "MG Vcmax", col = "Leaf")+
#   theme(axis.title.x=element_text(size=14, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif")) +
#   scale_x_continuous(limits = c(1, 90)) + 
#   scale_y_continuous(limits = c(1, 90))+
#   annotate(geom = "text", label = "r = 0.963", x = 65, y = 20)
# photovMG_vcmax
# ggsave("Figures/Photo_MG_vcmax.pdf")


#Jmax
# Ecophys vs MG
# cor(params_ecophys$Jmax, params_mg$Best.Jmax_25C) ## see if we can fix the 800000 and try again
# ecovMG_jmax <- ggplot(mapping = aes(x = params_ecophys$Jmax, y = params_mg$Best.Jmax_25C,
#                                     color = params_mg$Tree_id))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Plantecophys Jmax", y="MG Jmax", col = "Leaf")+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif")) +
#   scale_x_continuous(limits = c(1, 200)) + 
#   scale_y_continuous(limits = c(1, 200)) +
#   annotate(geom = "text", label = "RMSE = 23.27", x = 125, y = 50)
# ecovMG_jmax
# Note: K6706L1 jmax for plantecophys is 800000, and is not included in the graph

# Ecophys vs Photosynthesis
# rmse(params_ecophys$Jmax, params_photo$Best_Jmax_25C)
# ecovphoto_jmax <- ggplot(mapping = aes(x = params_ecophys$Jmax, y = params_photo$J_max, 
#                                        color = params_photo$ID))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Plantecophys Jmax", y="Photosynthesis Jmax", col = "Leaf")+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif")) +
#   scale_x_continuous(limits = c(1, 170)) + 
#   scale_y_continuous(limits = c(1, 170)) +
#   annotate(geom = "text", label = "RMSE = 12.74", x = 125, y = 50)
# ecovphoto_jmax
# Note: K6706L1 jmax for plantecophys is 800000, and is not included in the graph

#Photosynthesis vs MG
# cor(params_photo$Best_Jmax_25C, params_mg$Best.Jmax_25C)
# photovMG_jmax <- ggplot(mapping = aes(x = params_photo$Best_Jmax_25C, y = params_mg$Best.Jmax_25C,
#                                       color = params_photo$ID))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Photosynthesis Jmax", y="MG Jmax", col = "Leaf")+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif")) +
#   scale_x_continuous(limits = c(1, 150)) + 
#   scale_y_continuous(limits = c(1, 150)) +
#   annotate(geom = "text", label = "r = 0.996", x = 125, y = 50)
# photovMG_jmax
# ggsave("Figures/Photo_MG_jmax.pdf")

#TPU
# rmse(params_ecophys$TPU, params_mg$TPU_Best)
# ecovMG_tpu <- ggplot(mapping = aes(x = params_ecophys$TPU, y = params_mg$TPU_Best,
#                                    color = params_mg$Tree_id))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Plantecophys TPU", y="MG TPU", col = "Leaf")+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif")) +
#   scale_x_continuous(limits = c(1, 11)) + 
#   scale_y_continuous(limits = c(1, 11)) +
#   annotate(geom = "text", label = "RMSE = 1.03", x = 7, y = 3)
# ecovMG_tpu


# rmse(params_ecophys$TPU, params_photo$V_TPU)
# ecovphoto_tpu <- ggplot(mapping = aes(x = params_ecophys$TPU, y = params_photo$V_TPU,
#                                       color = params_photo$ID))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Plantecophys TPU", y="Photosynthesis TPU", col = "Leaf")+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif")) +
#   scale_x_continuous(limits = c(1, 10)) + 
#   scale_y_continuous(limits = c(1, 10)) +
#   annotate(geom = "text", label = "RMSE = 0.32", x = 7, y = 3)
# ecovphoto_tpu
# Note: TPU for several of the curves via the photosynthesis package are at 1000 (likely meaning it
# wasn't fit). These are not included in the chart


# cor(params_photo_tpu$V_TPU, params_mg$TPU_Best)
# photovMG_tpu <- ggplot(mapping = aes(x = params_photo_tpu$V_TPU, y = params_mg$TPU_Best, 
#                                      color = params_photo$ID))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Photosynthesis TPU", y="MG TPU", col = "Leaf")+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif")) +
#   scale_x_continuous(limits = c(1, 11)) + 
#   scale_y_continuous(limits = c(1, 11)) +
#   annotate(geom = "text", label = "RMSE = 1.47", x = 7, y = 3)
# photovMG_tpu





# Photosynthesis results visualization ------------------------------------

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
pho_nd_stat <- select(pho_nd_leaf, 'DAT', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'leaf_id', 'leaf_unique')
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

grp_pho_dat <- pho_stat %>%
    filter(method == "DAT") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_vcmax = length(vcmax),
              mean_vcmax = mean(vcmax)) %>% 
    mutate(method = "DAT")

grp_pho_trad <- pho_stat %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_vcmax = length(vcmax),
              mean_vcmax = mean(vcmax)) %>% 
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
              mean_jmax = mean(jmax)) %>% 
    mutate(method = "DAT")

max(grp_pho_jmax_dat$mean_jmax)

grp_pho_jmax_trad <- pho_stat %>%
    filter(method == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarise(n_jmax = length(jmax),
              mean_jmax = mean(jmax)) %>% 
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

nd_vcmax_box <- ggplot(grp_pho_nd_all, aes(x=method, y=mean_vcmax)) +
    geom_boxplot()+
    labs(x="Method", y = "Vcmax")+
    theme_classic()+
    theme(axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"),
          legend.position="none")
nd_vcmax_box
ggsave("Figures/photo_box_nodip_datvtrad_vcmax.pdf")

nd_jmax_box <- ggplot(grp_pho_nd_all, aes(x=method, y=mean_jmax)) +
    geom_boxplot()+
    labs(x="Method", y = "Jmax")+
    theme_classic()+
    theme(axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=15, family = "serif"),
          axis.text.y=element_text(size=15, family = "serif"),
          legend.position="none")
nd_jmax_box
ggsave("Figures/photo_box_nodip_datvtrad_jmax.pdf")

# Visualization of DAT vs Traditional in Photosynthesis package -----

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
ggsave("Figures/photo_box_datvtrad_vcmax.pdf")

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
ggsave("Figures/photo_box_datvtrad_jmax.pdf")

# # Stacked scatters for photosynthesis package
# filt_par_dummy <- mutate(.data = grp_pho_leaf,# makes a dummy variable to plot
#                          dummy = if_else(grp_pho_leaf$DAT == "Before_DAT", 0,1))
# 
# scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_vcmax,
#                                                           color = leaf_unique,
#                                                           label = leaf_unique)) +
#     geom_line() + 
#     geom_point() +
#     theme_classic() +
#     geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
#                     size          = 2.4,
#                     box.padding   = 0.25,
#                     point.padding = 0.25,
#                     segment.size  = 0.2,
#                     segment.linetype = 3,
#                     direction     = "y",
#                     nudge_x = -0.2)+ ## Playing around with this to help visualize
#     scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
#                        expand = expansion(mult=0.3)) +
#     labs(x="Method", y="Vcmax")+
#     theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
#           axis.title.x=element_text(size=18, family = "serif"),
#           axis.title.y=element_text(size=18, family = "serif"),
#           axis.text.x=element_text(size=15, family = "serif"),
#           axis.text.y=element_text(size=15, family = "serif"),
#           legend.text=element_text(size=9, family = "serif"),
#           legend.title=element_text(size=11, family = "serif"),
#           legend.position="none")
# scat_vcmax
# ggsave("Figures/photo_grpscat_datvtrad_vcmax.pdf")
# 
# scat_jmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_jmax,
#                                                          color = leaf_unique,
#                                                          label = leaf_unique)) +
#     geom_line() + 
#     geom_point() +
#     theme_classic() +
#     geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
#                     size          = 2.8,
#                     box.padding   = 0.25,
#                     point.padding = 0.25,
#                     segment.size  = 0.2,
#                     direction     = "y",
#                     nudge_x = -0.2)+ ## Playing around with this to help visualize
#     scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
#                        expand = expansion(mult=0.3)) +
#     labs(x="Method", y="Jmax")+
#     theme(aspect.ratio = 1.5,
#           axis.title.x=element_text(size=18, family = "serif"),
#           axis.title.y=element_text(size=18, family = "serif"),
#           axis.text.x=element_text(size=15, family = "serif"),
#           axis.text.y=element_text(size=15, family = "serif"),
#           legend.text=element_text(size=9, family = "serif"),
#           legend.title=element_text(size=11, family = "serif"),
#           legend.position="none")+
#     guides(color = guide_legend(title = "Leaf Identifier"))
# scat_jmax
# ggsave("Figures/photo_grpscat_datvtrad_jmax.pdf")

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
photo_leaf_vcmax
ggsave("Figures/photo_1to1_datvtrad_vcmax.pdf")

#Jmax
leaf_sub_jmax <- select(grp_pho_leaf, mean_jmax, DAT, leaf_unique)
leaf_wide_jmax <- reshape(leaf_sub_jmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_jmax)[2:3]=c("jmax_DAT", "jmax_Trad")
cor1 <- round(cor(leaf_wide_jmax$jmax_DAT, leaf_wide_jmax$jmax_Trad), 3)
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
    scale_y_continuous(limits = c(1, 130)) +
    annotate(geom = "text", label = paste0("r = ", cor1), x = 100, y = 30)
photo_leaf_jmax
ggsave("Figures/photo_1to1_datvtrad_jmax.pdf")

# MG data processing ------------------------------------------------

## Separate the concatenated tree ID column
# mg_split <- unlist(str_split(init_mg$Tree_id, "_", n=2))
# mg_sub <- subset(mg_split, mg_split != "Before_DAT" & mg_split != "Traditional")
# init_mg$leaf_id <- mg_sub
# mg_complete <- subset(init_mg, select = -c(Tree_id))
# 
# leaf_split <- unlist(str_split(mg_complete$leaf_id, "L", n=2))
# leaf_sub <- subset(leaf_split, leaf_split != "1" & leaf_split != "2" & leaf_split != "1-1"
#                    & leaf_split != 3 & leaf_split != "6" & leaf_split != "8" & leaf_split != "2-2")
# mg_complete$tree_id <- leaf_sub
# 
# #Add in the relative canopy position
# codes <- read.csv("~/Documents/GitHub/DAT_Proj/Inputs/unique_ids.csv") # need to keep this in for leaf_id
# canopy_pos <- read.csv("~/Documents/GitHub/DAT_Proj/Inputs/rel_canopy_pos.csv")
# names(codes)[1] ="leaf_id"
# codes_and_can <- left_join(codes, canopy_pos, by = "k67.id")
# names(codes_and_can)[3]="code4let"
# codes_and_can <- subset(codes_and_can, select = -code.y)
# 
# mg_complete <- left_join(mg_complete, codes_and_can, by = "leaf_id") %>% 
#   select(-c(k67.id,code4let))
# 
# 
# #This is the data with the back correction filtered out
# mg_leaf <- mg_complete %>%
#   mutate(leaf_unique = substring(leaf_id, 1, 7))
# 
# #group data for further analysis
# grp_leaf <- mg_leaf %>% group_by(DAT, leaf_unique) %>%
#   summarise(mean_vcmax=mean(Best_Vcmax_25C),
#             mean_jmax= mean(Best.Jmax_25C),
#             #mean_tpumax= mean(TPU_Best),
#             tree_id = tree_id,
#             leaf_unique=leaf_unique,
#             leaf_id=leaf_id,
#             rel_can_pos=rel_can_pos) %>%
#   as.data.frame()
# 
# summary(grp_leaf)
# 
# grp_tree <- mg_leaf %>% group_by(DAT,tree_id) %>% 
#   summarise(mean_vcmax=mean(Best_Vcmax_25C),
#             mean_jmax= mean(Best.Jmax_25C),
#             #mean_tpumax= mean(TPU_Best),
#             tree_id = tree_id,
#             leaf_unique=leaf_unique,
#             leaf_id=leaf_id,
#             rel_can_pos=rel_can_pos,
#             .groups = 'drop') %>%
#   as.data.frame()
# 
# 
# ### MG Visualizations -------------------------------
# #### Visualization of DAT vs Traditional
# 
# ## Boxplots
# lab_DATTrad <- c('DAT', 'Traditional')
# b4 <- ggplot(mg_leaf, aes(x=DAT, y=Best_Vcmax_25C)) +
#   geom_boxplot()+
#   labs(x="Method", y = "Vcmax")+
#   theme_classic()+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.position="none")+
#   scale_x_discrete(labels=lab_DATTrad)
# b4
# 
# b5 <- ggplot(mg_leaf, aes(x=DAT, y=Best.Jmax_25C)) +
#   geom_boxplot()+
#   labs(x="Method", y = "Jmax")+
#   theme_classic()+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.position="none")+
#   scale_x_discrete(labels=lab_DATTrad)
# b5
# 
# b6 <- ggplot(mg_leaf, aes(x=DAT, y=TPU_Best)) +
#   geom_boxplot()+
#   labs(x="Method", y = "TPU")+
#   theme_classic()+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.position="none")+
#   scale_x_discrete(labels=lab_DATTrad)
# b6
# 



# ## Stacked Scatters by leaf
# filt_par_dummy <- mutate(.data = grp_leaf,# makes a dummy variable to plot
#                          dummy = if_else(grp_leaf$DAT == "Before_DAT", 0,1))
# 
# scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_vcmax,
#                                                           color = leaf_unique,
#                                                           label = tree_id)) +
#   geom_line() + 
#   geom_point() +
#   theme_classic() +
#   geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
#                   size          = 2.4,
#                   box.padding   = 0.25,
#                   point.padding = 0.25,
#                   segment.size  = 0.2,
#                   segment.linetype = 3,
#                   direction     = "y",
#                   nudge_x = -0.2)+ ## Playing around with this to help visualize
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
# scat_vcmax
# 
# scat_jmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_jmax,
#                                                          color = leaf_unique,
#                                                          label = tree_id)) +
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
#   labs(x="Method", y="Jmax")+
#   theme(aspect.ratio = 1.5,
#         axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=9, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"),
#         legend.position="none")+
#   guides(color = guide_legend(title = "Leaf Identifier"))
# scat_jmax
# 
# # scat_tpu <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_tpumax,
# #                                                         color = leaf_unique,
# #                                                         label = tree_id)) +
# #   geom_line() + 
# #   geom_point() +
# #   theme_classic() +
# #   geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
# #                   size          = 2.8,
# #                   box.padding   = 0.25,
# #                   point.padding = 0.25,
# #                   segment.size  = 0.2,
# #                   direction     = "y",
# #                   nudge_x = -0.2)+ ## Playing around with this to help visualize
# #   scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
# #                      expand = expansion(mult=0.3)) +
# #   labs(x="Method", y="TPU")+
# #   theme(aspect.ratio = 1.5,
# #         axis.title.x=element_text(size=18, family = "serif"),
# #         axis.title.y=element_text(size=18, family = "serif"),
# #         axis.text.x=element_text(size=15, family = "serif"),
# #         axis.text.y=element_text(size=15, family = "serif"),
# #         legend.text=element_text(size=9, family = "serif"),
# #         legend.title=element_text(size=11, family = "serif"),
# #         legend.position ="none")+
# #   guides(color = guide_legend(title = "Leaf Identifier"))
# # scat_tpu



## 1:1 plots DAT vs Trad
# By leaf

#Vcmax
# leaf_sub_vcmax <- select(grp_leaf, mean_vcmax, DAT, leaf_unique, tree_id)
# leaf_wide_vcmax <- reshape(leaf_sub_vcmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
# names(leaf_wide_vcmax)[2:4]=c("vcmax_DAT", "tree_id", "vcmax_Trad")
# leaf_wide_vcmax <- subset(leaf_wide_vcmax, select = -tree_id.Traditional)
# mng_leaf_vcmax <- ggplot(data = leaf_wide_vcmax, mapping = aes(x = vcmax_Trad,
#                                                                y = vcmax_DAT,
#                                                                color = leaf_unique))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Traditional Vcmax", y="DAT Vcmax", col = "Unique Leaf")+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"))+
#     scale_x_continuous(limits = c(1, 80)) + 
#     scale_y_continuous(limits = c(1, 80))
# mng_leaf_vcmax
# 
# #Jmax
# leaf_sub_jmax <- select(grp_leaf, mean_jmax, DAT, leaf_unique, tree_id)
# leaf_wide_jmax <- reshape(leaf_sub_jmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
# names(leaf_wide_jmax)[2:4]=c("jmax_DAT", "tree_id", "jmax_Trad")
# leaf_wide_jmax <- subset(leaf_wide_jmax, select = -tree_id.Traditional)
# mng_leaf_jmax <- ggplot(data = leaf_wide_jmax, mapping = aes(x = jmax_Trad,
#                                                              y = jmax_DAT,
#                                                              color = leaf_unique))+
#   geom_point()+
#   geom_abline(intercept = 0, slope = 1)+
#   theme_classic()+
#   labs(x="Traditional Jmax", y="DAT Jmax", col = "Unique Leaf")+
#   theme(axis.title.x=element_text(size=18, family = "serif"),
#         axis.title.y=element_text(size=18, family = "serif"),
#         axis.text.x=element_text(size=15, family = "serif"),
#         axis.text.y=element_text(size=15, family = "serif"),
#         legend.text=element_text(size=7, family = "serif"),
#         legend.title=element_text(size=11, family = "serif"))+
#     scale_x_continuous(limits = c(1, 120)) + 
#     scale_y_continuous(limits = c(1, 120))
# mng_leaf_jmax

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



# MG stat analysis for vcmax --------------------------------

# #Subset variables of interest
# leaf_stat <- select(mg_leaf, 'DAT', 'Best_Vcmax_25C', 'Best.Jmax_25C', 'leaf_id', 'tree_id', 'leaf_unique', 'rel_can_pos')
# leaf_stat$DAT[leaf_stat$DAT == "Before_DAT"] <- "DAT"
# leaf_stat <- rename(leaf_stat,
#                     method = DAT,
#                     vcmax = Best_Vcmax_25C,
#                     jmax = Best.Jmax_25C,
#                     #tpu = TPU_Best,
#                     )
# #Describe factor levels. 0 is traditional, 1 is DAT
# leaf_stat$method <- factor(leaf_stat$method)
# 
# #displays grouped summary
# leaf_summary <- leaf_stat %>%
#   group_by(method) %>%
#   summarise(n_vcmax = length(vcmax),
#             mean_vcmax = round(mean(vcmax),3),
#             sd_vcmax = round(sd(vcmax),3),
#             se_vcmax = sd(vcmax)/sqrt(n())) %>% 
#   mutate(low_ci_vcmax = mean_vcmax - qt(1 - (0.05 / 2), n_vcmax - 1) * se_vcmax,
#          up_ci_vcmax = mean_vcmax + qt(1 - (0.05 / 2), n_vcmax - 1) * se_vcmax)
# print(leaf_summary)
# 
# #Visualize Vcmax by Method
# ci.mean(vcmax ~ method, data = leaf_stat)
# ci1<-ci.mean(vcmax~method, data=leaf_stat)
# plot(ci1,title.labels="Method")
# 
# #Histogram to visualize
# leaf_hist<-ggplot(leaf_stat, aes(x=vcmax)) + 
#   geom_histogram(color="black", fill="white", bins = 8)+
#   geom_vline(aes(xintercept=mean(vcmax)),
#              color="red", linetype="dashed", linewidth=0.5)+
#   theme_classic()
# leaf_hist
# #data are positively skewed
# skewness(leaf_stat$vcmax)
# #value close to 1, should be okay
# kurtosis(leaf_stat$vcmax)
# #It's a bit high
# 
# #Test the assumption of equal variances for each group for t-test with Levene's
# leveneTest(vcmax ~ method, data = leaf_stat)
# #Non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances
# 
# # Shapiro-Wilk normality test for Vcmax for the one-sample t-test
# with(leaf_stat, shapiro.test(vcmax))
# # Shapiro-Wilk normality test for the DAT measurement methodology
# with(leaf_stat, shapiro.test(vcmax[method == "DAT"]))
# # Shapiro-Wilk normality test for the Traditional measurement methodology
# with(leaf_stat, shapiro.test(vcmax[method == "Traditional"])) 
# #All are significant. Rejects null hypothesis that these data are not normally distributed.
# # Therefore the sample varies from the normal distribution.
# 
# grp_dat <- leaf_stat %>%
#   filter(method == "DAT") %>% 
#   group_by(leaf_unique) %>% 
#   summarise(n_vcmax = length(vcmax),
#             mean_vcmax = mean(vcmax),
#             sd_vcmax = sd(vcmax),
#             se_vcmax = sd(vcmax)/sqrt(n())) %>% 
#   mutate(method = "DAT")
#   
# grp_trad <- leaf_stat %>%
#   filter(method == "Traditional") %>% 
#   group_by(leaf_unique) %>% 
#   summarise(n_vcmax = length(vcmax),
#             mean_vcmax = mean(vcmax),
#             sd_vcmax = sd(vcmax),
#             se_vcmax = sd(vcmax)/sqrt(n())) %>% 
#   mutate(method = "Traditional")
# 
# grp_all <- rbind(grp_dat, grp_trad)
# 
# with(grp_all, shapiro.test(mean_vcmax))
# # Shapiro-Wilk normality test for the DAT measurement methodology
# with(grp_all, shapiro.test(mean_vcmax[method == "DAT"]))
# # Shapiro-Wilk normality test for the Traditional measurement methodology
# with(grp_all, shapiro.test(mean_vcmax[method == "Traditional"])) 
# #All are significant. Rejects null hypothesis that these data are not normally distributed.
# # Therefore the sample varies from the normal distribution.
# 
# wilcox.test(mean_vcmax ~ method, data = grp_all, paired = TRUE)
# #This result is significant
# 
# #Effect size for the independent sample t-test:
# cohen.ES(test = "t", size = "large") # To remind oneself
# cohen.d(mean_vcmax ~ method | Subject(leaf_unique), data=grp_all, paired = TRUE)
# #small effect size
# 
# #### Power Analysis
# 
# d <- cohen.d(mean_vcmax ~ method | Subject(leaf_unique), data=grp_all, paired = TRUE)
# d[["estimate"]]
# 
# #How many samples to achieve a certain power?
# pwr1 <- pwr.t.test(n = , d = d[["estimate"]], power = 0.7, sig.level = 0.05, type = "paired", alternative = "two.sided")
# plot(pwr1)
# 
# #What was the power of our study?
# pwr.t.test(n = 28, d = d[["estimate"]], power = , sig.level = 0.05, type = "paired", alternative = "two.sided")
# plot(pwr2)
# 
# 
# # MG stat analysis for jmax ------------------------------
# 
# leaf_jmax_summary <- leaf_stat %>%
#     group_by(method) %>%
#     summarise(n_jmax = length(jmax),
#               mean_jmax = round(mean(jmax),3),
#               sd_jmax = round(sd(jmax),3),
#               se_jmax = sd(jmax)/sqrt(n())) %>% 
#     mutate(low_ci_jmax = mean_jmax - qt(1 - (0.05 / 2), n_jmax - 1) * se_jmax,
#            up_ci_jmax = mean_jmax + qt(1 - (0.05 / 2), n_jmax - 1) * se_jmax)
# print(leaf_jmax_summary)
# 
# #Visualize Vcmax by Method
# ci.mean(jmax ~ method, data = leaf_stat)
# ci1<-ci.mean(jmax~method, data=leaf_stat)
# plot(ci1,title.labels="Method")
# 
# #Histogram to visualize
# leaf_jmax_hist<-ggplot(leaf_stat, aes(x=jmax)) + 
#     geom_histogram(color="black", fill="white", bins = 8)+
#     geom_vline(aes(xintercept=mean(vcmax)),
#                color="red", linetype="dashed", linewidth=0.5)+
#     theme_classic()
# leaf_jmax_hist
# #data are positively skewed
# skewness(leaf_stat$jmax)
# #value close to 1, should be okay
# kurtosis(leaf_stat$jmax)
# #It's a bit high
# 
# #Test the assumption of equal variances for each group for t-test with Levene's
# leveneTest(jmax ~ method, data = leaf_stat)
# #Non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances
# 
# # Shapiro-Wilk normality test for Vcmax for the one-sample t-test
# with(leaf_stat, shapiro.test(jmax))
# # Shapiro-Wilk normality test for the DAT measurement methodology
# with(leaf_stat, shapiro.test(jmax[method == "DAT"]))
# # Shapiro-Wilk normality test for the Traditional measurement methodology
# with(leaf_stat, shapiro.test(jmax[method == "Traditional"])) 
# #2/3 significant. Rejects null hypothesis that these data are not normally distributed.
# # Therefore the sample varies from the normal distribution.
# 
# grp_jmax_dat <- leaf_stat %>%
#     filter(method == "DAT") %>% 
#     group_by(leaf_unique) %>% 
#     summarise(n_jmax = length(jmax),
#               mean_jmax = mean(jmax),
#               sd_jmax = sd(jmax),
#               se_jmax = sd(jmax)/sqrt(n())) %>% 
#     mutate(method = "DAT")
# 
# grp_jmax_trad <- leaf_stat %>%
#     filter(method == "Traditional") %>% 
#     group_by(leaf_unique) %>% 
#     summarise(n_jmax = length(jmax),
#               mean_jmax = mean(jmax),
#               sd_jmax = sd(jmax),
#               se_jmax = sd(jmax)/sqrt(n())) %>% 
#     mutate(method = "Traditional")
# 
# grp_jmax_all <- rbind(grp_jmax_dat, grp_jmax_trad)
# 
# with(grp_jmax_all, shapiro.test(mean_jmax))
# # Shapiro-Wilk normality test for the DAT measurement methodology
# with(grp_jmax_all, shapiro.test(mean_jmax[method == "DAT"]))
# # Shapiro-Wilk normality test for the Traditional measurement methodology
# with(grp_jmax_all, shapiro.test(mean_jmax[method == "Traditional"])) 
# #2/3 significant. Rejects null hypothesis that these data are not normally distributed.
# # Therefore the sample varies from the normal distribution.
# 
# wilcox.test(mean_jmax ~ method, data = grp_jmax_all, paired = TRUE)
# #This result is significant
# 
# #Effect size for the independent sample t-test:
# cohen.ES(test = "t", size = "large") # To remind oneself
# cohen.d(mean_jmax ~ method | Subject(leaf_unique), data=grp_jmax_all, paired = TRUE)
# #small effect size
# 
# #### Power Analysis
# 
# d <- cohen.d(mean_jmax ~ method | Subject(leaf_unique), data=grp_jmax_all, paired = TRUE)
# d[["estimate"]]
# 
# #How many samples to achieve a certain power?
# pwr1 <- pwr.t.test(n = , d = d[["estimate"]], power = 0.7, sig.level = 0.05, type = "paired", alternative = "two.sided")
# plot(pwr1)
# 
# #What was the power of our study?
# pwr.t.test(n = 28, d = d[["estimate"]], power = , sig.level = 0.05, type = "paired", alternative = "two.sided")
# plot(pwr2)


### comparing photosynthesis and MG, barplots and ANOVA? DECIDE WHETHER TO KEEP THIS. IF SO NEED TO UNCOMMENT SOME OF THE ABOVE -----------------------

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
ggsave("Figures/mg_photo_method_vcmax.pdf")

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
ggsave("Figures/mg_photo_method_jmax.pdf")

#Analysis of variance
library(AICcmodavg)

#vcmax
lm_method <- lm(vcmax ~ method, data = all_results)
lm_fit <- lm(vcmax ~ fit_type, data = all_results)
lm_both <- lm(vcmax ~ method + fit_type, data = all_results)
lm_null <- lm(vcmax ~ 1, data = all_results)

mod_names <- c("Method", "Fit Type", "Method + Fit Type", "Intercept Only")

mod_table <- aictab(list(lm_method, lm_fit, lm_both, lm_null), modnames = mod_names)
mod_table #This suggests the best model is the intercept model! Fit type has minor support
summary(lm_null)
summary(lm_fit) # this explains a negligible amount of the variance in vcmax
summary(lm_method) # this explains a negligible amount of the variance in vcmax

#kruskal.test(vcmax ~ method, data = all_results) #non-parametric ANOVA
#chi-squared = 0.00042, df = 1, p-value = 0.9836. Not significant.

#RUN THIS!!!!!!!!!!
wilcox.test(vcmax ~ method, data = all_results, conf.int = TRUE)
wilcox.test(vcmax ~ fit_type, data = all_results, conf.int = TRUE)


car::vif(lm_null)
check_heteroscedasticity(lm_null)


#jmax
lm2_method <- lm(jmax ~ method, data = all_results)
lm2_fit <- lm(jmax ~ fit_type, data = all_results)
lm2_both <- lm(jmax ~ method + fit_type, data = all_results)
lm2_null <- lm(jmax ~ 1, data = all_results)

mod2_names <- c("Method", "Fit Type", "Method + Fit Type", "Intercept Only")

mod2_table <- aictab(list(lm2_method, lm2_fit, lm2_both, lm2_null), modnames = mod_names)
mod2_table # The best model is the method, followed by the method and fit type models
summary(lm2_method) # We're still only explaining 3% of the variation here...
summary(lm2_both) # only 2% here

kruskal.test(jmax ~ method, data = all_results) #non-parametric ANOVA
#chi-squared = 2.72, df = 1, p-value =0.0991. Not significant.

#Maybe should do wilcox.test(jmax ~ method, data = all_results, conf.int = TRUE)

wilcox.test(jmax ~ method, data = all_results, conf.int = TRUE)
wilcox.test(jmax ~ fit_type, data = all_results, conf.int = TRUE)


#Testing outliers
leveragePlots(lm2_method) #70, 106, 34, 35
outlierTest(lm2_method) #70

cutoff <- 4/122
plot(lm2_method, which=4, cook.levels=cutoff)
abline(h = 4/122,  lty = 4, col = "red")
#We see there's 6 points of influence above our Cook's threshold

# Saving cook's distance values as a separate object
cd <- cooks.distance(lm2_method)
cd.idx <- cd[cd > cutoff] #creating an index with high CD values
view(cd.idx)

#Now let's filter out Case #70, which were both present in all 3 outlier tests
all_res_filt <- all_results[-70,]

#create a new model with the new dataset
new_lm2_method <- lm(jmax ~ method, data = all_res_filt)

summary(lm2_method) #Now compare original model
summary(new_lm2_method) #Now compare original model to the new model
#not a huge difference.




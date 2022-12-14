### For statistics and analysis

library(tidyverse)
library(ggpubr)

wd <- "/Users/charlessouthwick/Documents/GitHub/DAT_Proj"

#------------------------------------------ 
#Change all "no_back" to "original" for clarity


#Change the directory where you want to save your data
setwd(paste0(wd, "/Figures"))


curves_df <- read.csv("~/Documents/GitHub/DAT_Proj/Results/curve_fitting_out.csv")

#All this is to separate the concatenated tree_name column ------------------------
curve_split <- unlist(str_split(curves_df$Tree_id, "_", n=2))
curve_sub <- subset(curve_split, curve_split != "Before_DAT" & curve_split != "Traditional")
curves_df$leaf_id <- curve_sub
curves_df_fixed <- subset(curves_df, select = -c(Tree_id))

leaf_split <- unlist(str_split(curves_df_fixed$leaf_id, "L", n=2))
leaf_sub <- subset(leaf_split, leaf_split != "1" & leaf_split != "2" & leaf_split != "1-1"
                   & leaf_split != 3 & leaf_split != "6" & leaf_split != "8" & leaf_split != "2-2")
curves_df_fixed$tree_id <- leaf_sub
curves_final <- curves_df_fixed

#This is the data without the back-correction filter
curves_no_back <- subset(curves_final, back_filt == "no_back")

#group means for DAT and Trad results -------------------------------
leaf2 <- mutate(curves_no_back, leaf_unique = substring(curves_no_back$leaf_id, 1, 7))

dat_all <- subset(curves_final, DAT=="Before_DAT")

b1 <- ggplot(dat_all, aes(x=back_filt, y=vcmax_Best_Model)) +
  geom_boxplot()+
  labs(x="Dataset", y = "Vcmax")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(legend.position="none")
b1

b2 <- ggplot(dat_all, aes(x=back_filt, y=Jmax_Best)) +
  geom_boxplot()+
  labs(x="Dataset", y = "Jmax")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(legend.position="none")
b2

b3 <- ggplot(dat_all, aes(x=back_filt, y=TPU_Best)) +
  geom_boxplot()+
  labs(x="Dataset", y = "TPU")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(legend.position="none")
b3

#Plotting DAT vs Trad on the original, "no_back" data -----------------------------
b4 <- ggplot(leaf2, aes(x=DAT, y=vcmax_Best_Model)) +
  geom_boxplot()+
  labs(x="Method", y = "Vcmax")+
  theme_classic()+
  theme(legend.position="none")
b4

b5 <- ggplot(leaf2, aes(x=DAT, y=Jmax_Best)) +
  geom_boxplot()+
  labs(x="Method", y = "Jmax")+
  theme_classic()+
  theme(legend.position="none")
b5

b6 <- ggplot(leaf2, aes(x=DAT, y=TPU_Best)) +
  geom_boxplot()+
  labs(x="Method", y = "TPU")+
  theme_classic()+
  theme(legend.position="none")
b6

#Just a scatter to understand the spread a bit
g1 <- ggplot(leaf2, aes(x = tree_id, y = vcmax_Best_Model)) +
  geom_point() +
  xlab("Tree")+
  ylab("Vcmax") +
  theme_classic() +
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"))
g1

# Linked point scatter, by leaf_unique ----------------------------------------
filt_par_dummy <- mutate(.data = grp_leaf,# makes a dummy variable to plot
                         dummy = if_else(grp_leaf$DAT == "Before_DAT", 0,1))

scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_vcmax,
                                                          color = leaf_unique)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional")) +
  labs(x="Method", y="Vcmax")+
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")+
  guides(color = guide_legend(title = "Leaf Identifier"))
scat_vcmax

scat_jmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_jmax,
                                                         color = leaf_unique)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional")) +
  labs(x="Method", y="Jmax")+
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")+
  guides(color = guide_legend(title = "Leaf Identifier"))
scat_jmax

scat_tpu <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_tpumax,
                                                          color = leaf_unique)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional")) +
  labs(x="Method", y="TPU")+
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))+
  guides(color = guide_legend(title = "Leaf Identifier"))
scat_tpu

# The "Maquelle Special", aka 1-to-1

## Maquelle will send code

#Testing for normality, first by the leaf ID ------------------------------------
#Vcmax
leafnorm_vcmax <- leaf2 %>%
  select(leaf_id,vcmax_Best_Model,DAT)%>%
  group_by(leaf_id) %>%
  summarise(dif_vcmax=diff(vcmax_Best_Model)) %>%
  ungroup()

shapiro.test(leafnorm_vcmax$dif_vcmax) #This test is super conservative
ks.test(leafnorm_vcmax$dif_vcmax, 'pnorm')
#If in the output p-value > 0.05, it implies that the distribution of the data
#are not significantly different from normal distribution.
#In other words, we can assume the normality.
hist(leafnorm_vcmax$dif_vcmax)
qqnorm(leafnorm_vcmax$dif_vcmax)
qqline(leafnorm_vcmax$dif_vcmax)

#Jmax
leafnorm_jmax <- leaf2 %>%
  select(leaf_id,Jmax_Best,DAT)%>%
  group_by(leaf_id) %>%
  summarise(dif_jmax=diff(Jmax_Best)) %>%
  ungroup()

shapiro.test(leafnorm_jmax$dif_jmax) #This test is super conservative
ks.test(leafnorm_jmax$dif_jmax, 'pnorm')
#If in the output p-value > 0.05, it implies that the distribution of the data
#are not significantly different from normal distribution.
#In other words, we can assume the normality.
hist(leafnorm_jmax$dif_jmax)
qqnorm(leafnorm_jmax$dif_jmax)
qqline(leafnorm_jmax$dif_jmax)

#TPU
leafnorm_tpu <- leaf2 %>%
  select(leaf_id,TPU_Best,DAT)%>%
  group_by(leaf_id) %>%
  summarise(dif_tpu=diff(TPU_Best)) %>%
  ungroup()

shapiro.test(leafnorm_tpu$dif_tpu) #This test is super conservative
ks.test(leafnorm_tpu$dif_tpu, 'pnorm')
#If in the output p-value > 0.05, it implies that the distribution of the data
#are not significantly different from normal distribution.
#In other words, we can assume the normality.
hist(leafnorm_tpu$dif_tpu)
qqnorm(leafnorm_tpu$dif_tpu)
qqline(leafnorm_tpu$dif_tpu)


#Try log transform
#Vcmax
leafnorm_vcmax$logdif_vcmax <- log(abs(leafnorm_vcmax$dif_vcmax))
shapiro.test(leafnorm_vcmax$logdif_vcmax) #This test is super conservative
ks.test(leafnorm_vcmax$logdif_vcmax, 'pnorm')
#If in the output, the p-value > 0.05 that imply that the distribution of the data
#are not significantly different from normal distribution.
#In other words, we can assume the normality.
hist(leafnorm_vcmax$logdif_vcmax)
qqnorm(leafnorm_vcmax$logdif_vcmax)
qqline(leafnorm_vcmax$logdif_vcmax)

#Jmax
leafnorm_jmax$logdif_jmax <- log(abs(leafnorm_jmax$dif_jmax))
shapiro.test(leafnorm_jmax$logdif_jmax) #This test is super conservative
ks.test(leafnorm_jmax$logdif_jmax, 'pnorm')
#If in the output, the p-value > 0.05 that imply that the distribution of the data
#are not significantly different from normal distribution.
#In other words, we can assume the normality.
hist(leafnorm_jmax$logdif_jmax)
qqnorm(leafnorm_jmax$logdif_jmax)
qqline(leafnorm_jmax$logdif_jmax)

#TPU
leafnorm_tpu$logdif_tpu <- log(abs(leafnorm_tpu$dif_tpu))
shapiro.test(leafnorm_tpu$logdif_tpu) #This test is super conservative
ks.test(leafnorm_tpu$logdif_tpu, 'pnorm')
#If in the output, the p-value > 0.05 that imply that the distribution of the data
#are not significantly different from normal distribution.
#In other words, we can assume the normality.
hist(leafnorm_tpu$logdif_tpu)
qqnorm(leafnorm_tpu$logdif_tpu)
qqline(leafnorm_tpu$logdif_tpu)

#Testing for normality, now by tree ID ------------------------------------
#Vcmax
treenorm_vcmax <- leaf2 %>%
  select(tree_id,vcmax_Best_Model,DAT)%>%
  group_by(tree_id) %>%
  summarise(dif_vcmax=diff(vcmax_Best_Model)) %>%
  ungroup()

shapiro.test(treenorm_vcmax$dif_vcmax) #This test is super conservative
ks.test(treenorm_vcmax$dif_vcmax, 'pnorm')
#If in the output p-value > 0.05, it implies that the distribution of the data
#are not significantly different from normal distribution.
#In other words, we can assume the normality.
hist(treenorm_vcmax$dif_vcmax)
qqnorm(treenorm_vcmax$dif_vcmax)
qqline(treenorm_vcmax$dif_vcmax)

#Jmax
treenorm_jmax <- leaf2 %>%
  select(tree_id,Jmax_Best,DAT)%>%
  group_by(tree_id) %>%
  summarise(dif_jmax=diff(Jmax_Best)) %>%
  ungroup()

shapiro.test(treenorm_jmax$dif_jmax) #This test is super conservative
ks.test(treenorm_jmax$dif_jmax, 'pnorm')
#If in the output p-value > 0.05, it implies that the distribution of the data
#are not significantly different from normal distribution.
#In other words, we can assume the normality.
hist(treenorm_jmax$dif_jmax)
qqnorm(treenorm_jmax$dif_jmax)
qqline(treenorm_jmax$dif_jmax)

#TPU
treenorm_tpu <- leaf2 %>%
  select(tree_id,TPU_Best,DAT)%>%
  group_by(tree_id) %>%
  summarise(dif_tpu=diff(TPU_Best)) %>%
  ungroup()

shapiro.test(treenorm_tpu$dif_tpu) #This test is super conservative
ks.test(treenorm_tpu$dif_tpu, 'pnorm')
#If in the output p-value > 0.05, it implies that the distribution of the data
#are not significantly different from normal distribution.
#In other words, we can assume the normality.
hist(treenorm_tpu$dif_tpu)
qqnorm(treenorm_tpu$dif_tpu)
qqline(treenorm_tpu$dif_tpu)

#Try log transform
#Vcmax
# treenorm_vcmax$logdif_vcmax <- log(abs(treenorm_vcmax$dif_vcmax))
# shapiro.test(treenorm_vcmax$logdif_vcmax) #This test is super conservative
# ks.test(treenorm_vcmax$logdif_vcmax, 'pnorm')
# #If in the output, the p-value > 0.05 that imply that the distribution of the data
# #are not significantly different from normal distribution.
# #In other words, we can assume the normality.
# hist(treenorm_vcmax$logdif_vcmax)
# qqnorm(treenorm_vcmax$logdif_vcmax)
# qqline(treenorm_vcmax$logdif_vcmax)

#Jmax
# treenorm_jmax$logdif_jmax <- log(abs(treenorm_jmax$dif_jmax))
# shapiro.test(treenorm_jmax$logdif_jmax) #This test is super conservative
# ks.test(treenorm_jmax$logdif_jmax, 'pnorm')
# #If in the output, the p-value > 0.05 that imply that the distribution of the data
# #are not significantly different from normal distribution.
# #In other words, we can assume the normality.
# hist(treenorm_jmax$logdif_jmax)
# qqnorm(treenorm_jmax$logdif_jmax)
# qqline(treenorm_jmax$logdif_jmax)

#TPU
treenorm_tpu$logdif_tpu <- log(abs(treenorm_tpu$dif_tpu))
shapiro.test(treenorm_tpu$logdif_tpu) #This test is super conservative
ks.test(treenorm_tpu$logdif_tpu, 'pnorm')
#If in the output, the p-value > 0.05 that imply that the distribution of the data
#are not significantly different from normal distribution.
#In other words, we can assume the normality.
hist(treenorm_tpu$logdif_tpu)
qqnorm(treenorm_tpu$logdif_tpu)
qqline(treenorm_tpu$logdif_tpu)


#T-testing ----------------------------------------------------------------

# Group data for t-tests
grp_leaf <- leaf2 %>% group_by(DAT, leaf_unique) %>%
  summarise(mean_vcmax=mean(vcmax_Best_Model),
            mean_jmax= mean(Jmax_Best),
            mean_tpumax= mean(TPU_Best)) %>%
  as.data.frame()

grp_tree <- leaf2 %>% group_by(DAT,tree_id) %>% 
  summarise(mean_vcmax=mean(vcmax_Best_Model),
            sd_vcmax = sd(vcmax_Best_Model),
            mean_jmax= mean(Jmax_Best),
            sd_jmax= sd(Jmax_Best),
            mean_tpumax= mean(TPU_Best),
            sd_tpumax= sd(TPU_Best),
            .groups = 'drop') %>%
  as.data.frame()

#DAT vs Trad t-tests -- tree level -- Will need to re-run with normal transformations!
dat_tree_df <- subset(grp_tree, DAT=="Before_DAT")
trad_tree_df <- subset(grp_tree, DAT=="Traditional")

res1<-t.test(dat_tree_df$mean_vcmax, trad_tree_df$mean_vcmax, paired=TRUE)
res1

res2<-t.test(dat_tree_df$mean_jmax, trad_tree_df$mean_jmax, paired=TRUE)
res2

res3<-t.test(dat_tree_df$mean_tpumax, trad_tree_df$mean_tpumax, paired=TRUE)
res3

#DAT vs Trad t-tests -- leaf level -- Will need to re-reun with normal transformations!
dat_leaf_df <- subset(grp_leaf, DAT=="Before_DAT")
trad_leaf_df <- subset(grp_leaf, DAT=="Traditional")

res4<-t.test(dat_leaf_df$mean_vcmax, trad_leaf_df$mean_vcmax, paired=TRUE)
res4

res5<-t.test(dat_leaf_df$mean_jmax, trad_leaf_df$mean_jmax, paired=TRUE)
res5

res6<-t.test(dat_leaf_df$mean_tpumax, trad_leaf_df$mean_tpumax, paired=TRUE)
res6

#Testing whether the back-corrected are different from the non-back-corrected -----------------
#T-tests
dat_noback <- subset(dat_all, back_filt=="no_back")
dat_filt <- subset(dat_all, back_filt=="back_filtered")
res7<-t.test(dat_noback$vcmax_Best_Model, dat_filt$vcmax_Best_Model, paired=TRUE) #Need same number
res7 #These are not significantly different!

res8<-t.test(dat_noback$Jmax_Best, dat_filt$Jmax_Best, paired=TRUE) #Need same number
res8 #These are not significantly different!

res9<-t.test(dat_noback$TPU_Best, dat_filt$TPU_Best, paired=TRUE) #Need same number
res9

## Leaf/Leaflet position? ------------------------------------------------


## Relative canopy height? ------------------------------------------------

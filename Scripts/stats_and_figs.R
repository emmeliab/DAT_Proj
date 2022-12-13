### For statistics and analysis

library(tidyverse)
library(ggpubr)

wd <- "/Users/charlessouthwick/Documents/GitHub/DAT_Proj"

#------------------------------------------ 
#Change all "no_back" to "original" for clarity


#Change the directory where you want to save your data
setwd(paste0(wd, "/Figures"))


curves_df <- read.csv("~/Documents/GitHub/DAT_Proj/Results/curve_fitting_out.csv")

#All this is to separate the concatenated tree_name column
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

#group means for DAT and Trad results
leaf2 <- mutate(curves_no_back, leaf_unique = substring(curves_no_back$leaf_id, 1, 7))

grp_leaf <- leaf2 %>% group_by(DAT, leaf_unique) %>%
  summarise(mean_vcmax=mean(vcmax_Best_Model),
            mean_jmax= mean(Jmax_Best),
            mean_tpumax= mean(TPU_Best)) %>%
  as.data.frame()

grp_tree <- curves_no_back %>% group_by(DAT,tree_id) %>% 
  summarise(mean_vcmax=mean(vcmax_Best_Model),
            sd_vcmax = sd(vcmax_Best_Model),
            mean_jmax= mean(Jmax_Best),
            sd_jmax= sd(Jmax_Best),
            mean_tpumax= mean(TPU_Best),
            sd_tpumax= sd(TPU_Best),
            .groups = 'drop') %>%
  as.data.frame()

#DAT vs Trad t-tests
dat_mean_df <- subset(grp_tree, DAT=="Before_DAT")
trad_mean_df <- subset(grp_tree, DAT=="Traditional")

res2<-t.test(dat_mean_df$mean_vcmax, trad_mean_df$mean_vcmax, paired=TRUE)
res2

res3<-t.test(dat_mean_df$mean_jmax, trad_mean_df$mean_jmax, paired=TRUE)
res3

res4<-t.test(dat_mean_df$mean_tpumax, trad_mean_df$mean_tpumax, paired=TRUE)
res4


#Testing whether the back-corrected are different from the non-back-filtered
dat_all <- subset(curves_final, DAT=="Before_DAT")

b2 <- ggplot(dat_all, aes(x=back_filt, y=vcmax_Best_Model)) +
  geom_boxplot()+
  labs(x="Dataset", y = "Vcmax")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(legend.position="none")
b2

b3 <- ggplot(dat_all, aes(x=back_filt, y=Jmax_Best)) +
  geom_boxplot()+
  labs(x="Dataset", y = "Jmax")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(legend.position="none")
b3

#T-tests
#back filtered vs non-back-filtered
dat_noback <- subset(dat_all, back_filt=="no_back")
dat_filt <- subset(dat_all, back_filt=="back_filtered")
res<-t.test(dat_noback$vcmax_Best_Model, dat_filt$vcmax_Best_Model, paired=TRUE) #Need same number
res #These are not significantly different!

res8<-t.test(dat_noback$Jmax_Best, dat_filt$Jmax_Best, paired=TRUE) #Need same number
res8 #These are not significantly different!

#Some plotting
b1 <- ggplot(curves_no_back, aes(x=curves_no_back$DAT, y=vcmax_Best_Model)) +
  geom_boxplot()+
  labs(x="Method", y = "Vcmax")+
  theme_classic()+
  theme(legend.position="none")
b1

b2 <- ggplot(curves_no_back, aes(x=curves_no_back$DAT, y=Jmax_Best)) +
  geom_boxplot()+
  labs(x="Method", y = "Jmax")+
  theme_classic()+
  theme(legend.position="none")
b2

g1 <- ggplot(data = curves_no_back, aes(x = tree_id, y = vcmax_Best_Model)) +
  geom_point() +
  xlab("Tree")+
  ylab("Vcmax") +
  theme_classic() +
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"))
g1

# Creates the linked point scatter


#GROUP THESE BY TREE!!
filt_par_dummy <- mutate(.data = grp_leaf,# makes a dummy variable to plot
                         dummy = if_else(grp_leaf$DAT == "Before_DAT", 0,1))

scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_vcmax,
                                                          color = leaf_unique)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  scale_color_viridis_d()+
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional")) +
  labs(x="Method", y="Vcmax")+
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))+
  guides(color = guide_legend(title = "Leaf Identifier"))
scat_vcmax

scat_tpu <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_tpumax,
                                                          color = substring(filt_par_dummy$leaf_unique, 1, 5))) +
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


#Testing for normality
mydata1 <- curves_no_back %>%
  select(tree_id,vcmax_Best_Model,DAT)%>%
  group_by(tree_id) %>%
  summarise(dif_vcmax=diff(vcmax_Best_Model)) %>%
  ungroup()

shapiro.test(mydata1$dif_vcmax)
ks.test(mydata1$dif_vcmax, 'pnorm')
#If in the output p-value > 0.05, it implies that the distribution of the data
#are not significantly different from normal distribution.
#In other words, we can assume the normality.
hist(mydata1$dif_vcmax) # Looks pretty normal
qqnorm(mydata1$dif_vcmax)
qqline(mydata1$dif_vcmax)
#it worked






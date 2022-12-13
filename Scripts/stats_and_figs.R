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
curv_nobk <- group_by(curves_no_back, DAT) %>%
  summarise(
    count = n(),
    mean = mean(vcmax_Best_Model, na.rm = TRUE),
    sd = sd(vcmax_Best_Model, na.rm = TRUE))

#Testing whether the back-corrected are different from the non-back-filtered
dat_all <- subset(curves_final, DAT=="Before_DAT")

b2 <- ggplot(dat_all, aes(x=back_filt, y=vcmax_Best_Model)) +
  geom_boxplot()+
  labs(x="Dataset", y = "Vcmax")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(legend.position="none")
b2

#T-tests
#back filtered vs non-back-filtered
dat_noback <- subset(dat_all, back_filt=="no_back")
dat_filt <- subset(dat_all, back_filt=="back_filtered")
res<-t.test(dat_noback$vcmax_Best_Model, dat_filt$vcmax_Best_Model, paired=TRUE) #Need same number
res #These are not significantly different!


#DAT vs Trad -- How to do this with different numbers of curves for DAT vs Trad???
dat_df <- subset(curves_no_back, DAT=="Before_DAT")
trad_df <- subset(curves_no_back, DAT=="Traditional")

# res2<-t.test(dat_grp$vcmax_Best_Model, trad_grp$vcmax_Best_Model, paired=TRUE) #Need same number
# res2
# 
# res3<-t.test(dat_df$Jmax_Best, trad_df$Jmax_Best, paired=TRUE) #Need same number
# res3
# 
# res4<-t.test(dat_df$TPU_Best, trad_df$TPU_Best, paired=TRUE) #Need same number
# res4
# 
# res5<-t.test(dat_df$Rd_Best, trad_df$Rd_Best, paired=TRUE) #Need same number
# res5

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
filt_par_dummy <- mutate(.data = curves_no_back,# makes a dummy variable to plot
                         dummy = if_else(curves_no_back$DAT == "Before_DAT", 0,1))

scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = vcmax_Best_Model,
                                                          color = leaf_id)) +
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






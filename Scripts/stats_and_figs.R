### For statistics and analysis

library(tidyverse)
library(ggpubr)

wd <- "/Users/charlessouthwick/Documents/GitHub/DAT_Proj"

#------------------------------------------ 
#Change the directory where you want to save your data
setwd(paste0(wd, "/Figures"))


curves_df <- read.csv("~/Documents/GitHub/DAT_Proj/Results/curve_fitting_out.csv")

#All this is to separate the concatenated tree_name column
curve_split <- unlist(str_split(curves_df$Tree_id, "_", n=2))
curve_sub <- subset(curve_split, curve_split != "Before_DAT" & curve_split != "Traditional")
curves_df$leaf_id <- curv_sub
curves_df_fixed <- subset(curves_df, select = -c(Tree_id))

leaf_split <- unlist(str_split(curves_df_fixed$leaf_id, "L", n=2))
leaf_sub <- subset(leaf_split, leaf_split != "1" & leaf_split != "2" & leaf_split != "1-1"
                   & leaf_split != 3 & leaf_split != "6" & leaf_split != "8" & leaf_split != "2-2")
curves_df_fixed$tree_id <- leaf_sub
curves_final <- curves_df_fixed

#This is the data without the back-correction filter
curves_no_back <- subset(curves_final, back_filt == "no_back")

#Some plotting
g1 <- ggplot(data = curves_no_back, aes(x = tree_id, y = vcmax_Best_Model)) +
  geom_point() +
  xlab("Tree")+
  ylab("Vcmax") +
  stat_regline_equation(label.y = 0.8, aes(label = ..rr.label..)) +
  theme_classic() +
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"))
g1

b1 <- ggplot(curves_no_back, aes(x=curves_no_back$DAT, y=vcmax_Best_Model)) +
  geom_boxplot()+
  labs(x="Method", y = "Vcmax")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(legend.position="none")
b1

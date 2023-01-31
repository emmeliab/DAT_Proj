# Statistical Analysis of A/Ci curve data

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(hydroGOF)

wd <- "C://Users/emmel/Desktop/DAT_proj/"
setwd(wd)

## Read in the datasets
params_ecophys <- read.csv(file = paste0(wd, "Results/params_ecophys.csv"), sep = ",", 
                           header = TRUE) %>% 
  filter(method == "dat")
params_ecophys <- add_row(.data = params_ecophys, unique = "K6702L1", Vcmax = 15.361861, 
                          Jmax = 27.750436, Rd = -1.010354,
                          TPU = 1.78906, Vcmax_SE = NA, Jmax_SE = NA, Rd_SE = NA, TPU_SE = NA, 
                          unique.1 = "K6702L1", method = "dat") %>% 
  arrange(unique) # since this curve did not run with the rest of them
#params_ecophys[7] <- c("K6706L1", Vcmax = 23.684034, 
#Jmax = 27.750436, Rd = -1.010354,
#TPU = 1.78906, Vcmax_SE = NA, Jmax_SE = NA, Rd_SE = NA, TPU_SE = NA, 
#unique.1 = "K6702L1", method = "dat")
params_photo <- read.csv(file = paste0(wd, "Results/dat_fit_ex_photo_pars.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000) ## TPU values at 1000 are coded as NA
params_mg <- read.csv(file = paste0(wd, "Results/curve_fitting_MG_out.csv"), sep = ",",
                      header = TRUE) %>% 
  filter(DAT == "Before_DAT", back_filt == "back_filtered")


# Comparing fits across photosynthesis, plantecophys, and MG code -----------
## Note: for comparisons of plantecophys and MG results, I used Vcmax and Jmax corrected to 25 degrees


# Vcmax
rmse(params_ecophys$Vcmax, params_photo$V_cmax)
ecovphoto_vcmax <- ggplot(mapping = aes(x = params_ecophys$Vcmax, y = params_photo$V_cmax, 
                                        color = params_photo$ID)) +
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys Vcmax", y="Photosynthesis Vcmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 150)) + 
  scale_y_continuous(limits = c(1, 150)) +
  annotate(geom = "text", label = "RMSE = 33.32", x = 125, y = 50)
ecovphoto_vcmax


rmse(params_ecophys$Vcmax, params_mg$Best_Vcmax_25C)
ecovMG_vcmax <- ggplot(mapping = aes(x = params_ecophys$Vcmax, y = params_mg$Best_Vcmax_25C, 
                                     color = params_mg$Tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys Vcmax", y="MG Vcmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 170)) + 
  scale_y_continuous(limits = c(1, 170)) +
  annotate(geom = "text", label = "RMSE = 16.63", x = 125, y = 50)
ecovMG_vcmax


rmse(params_photo$V_cmax, params_mg$vcmax_Best_Model)
photovMG_vcmax <- ggplot(mapping = aes(x = params_photo$V_cmax, y = params_mg$vcmax_Best_Model,
                                       color = params_photo$ID))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x = "Photosynthesis Vcmax", y = "MG Vcmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 170)) + 
  scale_y_continuous(limits = c(1, 170)) +
  annotate(geom = "text", label = "RMSE = 22.06", x = 125, y = 50)
photovMG_vcmax


#Jmax
rmse(params_ecophys$Jmax, params_mg$Best.Jmax_25C) ## see if we can fix the 800000 and try again
ecovMG_jmax <- ggplot(mapping = aes(x = params_ecophys$Jmax, y = params_mg$Best.Jmax_25C,
                                    color = params_mg$Tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys Jmax", y="MG Jmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 200)) + 
  scale_y_continuous(limits = c(1, 200)) +
  annotate(geom = "text", label = "RMSE = 23.27", x = 125, y = 50)
ecovMG_jmax
# Note: K6706L1 jmax for plantecophys is 800000, and is not included in the graph


rmse(params_ecophys$Jmax, params_photo$J_max)
ecovphoto_jmax <- ggplot(mapping = aes(x = params_ecophys$Jmax, y = params_photo$J_max, 
                                       color = params_photo$ID))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys Jmax", y="Photosynthesis Jmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 170)) + 
  scale_y_continuous(limits = c(1, 170)) +
  annotate(geom = "text", label = "RMSE = 12.74", x = 125, y = 50)
ecovphoto_jmax
# Note: K6706L1 jmax for plantecophys is 800000, and is not included in the graph


rmse(params_photo$J_max, params_mg$Jmax_Best)
photovMG_jmax <- ggplot(mapping = aes(x = params_photo$J_max, y = params_mg$Jmax_Best,
                                      color = params_photo$ID))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Photosynthesis Jmax", y="MG Jmax", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 200)) + 
  scale_y_continuous(limits = c(1, 200)) +
  annotate(geom = "text", label = "RMSE = 37.64", x = 125, y = 50)
photovMG_jmax


#TPU
rmse(params_ecophys$TPU, params_mg$TPU_Best)
ecovMG_tpu <- ggplot(mapping = aes(x = params_ecophys$TPU, y = params_mg$TPU_Best,
                                   color = params_mg$Tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys TPU", y="MG TPU", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 11)) + 
  scale_y_continuous(limits = c(1, 11)) +
  annotate(geom = "text", label = "RMSE = 1.03", x = 7, y = 3)
ecovMG_tpu


rmse(params_ecophys$TPU, params_photo$V_TPU)
ecovphoto_tpu <- ggplot(mapping = aes(x = params_ecophys$TPU, y = params_photo$V_TPU,
                                      color = params_photo$ID))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Plantecophys TPU", y="Photosynthesis TPU", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 10)) + 
  scale_y_continuous(limits = c(1, 10)) +
  annotate(geom = "text", label = "RMSE = 0.32", x = 7, y = 3)
ecovphoto_tpu
# Note: TPU for several of the curves via the photosynthesis package are at 1000 (likely meaning it
# wasn't fit). These are not included in the chart


rmse(params_photo$V_TPU, params_mg$TPU_Best)
photovMG_tpu <- ggplot(mapping = aes(x = params_photo$V_TPU, y = params_mg$TPU_Best, 
                                     color = params_photo$ID))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Photosynthesis TPU", y="MG TPU", col = "Leaf")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif")) +
  scale_x_continuous(limits = c(1, 11)) + 
  scale_y_continuous(limits = c(1, 11)) +
  annotate(geom = "text", label = "RMSE = 1.47", x = 7, y = 3)
photovMG_tpu






# Plantecophys result visualization (fix code to make sense with new variables)---------------------------------------
## Boxplots

### Vcmax
box_both_vcmax <- filt_par_species %>%
  ggplot() +
  geom_boxplot(aes(x = method, y = Vcmax))
box_both_vcmax


### Jmax
box_both_jmax <- filt_par_species %>% 
  ggplot() +
  geom_boxplot(aes(x = method, y = Jmax))
box_both_jmax



# Stacked Scatters
filt_par_dummy <- mutate(.data = filt_par_species,# makes a dummy variable to plot
                         dummy = if_else(filt_par_species$method == "dat", 0,1))



## Vcmax
scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = Vcmax,
                                                          color = unique)) + 
  geom_line() + 
  geom_point() +
  theme_light() +
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Trad")) +
  xlab("Method")
scat_vcmax


## Jmax
scat_jmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = Jmax,
                                                         color = unique)) + 
  geom_line() + 
  geom_point() +
  theme_light() +
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Trad")) +
  xlab("Method")
scat_jmax


# Testing Assumptions of Plantecophys results -------------------------------------------------------


#### Delete this after fixing code v
# ## read in the data
# par_species <- read.csv(file = paste0(wd, "/Results/params_ecophys.csv"), header = TRUE, sep = ",")
# 
# ## Filter some outliers
# which(par_species$Vcmax > 100)
# which(par_species$Vcmax < 0)
# filt_par_species <- par_species %>% 
#   filter(Vcmax < 100 & Vcmax > 0)
# head(filt_par_species)
###### ^

## Summary Stats
table(filt_par_species$method) ## number of each type of curve

group_by(filt_par_species, method) %>%
  summarise(
    count = n(),
    mean = mean(Vcmax, na.rm = TRUE),
    sd = sd(Vcmax, na.rm = TRUE))


# Q-Q plots
ggqqplot(filt_par_species$Vcmax)
ggqqplot(filt_par_species$Jmax)



# Shapiro-Wilk test
shapiro.test(filt_par_species$Vcmax)
shapiro.test(filt_par_species$Jmax)



#I think a wilcox signed-rank test uses non-parametric, paired data.
wilcox.test(filt_par_species$Vcmax ~ filt_par_species$method, paired = TRUE)






# Photosynthesis results visualization ------------------------------------



# MG results visualization ------------------------------------------------

## Separate the concatenated tree ID column
curve_split <- unlist(str_split(curves_df$Tree_id, "_", n=2))
curve_sub <- subset(curve_split, curve_split != "Before_DAT" & curve_split != "Traditional")
curves_df$leaf_id <- curve_sub
curves_df_fixed <- subset(curves_df, select = -c(Tree_id))

leaf_split <- unlist(str_split(curves_df_fixed$leaf_id, "L", n=2))
leaf_sub <- subset(leaf_split, leaf_split != "1" & leaf_split != "2" & leaf_split != "1-1"
                   & leaf_split != 3 & leaf_split != "6" & leaf_split != "8" & leaf_split != "2-2")
curves_df_fixed$tree_id <- leaf_sub
curves_final2 <- curves_df_fixed

#Add in the relative canopy position and four letter code
codes <- read.csv("~/Documents/GitHub/DAT_Proj/Inputs/unique_ids.csv")
canopy_pos <- read.csv("~/Documents/GitHub/DAT_Proj/Inputs/rel_canopy_pos.csv")
names(codes)[1] ="leaf_id"
codes_and_can <- left_join(codes, canopy_pos, by = "k67.id")
names(codes_and_can)[3]="code4let"
codes_and_can <- subset(codes_and_can, select = -code.y)




## Boxplots of filtered vs. nonfiltered data
dat_all <- subset(curves_final, DAT=="Before_DAT")

label_backfilt <- c('Back Filtered', 'Original')
b1 <- ggplot(dat_all, aes(x=back_filt, y=vcmax_Best_Model)) +
  geom_boxplot()+
  labs(x="Dataset", y = "Vcmax")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.position="none")+
  scale_x_discrete(labels=label_backfilt)
b1

b2 <- ggplot(dat_all, aes(x=back_filt, y=Jmax_Best)) +
  geom_boxplot()+
  labs(x="Dataset", y = "Jmax")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.position="none")+
  scale_x_discrete(labels=label_backfilt)
b2

b3 <- ggplot(dat_all, aes(x=back_filt, y=TPU_Best)) +
  geom_boxplot()+
  labs(x="Dataset", y = "TPU")+
  scale_fill_manual(values=c("#E69F00","#7bccc4"))+
  theme_classic()+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.position="none")+
  scale_x_discrete(labels=label_backfilt)
b3




#### Visualization of DAT vs Traditional

#This is the data without the back-correction filter

## Boxplots
curves_no_back <- subset(curves_final, back_filt == "no_back")
leaf2 <- mutate(curves_no_back, leaf_unique = substring(curves_no_back$leaf_id, 1, 7))

lab_DATTrad <- c('DAT', 'Traditional')
b4 <- ggplot(leaf2, aes(x=DAT, y=vcmax_Best_Model)) +
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

b5 <- ggplot(leaf2, aes(x=DAT, y=Jmax_Best)) +
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

b6 <- ggplot(leaf2, aes(x=DAT, y=TPU_Best)) +
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

# Group data for further analysis
grp_leaf <- leaf2 %>% group_by(DAT, leaf_unique) %>%
  summarise(mean_vcmax=mean(vcmax_Best_Model),
            mean_jmax= mean(Jmax_Best),
            mean_tpumax= mean(TPU_Best),
            code4let=code4let,
            rel_can_pos=rel_can_pos) %>%
  as.data.frame()

grp_tree <- leaf2 %>% group_by(DAT,tree_id) %>% 
  summarise(mean_vcmax=mean(vcmax_Best_Model),
            mean_jmax= mean(Jmax_Best),
            mean_tpumax= mean(TPU_Best),
            code4let=code4let,
            rel_can_pos=rel_can_pos,
            .groups = 'drop') %>%
  as.data.frame()


## Stacked Scatters by leaf
filt_par_dummy <- mutate(.data = grp_leaf,# makes a dummy variable to plot
                         dummy = if_else(grp_leaf$DAT == "Before_DAT", 0,1))

scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_vcmax,
                                                          color = leaf_unique,
                                                          label = code4let)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
                  size          = 2.4,
                  box.padding   = 0.25,
                  point.padding = 0.25,
                  segment.size  = 0.2,
                  segment.linetype = 3,
                  direction     = "y",
                  nudge_x = -0.2)+ ## Playing around with this to help visualize
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                     expand = expansion(mult=0.3)) +
  labs(x="Method", y="Vcmax")+
  theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
        axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")
scat_vcmax

scat_jmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_jmax,
                                                         color = leaf_unique,
                                                         label = code4let)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
                  size          = 2.8,
                  box.padding   = 0.25,
                  point.padding = 0.25,
                  segment.size  = 0.2,
                  direction     = "y",
                  nudge_x = -0.2)+ ## Playing around with this to help visualize
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                     expand = expansion(mult=0.3)) +
  labs(x="Method", y="Jmax")+
  theme(aspect.ratio = 1.5,
        axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")+
  guides(color = guide_legend(title = "Leaf Identifier"))
scat_jmax

scat_tpu <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_tpumax,
                                                        color = leaf_unique,
                                                        label = code4let)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  geom_text_repel(data          = subset(filt_par_dummy, DAT == "Before_DAT"),
                  size          = 2.8,
                  box.padding   = 0.25,
                  point.padding = 0.25,
                  segment.size  = 0.2,
                  direction     = "y",
                  nudge_x = -0.2)+ ## Playing around with this to help visualize
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                     expand = expansion(mult=0.3)) +
  labs(x="Method", y="TPU")+
  theme(aspect.ratio = 1.5,
        axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position ="none")+
  guides(color = guide_legend(title = "Leaf Identifier"))
scat_tpu


## Stacked scatter by tree
filt_par_dummy2 <- mutate(.data = grp_tree,# makes a dummy variable to plot
                          dummy = if_else(grp_tree$DAT == "Before_DAT", 0,1))

scat_vcmax2 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_vcmax,
                                                            color = tree_id,
                                                            label = code4let)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
                  size          = 2.8,
                  box.padding   = 0.25,
                  point.padding = 0.25,
                  segment.size  = 0.2,
                  direction     = "y",
                  nudge_x = -0.2)+
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                     expand = expansion(mult=0.3)) +
  labs(x="Method", y="Vcmax")+
  theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
        axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")
scat_vcmax2

scat_jmax2 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_jmax,
                                                           color = tree_id,
                                                           label = code4let)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
                  size          = 2.8,
                  box.padding   = 0.25,
                  point.padding = 0.25,
                  segment.size  = 0.2,
                  direction     = "y",
                  nudge_x = -0.2)+
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                     expand = expansion(mult=0.3)) +
  labs(x="Method", y="Jmax")+
  theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
        axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")
scat_jmax2

scat_tpu2 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_tpumax,
                                                          color = tree_id,
                                                          label = code4let)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
                  size          = 2.8,
                  box.padding   = 0.25,
                  point.padding = 0.25,
                  segment.size  = 0.2,
                  direction     = "y",
                  nudge_x = -0.2)+
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                     expand = expansion(mult=0.3)) +
  labs(x="Method", y="TPU")+
  theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
        axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")
scat_tpu2

#grouped scatters labelled with relative canopy heights
scat_vcmax3 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_vcmax,
                                                            color = tree_id,
                                                            label = rel_can_pos)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
                  size          = 2.8,
                  box.padding   = 0.25,
                  point.padding = 0.25,
                  segment.size  = 0.2,
                  direction     = "y",
                  nudge_x = -0.2)+
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                     expand = expansion(mult=0.3)) +
  labs(x="Method", y="Vcmax")+
  theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
        axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")
scat_vcmax3

scat_jmax3 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_jmax,
                                                           color = tree_id,
                                                           label = rel_can_pos)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
                  size          = 2.8,
                  box.padding   = 0.25,
                  point.padding = 0.25,
                  segment.size  = 0.2,
                  direction     = "y",
                  nudge_x = -0.2)+
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                     expand = expansion(mult=0.3)) +
  labs(x="Method", y="Jmax")+
  theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
        axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")
scat_jmax3

scat_tpu3 <- ggplot(data = filt_par_dummy2, mapping = aes(x = dummy, y = mean_tpumax,
                                                          color = tree_id,
                                                          label = rel_can_pos)) +
  geom_line() + 
  geom_point() +
  theme_classic() +
  geom_text_repel(data          = subset(filt_par_dummy2, DAT == "Before_DAT"),
                  size          = 2.8,
                  box.padding   = 0.25,
                  point.padding = 0.25,
                  segment.size  = 0.2,
                  direction     = "y",
                  nudge_x = -0.2)+
  scale_x_continuous(breaks = c(0,1), labels = c("DAT", "Traditional"),
                     expand = expansion(mult=0.3)) +
  labs(x="Method", y="TPU")+
  theme(aspect.ratio = 1.5, #Trying to adjust the sizing to look a bit better
        axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=9, family = "serif"),
        legend.title=element_text(size=11, family = "serif"),
        legend.position="none")
scat_tpu3




## 1:1 plots DAT vs Trad
# By leaf

#Vcmax
leaf_sub_vcmax <- select(grp_leaf, mean_vcmax, DAT, leaf_unique, code4let)
leaf_wide_vcmax <- reshape(leaf_sub_vcmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_vcmax)[2:4]=c("vcmax_DAT", "code4let", "vcmax_Trad")
leaf_wide_vcmax <- subset(leaf_wide_vcmax, select = -code4let.Traditional)
mng_leaf_vcmax <- ggplot(data = leaf_wide_vcmax, mapping = aes(x = vcmax_Trad,
                                                               y = vcmax_DAT,
                                                               color = code4let))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional Vcmax", y="DAT Vcmax", col = "Species Code")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_leaf_vcmax

#Jmax
leaf_sub_jmax <- select(grp_leaf, mean_jmax, DAT, leaf_unique, code4let)
leaf_wide_jmax <- reshape(leaf_sub_jmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_jmax)[2:4]=c("jmax_DAT", "code4let", "jmax_Trad")
leaf_wide_jmax <- subset(leaf_wide_jmax, select = -code4let.Traditional)
mng_leaf_jmax <- ggplot(data = leaf_wide_jmax, mapping = aes(x = jmax_Trad,
                                                             y = jmax_DAT,
                                                             color = code4let))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional Jmax", y="DAT Jmax", col = "Species Code")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_leaf_jmax

#TPU
leaf_sub_tpu <- select(grp_leaf, mean_tpumax, DAT, leaf_unique, code4let)
leaf_wide_tpu <- reshape(leaf_sub_tpu, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_tpu)[2:4]=c("tpu_DAT", "code4let", "tpu_Trad")
leaf_wide_tpu <- subset(leaf_wide_tpu, select = -code4let.Traditional)
mng_leaf_tpu <- ggplot(data = leaf_wide_tpu, mapping = aes(x = tpu_Trad,
                                                           y = tpu_DAT,
                                                           color = code4let))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional TPU", y="DAT TPU", col = "Species Code")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_leaf_tpu



# by tree

#Vcmax
tree_sub_vcmax <- select(grp_tree, mean_vcmax, DAT, tree_id, code4let)
tree_wide_vcmax <- reshape(tree_sub_vcmax, idvar = "tree_id", timevar = "DAT", direction = "wide")
names(tree_wide_vcmax)[2:4]=c("vcmax_DAT", "code4let", "vcmax_Trad")
tree_wide_vcmax <- subset(tree_wide_vcmax, select = -code4let.Traditional)
mng_tree_vcmax <- ggplot(data = tree_wide_vcmax, mapping = aes(x = vcmax_Trad,
                                                               y = vcmax_DAT,
                                                               color = code4let))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional Vcmax", y="DAT Vcmax", col = "Species Code")+
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_tree_vcmax

#Jmax
tree_sub_jmax <- select(grp_tree, mean_jmax, DAT, tree_id, code4let)
tree_wide_jmax <- reshape(tree_sub_jmax, idvar = "tree_id", timevar = "DAT", direction = "wide")
names(tree_wide_jmax)[2:4]=c("jmax_DAT", "code4let", "jmax_Trad")
tree_wide_jmax <- subset(tree_wide_jmax, select = -code4let.Traditional)
mng_tree_jmax <- ggplot(data = tree_wide_jmax, mapping = aes(x = jmax_Trad,
                                                             y = jmax_DAT,
                                                             color = code4let))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional Jmax", y="DAT Jmax", col = "Species Code")+
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_tree_jmax

#TPU
tree_sub_tpu <- select(grp_tree, mean_tpumax, DAT, tree_id, code4let)
tree_wide_tpu <- reshape(tree_sub_tpu, idvar = "tree_id", timevar = "DAT", direction = "wide")
names(tree_wide_tpu)[2:4]=c("tpu_DAT", "code4let", "tpu_Trad")
tree_wide_tpu <- subset(tree_wide_vcmax, select = -code4let.Traditional)
mng_tree_tpu <- ggplot(data = tree_wide_tpu, mapping = aes(x = tpu_Trad,
                                                           y = tpu_DAT,
                                                           color = code4let))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional TPU", y="DAT TPU", col = "Species Code")+
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_tree_tpu



# Stats for MG code --------------------------------

#Testing whether the back-corrected are different from the non-back-corrected
#T-tests
dat_noback <- subset(dat_all, back_filt=="no_back")
dat_filt <- subset(dat_all, back_filt=="back_filtered")
res7<-t.test(dat_noback$vcmax_Best_Model, dat_filt$vcmax_Best_Model, paired=TRUE) #Need same number
res7 #These are not significantly different!

res8<-t.test(dat_noback$Jmax_Best, dat_filt$Jmax_Best, paired=TRUE) #Need same number
res8 #These are not significantly different!

res9<-t.test(dat_noback$TPU_Best, dat_filt$TPU_Best, paired=TRUE) #Need same number
res9


## Test for Normality by leaf
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


## test for normality by tree
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



## T-tests
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






# 
# setwd(paste0(wd, "/Figures"))
# 
# curves_df <- read.csv("~/Documents/GitHub/DAT_Proj/Results/curve_fitting_MG_out.csv")
# 
# 
# 
# 
# curves_final <- left_join(curves_final2, codes_and_can, by = "leaf_id")
# 









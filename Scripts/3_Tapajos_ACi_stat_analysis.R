# Statistical Analysis of A/Ci curve data

library(tidyverse)
library(ggpubr)
library(ggrepel)
library(hydroGOF)

wd <- "/Users/charlessouthwick/Documents/GitHub/DAT_Proj/"
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
init_mg <- read.csv(file = paste0(wd, "Results/curve_fitting_MG_out.csv"), sep = ",",
                      header = TRUE)
init_mg$back_filt[init_mg$DAT == "Traditional"] <- "irrelevant"
params_mg <- init_mg %>% 
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



# MG data processing ------------------------------------------------

## Separate the concatenated tree ID column
mg_split <- unlist(str_split(init_mg$Tree_id, "_", n=2))
mg_sub <- subset(mg_split, mg_split != "Before_DAT" & mg_split != "Traditional")
init_mg$leaf_id <- mg_sub
mg_complete <- subset(init_mg, select = -c(Tree_id))

leaf_split <- unlist(str_split(mg_complete$leaf_id, "L", n=2))
leaf_sub <- subset(leaf_split, leaf_split != "1" & leaf_split != "2" & leaf_split != "1-1"
                   & leaf_split != 3 & leaf_split != "6" & leaf_split != "8" & leaf_split != "2-2")
mg_complete$tree_id <- leaf_sub

#Add in the relative canopy position
codes <- read.csv("~/Documents/GitHub/DAT_Proj/Inputs/unique_ids.csv") # need to keep this in for leaf_id
canopy_pos <- read.csv("~/Documents/GitHub/DAT_Proj/Inputs/rel_canopy_pos.csv")
names(codes)[1] ="leaf_id"
codes_and_can <- left_join(codes, canopy_pos, by = "k67.id")
names(codes_and_can)[3]="code4let"
codes_and_can <- subset(codes_and_can, select = -code.y)

mg_complete <- left_join(mg_complete, codes_and_can, by = "leaf_id") %>% 
  select(-c(k67.id,code4let))

# For back-filter vs no_back analysis
mg_all_dat <- subset(mg_complete, DAT=="Before_DAT")

#This is the data with the back correction filtered out
#mg_no_back <- subset(mg_complete, back_filt == "back_filtered")
mg_leaf <- mg_complete %>%
  subset(back_filt != "no_back") %>%
  mutate(leaf_unique = substring(leaf_id, 1, 7))

#group data for further analysis
grp_leaf <- mg_leaf %>% group_by(DAT, leaf_unique) %>%
  summarise(mean_vcmax=mean(Best_Vcmax_25C),
            mean_jmax= mean(Best.Jmax_25C),
            mean_tpumax= mean(TPU_Best),
            tree_id = tree_id,
            leaf_unique=leaf_unique,
            leaf_id=leaf_id,
            rel_can_pos=rel_can_pos) %>%
  as.data.frame()
summary(grp_leaf)

grp_tree <- mg_leaf %>% group_by(DAT,tree_id) %>% 
  summarise(mean_vcmax=mean(Best_Vcmax_25C),
            mean_jmax= mean(Best.Jmax_25C),
            mean_tpumax= mean(TPU_Best),
            tree_id = tree_id,
            leaf_unique=leaf_unique,
            leaf_id=leaf_id,
            rel_can_pos=rel_can_pos,
            .groups = 'drop') %>%
  as.data.frame()


### MG Visualizations -------------------------------
## Boxplots of filtered vs. nonfiltered data
label_backfilt <- c('Back Filtered', 'Original')
b1 <- ggplot(mg_all_dat, aes(x=back_filt, y=Best_Vcmax_25C)) +
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

b2 <- ggplot(mg_all_dat, aes(x=back_filt, y=Best.Jmax_25C)) +
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

b3 <- ggplot(mg_all_dat, aes(x=back_filt, y=TPU_Best)) +
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

## Boxplots
lab_DATTrad <- c('DAT', 'Traditional')
b4 <- ggplot(mg_leaf, aes(x=DAT, y=Best_Vcmax_25C)) +
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

b5 <- ggplot(mg_leaf, aes(x=DAT, y=Best.Jmax_25C)) +
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

b6 <- ggplot(mg_leaf, aes(x=DAT, y=TPU_Best)) +
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
g1 <- ggplot(mg_leaf, aes(x = tree_id, y = Best_Vcmax_25C)) +
  geom_point() +
  xlab("Tree")+
  ylab("Vcmax") +
  theme_classic() +
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"))
g1


## Stacked Scatters by leaf
filt_par_dummy <- mutate(.data = grp_leaf,# makes a dummy variable to plot
                         dummy = if_else(grp_leaf$DAT == "Before_DAT", 0,1))

scat_vcmax <- ggplot(data = filt_par_dummy, mapping = aes(x = dummy, y = mean_vcmax,
                                                          color = leaf_unique,
                                                          label = tree_id)) +
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
                                                         label = tree_id)) +
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
                                                        label = tree_id)) +
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
                                                            label = tree_id)) +
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
                                                           label = tree_id)) +
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
                                                          label = tree_id)) +
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
leaf_sub_vcmax <- select(grp_leaf, mean_vcmax, DAT, leaf_unique, tree_id)
leaf_wide_vcmax <- reshape(leaf_sub_vcmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_vcmax)[2:4]=c("vcmax_DAT", "tree_id", "vcmax_Trad")
leaf_wide_vcmax <- subset(leaf_wide_vcmax, select = -tree_id.Traditional)
mng_leaf_vcmax <- ggplot(data = leaf_wide_vcmax, mapping = aes(x = vcmax_Trad,
                                                               y = vcmax_DAT,
                                                               color = tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional Vcmax", y="DAT Vcmax", col = "Tree ID")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_leaf_vcmax

#Jmax
leaf_sub_jmax <- select(grp_leaf, mean_jmax, DAT, leaf_unique, tree_id)
leaf_wide_jmax <- reshape(leaf_sub_jmax, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_jmax)[2:4]=c("jmax_DAT", "tree_id", "jmax_Trad")
leaf_wide_jmax <- subset(leaf_wide_jmax, select = -tree_id.Traditional)
mng_leaf_jmax <- ggplot(data = leaf_wide_jmax, mapping = aes(x = jmax_Trad,
                                                             y = jmax_DAT,
                                                             color = tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional Jmax", y="DAT Jmax", col = "Tree ID")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_leaf_jmax

#TPU
leaf_sub_tpu <- select(grp_leaf, mean_tpumax, DAT, leaf_unique, tree_id)
leaf_wide_tpu <- reshape(leaf_sub_tpu, idvar = "leaf_unique", timevar = "DAT", direction = "wide")
names(leaf_wide_tpu)[2:4]=c("tpu_DAT", "tree_id", "tpu_Trad")
leaf_wide_tpu <- subset(leaf_wide_tpu, select = -tree_id.Traditional)
mng_leaf_tpu <- ggplot(data = leaf_wide_tpu, mapping = aes(x = tpu_Trad,
                                                           y = tpu_DAT,
                                                           color = tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional TPU", y="DAT TPU", col = "Tree ID")+
  theme(axis.title.x=element_text(size=18, family = "serif"),
        axis.title.y=element_text(size=18, family = "serif"),
        axis.text.x=element_text(size=15, family = "serif"),
        axis.text.y=element_text(size=15, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_leaf_tpu


# by tree

#Vcmax
tree_sub_vcmax <- select(grp_tree, mean_vcmax, DAT, tree_id)
tree_wide_vcmax <- reshape(tree_sub_vcmax, idvar = "tree_id", timevar = "DAT", direction = "wide")
names(tree_wide_vcmax)[2:4]=c("vcmax_DAT", "tree_id", "vcmax_Trad")
tree_wide_vcmax <- subset(tree_wide_vcmax, select = -tree_id.Traditional)
mng_tree_vcmax <- ggplot(data = tree_wide_vcmax, mapping = aes(x = vcmax_Trad,
                                                               y = vcmax_DAT,
                                                               color = tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional Vcmax", y="DAT Vcmax", col = "Tree ID")+
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_tree_vcmax

#Jmax
tree_sub_jmax <- select(grp_tree, mean_jmax, DAT, tree_id)
tree_wide_jmax <- reshape(tree_sub_jmax, idvar = "tree_id", timevar = "DAT", direction = "wide")
names(tree_wide_jmax)[2:4]=c("jmax_DAT", "tree_id", "jmax_Trad")
tree_wide_jmax <- subset(tree_wide_jmax, select = -tree_id.Traditional)
mng_tree_jmax <- ggplot(data = tree_wide_jmax, mapping = aes(x = jmax_Trad,
                                                             y = jmax_DAT,
                                                             color = tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional Jmax", y="DAT Jmax", col = "Tree ID")+
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_tree_jmax

#TPU
tree_sub_tpu <- select(grp_tree, mean_tpumax, DAT, tree_id)
tree_wide_tpu <- reshape(tree_sub_tpu, idvar = "tree_id", timevar = "DAT", direction = "wide")
names(tree_wide_tpu)[2:4]=c("tpu_DAT", "tree_id", "tpu_Trad")
tree_wide_tpu <- subset(tree_wide_vcmax, select = -tree_id.Traditional)
mng_tree_tpu <- ggplot(data = tree_wide_tpu, mapping = aes(x = tpu_Trad,
                                                           y = tpu_DAT,
                                                           color = tree_id))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  theme_classic()+
  labs(x="Traditional TPU", y="DAT TPU", col = "Tree ID")+
  theme(axis.title.x=element_text(size=11, family = "serif"),
        axis.title.y=element_text(size=11, family = "serif"),
        axis.text.x=element_text(size=11, family = "serif"),
        axis.text.y=element_text(size=11, family = "serif"),
        legend.text=element_text(size=7, family = "serif"),
        legend.title=element_text(size=11, family = "serif"))
mng_tree_tpu


# Stats for MG --------------------------------

library(Publish) 
library(moments) 
library(vcd)
library(effsize)
library(car)
library(pwr)

#Subset variables of interest
leaf_stat <- select(mg_leaf, 'DAT', 'Best_Vcmax_25C', 'Best.Jmax_25C', 'TPU_Best', 'leaf_id', 'tree_id', 'leaf_unique', 'rel_can_pos')
leaf_stat$DAT[leaf_stat$DAT == "Before_DAT"] <- "DAT"
leaf_stat <- rename(leaf_stat,
                    method = DAT,
                    vcmax = Best_Vcmax_25C,
                    jmax = Best.Jmax_25C,
                    tpu = TPU_Best)
#Describe factor levels. 0 is traditional, 1 is DAT
leaf_stat$method <- factor(leaf_stat$method)

#displays grouped summary
leaf_summary <- leaf_stat %>%
  group_by(method) %>%
  summarise(n_vcmax = length(vcmax),
            mean_vcmax = round(mean(vcmax),3),
            sd_vcmax = round(sd(vcmax),3),
            se_vcmax = sd(vcmax)/sqrt(n())) %>% 
  mutate(low_ci_vcmax = mean_vcmax - qt(1 - (0.05 / 2), n_vcmax - 1) * se_vcmax,
         up_ci_vcmax = mean_vcmax + qt(1 - (0.05 / 2), n_vcmax - 1) * se_vcmax)
print(leaf_summary)

#Visualize Vcmax by Method
ci.mean(vcmax ~ method, data = leaf_stat)
ci1<-ci.mean(vcmax~method, data=leaf_stat)
plot(ci1,title.labels="Method")

#Histogram to visualize
leaf_hist<-ggplot(leaf_stat, aes(x=vcmax)) + 
  geom_histogram(color="black", fill="white", bins = 8)+
  geom_vline(aes(xintercept=mean(vcmax)),
             color="red", linetype="dashed", linewidth=0.5)+
  theme_classic()
leaf_hist
#data are positively skewed
skewness(leaf_stat$vcmax)
#value close to 1, should be okay
kurtosis(leaf_stat$vcmax)
#It's a bit high

#Test the assumption of equal variances for each group for t-test with Levene's
leveneTest(vcmax ~ method, data = leaf_stat)
#Non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances

# Shapiro-Wilk normality test for Vcmax for the one-sample t-test
with(leaf_stat, shapiro.test(vcmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(leaf_stat, shapiro.test(vcmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(leaf_stat, shapiro.test(vcmax[method == "Traditional"])) 
#All are significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

grp_dat <- leaf_stat %>%
  filter(method == "DAT") %>% 
  group_by(tree_id) %>% 
  summarise(n_vcmax = length(vcmax),
            mean_vcmax = mean(vcmax),
            sd_vcmax = sd(vcmax),
            se_vcmax = sd(vcmax)/sqrt(n())) %>% 
  mutate(method = "DAT")
  
grp_trad <- leaf_stat %>%
  filter(method == "Traditional") %>% 
  group_by(tree_id) %>% 
  summarise(n_vcmax = length(vcmax),
            mean_vcmax = mean(vcmax),
            sd_vcmax = sd(vcmax),
            se_vcmax = sd(vcmax)/sqrt(n())) %>% 
  mutate(method = "Traditional")

grp_all <- rbind(grp_dat, grp_trad)

with(grp_all, shapiro.test(mean_vcmax))
# Shapiro-Wilk normality test for the DAT measurement methodology
with(grp_all, shapiro.test(mean_vcmax[method == "DAT"]))
# Shapiro-Wilk normality test for the Traditional measurement methodology
with(grp_all, shapiro.test(mean_vcmax[method == "Traditional"])) 
#All are significant. Rejects null hypothesis that these data are not normally distributed.
# Therefore the sample varies from the normal distribution.

wilcox.test(mean_vcmax ~ method, data = grp_all, paired = TRUE)
#This result is significant

#Effect size for the independent sample t-test:
cohen.ES(test = "t", size = "large") # To remind oneself
cohen.d(mean_vcmax ~ method | Subject(tree_id), data=grp_all, paired = TRUE)
#small effect size

#### Power Analysis

d <- cohen.d(mean_vcmax ~ method | Subject(tree_id), data=grp_all, paired = TRUE)
d[["estimate"]]

#How many samples to achieve a certain power?
pwr.t.test(n = , d = d[["estimate"]], power = 0.8, sig.level = 0.05, type = "paired", alternative = "two.sided")

#What was the power of our study?
pwr.t.test(n = 13, d = d[["estimate"]], power = , sig.level = 0.05, type = "paired", alternative = "two.sided")

#Try log transform!!!!!

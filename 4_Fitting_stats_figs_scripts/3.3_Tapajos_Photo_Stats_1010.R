# Computing statistical analysis on the A/Ci curve fits from 2_Tapajos_Fit_ACi.R
# Using the 'photosynthesis' package

######### Check to see which of these packages we actually need
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(hydroGOF)
library(Publish)
library(moments) 
library(vcd)
library(car)
library(rstatix) #For wilcox_effsize function
library(gridExtra)
library(reshape2)
library(here)


# Load and Subest the datasets -------------------------------------------

### Excluding K6709L3 as the SS curve was not a matched pair with the DAT curve

pho_dat <- read.csv(file = here("5_Results/DAT_photo_pars_crct_noTPU.csv"), 
                    sep = ",", 
                    header = TRUE, na.strings = 1000) %>% 
    ## TPU values at 1000 are coded as NA
    subset(ID != "K6709L3")

pho_SS <- read.csv(file = here("5_Results/SS_photo_pars_crct_noTPU.csv"),
                   sep = ",", 
                   header = TRUE, na.strings = 1000) %>% 
    subset(ID != "K6709L3")


pho_dat_tpu <- read.csv(file = here("5_Results/dat_fits_photo_pars_filt_correct_with_TPU.csv"),
                        sep = ",", 
                        header = TRUE, na.strings = 1000) %>%  
    ## TPU values at 1000 are coded as NA
    subset(ID != "K6709L3")

pho_SS_tpu <- read.csv(file = here("5_Results/trad_fits_photo_pars_correct_with_TPU.csv"),
                       sep = ",", 
                       header = TRUE, na.strings = 1000) %>% 
    subset(ID != "K6709L3")


### Append the datasets together (separate for TPU and noTPU)
pho_both <- rbind(pho_dat, pho_SS)
pho_both_tpu <- rbind(pho_dat_tpu, pho_SS_tpu)


# Subset out the No overshoot data

pho_nd <- pho_dat %>% 
    subset(ID != "K6707L2" & ID != "K6707L2-2" & ID != "K6707L1" & ID != "K6707L1-1" & ID != "K6709L6" & ID != "K6714L2" & ID != "K6714L1" & ID != "K6702L1" & ID != "K6706L2" & ID != "K6706L1" & ID != "K6709L2")

pho_nd_tpu <- pho_dat_tpu %>% 
    subset(ID != "K6707L2" & ID != "K6707L2-2" & ID != "K6707L1" & ID != "K6707L1-1" & ID != "K6709L6" & ID != "K6714L2" & ID != "K6714L1" & ID != "K6702L1" & ID != "K6706L2" & ID != "K6706L1" & ID != "K6709L2")

### Append those datasets together
pho_nd_both <- rbind(pho_nd, pho_SS)
pho_nd_both_tpu <- rbind(pho_SS_tpu, pho_nd_tpu)


# Initial Processing/Grouping ----------------------------------

# Full data, no TPU

pho_leaf <- mutate(pho_both, leaf_unique = substring(ID, 1, 7))
##### CHECK! MAY NOT NEED 'V_TPU'
pho_stat <- select(pho_leaf, 
                   'curv_meth', 'Best_Vcmax_25C', 'Best_Jmax_25C',
                   'V_TPU', 'ID', 
                   'leaf_unique', "V_cmax_se", "J_se", "V_TPU_se")
pho_stat <- rename(pho_stat,
                   vcmax = Best_Vcmax_25C,
                   jmax = Best_Jmax_25C,
                   tpu = V_TPU,
                   leaf_id = ID) %>% 
    mutate(fit_type = "no_tpu")
# Describe factor levels: 0 is SS, 1 is DAT
pho_stat$curv_meth <- factor(pho_stat$curv_meth) ########## Is this necessary?


#Full data, with TPU
pho_leaf_tpu <- mutate(pho_both_tpu, leaf_unique = substring(ID, 1, 7))
pho_stat_tpu <- select(pho_leaf_tpu, 'curv_meth', 
                       'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU', 
                       'ID', 'leaf_unique', 'V_cmax_se', 'J_se', 'V_TPU_se')
pho_stat_tpu <- rename(pho_stat_tpu,
                       vcmax = Best_Vcmax_25C,
                       jmax = Best_Jmax_25C,
                       tpu = V_TPU,
                       leaf_id = ID) %>% 
    mutate(fit_type = "tpu")
# Describe factor levels: 0 is SS, 1 is DAT
pho_stat_tpu$curv_meth <- factor(pho_stat_tpu$curv_meth)

all_results <- rbind(pho_stat, pho_stat_tpu)




#No overshoot, no TPU
pho_nd_leaf <- mutate(pho_nd_both, leaf_unique = substring(ID, 1, 7))
pho_nd_stat <- select(pho_nd_leaf, 'curv_meth', 
                      'Best_Vcmax_25C', 'Best_Jmax_25C',
                      'V_TPU', 'ID', 
                      'leaf_unique', 'V_cmax_se', 'J_se', 'V_TPU_se')
pho_nd_stat <- rename(pho_nd_stat,
                      vcmax = Best_Vcmax_25C,
                      jmax = Best_Jmax_25C,
                      tpu = V_TPU,
                      leaf_id = ID) %>%
    mutate(fit_type = "no_tpu")
#Describe factor levels: 0 is SS, 1 is DAT
pho_nd_stat$curv_meth <- factor(pho_nd_stat$curv_meth)



#No overshoot, with TPU
pho_nd_leaf_tpu <- mutate(pho_nd_both_tpu, leaf_unique = substring(ID, 1, 7))
pho_nd_stat_tpu <- select(pho_nd_leaf_tpu, 'curv_meth', 'Best_Vcmax_25C', 
                          'Best_Jmax_25C', 'V_TPU',
                          'ID', 'leaf_unique', "V_cmax_se", "J_se", 'V_TPU_se')
pho_nd_stat_tpu <- rename(pho_nd_stat_tpu,
                          vcmax = Best_Vcmax_25C,
                          jmax = Best_Jmax_25C,
                          tpu = V_TPU,
                          leaf_id = ID) %>% 
    mutate(fit_type = "tpu")
#Describe factor levels: 0 is SS, 1 is DAT
pho_nd_stat_tpu$curv_meth <- factor(pho_nd_stat_tpu$curv_meth)



# Group data for further analysis -------------------------------------------

### Create summary table for WITH TPU results
tpu_results_grp <- pho_stat_tpu %>%
    group_by(curv_meth, leaf_unique) %>% 
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>% 
    mutate(fit_type = "tpu")

### Create summary table for WITHOUT TPU results
notpu_results_grp <- pho_stat %>%
    group_by(curv_meth, leaf_unique) %>% 
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>%
    mutate(fit_type = "no_tpu")

### Merge the datasets together
all_results2 <- rbind(tpu_results_grp, notpu_results_grp)
all_results2$fit_type <- factor(all_results2$fit_type)


# Creating a table with Vcmax and Jmax differences by species ----------------------------

########## Check how necessary this is
treename <- as.numeric(substring(all_results2$leaf_unique, 4, 5))

all_results2$treename <- treename



## Filter by TPU ---
all_res_tpu <- all_results2 %>% filter(fit_type == "tpu")


## Group by tree and calculate the mean and sd for Vcmax and Jmax
all_res_tpu_summ <- all_res_tpu %>% 
    group_by(curv_meth, leaf_unique) %>% 
    summarize(vcmax = mean(vcmax),
              jmax = mean(jmax),
              treename = treename)

dat_res_tpu_summ <- filter(all_res_tpu_summ, curv_meth == "DAT") %>% 
    rename(dat_vcmax = vcmax,
           dat_jmax = jmax)


# This is Steady State vcmax and jmax for each leaf
ss_res_tpu_summ <- filter(all_res_tpu_summ, curv_meth == "SS") %>% 
    rename(ss_vcmax = vcmax,
           ss_jmax = jmax)


## Group datasets, calculate diff and SE
all_diff_tpu2 <- full_join(dat_res_tpu_summ, ss_res_tpu_summ, 
                          by = c("leaf_unique", "treename")) %>% 
    mutate(vc_diff = ss_vcmax - dat_vcmax,
           j_diff = ss_jmax - dat_jmax) 

all_diff_tpu <- ungroup(all_diff_tpu2) %>%
    group_by(treename) %>% 
    mutate(vc_diff_se = sd(vc_diff)/sqrt(nrow(dat_res_tpu_summ))) %>% 
    mutate(j_diff_se = sd(j_diff)/sqrt(nrow(dat_res_tpu_summ))) %>% 
    summarize(vc_diff = mean(vc_diff),
              j_diff = mean(j_diff),
              vc_diff_se = mean(vc_diff_se),
              j_diff_se = mean(j_diff_se))


rel_can_pos <- c(11, 5.1, 6, 4, 2, 3, 1, 10, 12, 8, 7, 9, 5.2) ###### Again, check how necessary this is
all_diff_tpu$rel_can_pos <- rel_can_pos
all_diff_tpu <- all_diff_tpu[order(all_diff_tpu$rel_can_pos, decreasing = TRUE),]

codebook <- read_csv(here("5_Results/id_codebook.csv")) %>%
    arrange(desc(rel_can_pos)) %>%
    select(-c(overshoot, treeid, rel_can_pos))

all_diff_tpu_codes <- left_join(all_diff_tpu, codebook, by = "treename")

write.csv(all_diff_tpu_codes, here("5_Results/species_diffs_summary_TPU.csv"))



### Filter for no-TPU data ---
all_res_notpu <- all_results2 %>% 
    filter(fit_type == "no_tpu")


## Group by tree and calculate the mean and sd for Vcmax and Jmax
all_res_notpu_summ <- all_res_notpu %>% 
    group_by(curv_meth, leaf_unique) %>% 
    summarize(vcmax = mean(vcmax),
              jmax = mean(jmax),
              treename = treename)

dat_res_notpu_summ <- filter(all_res_notpu_summ, curv_meth == "DAT") %>% 
    rename(dat_vcmax = vcmax,
           dat_jmax = jmax)

ss_res_notpu_summ <- filter(all_res_notpu_summ, curv_meth == "SS") %>% 
    rename(ss_vcmax = vcmax,
           ss_jmax = jmax)

## Group datasets, calculate diff and SE
all_diff_notpu2 <- full_join(dat_res_notpu_summ, ss_res_notpu_summ, 
                             by = c("leaf_unique", "treename")) %>% 
    mutate(vc_diff = ss_vcmax - dat_vcmax) %>% 
    mutate(j_diff = ss_jmax - dat_jmax)

all_diff_notpu <- ungroup(all_diff_notpu2) %>% 
    group_by(treename) %>% 
    mutate(vc_diff_se = sd(vc_diff) / sqrt(nrow(dat_res_notpu_summ))) %>% 
    mutate(j_diff_se = sd(j_diff) / sqrt(nrow(dat_res_notpu_summ))) %>% 
    group_by(treename) %>% 
    summarize(vc_diff = mean(vc_diff),
              j_diff = mean(j_diff),
              vc_diff_se = mean(vc_diff_se),
              j_diff_se = mean(j_diff_se))


all_diff_notpu$rel_can_pos <- rel_can_pos
all_diff_notpu <- all_diff_notpu[order(all_diff_notpu$rel_can_pos, decreasing = TRUE),]
all_diff_notpu_codes <- left_join(all_diff_notpu, codebook, by = "treename")


write.csv(all_diff_notpu_codes, here("5_Results/species_diffs_summary_noTPU.csv"))


# Comparing mean differences --------------------------
all_res_dat_tpu <- all_results2 %>%
    filter(curv_meth == "DAT") %>%
    filter(fit_type == "tpu") %>%
    mutate(dat_vcmax = vcmax,
           dat_jmax = jmax) %>%
    select(-c(vcmax, jmax))
all_res_SS_tpu <- all_results2 %>%
    filter(curv_meth == "SS") %>%
    filter(fit_type == "tpu") %>%
    mutate(SS_vcmax = vcmax,
           SS_jmax = jmax) %>%
    select(-c(vcmax, jmax))

res_tpu_summ <- cbind(all_res_dat_tpu, all_res_SS_tpu) %>% 
    select(-c(1, 7, 8, 9, 10)) %>%
    rename(leaf_unique = leaf_unique...2,
           fit_type = fit_type...3,
           treename = treename...4)

res_tpu_summ$vc_diff <- res_tpu_summ$SS_vcmax - res_tpu_summ$dat_vcmax
res_tpu_summ$j_diff <- res_tpu_summ$SS_jmax - res_tpu_summ$dat_jmax

### WITHOUT TPU
all_res_dat_notpu <- all_results2 %>%
    filter(curv_meth == "DAT") %>%
    filter(fit_type == "no_tpu") %>%
    mutate(dat_vcmax = vcmax,
           dat_jmax = jmax) %>%
    select(-c(vcmax, jmax))
all_res_SS_notpu <- all_results2 %>%
    filter(curv_meth == "SS") %>%
    filter(fit_type == "no_tpu") %>%
    mutate(SS_vcmax = vcmax,
           SS_jmax = jmax) %>%
    select(-c(vcmax, jmax))

res_notpu_summ <- cbind(all_res_dat_notpu, all_res_SS_notpu) %>%
    select(-c(1, 7, 8, 9, 10)) %>%
    rename(leaf_unique = leaf_unique...2,
           fit_type = fit_type...3,
           treename = treename...4)

res_notpu_summ$vc_diff <- res_notpu_summ$SS_vcmax - res_notpu_summ$dat_vcmax
res_notpu_summ$j_diff <- res_notpu_summ$SS_jmax - res_notpu_summ$dat_jmax


# mean differences on leaf basis
mean(all_diff_tpu2$vc_diff)
mean(all_diff_tpu2$j_diff)
mean(all_diff_notpu2$vc_diff)
mean(all_diff_notpu2$j_diff)

# mean differences on tree basis
mean(all_diff_tpu_codes$vc_diff)
mean(all_diff_tpu_codes$j_diff)
mean(all_diff_notpu_codes$vc_diff)
mean(all_diff_notpu_codes$j_diff)










#No-overshoot data, no TPU
grp_pho_nd_dat <- pho_nd_stat %>%
    filter(curv_meth == "DAT") %>%
    group_by(leaf_unique) %>%
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>% 
    mutate(curv_meth = "DAT",
           fit_type = "no_tpu")

grp_pho_nd_SS <- pho_nd_stat %>%
    filter(curv_meth == "SS") %>% 
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
    mutate(curv_meth = "SS",
           fit_type = "no_tpu")

grp_pho_nd_all <- rbind(grp_pho_nd_dat, grp_pho_nd_SS)


#No-overshoot data, with TPU
grp_pho_nd_dat_tpu <- pho_nd_stat_tpu %>%
    filter(curv_meth == "DAT") %>%
    group_by(leaf_unique) %>%
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>% 
    mutate(curv_meth = "DAT",
           fit_type = "tpu")

grp_pho_nd_SS_tpu <- pho_nd_stat_tpu %>%
    filter(curv_meth == "SS") %>% 
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
    mutate(curv_meth = "SS",
           fit_type = "tpu")

grp_pho_nd_all_tpu <- rbind(grp_pho_nd_dat_tpu, grp_pho_nd_SS_tpu)


nd_complete <- rbind(grp_pho_nd_all, grp_pho_nd_all_tpu)
nd_complete$curv_meth <- factor(nd_complete$curv_meth)
nd_complete$fit_type <- factor(nd_complete$fit_type)

#Just for the TPU analysis

grp_narm_pho_tpu <- pho_stat_tpu %>%
    na.omit(.$tpu) #This leaves 6 SS and 22 DAT

grp_tpu_6SS <- grp_narm_pho_tpu%>%
    filter(curv_meth == "SS") %>% 
    group_by(leaf_unique) %>% 
    summarize(vcmax = vcmax,
              jmax = jmax,
              tpu = tpu,
              curv_meth = "SS")

grp_tpu_6dat <- grp_narm_pho_tpu %>%
    filter(curv_meth == "DAT") %>%
    filter(leaf_unique %in% c("K6706L1", "K6707L1", "K6707L2", "K6709L6", "K6714L1", "K6714L2")) %>%
    group_by(leaf_unique) %>% 
    summarize(vcmax = mean(vcmax),
              jmax = mean(jmax),
              tpu = mean(tpu),
              curv_meth = "DAT")

tpu_just6_all <- rbind(grp_tpu_6SS, grp_tpu_6dat)
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


all_results2 %>% group_by(fit_type, curv_meth) %>% get_summary_stats(vcmax, jmax, show = c("mean", "median", "sd", "min", "max"))


#No TPU

all_results2 %>% filter(fit_type == "no_tpu") %>%
ggplot(aes(x=vcmax)) + 
    geom_histogram()

all_results2 %>% filter(fit_type == "no_tpu") %>%
ggplot(aes(x=jmax)) + 
    geom_histogram()


#Levene's for homogeneity of variance
all_results2 %>% filter(fit_type == "no_tpu") %>% leveneTest(vcmax ~ curv_meth, data = .)
all_results2 %>% filter(fit_type == "no_tpu") %>% leveneTest(jmax ~ curv_meth, data = .)
# Both non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances

#Shapiro-Wilk test for normality
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax))
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax[curv_meth == "DAT"]))
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax[curv_meth == "SS"]))

all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax))
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "SS"]))
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
all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax[curv_meth == "SS"]))

all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax))
all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "SS"]))
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


nd_complete %>% group_by(fit_type, curv_meth) %>% get_summary_stats(vcmax, jmax, show = c("mean", "median", "sd", "min", "max"))


#Levene's test
nd_complete %>% filter(fit_type == "tpu") %>% leveneTest(vcmax ~ curv_meth, data = .)
nd_complete %>% filter(fit_type == "tpu") %>% leveneTest(jmax ~ curv_meth, data = .)
nd_complete %>% filter(fit_type == "no_tpu") %>% leveneTest(vcmax ~ curv_meth, data = .)
nd_complete %>% filter(fit_type == "no_tpu") %>% leveneTest(jmax ~ curv_meth, data = .)
#Met assumption of homogeneity of variances

#Shapiro Wilk Test
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax[curv_meth == "SS"]))

nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "SS"]))

nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax[curv_meth == "SS"]))

nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "SS"]))
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
    mutate(differences = SS - DAT)

all_res_wide_vcmax %>% filter(fit_type == "no_tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#That should be fine

all_res_wide_vcmax %>% filter(fit_type == "tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#okay

#Jmax
all_res_wide_jmax <- all_results2 %>%
    dcast(., leaf_unique + fit_type ~ curv_meth, value.var="jmax") %>%
    mutate(differences = SS - DAT)

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
    mutate(differences = SS - DAT)

nd_comp_wide_vcmax %>% filter(fit_type == "no_tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#That should be fine

nd_comp_wide_vcmax %>% filter(fit_type == "tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#okay

#Jmax
nd_comp_wide_jmax <- nd_complete %>%
    dcast(., leaf_unique + fit_type ~ curv_meth, value.var="jmax") %>%
    mutate(differences = SS - DAT)

nd_comp_wide_jmax %>% filter(fit_type == "no_tpu") %>%
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#This is fine

nd_comp_wide_jmax %>% filter(fit_type == "tpu") %>%
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#This is fine


#So it's just the Jmax where we don't totally meet the Wilcoxon assumptions. We will run both Sign and Wilcoxon for both of these, and make a note of this.



#Stat analysis, clean ------------------------------------------

#all results, Wilcoxon signed rank test on paired samples (and a few sign tests)

#Vcmax
all_results2 %>%
    group_by(fit_type) %>%
    wilcox_test(data =., vcmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    add_significance()

all_results2 %>%
    group_by(fit_type) %>% wilcox_effsize(data = ., vcmax ~ curv_meth, paired = TRUE)


# #try without Tachi to see what changes
# all_results2 %>%
#     group_by(fit_type) %>%
#     filter(leaf_unique != "K6707L1" & leaf_unique != "K6707L2") %>% 
#     wilcox_test(data =., vcmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
#     add_significance()
# 
# 
# all_results2 %>%
#     group_by(fit_type) %>% filter(leaf_unique != "K6707L1" & leaf_unique != "K6707L2") %>% wilcox_effsize(data = ., vcmax ~ curv_meth, paired = TRUE)
# #Not significant without Tachi!


all_results2 %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>% 
    add_significance()

all_results2 %>%
    group_by(curv_meth) %>% wilcox_effsize(data = ., vcmax ~ fit_type, paired = TRUE)


all_results2 %>%
    group_by(curv_meth) %>%
    filter(leaf_unique != "K6707L1" & leaf_unique != "K6707L2") %>% 
    wilcox_test(data =., vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>% 
    add_significance()





#Jmax
all_results2 %>%
    group_by(fit_type) %>%
    wilcox_test(data =., jmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    add_significance()

#Note we're running the sign test here in addition!
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






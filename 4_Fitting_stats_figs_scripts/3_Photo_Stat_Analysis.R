# Computing statistical analysis on the A/Ci curve fits from 2_Tapajos_Fit_ACi.R
# Using the 'photosynthesis' package

library(tidyverse) 
library(ggpubr) 
library(car)
library(rstatix)
library(reshape2) 
library(here) 

###

# Load and Subset the datasets -------------------------------------------

### Tree ids
ids <- read.csv(here("3_Clean_data/id_codebook.csv")) #%>% 
    #rename(treeid = ï..treeid)

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


pho_dat_tpu <- read.csv(file = here("5_Results/DAT_photo_pars_crct_TPU.csv"),
                        sep = ",", 
                        header = TRUE, na.strings = 1000) %>%  
    ## TPU values at 1000 are coded as NA
    subset(ID != "K6709L3")



pho_SS_tpu <- read.csv(file = here("5_Results/SS_photo_pars_crct_TPU.csv"),
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
pho_nd_both_tpu <- rbind(pho_nd_tpu, pho_SS_tpu)


###

# Initial Processing/Grouping ----------------------------------

# Full data, no TPU

pho_leaf <- mutate(pho_both, leaf_unique = substring(ID, 1, 7))
pho_stat <- select(pho_leaf, 
                   'curv_meth', 'Best_Vcmax_25C', 'Best_Jmax_25C',
                   'V_TPU', 'ID', 'treeid',
                   'leaf_unique', "V_cmax_se", "J_se", "V_TPU_se") %>% 
    rename(.,
           vcmax = Best_Vcmax_25C,
           jmax = Best_Jmax_25C,
           tpu = V_TPU,
           leaf_id = ID) %>% 
    mutate(fit_type = "no_tpu")


write.csv(pho_stat, file = here("5_Results/pho_stat.csv"))

# Full data, with TPU
pho_leaf_tpu <- mutate(pho_both_tpu, leaf_unique = substring(ID, 1, 7))
pho_stat_tpu <- select(pho_leaf_tpu,
                       'curv_meth', 'Best_Vcmax_25C', 'Best_Jmax_25C', 
                       'V_TPU', 'ID', 'treeid',
                       'leaf_unique', 'V_cmax_se', 'J_se', 'V_TPU_se') %>% 
    rename(.,
           vcmax = Best_Vcmax_25C,
           jmax = Best_Jmax_25C,
           tpu = V_TPU,
           leaf_id = ID) %>% 
    mutate(fit_type = "tpu")

write.csv(pho_stat_tpu, file = here("5_Results/pho_stat_tpu.csv"))



### Bind the results together
all_results <- rbind(pho_stat, pho_stat_tpu)




# No overshoot, no TPU
pho_nd_leaf <- mutate(pho_nd_both, leaf_unique = substring(ID, 1, 7))
pho_nd_stat <- select(pho_nd_leaf, 'curv_meth', 
                      'Best_Vcmax_25C', 'Best_Jmax_25C',
                      'V_TPU', 'ID', 'treeid',
                      'leaf_unique', 'V_cmax_se', 'J_se', 'V_TPU_se') %>% 
    rename(.,
           vcmax = Best_Vcmax_25C,
           jmax = Best_Jmax_25C,
           tpu = V_TPU,
           leaf_id = ID) %>%
    mutate(fit_type = "no_tpu")

write.csv(pho_nd_stat, file = here("5_Results/pho_nd_stat.csv"))



# No overshoot, with TPU
pho_nd_leaf_tpu <- mutate(pho_nd_both_tpu, leaf_unique = substring(ID, 1, 7))
pho_nd_stat_tpu <- select(pho_nd_leaf_tpu, 'curv_meth', 'Best_Vcmax_25C', 
                          'Best_Jmax_25C', 'V_TPU', 'treeid',
                          'ID', 'leaf_unique', "V_cmax_se", "J_se", 'V_TPU_se') %>% 
    rename(.,
           vcmax = Best_Vcmax_25C,
           jmax = Best_Jmax_25C,
           tpu = V_TPU,
           leaf_id = ID) %>% 
    mutate(fit_type = "tpu")

write.csv(pho_nd_stat_tpu, file = here("5_Results/pho_nd_stat_tpu.csv"))


###

# Calculate mean Vcmax and Jmax values: All Data ------------------------------------

### Create summary table for WITH TPU results
tpu_results_grp <- pho_stat_tpu %>%
    group_by(curv_meth, leaf_unique) %>% 
    summarise(vcmax = mean(vcmax),
            jmax = mean(jmax)) %>% 
    mutate(fit_type = "tpu",
           treeid = substring(leaf_unique, 1, 5))


### Create summary table for WITHOUT TPU results
notpu_results_grp <- pho_stat %>%
    group_by(curv_meth, leaf_unique) %>% 
    summarise(vcmax = mean(vcmax),
            jmax = mean(jmax)) %>%
    mutate(fit_type = "no_tpu",
           treeid = substring(leaf_unique, 1, 5))

### Merge the summary tables together
all_avg_lf_res <- rbind(tpu_results_grp, notpu_results_grp)
all_avg_lf_res$curv_meth <- factor(all_avg_lf_res$curv_meth)
all_avg_lf_res$fit_type <- factor(all_avg_lf_res$fit_type)



# WITH TPU results

### Pull out the DAT and SS results to concatenate width-wise to calculate differences

## Pull out the DAT results
dat_res_tpu_summ <- filter(tpu_results_grp, curv_meth == "DAT") %>% 
    rename(dat_vcmax = vcmax,
           dat_jmax = jmax)

## Pull out the SS results
ss_res_tpu_summ <- filter(tpu_results_grp, curv_meth == "SS") %>% 
    rename(ss_vcmax = vcmax,
           ss_jmax = jmax)


## Join datasets, calculate diff and SE
diff_tpu_lf <- full_join(dat_res_tpu_summ, ss_res_tpu_summ, 
                          by = c("leaf_unique", "treeid", "fit_type")) %>% 
    mutate(vc_diff = ss_vcmax - dat_vcmax,
           j_diff = ss_jmax - dat_jmax) 


diff_tpu_tree <- ungroup(diff_tpu_lf) %>%
    group_by(treeid) %>% 
    mutate(vc_diff_se = sd(vc_diff)/sqrt(nrow(dat_res_tpu_summ))) %>% 
    mutate(j_diff_se = sd(j_diff)/sqrt(nrow(dat_res_tpu_summ))) %>% 
    summarize(vc_diff = mean(vc_diff),
              j_diff = mean(j_diff),
              vc_diff_se = mean(vc_diff_se),
              j_diff_se = mean(j_diff_se))


### Adding back in the relative canopy position and sorting by height
rel_can_pos <- select(ids, treeid, rel_can_pos)
diff_tpu_tree <- left_join(diff_tpu_tree, rel_can_pos, by = "treeid")
diff_tpu_tree <- diff_tpu_tree[order(diff_tpu_tree$rel_can_pos, decreasing = TRUE),]

diff_tpu_tree_codes <- left_join(diff_tpu_tree, ids, by = "treeid")


### Write .csv file of WITH TPU differences
write.csv(diff_tpu_tree_codes, here("5_Results/tree_diffs_summary_TPU.csv"))
write.csv(diff_tpu_lf, here("3_Clean_data/lf_diffs_summ_TPU.csv"))



# WITHOUT TPU

### Pull out the DAT and SS results to concatenate width-wise to calculate differences

### Pull out the DAT results
dat_res_notpu_summ <- filter(notpu_results_grp, curv_meth == "DAT") %>% 
    rename(dat_vcmax = vcmax,
           dat_jmax = jmax)

### Pull out the SS results
ss_res_notpu_summ <- filter(notpu_results_grp, curv_meth == "SS") %>% 
    rename(ss_vcmax = vcmax,
           ss_jmax = jmax)


### Join datasets, calculate diff and SE
diff_notpu_lf <- full_join(dat_res_notpu_summ, ss_res_notpu_summ, 
                             by = c("leaf_unique", "treeid", "fit_type")) %>% 
    mutate(vc_diff = ss_vcmax - dat_vcmax) %>% 
    mutate(j_diff = ss_jmax - dat_jmax)

diff_notpu_tree <- ungroup(diff_notpu_lf) %>% 
    group_by(treeid) %>% 
    mutate(vc_diff_se = sd(vc_diff) / sqrt(nrow(dat_res_notpu_summ))) %>% 
    mutate(j_diff_se = sd(j_diff) / sqrt(nrow(dat_res_notpu_summ))) %>% 
    group_by(treeid) %>% 
    summarize(vc_diff = mean(vc_diff),
              j_diff = mean(j_diff),
              vc_diff_se = mean(vc_diff_se),
              j_diff_se = mean(j_diff_se))



### Adding back in the relative canopy position and sorting by height
rel_can_pos <- select(ids, treeid, rel_can_pos)
diff_notpu_tree <- left_join(diff_notpu_tree, rel_can_pos, by = "treeid")
diff_notpu_tree <- diff_notpu_tree[order(diff_notpu_tree$rel_can_pos, decreasing = TRUE),]
diff_notpu_tree_codes <- left_join(diff_notpu_tree, ids, by = "treeid")


### Write .csv file of WITHOUT TPU differences
write.csv(diff_notpu_tree_codes, here("5_Results/tree_diffs_summary_noTPU.csv"))
write.csv(diff_notpu_lf, here("5_Results/lf_diffs_summ_noTPU.csv"))


# Calculate mean Vcmax and Jmax values: No Overshoot Data -----------------------------

### Create summary table for WITH TPU results
tpu_nd_results_grp <- pho_nd_stat_tpu %>%
    group_by(curv_meth, leaf_unique) %>% 
    reframe(vcmax = mean(vcmax),
            jmax = mean(jmax)) %>% 
    mutate(fit_type = "tpu",
           treeid = substring(leaf_unique, 1, 5))


### Create summary table for WITHOUT TPU results
notpu_nd_results_grp <- pho_nd_stat %>%
    group_by(curv_meth, leaf_unique) %>% 
    reframe(vcmax = mean(vcmax),
            jmax = mean(jmax)) %>%
    mutate(fit_type = "no_tpu",
           treeid = substring(leaf_unique, 1, 5))

### Merge the summary tables together
nd_avg_lf_res <- rbind(tpu_nd_results_grp, notpu_nd_results_grp)



# WITH TPU results

### Pull out the DAT and SS results to concatenate width-wise to calculate differences

## Pull out the DAT results
dat_nd_res_tpu_summ <- filter(tpu_nd_results_grp, curv_meth == "DAT") %>% 
    rename(dat_vcmax = vcmax,
           dat_jmax = jmax)

## Pull out the SS results
ss_nd_res_tpu_summ <- filter(tpu_nd_results_grp, curv_meth == "SS") %>% 
    rename(ss_vcmax = vcmax,
           ss_jmax = jmax)


## Join datasets, calculate diff and SE
nd_diff_tpu_lf <- full_join(dat_nd_res_tpu_summ, ss_nd_res_tpu_summ, 
                           by = c("leaf_unique", "treeid", "fit_type")) %>% 
    mutate(vc_diff = ss_vcmax - dat_vcmax,
           j_diff = ss_jmax - dat_jmax) 


nd_diff_tpu_tree <- ungroup(nd_diff_tpu_lf) %>%
    group_by(treeid) %>% 
    mutate(vc_diff_se = sd(vc_diff)/sqrt(nrow(dat_nd_res_tpu_summ))) %>% 
    mutate(j_diff_se = sd(j_diff)/sqrt(nrow(dat_nd_res_tpu_summ))) %>% 
    summarize(vc_diff = mean(vc_diff),
              j_diff = mean(j_diff),
              vc_diff_se = mean(vc_diff_se),
              j_diff_se = mean(j_diff_se))


### Adding back in the relative canopy position and sorting by height
rel_can_pos <- select(ids, treeid, rel_can_pos)
nd_diff_tpu_tree <- left_join(nd_diff_tpu_tree, rel_can_pos, by = "treeid")
nd_diff_tpu_tree <- nd_diff_tpu_tree[order(nd_diff_tpu_tree$rel_can_pos, 
                                           decreasing = TRUE),]

nd_diff_tpu_tree_codes <- left_join(nd_diff_tpu_tree, ids, by = "treeid")


### Write .csv file of WITH TPU differences
write.csv(nd_diff_tpu_tree_codes, here("5_Results/tree_nOS_diffs_summary_TPU.csv"))




# WITHOUT TPU

### Pull out the DAT and SS results to concatenate width-wise to calculate differences

### Pull out the DAT results
dat_nd_res_notpu_summ <- filter(notpu_nd_results_grp, curv_meth == "DAT") %>% 
    rename(dat_vcmax = vcmax,
           dat_jmax = jmax)

### Pull out the SS results
ss_nd_res_notpu_summ <- filter(notpu_nd_results_grp, curv_meth == "SS") %>% 
    rename(ss_vcmax = vcmax,
           ss_jmax = jmax)


### Join datasets, calculate diff and SE
nd_diff_notpu_lf <- full_join(dat_nd_res_notpu_summ, ss_nd_res_notpu_summ, 
                             by = c("leaf_unique", "treeid", "fit_type")) %>% 
    mutate(vc_diff = ss_vcmax - dat_vcmax) %>% 
    mutate(j_diff = ss_jmax - dat_jmax)

nd_diff_notpu_tree <- ungroup(nd_diff_notpu_lf) %>% 
    group_by(treeid) %>% 
    mutate(vc_diff_se = sd(vc_diff) / sqrt(nrow(dat_nd_res_notpu_summ))) %>% 
    mutate(j_diff_se = sd(j_diff) / sqrt(nrow(dat_nd_res_notpu_summ))) %>% 
    group_by(treeid) %>% 
    summarize(vc_diff = mean(vc_diff),
              j_diff = mean(j_diff),
              vc_diff_se = mean(vc_diff_se),
              j_diff_se = mean(j_diff_se))




### Adding back in the relative canopy position and sorting by height
rel_can_pos <- select(ids, treeid, rel_can_pos)
nd_diff_notpu_tree <- left_join(nd_diff_notpu_tree, rel_can_pos, by = "treeid")
nd_diff_notpu_tree <- nd_diff_notpu_tree[order(nd_diff_notpu_tree$rel_can_pos, decreasing = TRUE),]
nd_diff_notpu_tree_codes <- left_join(nd_diff_notpu_tree, ids, by = "treeid")


### Write .csv file of WITHOUT TPU differences
write.csv(nd_diff_notpu_tree_codes, here("5_Results/tree_nOS_diffs_summary_noTPU.csv"))




# Taking out the SS leaves which had OS DAT pairs and calculating the mean

## WITHOUT TPU 
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
    reframe(vcmax = mean(vcmax),
            jmax = mean(jmax), 
            curv_meth = curv_meth, 
            fit_type = fit_type,
            treeid = treeid) 

### Combine with the DAT noOS results
grp_pho_nd_all <- dat_nd_res_notpu_summ %>% 
    rename(vcmax = dat_vcmax, jmax = dat_jmax) %>% 
    rbind(., grp_pho_nd_SS)


# Taking out the SS leaves which had OS DAT pairs and calculating the mean

# WITH TPU
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
            jmax = mean(jmax),
            treeid = substring(leaf_unique, 1, 5)) %>% 
    mutate(curv_meth = "SS", 
           fit_type = "tpu")

### Combine with the DAT noOS results
grp_pho_nd_all_tpu <- dat_nd_res_tpu_summ %>% 
    rename(vcmax = dat_vcmax, jmax = dat_jmax) %>% 
    rbind(., grp_pho_nd_SS_tpu)



nd_complete <- rbind(grp_pho_nd_all, grp_pho_nd_all_tpu)
nd_complete$curv_meth <- factor(nd_complete$curv_meth)
nd_complete$fit_type <- factor(nd_complete$fit_type)






# Subset curves for the SS TPU v DAT TPU analysis -------------------------

### Pulling out the pairs of curves for which TPU was fit in both
### This leaves 6 SS and 22 DAT

grp_tpu_6SS <- pho_stat_tpu %>% 
    na.omit(tpu) %>%
    filter(curv_meth == "SS") %>% 
    group_by(leaf_unique) %>% 
    reframe(vcmax = vcmax,
            jmax = jmax,
            tpu = tpu,
            curv_meth = "SS")

grp_tpu_6dat <- pho_stat_tpu %>% 
    na.omit(tpu) %>%
    filter(curv_meth == "DAT") %>%
    filter(leaf_unique %in% c("K6706L1", 
                              "K6707L1", 
                              "K6707L2", 
                              "K6709L6",
                              "K6714L1",
                              "K6714L2")) %>%
    group_by(leaf_unique) %>% 
    reframe(vcmax = mean(vcmax),
            jmax = mean(jmax),
            tpu = mean(tpu),
            curv_meth = "DAT")  


tpu_just6_all <- rbind(grp_tpu_6SS, grp_tpu_6dat)
tpu_just6_all$curv_meth <- factor(tpu_just6_all$curv_meth)


# Produce Summary statistics of results ----------------------

### All data
all_avg_lf_res %>% 
    group_by(fit_type, curv_meth) %>%
    get_summary_stats(vcmax, jmax, show = c("mean", "median", "sd", "min", "max")) %>% as.data.frame()


### No-overshoot data
nd_complete %>% 
    group_by(fit_type, curv_meth) %>% 
    get_summary_stats(vcmax, jmax, show = c("mean", "median", "sd", "min", "max")) %>% as.data.frame()


### Just the TPU data
tpu_just6_all %>% 
    group_by(curv_meth)%>%
    get_summary_stats(vcmax, jmax, show = c("mean", "median", "sd", "min", "max"))



# Exploring random effect models ---------------------------------
library(lmerTest)
library(lme4)
library(performance)

# Mixed model with random effect for intercept
# Our leaves are nested within trees. Trees will be the random effect.
# These models have the difference between DAT/SS vcmax or jmax as the independent variable.

#Set up datasets ------------------------
#Convert to factor. n = 27 pairs
diff_notpu_lf$treeid <- as.factor(diff_notpu_lf$treeid)
diff_tpu_lf$treeid <- as.factor(diff_tpu_lf$treeid)

#Set up a dataframe to exclude K6709L6 (for testing; the one with the big diff). n = 25 pairs
diff_notpu_nol6_lf <- diff_notpu_lf %>% filter(leaf_unique != 'K6709L6')
diff_tpu_nol6_lf <- diff_tpu_lf %>% filter(leaf_unique != 'K6709L6')

# correct nd_diff datasets. n = 19 pairs
nd_diff_notpu_lf$treeid <- as.factor(nd_diff_notpu_lf$treeid)
nd_diff_tpu_lf$treeid <- as.factor(nd_diff_tpu_lf$treeid)

nd_diff_notpu_lf <- nd_diff_notpu_lf %>% filter(!is.na(vc_diff))
nd_diff_tpu_lf <- nd_diff_tpu_lf %>% filter(!is.na(vc_diff))


# Define models ----------------------
#Models with full dataset

## Vcmax without TPU
mod_notpu_v <- lmer(vc_diff ~ 1 + (1|treeid),
                    data = diff_notpu_lf)
no_tpu_res <- filter(all_results, fit_type == "no_tpu")
mod_notpu_v2 <- lmer(vcmax ~ curv_meth + (1|treeid) + (1|treeid:leaf_unique), data = no_tpu_res)

summary(mod_notpu_v)
### Slightly sig diff, estimated diff 2.1 +- 0.92

## Jmax without TPU
mod_notpu_j <- lmer(j_diff ~ 1 + (1|treeid),
                    data = diff_notpu_lf)

summary(mod_notpu_j)
### Sig diff, estimated diff 10.8 +- 3.0 


## Vcmax with TPU
mod_tpu_v <- lmer(vc_diff ~ 1 + (1|treeid),
                    data = diff_tpu_lf)
tpu_res <- filter(all_results, fit_type == "tpu")
mod_tpu_v2 <- lmer(vcmax ~ curv_meth + (1|treeid) + (1|treeid:leaf_unique), data = tpu_res)
#This leads to a singularity


summary(mod_tpu_v)
### Not sig, estimate diff 1.00 +- 0.67; gives singularity

#Trying to avoid the singularity
# nlme_tpu_v <- nlmer(vc_diff ~ 1 + (1|treeid),
#                   data = diff_tpu_lf)
# #Still no luck...
# nlme_tpu_1 <- nlme::nlme(vc_diff ~ 1,
#                          random = treeid ~ 1,
#                          data = diff_tpu_lf)



## Jmax with TPU
mod_tpu_j <- lmer(j_diff ~ 1 + (1|treeid),
                  data = diff_tpu_lf)

summary(mod_tpu_j)
### Sig, estimated differences 7.9 +- 2.46



# Models without K6709L6 -- Testing the effects of this potentially influential point.

## Vcmax no TPU
mod_notpu_nol6_v <- lmer(vc_diff ~ 1 + (1|treeid),
                    data = diff_notpu_nol6_lf)
summary(mod_notpu_nol6_v)
### almost sig, estimated diff 1.9 +- 0.06

## Jmax no TPU
mod_notpu_nol6_j <- lmer(j_diff ~ 1 + (1|treeid),
                    data = diff_notpu_nol6_lf)
summary(mod_notpu_nol6_j)
### sig, estimated diff 10.2 +-2.8

## Vcmax with TPU
mod_tpu_nol6_v <- lmer(vc_diff ~ 1 + (1|treeid),
                  data = diff_tpu_nol6_lf) #This leads to a singularity
summary(mod_tpu_nol6_v)
### sig, estimated diff 0.5 +- 0.43


## Jmax with TPU
mod_tpu_nol6_j <- lmer(j_diff ~ 1 + (1|treeid),
                  data = diff_tpu_nol6_lf)
summary(mod_tpu_nol6_j)
### sig, estimated diff 6.1 +- 1.2


# Models with no-overshoot data. Just a quick check
mod_nd_notpu_v <- lmer(vc_diff ~ 1 + (1|treeid),
                    data = nd_diff_notpu_lf) #This leads to a singularity
mod_nd_notpu_j <- lmer(j_diff ~ 1 + (1|treeid),
                    data = nd_diff_notpu_lf)
summary(mod_nd_notpu_v)
summary(mod_nd_notpu_j)


mod_nd_tpu_v <- lmer(vc_diff ~ 1 + (1|treeid),
                  data = nd_diff_tpu_lf)
mod_nd_tpu_j <- lmer(j_diff ~ 1 + (1|treeid),
                  data = nd_diff_tpu_lf)
summary(mod_nd_tpu_v)
summary(mod_nd_tpu_j)


# Residuals for the full dataset

plot(residuals(mod_notpu_v) ~ diff_notpu_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red")

plot(residuals(mod_notpu_j) ~ diff_notpu_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red")

plot(residuals(mod_tpu_v) ~ diff_tpu_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red")

plot(residuals(mod_tpu_j) ~ diff_tpu_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red")


#Summary of the models ------------------
#Complete dataset summary
summary(mod_notpu_v)
summary(mod_notpu_j)

summary(mod_tpu_v)
summary(mod_tpu_j)

#No Tachi dataset summary
summary(mod_notpu_nol6_v)
summary(mod_notpu_nol6_j)

summary(mod_tpu_nol6_v)
summary(mod_tpu_nol6_j)

#No overshoot dataset summary
summary(mod_nd_notpu_v)
summary(mod_nd_notpu_j)

summary(mod_nd_tpu_v)
summary(mod_nd_tpu_j)

# Pulls out key coefficients of interest
set.seed(304)
mod_coefs <- list(
    "mod_tpu_v",
    "mod_tpu_j",
    "mod_notpu_v",
    "mod_notpu_j",
    "mod_nd_tpu_v",
    "mod_nd_tpu_j",
    "mod_nd_notpu_v",
    "mod_nd_notpu_j",
    "mod_tpu_nol6_v",
    "mod_tpu_nol6_j",
    "mod_notpu_nol6_v",
    "mod_notpu_nol6_j"
) %>%
    map_dfr(~ {
        model_name <- .x
        model <- get(model_name)
        coeff <- lmerTest:::get_coefmat(model) %>% 
            as.data.frame() %>% 
            rownames_to_column()
         confints <- confint(model, method = 'boot', oldNames = FALSE) %>%
             as.data.frame() %>%
             rownames_to_column() %>% 
             filter(rowname == '(Intercept)') %>% #These are the fixed effect CIs
             select(-rowname)
        icc_value <- icc(model) %>%
            as.data.frame() %>% 
            rownames_to_column()
        std.d.ranef <- as.data.frame(VarCorr(model)) %>%
            rename(Ranef.Var = vcov, Ranef.StdDev = sdcor) %>% 
            rownames_to_column() %>% 
            filter(rowname == 1) %>% 
            select(-c(var1,var2, rowname))
        
        tibble(
            model = model_name,
            coeff = list(coeff),
            conf = confints,
            icc = icc_value[2],
            stdranef = std.d.ranef
        )
    }) %>%
    unnest_wider(coeff)

mod_coefs$Estimate <- round(mod_coefs$Estimate, 1)
mod_coefs$'Std. Error' <- round(mod_coefs$'Std. Error', 1)
mod_coefs$df <- round(mod_coefs$df, 1)
mod_coefs$'t value' <- round(mod_coefs$'t value', 2)
mod_coefs$'Pr(>|t|)' <- round(mod_coefs$'Pr(>|t|)', 4)
mod_coefs$conf$`2.5 %` <- round(mod_coefs$conf$`2.5 %`, 2)
mod_coefs$conf$`97.5 %` <- round(mod_coefs$conf$`97.5 %`, 2)
mod_coefs$icc$ICC_adjusted <- round(mod_coefs$icc$ICC_adjusted, 2)
mod_coefs$stdranef$Ranef.Var <- round(mod_coefs$stdranef$Ranef.Var, 1)
mod_coefs$stdranef$Ranef.StdDev <- round(mod_coefs$stdranef$Ranef.StdDev, 1)

print(mod_coefs)


#
write.csv(x = mod_coefs, 
          file = here("5_Results/model_output.csv"),
          row.names = FALSE)

#Compute bootstrapped confidence intervals
# This allows for additional confidence intervals to the ones pulled out in the code above.

# Complete dataset
#sd_(Intercept)|treeid: CI for random intercept
# sigma: CI for random residuals
#(Intercept): CI for the fixed effects

set.seed(304)
confint(mod_notpu_v, method = 'boot', oldNames = FALSE)
set.seed(304)
confint(mod_notpu_j, method = 'boot', oldNames = FALSE)
set.seed(304)
confint(mod_tpu_v, method = 'boot', oldNames = FALSE)
set.seed(304)
confint(mod_tpu_j, method = 'boot', oldNames = FALSE)

#No Tachi dataset
set.seed(304)
confint(mod_notpu_nol6_v, method = 'boot', oldNames = FALSE)
set.seed(304)
confint(mod_notpu_nol6_j, method = 'boot', oldNames = FALSE)
set.seed(304)
confint(mod_tpu_nol6_v, method = 'boot', oldNames = FALSE)
set.seed(304)
confint(mod_tpu_nol6_j, method = 'boot', oldNames = FALSE)

# No overshoot dataset
set.seed(304)
confint(mod_nd_notpu_v, method = 'boot', oldNames = FALSE)
set.seed(304)
confint(mod_nd_notpu_j, method = 'boot', oldNames = FALSE)
set.seed(304)
confint(mod_nd_tpu_v, method = 'boot', oldNames = FALSE)
set.seed(304)
confint(mod_nd_tpu_j, method = 'boot', oldNames = FALSE)


# Look at random effect coefficients
# Full dataset
ranef(mod_notpu_v)
ranef(mod_notpu_j)

ranef(mod_tpu_v) #This doesn't work because it's a singular model
ranef(mod_tpu_j)

# No K6709L6 dataset
ranef(mod_notpu_nol6_v)
ranef(mod_notpu_nol6_j)

ranef(mod_tpu_nol6_v) #This doesn't work because it's a singular model
ranef(mod_tpu_nol6_j)

# No overshoot dataset
ranef(mod_nd_notpu_v) #This doesn't work because it's a singular model
ranef(mod_nd_notpu_j)

ranef(mod_nd_tpu_v)
ranef(mod_nd_tpu_j)



# Checking Assumptions of paired t-tests -----------------------------------------

# Levene's for homogeneity of variance is not necessary because the paired t-test cares only about the distribution of the difference between the two variables. The variances of the two variables separately is irrelevant.

# Checking the symmetric distribution of the differences -----------------------------------------

# Includes Shapiro Wilk normality test for paired differences; no overshoot, WITHOUT TPU

gghistogram(diff_notpu_lf, x = "vc_diff", bins = 10, add_density = TRUE) + geom_vline(xintercept = median(diff_notpu_lf$vc_diff), color = 'red')
#That should be fine

diff_notpu_lf %>%
    with(., shapiro.test(vc_diff))

gghistogram(diff_tpu_lf, x = "vc_diff", bins = 10, add_density = TRUE) + geom_vline(xintercept = median(diff_tpu_lf$vc_diff), color = 'red')

diff_tpu_lf %>%
    with(., shapiro.test(vc_diff))

#Vcmax differences are approximately normally distributed.

### Jmax
gghistogram(diff_notpu_lf, x = "j_diff", bins = 10, add_density = TRUE) + geom_vline(xintercept = median(diff_notpu_lf$j_diff), color = 'red')

diff_notpu_lf %>%
    with(., shapiro.test(j_diff))

gghistogram(diff_tpu_lf, x = "j_diff", bins = 10, add_density = TRUE) + geom_vline(xintercept = median(diff_tpu_lf$j_diff), color = 'red')

diff_tpu_lf %>%
    with(., shapiro.test(j_diff))

#### Jmax not great! Differences not normally distributed; far worse than Vcmax differences.

#Testing log transformation
diff_notpu_lf$log_j_diff <- log(diff_notpu_lf$j_diff)
diff_tpu_lf$log_j_diff <- log(diff_tpu_lf$j_diff)



# No Overshoot data

### Vcmax
gghistogram(nd_diff_notpu_lf, x = "vc_diff", bins = 10, add_density = TRUE)
#That should be fine

nd_diff_notpu_lf %>%
    with(., shapiro.test(vc_diff))

gghistogram(nd_diff_tpu_lf, x = "vc_diff", bins = 10, add_density = TRUE)
#okay

nd_diff_tpu_lf %>%
    with(., shapiro.test(vc_diff))
# All no-overshoot Vcmax looks fine.


### Jmax
gghistogram(nd_diff_notpu_lf, x = "j_diff", bins = 10, add_density = TRUE)

nd_diff_notpu_lf %>%
    with(., shapiro.test(j_diff))

gghistogram(nd_diff_tpu_lf, x = "j_diff", bins = 10, add_density = TRUE)

nd_diff_tpu_lf %>%
    with(., shapiro.test(j_diff))
# no-overshoot Jmax is only okay.

# To use wilcoxon paired, we assume the differences between paired samples are distributed symmetrically about the median.

# It's just Jmax where we don't totally meet the Wilcoxon assumptions. 
# We will run both Sign test and Wilcoxon test on Jmax to see if their results agree.



#Wilcoxon tests for the data grouped on a tree level. ---------------------------

all_avg_tr_res <- all_avg_lf_res %>% 
    group_by(fit_type, curv_meth, treeid) %>% 
    summarize(vcmax = mean(vcmax),
              jmax = mean(jmax))

### Vcmax Wilcoxon by curve method 
w_vc_cm <- all_avg_tr_res %>%
    group_by(fit_type) %>%
    wilcox_test(data =., vcmax ~ curv_meth, ref.group = 'SS', paired = TRUE, detailed = TRUE) %>%
    add_significance()
wes_vc_cm <- all_avg_tr_res %>%
    group_by(fit_type) %>%
    wilcox_effsize(data = ., vcmax ~ curv_meth, ref.group = 'SS', paired = TRUE)
w_vc_cm_full <- left_join(w_vc_cm, wes_vc_cm)
w_vc_cm_full

### Jmax Wilcoxon by curve method
w_j_cm <- all_avg_tr_res %>%
    group_by(fit_type) %>%
    wilcox_test(data =., jmax ~ curv_meth, ref.group = 'SS', paired = TRUE, detailed = TRUE) %>%
    add_significance()
wes_j_cm <- all_avg_tr_res %>%
    group_by(fit_type) %>%
    wilcox_effsize(data = ., jmax ~ curv_meth, ref.group = 'SS', paired = TRUE)
w_j_cm_full <- left_join(w_j_cm, wes_j_cm)
w_j_cm_full

### Vcmax Wilcoxon by TPU v. no TPU
w_vc_ft <- all_avg_tr_res %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>% 
    add_significance()
wes_vc_ft <- all_avg_tr_res %>%
    group_by(curv_meth) %>%
    wilcox_effsize(data = ., vcmax ~ fit_type, paired = TRUE)
w_vc_ft_full <- left_join(w_vc_ft, wes_vc_ft)
w_vc_ft_full

### Jmax Wilcoxon by TPU v. no TPU
w_j_ft <- all_avg_tr_res %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., jmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    add_significance()
wes_j_ft <- all_avg_tr_res %>%
    group_by(curv_meth) %>%
    wilcox_effsize(data = ., jmax ~ fit_type, paired = TRUE)
w_j_ft_full <- left_join(w_j_ft, wes_j_ft)
w_j_ft_full


#Wilcoxon tests for the data grouped on a tree level, without MAEL Leaf 6 (testing for influence). ---------------------------

all_nol6_tr_res <- all_avg_lf_res %>% 
    filter(leaf_unique != 'K6709L6') %>% 
    group_by(fit_type, curv_meth, treeid) %>% 
    summarize(vcmax = mean(vcmax),
              jmax = mean(jmax))

### Vcmax Wilcoxon by curve method 
w_nol6_vc_cm <- all_nol6_tr_res %>%
    group_by(fit_type) %>%
    wilcox_test(data =., vcmax ~ curv_meth, ref.group = 'SS', paired = TRUE, detailed = TRUE) %>%
    add_significance()
wes_nol6_vc_cm <- all_nol6_tr_res %>%
    group_by(fit_type) %>%
    wilcox_effsize(data = ., vcmax ~ curv_meth, ref.group = 'SS', paired = TRUE)
w_nol6_vc_cm_full <- left_join(w_nol6_vc_cm, wes_nol6_vc_cm)
w_nol6_vc_cm_full

### Jmax Wilcoxon by curve method
w_nol6_j_cm <- all_nol6_tr_res %>%
    group_by(fit_type) %>%
    wilcox_test(data =., jmax ~ curv_meth, ref.group = 'SS', paired = TRUE, detailed = TRUE) %>%
    add_significance()
wes_nol6_j_cm <- all_nol6_tr_res %>%
    group_by(fit_type) %>%
    wilcox_effsize(data = ., jmax ~ curv_meth, ref.group = 'SS', paired = TRUE)
w_nol6_j_cm_full <- left_join(w_nol6_j_cm, wes_nol6_j_cm)
w_nol6_j_cm_full

### Vcmax Wilcoxon by TPU v. no TPU
w_nol6_vc_ft <- all_nol6_tr_res %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>% 
    add_significance()
wes_nol6_vc_ft <- all_nol6_tr_res %>%
    group_by(curv_meth) %>%
    wilcox_effsize(data = ., vcmax ~ fit_type, paired = TRUE)
w_nol6_vc_ft_full <- left_join(w_nol6_vc_ft, wes_nol6_vc_ft)
w_nol6_vc_ft_full

### Jmax Wilcoxon by TPU v. no TPU
w_nol6_j_ft <- all_nol6_tr_res %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., jmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    add_significance()
wes_nol6_j_ft <- all_nol6_tr_res %>%
    group_by(curv_meth) %>%
    wilcox_effsize(data = ., jmax ~ fit_type, paired = TRUE)
w_nol6_j_ft_full <- left_join(w_nol6_j_ft, wes_nol6_j_ft)
w_nol6_j_ft_full



# Wilcoxon Tests, no Overshoot subset ------------------------

### Note that all the SS curves for which TPU was fit also had overshoot, so the SS TPU and noTPU datasets are actually the same.
# Therefore, no steady-state-specific TPU vs noTPU comparisons are conducted here.

nd_tr_res <- nd_complete %>% 
    group_by(fit_type, curv_meth, treeid) %>% 
    summarize(vcmax = mean(vcmax),
              jmax = mean(jmax))

### Vcmax Wilcoxon by curve method
w_nd_vc_cm <- nd_tr_res %>%
    group_by(fit_type) %>%
    wilcox_test(data =., vcmax ~ curv_meth, ref.group = 'SS', paired = TRUE, detailed = TRUE) %>%
    add_significance()
wes_nd_vc_cm <- nd_tr_res %>%
    group_by(fit_type) %>%
    wilcox_effsize(data = ., vcmax ~ curv_meth, ref.group = 'SS', paired = TRUE)
w_nd_vc_cm_full <- left_join(w_nd_vc_cm, wes_nd_vc_cm)
w_nd_vc_cm_full

### Jmax Wilcoxon by curve method
w_nd_j_cm <- nd_tr_res %>%
    group_by(fit_type) %>%
    wilcox_test(data =., jmax ~ curv_meth, ref.group = 'SS', paired = TRUE, detailed = TRUE) %>%
    add_significance()
wes_nd_j_cm <-nd_tr_res %>%
    group_by(fit_type) %>%
    wilcox_effsize(data = ., jmax ~ curv_meth, ref.group = 'SS', paired = TRUE)
w_nd_j_cm_full <- left_join(w_nd_j_cm, wes_nd_j_cm)
w_nd_j_cm_full

### Vcmax Wilcoxon by TPU v. no TPU
# The SS TPU and noTPU datasets are actually the same.
w_nd_vc_ft <- nd_tr_res %>%
    filter(curv_meth == 'DAT') %>%
    group_by(curv_meth) %>% 
    wilcox_test(data =., vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    add_significance()
wes_nd_vc_ft <- nd_tr_res %>%
    filter(curv_meth == 'DAT') %>%
    group_by(curv_meth) %>%
    wilcox_effsize(data = ., vcmax ~ fit_type, paired = TRUE)
w_nd_vc_ft_full <- left_join(w_nd_vc_ft, wes_nd_vc_ft)
w_nd_vc_ft_full

### Jmax Wilcoxon by TPU v. no TPU
w_nd_j_ft <- nd_tr_res %>%
    filter(curv_meth == 'DAT') %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., jmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    add_significance()
wes_nd_j_ft <- nd_tr_res %>%
    filter(curv_meth == 'DAT') %>%
    group_by(curv_meth) %>%
    wilcox_effsize(data = ., jmax ~ fit_type, paired = TRUE)
w_nd_j_ft_full <- left_join(w_nd_j_ft, wes_nd_j_ft)
w_nd_j_ft_full


# TPU comparisons, Wilcoxon ---------------------
# Note this is a very small sample size and should be interpreted cautiously!
tpu_tr_res <- tpu_just6_all %>% 
    mutate(treeid = substring(leaf_unique, 1, 5)) %>% 
    group_by(curv_meth, treeid) %>% 
    summarize(vcmax = mean(vcmax),
              jmax = mean(jmax),
              tpu = mean(tpu)) %>%
    ungroup()

#tpu_tr_res$curv_meth <- as.factor(tpu_tr_res$curv_meth)

w_tpu_cm <- tpu_tr_res %>%
    wilcox_test(data = ., tpu ~ curv_meth, ref.group = 'SS', paired = TRUE, detailed = TRUE) %>% 
    add_significance()
wes_tpu_cm <- tpu_tr_res %>% 
    wilcox_effsize(data = ., tpu ~ curv_meth, ref.group = 'SS', paired = TRUE)
w_tpu_full <- left_join(w_tpu_cm, wes_tpu_cm) %>%
    mutate(fit_type = 'tpu') %>% 
    select(fit_type, everything())
w_tpu_full$fit_type <- factor(w_tpu_full$fit_type)
w_tpu_full


# List of models
wil_cm <- list(
    w_vc_cm_full,
    w_j_cm_full,
    w_tpu_full) %>%
    do.call(rbind, .) %>% 
    mutate(dataset = 'all_data') %>% 
    arrange(desc(fit_type))

wil_nd_cm <- list(
    w_nd_vc_cm_full,
    w_nd_j_cm_full) %>%
    do.call(rbind, .) %>%
    mutate(dataset = 'nd') %>% 
    arrange(desc(fit_type))

wil_nolf6_cm <- list(
    w_nol6_vc_cm_full,
    w_nol6_j_cm_full) %>%
    do.call(rbind, .) %>% 
    mutate(dataset = 'nolf6') %>% 
    arrange(desc(fit_type))

wilcox_cm_tab <- rbind(wil_cm, wil_nd_cm, wil_nolf6_cm) %>% 
    select(dataset, fit_type, .y., n1, n2, group1, group2, estimate, statistic, p, p.signif, effsize, conf.low, conf.high, magnitude)
wilcox_cm_tab

write.csv(x = wilcox_cm_tab, 
          file = here("5_Results/wilcox_table.csv"),
          row.names = FALSE)

wilcox_ft_tab <- list(
    w_vc_ft_full,
    w_j_ft_full,
    w_nol6_vc_ft_full,
    w_nol6_j_ft_full,
    w_nd_vc_ft_full,
    w_nd_j_ft_full
) %>% do.call(rbind, .)
wilcox_ft_tab


# Exploring rogme package (CDS added 6/27/24) -------------------------
#Need to use this to install the package:
#install.packages("remotes")
#remotes::install_github("GRousselet/rogme")
# library(rogme)
# 
# #If the goal is to detect differences anywhere in the distributions, a systematic approach consists of quantifying differences at multiple quantiles. First, for each participant (tree) and each condition (curv_meth), the sample deciles are computed over trials (leaves). Second, for each participant, condition 2 deciles are subtracted from condition 1 deciles - we’re dealing with a within-subject (repeated-measure) design. Third, for each decile, the distribution of differences is subjected to a one-sample test. Fourth, a correction for multiple comparisons is applied across the 9 one-sample tests. We call this procedure a hierarchical shift function. 
# 
# all_avg_lf_res$treeid <- factor(all_avg_lf_res$treeid)
# all_avg_lf_res$curv_meth <- factor(all_avg_lf_res$curv_meth)
# all_avg_lf_res$fit_type <- factor(all_avg_lf_res$fit_type)
# 
# tpu_res <- all_avg_lf_res %>% filter(fit_type == 'tpu')
# notpu_res <- all_avg_lf_res %>% filter(fit_type == 'no_tpu')
# 
# dat_res <- all_avg_lf_res %>% filter(curv_meth == 'DAT')
# ss_res <- all_avg_lf_res %>% filter(curv_meth == 'SS')
# 
# np <- length(unique(tpu_res$treeid)) #Number of 'participants' (trees)
# 
# #TPU vcmax hierarchical shift function (COMPLETE)
# set.seed(304)
# sf_v1 <- shiftdhd_pbci(tpu_res, formula = vcmax ~ curv_meth + treeid, nboot = 500)
# p_v1 <- plot_sf(sf_v1, plot_theme = 1)[[1]] + 
#     theme(axis.text = element_text(size = 16, colour="black"))
# p_v1
# 
# hsf_v1 <- hsf(tpu_res, vcmax ~ curv_meth + treeid) 
# #Plot hierarchical shift function
# p_hsf_v1 <- plot_hsf(hsf_v1, viridis_option = "D", ind_line_size = 0.8, gp_line_colour = "maroon3", gp_point_colour = "maroon3", gp_line_size = 1.2)  + ylim(-20,20) + ggtitle("Vcmax, TPU, DAT - SS")
# p_hsf_v1
# 
# hsf_v1$pvalues
# hsf_v1$adjusted_pvalues
# 
# #stochastic dominance
# nq_v1 <- length(hsf_v1$quantiles)
# pdmt0_v1 <- apply(hsf_v1$individual_sf > 0, 2, sum)
# print(paste0('In ',sum(pdmt0_v1 == nq_v1),' trees (',round(100 * sum(pdmt0_v1 == nq_v1) / np, digits = 1),'%), all quantile differences are more than zero at all points'))
# 
# pdlt0_v1 <- apply(hsf_v1$individual_sf < 0, 2, sum)
# print(paste0('In ',sum(pdlt0_v1 == nq_v1),' trees (',round(100 * sum(pdlt0_v1 == nq_v1) / np, digits = 1),'%), all quantile differences are less than zero at all points'))
# 
# #percentile bootstrap hierarchical shift function
# set.seed(304)
# hsf_pb_v1 <- hsf_pb(tpu_res, vcmax ~ curv_meth + treeid)
# 
# plot_hsf_pb(hsf_pb_v1, interv = "hdi")
# plot_hsf_pb_dist(hsf_pb_v1, point_interv = "median_ci", interval_width = .95, 
#                  int_colour = "blue", fill_colour = "grey")
# 
# 
# 
# 
# #TPU jmax hierarchical shift function (COMPLETE)
# set.seed(304)
# 
# hsf_j1 <- hsf(tpu_res, jmax ~ curv_meth + treeid)
# #Plot hierarchical shift function
# p_hsf_j1 <- plot_hsf(hsf_j1, viridis_option = "D", ind_line_size = 0.8, gp_line_colour = "maroon3", gp_point_colour = "maroon3", gp_line_size = 1.2) + ylim(-50,20) + ggtitle("Jmax, TPU, DAT - SS")
# p_hsf_j1
# 
# hsf_j1$pvalues
# hsf_j1$adjusted_pvalues
# 
# #stochastic dominance
# nq_j1 <- length(hsf_j1$quantiles)
# pdmt0_j1 <- apply(hsf_j1$individual_sf > 0, 2, sum)
# print(paste0('In ',sum(pdmt0_j1 == nq_j1),' trees (',round(100 * sum(pdmt0_j1 == nq_j1) / np, digits = 1),'%), all quantile differences are more than zero at all points'))
# 
# pdlt0_j1 <- apply(hsf_j1$individual_sf < 0, 2, sum)
# print(paste0('In ',sum(pdlt0_j1 == nq_j1),' trees (',round(100 * sum(pdlt0_j1 == nq_j1) / np, digits = 1),'%), all quantile differences are less than zero at all points'))
# 
# #percentile bootstrap hierarchical shift function
# set.seed(304)
# hsf_pb_j1 <- hsf_pb(tpu_res, jmax ~ curv_meth + treeid)
# 
# plot_hsf_pb(hsf_pb_j1, interv = "hdi")
# plot_hsf_pb_dist(hsf_pb_j1, point_interv = "median_ci", interval_width = .95, 
#                  int_colour = "blue", fill_colour = "grey")
# 
# 
# 
# 
# #No TPU vcmax hierarchical shift function
# set.seed(304)
# hsf_v2 <- hsf(notpu_res, vcmax ~ curv_meth + treeid)
# #Plot hierarchical shift function
# p_hsf_v2 <- plot_hsf(hsf_v2, viridis_option = "D", ind_line_size = 0.8, gp_line_colour = "maroon3", gp_point_colour = "maroon3", gp_line_size = 1.2)+ ylim(-20,20) + ggtitle("Vcmax, no TPU, DAT - SS")
# p_hsf_v2
# 
# hsf_v2$pvalues
# hsf_v2$adjusted_pvalues
# 
# #stochastic dominance
# nq_v2 <- length(hsf_v2$quantiles)
# pdmt0_v2 <- apply(hsf_v2$individual_sf > 0, 2, sum)
# print(paste0('In ',sum(pdmt0_v2 == nq_v2),' trees (',round(100 * sum(pdmt0_v2 == nq_v2) / np, digits = 1),'%), all quantile differences are more than zero at all points'))
# 
# pdlt0_v2 <- apply(hsf_v2$individual_sf < 0, 2, sum)
# print(paste0('In ',sum(pdlt0_v2 == nq_v2),' trees (',round(100 * sum(pdlt0_v2 == nq_v2) / np, digits = 1),'%), all quantile differences are less than zero at all points'))
# 
# #No TPU jmax hierarchical shift function
# set.seed(304)
# hsf_j2 <- hsf(notpu_res, jmax ~ curv_meth + treeid)
# #Plot hierarchical shift function
# p_hsf_j2 <- plot_hsf(hsf_j2, viridis_option = "D", ind_line_size = 0.8, gp_line_colour = "maroon3", gp_point_colour = "maroon3", gp_line_size = 1.2) + ylim(-50,20) + ggtitle("Jmax, no TPU, DAT - SS")
# p_hsf_j2
# 
# hsf_j2$pvalues
# hsf_j2$adjusted_pvalues
# 
# #stochastic dominance
# nq_j2 <- length(hsf_j2$quantiles)
# pdmt0_j2 <- apply(hsf_j2$individual_sf > 0, 2, sum)
# print(paste0('In ',sum(pdmt0_j2 == nq_j2),' trees (',round(100 * sum(pdmt0_j2 == nq_j2) / np, digits = 1),'%), all quantile differences are more than zero at all points'))
# 
# pdlt0_j2 <- apply(hsf_j2$individual_sf < 0, 2, sum)
# print(paste0('In ',sum(pdlt0_j2 == nq_j2),' trees (',round(100 * sum(pdlt0_j2 == nq_j2) / np, digits = 1),'%), all quantile differences are less than zero at all points'))
# 
# #DAT vcmax by fit type hierarchical shift function
# set.seed(304)
# hsf_v3 <- hsf(dat_res, vcmax ~ fit_type + treeid)
# #Plot hierarchical shift function
# p_hsf_v3 <- plot_hsf(hsf_v3, viridis_option = "D", ind_line_size = 0.8, gp_line_colour = "maroon3", gp_point_colour = "maroon3", gp_line_size = 1.2) + ylim(-20,20) + ggtitle("Vcmax, DAT, no TPU - TPU")
# p_hsf_v3
# 
# hsf_v3$pvalues
# hsf_v3$adjusted_pvalues
# 
# #stochastic dominance
# nq_v3 <- length(hsf_v3$quantiles)
# pdmt0_v3 <- apply(hsf_v3$individual_sf > 0, 2, sum)
# print(paste0('In ',sum(pdmt0_v3 == nq_v3),' trees (',round(100 * sum(pdmt0_v3 == nq_v3) / np, digits = 1),'%), all quantile differences are more than zero at all points'))
# 
# pdlt0_v3 <- apply(hsf_v3$individual_sf < 0, 2, sum)
# print(paste0('In ',sum(pdlt0_v3 == nq_v3),' trees (',round(100 * sum(pdlt0_v3 == nq_v3) / np, digits = 1),'%), all quantile differences are less than zero at all points'))
# 
# #DAT Jmax by fit type hierarchical shift function
# set.seed(304)
# hsf_j3 <- hsf(dat_res, jmax ~ fit_type + treeid)
# #Plot hierarchical shift function
# p_hsf_j3 <- plot_hsf(hsf_j3, viridis_option = "D", ind_line_size = 0.8, gp_line_colour = "maroon3", gp_point_colour = "maroon3", gp_line_size = 1.2) + ylim(-50,20) + ggtitle("Jmax, DAT, no TPU - TPU")
# p_hsf_j3
# 
# hsf_j3$pvalues
# hsf_j3$adjusted_pvalues
# 
# #stochastic dominance
# nq_j3 <- length(hsf_j3$quantiles)
# pdmt0_j3 <- apply(hsf_j3$individual_sf > 0, 2, sum)
# print(paste0('In ',sum(pdmt0_j3 == nq_j3),' trees (',round(100 * sum(pdmt0_j3 == nq_j3) / np, digits = 1),'%), all quantile differences are more than zero at all points'))
# 
# pdlt0_j3 <- apply(hsf_j3$individual_sf < 0, 2, sum)
# print(paste0('In ',sum(pdlt0_j3 == nq_j3),' trees (',round(100 * sum(pdlt0_j3 == nq_j3) / np, digits = 1),'%), all quantile differences are less than zero at all points'))
# 
# #SS Vcmax by fit type hierarchical shift function
# set.seed(304)
# hsf_v4 <- hsf(ss_res, vcmax ~ fit_type + treeid)
# #Plot hierarchical shift function
# p_hsf_v4 <- plot_hsf(hsf_v4, viridis_option = "D", ind_line_size = 0.8, gp_line_colour = "maroon3", gp_point_colour = "maroon3", gp_line_size = 1.2) + ylim(-20,20)+ ggtitle("Vcmax, SS, no TPU - TPU")
# p_hsf_v4
# 
# hsf_v4$pvalues
# hsf_v4$adjusted_pvalues
# 
# #stochastic dominance
# nq_v4 <- length(hsf_v4$quantiles)
# pdmt0_v4 <- apply(hsf_v4$individual_sf > 0, 2, sum)
# print(paste0('In ',sum(pdmt0_v4 == nq_v4),' trees (',round(100 * sum(pdmt0_v4 == nq_v4) / np, digits = 1),'%), all quantile differences are more than zero at all points'))
# 
# pdlt0_v4 <- apply(hsf_v4$individual_sf < 0, 2, sum)
# print(paste0('In ',sum(pdlt0_v4 == nq_v4),' trees (',round(100 * sum(pdlt0_v4 == nq_v4) / np, digits = 1),'%), all quantile differences are less than zero at all points'))
# 
# #DAT Jmax by fit type hierarchical shift function
# set.seed(304)
# hsf_j4 <- hsf(ss_res, jmax ~ fit_type + treeid)
# #Plot hierarchical shift function
# p_hsf_j4 <- plot_hsf(hsf_j4, viridis_option = "D", ind_line_size = 0.8, gp_line_colour = "maroon3", gp_point_colour = "maroon3", gp_line_size = 1.2) + ylim(-50,20) + ggtitle("Jmax, SS, no TPU - TPU")
# p_hsf_j4
# 
# hsf_j4$pvalues
# hsf_j4$adjusted_pvalues
# 
# #stochastic dominance
# nq_j4 <- length(hsf_j4$quantiles)
# pdmt0_j4 <- apply(hsf_j4$individual_sf > 0, 2, sum)
# print(paste0('In ',sum(pdmt0_j4 == nq_j4),' trees (',round(100 * sum(pdmt0_j4 == nq_j4) / np, digits = 1),'%), all quantile differences are more than zero at all points'))
# 
# pdlt0_j4 <- apply(hsf_j4$individual_sf < 0, 2, sum)
# print(paste0('In ',sum(pdlt0_j4 == nq_j4),' trees (',round(100 * sum(pdlt0_j4 == nq_j4) / np, digits = 1),'%), all quantile differences are less than zero at all points'))
# 
# library(patchwork)
# 
# curv_meth_grid <- p_hsf_v1 + p_hsf_v2 + p_hsf_j1 + p_hsf_j2
# curv_meth_grid
# 
# fit_type_grid <- p_hsf_v3 + p_hsf_v4 + p_hsf_j3 + p_hsf_j4
# fit_type_grid

# Run Wilcoxon and Sign tests ------------------------------------------


# All data, Wilcoxon signed rank test on paired samples (and a few sign tests)

# ### Vcmax Wilcoxon by curve method 
# all_avg_lf_res %>%
#     group_by(fit_type) %>%
#     wilcox_test(data =., vcmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
#     add_significance()
# 
# ## Effect size
# all_avg_lf_res %>%
#     group_by(fit_type) %>% 
#     wilcox_effsize(data = ., vcmax ~ curv_meth, paired = TRUE)
# 
# 
# ### Vcmax Wilcoxon by TPU v. no TPU
# all_avg_lf_res %>%
#     group_by(curv_meth) %>%
#     wilcox_test(data =., vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>% 
#     add_significance()
# 
# all_avg_lf_res %>%
#     group_by(curv_meth) %>% 
#     wilcox_effsize(data = ., vcmax ~ fit_type, paired = TRUE)
# 
# 
# ### Jmax Wilcoxon by curve method
# all_avg_lf_res %>%
#     group_by(fit_type) %>%
#     wilcox_test(data =., jmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
#     add_significance()
# 
# # Note we're running the sign test here in addition!
# all_avg_lf_res %>%
#     group_by(fit_type) %>%
#     sign_test(data =., jmax ~ curv_meth, detailed = TRUE) %>%
#     add_significance()
# 
# all_avg_lf_res %>%
#     group_by(fit_type) %>% 
#     wilcox_effsize(data = ., jmax ~ curv_meth, paired = TRUE)
# 
# 
# 
# ### Jmax Wilcoxon by TPU v. no TPU
# all_avg_lf_res %>%
#     group_by(curv_meth) %>%
#     wilcox_test(data =., jmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
#     add_significance()
# 
# # Sign test here. Note we have a different result.
# all_avg_lf_res %>%
#     group_by(curv_meth) %>%
#     sign_test(data =., jmax ~ fit_type, detailed = TRUE) %>%
#     add_significance()
# 
# all_avg_lf_res %>%
#     group_by(curv_meth) %>% 
#     wilcox_effsize(data = ., jmax ~ fit_type, paired = TRUE)
# 
# 
# # Wilcoxon Tests, no Overshoot subset
# 
# ### Vcmax Wilcoxon by curve method
# nd_complete %>%
#     group_by(fit_type) %>%
#     wilcox_test(data =., vcmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
#     add_significance()
# 
# nd_complete %>%
#     group_by(fit_type) %>% 
#     wilcox_effsize(data = ., vcmax ~ curv_meth, paired = TRUE)
# 
# 
# ### Vcmax Wilcoxon by TPU v. no TPU
# 
# ### Note that all the SS curves for which TPU was fit had overshoot, so the datasets are same
# nd_complete %>%
#     group_by(curv_meth) %>%
#     wilcox_test(data =., vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
#     add_significance()
# 
# 
# ### Jmax Wilcoxon by curve method
# nd_complete %>%
#     group_by(fit_type) %>%
#     wilcox_test(data =., jmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
#     add_significance()
# 
# nd_complete %>%
#     group_by(fit_type) %>%
#     wilcox_effsize(data = ., jmax ~ curv_meth, paired = TRUE)
# 
# 
# ### Jmax Wilcoxon by TPU v. no TPU
# 
# ### Note that all the SS curves for which TPU was fit had overshoot, so the datasets are same
# nd_complete %>%
#     group_by(curv_meth) %>%
#     wilcox_test(data =., jmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
#     add_significance()
# 
# #####
# 
# 
# # TPU Wilcoxon
# 
# tpu_just6_all %>% 
#     wilcox_test(data = ., tpu ~ curv_meth, paired = TRUE, detailed = TRUE) %>% 
#     add_significance()
# 
# tpu_just6_all %>% 
#     wilcox_effsize(data = ., tpu ~ curv_meth, paired = TRUE)

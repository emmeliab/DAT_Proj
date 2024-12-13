# Computing statistical analysis on the A/Ci curve fits from 2_Tapajos_Fit_ACi.R
# Code is associated with the article at DOI: 10.1093/treephys/tpae153
# Licence information:
# Questions can be directed to Loren Albert (corresponding author) at loren.albert@oregonstate.edu, Emmelia Braun (first-author) at emmelia.braun@oregonstate.edu, or Charles Southwick (first author) at charles.southwick@oregonstate.edu


library(tidyverse) 
library(ggpubr) 
library(car)
library(rstatix)
library(reshape2) 
library(here) 
library(nlme)
library(performance)

set.seed(304)

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
            jmax = mean(jmax),
            tpu = mean(tpu, na.rm = TRUE)) %>% 
    mutate(fit_type = "tpu",
           treeid = substring(leaf_unique, 1, 5))


### Create summary table for WITHOUT TPU results
notpu_results_grp <- pho_stat %>%
    group_by(curv_meth, leaf_unique) %>% 
    summarise(vcmax = mean(vcmax),
            jmax = mean(jmax),
            tpu = mean(tpu)) %>%
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
write.csv(diff_tpu_lf, here("5_Results/lf_diffs_summ_TPU.csv"))



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
            jmax = mean(jmax),
            tpu = mean(tpu, na.rm = TRUE)) %>% 
    mutate(fit_type = "tpu",
           treeid = substring(leaf_unique, 1, 5))


### Create summary table for WITHOUT TPU results
notpu_nd_results_grp <- pho_nd_stat %>%
    group_by(curv_meth, leaf_unique) %>% 
    reframe(vcmax = mean(vcmax),
            jmax = mean(jmax),
            tpu = mean(tpu)) %>%
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
            tpu = mean(tpu),
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
            tpu = mean(tpu, na.rm = TRUE),
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
    get_summary_stats(vcmax, jmax, tpu, show = c("mean", "median", "sd", "min", "max")) %>% 
    as.data.frame()


### No-overshoot data
nd_complete %>% 
    group_by(fit_type, curv_meth) %>% 
    get_summary_stats(vcmax, jmax, tpu, show = c("mean", "median", "sd", "min", "max")) %>%
    as.data.frame()


### Just the TPU data
tpu_just6_all %>% 
    group_by(curv_meth)%>%
    get_summary_stats(vcmax, jmax, show = c("mean", "median", "sd", "min", "max"))



# Look at data structure ---------------------------------
# library(lmerTest)
# library(lme4)
# library(performance)

## With TPU
plot(factor(diff_tpu_lf$treeid), diff_tpu_lf$vc_diff)
plot(factor(pho_both_tpu$curv_meth), pho_both_tpu$V_cmax)
hist(diff_tpu_lf$vc_diff)

plot(factor(diff_tpu_lf$treeid), diff_tpu_lf$j_diff)
plot(factor(pho_both_tpu$curv_meth), pho_both_tpu$J_max)
hist(diff_tpu_lf$j_diff)


## Without TPU
plot(factor(diff_notpu_lf$treeid), diff_notpu_lf$vc_diff)
plot(factor(pho_both$curv_meth), pho_both$V_cmax)
hist(diff_notpu_lf$vc_diff)

plot(factor(diff_notpu_lf$treeid), diff_notpu_lf$j_diff)
plot(factor(pho_both$curv_meth), pho_both$J_max)
hist(diff_notpu_lf$j_diff)






# Set up datasets for mixed models ------------------------

## Convert tree id to factor. n = 27 pairs
diff_notpu_lf$treeid <- as.factor(diff_notpu_lf$treeid)
diff_tpu_lf$treeid <- as.factor(diff_tpu_lf$treeid)

## Set up a dataframe to exclude K6709L6 (for testing; the one with the big diff). n = 26 pairs
diff_notpu_nol6_lf <- diff_notpu_lf %>% filter(leaf_unique != 'K6709L6')
diff_tpu_nol6_lf <- diff_tpu_lf %>% filter(leaf_unique != 'K6709L6')

## Convert tree id in no Overshoot subset to factor. n = 19 pairs
nd_diff_notpu_lf$treeid <- as.factor(nd_diff_notpu_lf$treeid)
nd_diff_tpu_lf$treeid <- as.factor(nd_diff_tpu_lf$treeid)

nd_diff_notpu_lf <- nd_diff_notpu_lf %>% filter(!is.na(vc_diff))
nd_diff_tpu_lf <- nd_diff_tpu_lf %>% filter(!is.na(vc_diff))

###

# Define linear mixed effects models ----------------------

### Mixed model with random effect for intercept
### Our leaves are nested within trees. Trees will be the random effect.
### These models have the difference (SS - DAT) Vcmax or Jmax as the response variable.

mod_list <- list()


# Models with full dataset

## Vcmax without TPU
mod_list[["mod_notpu_v"]] <- mod_notpu_v <- lme(vc_diff ~ 1,
                   random = ~ 1 | treeid,
                   weights = varIdent(form = ~ 1 | treeid),
                   data = diff_notpu_lf)


# mod_notpu_v <- lmer(vc_diff ~ 1 + (1|treeid),
#                     data = diff_notpu_lf)
# 
# summary(mod_notpu_v)
# ### Slightly sig diff, estimated diff 2.08 +- 0.92


## Jmax without TPU
mod_list[["mod_notpu_j"]] <- mod_notpu_j <- lme(j_diff ~ 1,
                   random = ~ 1 | treeid,
                   weights = varIdent(form = ~ 1 | treeid),
                   data = diff_notpu_lf)

# mod_notpu_j <- lmer(j_diff ~ 1 + (1|treeid),
#                     data = diff_notpu_lf)
# 
# mod_notpu_j2 <- lmer(jmax ~ curv_meth + (1|treeid) + (1|treeid:leaf_unique), data = pho_stat)
# 
# summary(mod_notpu_j)
# ### Sig diff, estimated diff 10.88 +- 3.0 



## Vcmax with TPU
mod_list[["mod_tpu_v"]] <- mod_tpu_v <- lme(vc_diff ~ 1,
                 random = ~ 1 | treeid,
                 weights = varIdent(form = ~ 1 | treeid),
                 data = diff_tpu_lf)


# mod_tpu_v <- lmer(vc_diff ~ 1 + (1|treeid),
#                     data = diff_tpu_lf)
# ### This leads to a singularity
# 
# summary(mod_tpu_v)
# ### Not sig, estimate diff 1.00 +- 0.67; gives singularity



## Jmax with TPU
mod_list[["mod_tpu_j"]] <- mod_tpu_j <- lme(j_diff ~ 1,
                 random = ~ 1 | treeid,
                 weights = varIdent(form = ~ 1 | treeid),
                 data = diff_tpu_lf)

# mod_tpu_j <- lmer(j_diff ~ 1 + (1|treeid),
#                   data = diff_tpu_lf)
# 
# mod_tpu_j2 <- lmer(jmax ~ curv_meth + (1|treeid) + (1|treeid:leaf_unique), data = pho_stat_tpu)
# 
# 
# summary(mod_tpu_j)
# ### Sig, estimated differences 7.9 +- 2.46




## TPU SS - TPU DAT comparison
# tpu_comparison <- pho_stat_tpu %>%
#     group_by(curv_meth, leaf_unique) %>%
#     summarise(tpu = mean(tpu)) %>%
#     mutate(fit_type = "tpu",
#            treeid = substring(leaf_unique, 1, 5))
# 
# ### Pull out the DAT results
# dat_res_tpu_comparison <- filter(tpu_comparison, curv_meth == "DAT") %>%
#     rename(dat_tpu = tpu)
# 
# ### Pull out the SS results
# ss_res_tpu_comparison <- filter(tpu_comparison, curv_meth == "SS") %>%
#     rename(ss_tpu = tpu)
# 
# ### Join datasets, calculate diff and SE
# diff_tpu_comparison <- full_join(dat_res_tpu_comparison, ss_res_tpu_comparison,
#                          by = c("leaf_unique", "treeid", "fit_type")) %>%
#     mutate(tpu_diff = ss_tpu - dat_tpu) %>%
#     filter(!is.na(tpu_diff))

### Pivot the dataframe wider for comparisons
diff_tpu_comparison <- pivot_wider(tpu_just6_all, 
                                   id_cols = leaf_unique, 
                                   names_from = curv_meth, values_from = tpu) %>% 
    mutate(tpu_diff = SS - DAT, treeid = substring(.$leaf_unique, 1, 5))

plot(factor(diff_tpu_comparison$treeid), diff_tpu_comparison$tpu_diff)


### Should we go with the unequal variance model for TPU??
# 
# tpu_only_mod <- lmer(tpu_diff ~ 1 + (1|treeid),
#                          data = diff_tpu_comparison)



mod_list[["tpu_only_mod"]] <- tpu_only_mod <- lme(tpu_diff ~ 1,
                 random = ~ 1 | treeid,
                 # weights = varIdent(form = ~ 1 | treeid), ## small sample size, so not doing unequal variances
                 data = diff_tpu_comparison)




# Models without K6709L6 -- Testing the effects of this potentially influential point.

## Vcmax no TPU
mod_list[["mod_notpu_nol6_v"]] <- mod_notpu_nol6_v <- lme(vc_diff ~ 1,
                 random = ~ 1 | treeid,
                 weights = varIdent(form = ~ 1 | treeid),
                 control = lmeControl(opt = "optim"),
                 data = diff_notpu_nol6_lf)

# mod_notpu_nol6_v <- lmer(vc_diff ~ 1 + (1|treeid),
#                     data = diff_notpu_nol6_lf)


## Jmax no TPU
mod_list[["mod_notpu_nol6_j"]] <- mod_notpu_nol6_j <- lme(j_diff ~ 1,
                                         random = ~ 1 | treeid,
                                         weights = varIdent(form = ~ 1 | treeid),
                                         control = lmeControl(opt = "optim"),
                                         data = diff_notpu_nol6_lf)
# mod_notpu_nol6_j <- lmer(j_diff ~ 1 + (1|treeid),
#                     data = diff_notpu_nol6_lf)


## Vcmax with TPU
mod_list[["mod_tpu_nol6_v"]] <- mod_tpu_nol6_v <- lme(vc_diff ~ 1,
                        random = ~ 1 | treeid,
                        weights = varIdent(form = ~ 1 | treeid),
                        #control = lmeControl(opt = "optim"),
                        data = diff_tpu_nol6_lf)
# mod_tpu_nol6_v <- lmer(vc_diff ~ 1 + (1|treeid),
#                   data = diff_tpu_nol6_lf) #This leads to a singularity


## Jmax with TPU
mod_list[["mod_tpu_nol6_j"]] <- mod_tpu_nol6_j <- lme(j_diff ~ 1,
                        random = ~ 1 | treeid,
                        weights = varIdent(form = ~ 1 | treeid),
                        #control = lmeControl(opt = "optim"),
                        data = diff_tpu_nol6_lf)
# mod_tpu_nol6_j <- lmer(j_diff ~ 1 + (1|treeid),
#                   data = diff_tpu_nol6_lf)



# Models with no-overshoot data. Just a quick check
mod_list[["mod_nd_notpu_v"]] <- mod_nd_notpu_v <- lme(vc_diff ~ 1,
                        random = ~ 1 | treeid,
                        weights = varIdent(form = ~ 1 | treeid),
                        #control = lmeControl(opt = "optim"),
                        data = nd_diff_notpu_lf)

mod_list[["mod_nd_notpu_j"]] <- mod_nd_notpu_j <- lme(j_diff ~ 1,
                      random = ~ 1 | treeid,
                      weights = varIdent(form = ~ 1 | treeid),
                      #control = lmeControl(opt = "optim"),
                      data = nd_diff_notpu_lf)



# mod_nd_notpu_v <- lmer(vc_diff ~ 1 + (1|treeid),
#                     data = nd_diff_notpu_lf) #This leads to a singularity
# mod_nd_notpu_j <- lmer(j_diff ~ 1 + (1|treeid),
#                     data = nd_diff_notpu_lf)


mod_list[["mod_nd_tpu_v"]] <- mod_nd_tpu_v <- lme(vc_diff ~ 1,
                      random = ~ 1 | treeid,
                      weights = varIdent(form = ~ 1 | treeid),
                      #control = lmeControl(opt = "optim"),
                      data = nd_diff_tpu_lf)

mod_list[["mod_nd_tpu_j"]] <- mod_nd_tpu_j <- lme(j_diff ~ 1,
                      random = ~ 1 | treeid,
                      weights = varIdent(form = ~ 1 | treeid),
                      #control = lmeControl(opt = "optim"),
                      data = nd_diff_tpu_lf)

# 
# mod_nd_tpu_v <- lmer(vc_diff ~ 1 + (1|treeid),
#                   data = nd_diff_tpu_lf)
# mod_nd_tpu_j <- lmer(j_diff ~ 1 + (1|treeid),
#                   data = nd_diff_tpu_lf)



###


# Check the model assumptions -----------------------------------------------

# Residuals for the full dataset

## Vcmax without TPU
plot(residuals(mod_notpu_v, type = "normalized") ~ fitted(mod_notpu_v)) # ok
plot(residuals(mod_notpu_v, type = "normalized") ~ factor(diff_notpu_lf$treeid))
abline(h = 0, 
       lty = 2, 
       col = "red") ## ok
hist(residuals(mod_notpu_v, type = "normalized")) # ok


plot(mod_notpu_j)
plot(residuals(mod_notpu_j, type = "normalized") ~ diff_notpu_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red") # ok
hist(residuals(mod_notpu_j, type = "normalized")) # ok


plot(mod_tpu_v) ### ok
plot(residuals(mod_tpu_v, type = "normalized") ~ diff_tpu_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red")# ok
hist(residuals(mod_tpu_v, type = "normalized")) # ok


plot(mod_tpu_j) ## not awesome
plot(residuals(mod_tpu_j, type = "normalized") ~ diff_tpu_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red") # not awesome
hist(residuals(mod_tpu_j, type = "normalized")) # skewed




## Residuals for models without K6709L6

plot(mod_notpu_nol6_v) # not great
plot(residuals(mod_notpu_nol6_v, type = "normalized") ~ diff_notpu_nol6_lf$treeid) #K6707 is a bit weird
abline(h = 0, 
       lty = 2, 
       col = "red") # ok
hist(residuals(mod_notpu_nol6_v, type = "normalized")) # ok


plot(mod_notpu_nol6_j) # ok
plot(residuals(mod_notpu_nol6_j, type = "normalized") ~ diff_notpu_nol6_lf$treeid) # K6707 is a bit weird
abline(h = 0, 
       lty = 2, 
       col = "red") # ok
hist(residuals(mod_notpu_nol6_j, type = "normalized")) # okayish


plot(mod_tpu_nol6_v) ### ok
plot(residuals(mod_tpu_nol6_v, type = "normalized") ~ diff_tpu_nol6_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red") # ok
hist(residuals(mod_tpu_nol6_v, type = "normalized")) # okayish


plot(mod_tpu_nol6_j) ## not great
plot(residuals(mod_tpu_nol6_j, type = "normalized") ~ diff_tpu_nol6_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red") # not great
hist(residuals(mod_tpu_nol6_j, type = "normalized")) # meh




## Residuals for no-overshoot subset

plot(mod_nd_notpu_v) ### ok
plot(residuals(mod_nd_notpu_v, type = "normalized") ~ nd_diff_tpu_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red") # ok
hist(residuals(mod_nd_notpu_v, type = "normalized")) # okayish


plot(mod_nd_notpu_j) ### okayish
plot(residuals(mod_nd_notpu_j, type = "normalized") ~ nd_diff_tpu_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red") # ok
hist(residuals(mod_nd_notpu_j)) # okayish


plot(mod_nd_tpu_v) ### ok
plot(residuals(mod_nd_tpu_v, type = "normalized") ~ nd_diff_tpu_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red") # ok
hist(residuals(mod_nd_tpu_v)) # ok


plot(mod_nd_tpu_j) ## non-linear
plot(residuals(mod_nd_tpu_j, type = "normalized") ~ nd_diff_tpu_lf$treeid)
abline(h = 0, 
       lty = 2, 
       col = "red") # okayish
hist(residuals(mod_nd_tpu_j, type = "normalized")) # ok



## Residuals for TPU vs TPU comparisons (Note small sample sizes!)

plot(tpu_only_mod) ## ok
# plot(residuals(tpu_only_mod, type = "normalized") ~ diff_tpu_comparison$treeid) ### doesn't work??
abline(h = 0, 
       lty = 2, 
       col = "red")
hist(residuals(tpu_only_mod)) # not great, but also tiny sample size



###

# Summary of the models ------------------

# Complete dataset summary
summary(mod_notpu_v)
summary(mod_notpu_j)

summary(mod_tpu_v)
summary(mod_tpu_j)

# No K6709L6 dataset summary
summary(mod_notpu_nol6_v)
summary(mod_notpu_nol6_j)

summary(mod_tpu_nol6_v)
summary(mod_tpu_nol6_j)

# No overshoot dataset summary
summary(mod_nd_notpu_v)
summary(mod_nd_notpu_j)

summary(mod_nd_tpu_v)
summary(mod_nd_tpu_j)

# Only TPU comparison dataset summary
summary(tpu_only_mod)

###


# Bootstrap model results ---------------------------------------------------

## Load the source code for bootstrapping
source(here("4_Fitting_stats_figs_scripts/bootstrapping_estimates.R"))

## Define the number of iterations
iter <- 500

## Create empty dataframes for the bootstrapped estimates and null estimates
estims_boot <- matrix(nrow = iter, ncol = length(mod_list))
null_boot <- matrix(nrow = iter, ncol = length(mod_list))
colnames(estims_boot) <- names(mod_list)
colnames(null_boot) <- names(mod_list)

## Apply the bootstrapped estimates to all models in the list
for(mod in 1:length(mod_list)){
    if (stringr::str_detect(names(mod_list[mod]), "v")){
        resp_var <- "vc_diff"
    } else if (stringr::str_detect(names(mod_list[mod]), "j")) {
        resp_var <- "j_diff"
    } else {
        resp_var <- "tpu_diff"
    }
    estims_boot[,names(mod_list)[mod]] <- sapply(
        1:iter,
        boot,
        mfit = mod_list[[mod]],
        resp_var = resp_var,
        calc_null = FALSE
    )
    null_boot[,names(mod_list)[mod]] <- sapply(
        1:iter,
        boot,
        mfit = mod_list[[mod]],
        resp_var = resp_var,
        calc_null = TRUE
    )
}

estims_boot <- as.data.frame(estims_boot)
null_boot <- as.data.frame(null_boot)


write.csv(estims_boot, file = here("5_Results/estimates_boot.csv"), row.names = FALSE)
write.csv(null_boot, file = here("5_Results/null_boot.csv"), row.names = FALSE)





# plot density of bootstrapped estimates against the
# theoretical sampling distribution

for(col in 1:ncol(estims_boot)){
    mod <- get(colnames(estims_boot)[col])
    
    hist(estims_boot[,col], freq = F, xlab = colnames(estims_boot)[col])
    
    # theoretical density
    m <-mod$coefficients$fixed
    se <- sqrt(mod$varFix[1,1])
    
    lines(
        x = seq(min(estims_boot[,col]) - 1, max(estims_boot[,col]) + 1, length.out = 100),
        y = dnorm(seq(min(estims_boot[,col]) - 1, max(estims_boot[,col]) + 1, length.out = 100),
                  mean = m, sd = se),
        col = "blue"
    )
    
}



mod_coefs <- names(mod_list) %>%
    map_dfr(~ {
        model_name <- .x
        model <- get(model_name)
        print(paste("Processing model:", model_name))
        model_summ <- summary(model)
        og_coeff <- model_summ$tTable %>% 
            as.data.frame() 
           # rownames_to_column()
        int_boot <- mean(estims_boot[,model_name])
        p_val_boot <- mean(abs(mean(estims_boot[,model_name])) <= null_boot[,model_name] ) * 2
        
        confints_boot <- t(quantile(estims_boot[,model_name], 
                                    probs = c(0.025, 0.975))) %>% 
            as.data.frame()
            # rownames_to_column()
            # filter(rowname == '(Intercept)') %>% #These are the fixed effect CIs
            # select(-rowname)
        icc_value <- icc(model) %>%
            as.data.frame() %>% 
            rownames_to_column()
        #std.d.ranef <- VarCorr(model) %>%
        std.d.ranef <- VarCorr(model)["(Intercept)",] %>%
            t() %>% 
            as.data.frame() %>% 
            rename(Ranef.Var = Variance, Ranef.StdDev = StdDev)
            # rownames_to_column() %>%
            # filter("rowname" == 1) %>%
            # select(-c(var1,var2, rowname))

        tibble(
            model = model_name,
            coeff = list(og_coeff),
            int_boot = int_boot,
            p_val_boot = p_val_boot,
            'ConfInt.Boot.2.5%' = confints_boot[,1],
            'ConfInt.Boot.97.5%' = confints_boot[,2],
            icc = icc_value[,2],
           # stdrandef = std.d.ranef
           Ranef.Var = std.d.ranef[,1],
           Ranef.StdDev = std.d.ranef[,2]
        )
    }) %>% 
    unnest_wider(coeff)
    
mod_coefs$Ranef.Var <- as.numeric(mod_coefs$Ranef.Var)
mod_coefs$Ranef.StdDev <- as.numeric(mod_coefs$Ranef.StdDev)

mod_coefs <- mod_coefs %>%
    mutate(across(where(is.numeric), ~ signif(., 3)))



write.csv(mod_coefs, here("5_Results/boot_res.csv"))






# # bootstrapped estimate
# mean(estims_boot)
# 
# # bootstrapped CI
# quantile(estims_boot, probs = c(0.025, 0.975))



# Checking Assumptions of paired t-tests -----------------------------------------

# Levene's for homogeneity of variance is not necessary because the paired t-test cares only about the distribution of the difference between the two variables. The variances of the two variables separately is irrelevant.

# Checking the symmetric distribution of the differences -----------------------------------------

# Includes Shapiro Wilk normality test for paired differences; no overshoot, WITHOUT TPU

## Vcmax WITHOUT TPU
gghistogram(diff_notpu_lf, x = "vc_diff", bins = 10, add_density = TRUE) + geom_vline(xintercept = median(diff_notpu_lf$vc_diff), color = 'red')
#That should be fine

diff_notpu_lf %>%
    with(., shapiro.test(vc_diff))

## Vcmax WITH TPU
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


# Sign test here.
all_avg_tr_res %>%
    group_by(fit_type) %>%
    sign_test(data =., vcmax ~ curv_meth, detailed = TRUE, ref.group = 'SS') %>%
    add_significance()



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


# Sign test here.
all_avg_tr_res %>%
    group_by(fit_type) %>%
    sign_test(data =., jmax ~ curv_meth, detailed = TRUE, ref.group = 'SS') %>%
    add_significance()


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


# Sign test here.
all_avg_tr_res %>%
    group_by(curv_meth) %>%
    sign_test(data =., vcmax ~ fit_type, detailed = TRUE) %>%
    add_significance()



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


# Sign test here.
all_avg_tr_res %>%
    group_by(curv_meth) %>%
    sign_test(data =., jmax ~ fit_type, detailed = TRUE) %>%
    add_significance()


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


# Sign test here.
all_nol6_tr_res %>%
    group_by(fit_type) %>%
    sign_test(data =., vcmax ~ curv_meth, detailed = TRUE, ref.group = 'SS') %>%
    add_significance()




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


# Sign test here.
all_nol6_tr_res %>%
    group_by(fit_type) %>%
    sign_test(data =., jmax ~ curv_meth, detailed = TRUE, ref.group = 'SS') %>%
    add_significance()



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

# Sign test here.
all_nol6_tr_res %>%
    group_by(curv_meth) %>%
    sign_test(data =., vcmax ~ fit_type, detailed = TRUE) %>%
    add_significance()


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


# Sign test here.
all_nol6_tr_res %>%
    group_by(curv_meth) %>%
    sign_test(data =., jmax ~ fit_type, detailed = TRUE) %>%
    add_significance()



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

# Sign test here.
nd_tr_res %>%
    group_by(fit_type) %>%
    sign_test(data =., vcmax ~ curv_meth, detailed = TRUE, ref.group = 'SS') %>%
    add_significance()




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


# Sign test here.
nd_tr_res %>%
    group_by(fit_type) %>%
    sign_test(data =., jmax ~ curv_meth, detailed = TRUE, ref.group = 'SS') %>%
    add_significance()



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


# Sign test here.
nd_tr_res %>%
    group_by(curv_meth) %>%
    sign_test(data =., vcmax ~ fit_type, detailed = TRUE) %>%
    add_significance()




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

# Sign test here.
nd_tr_res %>%
    group_by(curv_meth) %>%
    sign_test(data =., vcmax ~ fit_type, detailed = TRUE) %>%
    add_significance()





# Combine Wilcox models into a table --------------------------------------


# List of models
wil_cm <- list(
    w_vc_cm_full,
    w_j_cm_full) %>%
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


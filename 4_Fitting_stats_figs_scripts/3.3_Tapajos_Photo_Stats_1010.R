# Computing statistical analysis on the A/Ci curve fits from 2_Tapajos_Fit_ACi.R
# Using the 'photosynthesis' package

######### Check to see which of these packages we actually need
library(tidyverse) 
library(ggpubr) 
library(car)
library(rstatix)
library(reshape2) 
library(here) 


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
# Describe factor levels: 0 is SS, 1 is DAT
# pho_stat$curv_meth <- factor(pho_stat$curv_meth) ########## Is this necessary?

write.csv(pho_stat, file = here("3_Clean_data/pho_stat.csv"))

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
# Describe factor levels: 0 is SS, 1 is DAT
# pho_stat_tpu$curv_meth <- factor(pho_stat_tpu$curv_meth)

write.csv(pho_stat_tpu, file = here("3_Clean_data/pho_stat_tpu.csv"))



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
#Describe factor levels: 0 is SS, 1 is DAT
# pho_nd_stat$curv_meth <- factor(pho_nd_stat$curv_meth)

write.csv(pho_nd_stat, file = here("3_Clean_data/pho_nd_stat.csv"))



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
#Describe factor levels: 0 is SS, 1 is DAT
# pho_nd_stat_tpu$curv_meth <- factor(pho_nd_stat_tpu$curv_meth)

write.csv(pho_nd_stat_tpu, file = here("3_Clean_data/pho_nd_stat_tpu.csv"))

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
#all_avg_lf_res$fit_type <- factor(all_avg_lf_res$fit_type)


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
write.csv(diff_notpu_lf, here("3_Clean_data/lf_diffs_summ_noTPU.csv"))


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
#all_avg_lf_res$fit_type <- factor(all_avg_lf_res$fit_type)


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
nd_diff_tpu_tree <- nd_diff_tpu_tree[order(nd_diff_tpu_tree$rel_can_pos, decreasing = TRUE),]

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
# nd_complete$curv_meth <- factor(nd_complete$curv_meth)
# nd_complete$fit_type <- factor(nd_complete$fit_type)






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
#tpu_just6_all$curv_meth <- factor(tpu_just6_all$curv_meth)


# Produce Summary statistics of results ----------------------

### All data
all_avg_lf_res %>% 
    group_by(curv_meth, fit_type) %>%
    get_summary_stats(vcmax, jmax, show = c("mean", "median", "sd", "min", "max"))


### No-overshoot data
nd_complete %>% 
    group_by(curv_meth, fit_type) %>% 
    get_summary_stats(vcmax, jmax, show = c("mean", "median", "sd", "min", "max"))


### Just the TPU data
tpu_just6_all %>% 
    group_by(curv_meth)%>%
    get_summary_stats(vcmax, jmax, show = c("mean", "median", "sd", "min", "max"))


# Checking Assumptions of t-tests -----------------------------------------


### Looking at histograms for normality: WITHOUT TPU, all data

### Vcmax
all_avg_lf_res %>% 
    filter(fit_type == "no_tpu") %>%
    ggplot(aes(x=vcmax)) + 
    geom_histogram()

### Jmax
all_avg_lf_res %>% 
    filter(fit_type == "no_tpu") %>%
    ggplot(aes(x=jmax)) + 
    geom_histogram()



# Levene's for homogeneity of variance; WITHOUT TPU, all data

all_avg_lf_res %>%
    filter(fit_type == "no_tpu") %>% 
    leveneTest(vcmax ~ curv_meth, data = .)

all_avg_lf_res %>% 
    filter(fit_type == "no_tpu") %>% 
    leveneTest(jmax ~ curv_meth, data = .)

#### Both non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances



# Shapiro-Wilk test for normality; WITHOUT TPU, all data

# Vcmax
all_avg_lf_res %>%
    filter(fit_type == "no_tpu") %>%
    with(., shapiro.test(vcmax))

all_avg_lf_res %>% 
    filter(fit_type == "no_tpu") %>%
    with(., shapiro.test(vcmax[curv_meth == "DAT"]))

all_avg_lf_res %>%
    filter(fit_type == "no_tpu") %>% 
    with(., shapiro.test(vcmax[curv_meth == "SS"]))

# Jmax
all_avg_lf_res %>% 
    filter(fit_type == "no_tpu") %>% 
    with(., shapiro.test(jmax))

all_avg_lf_res %>%
    filter(fit_type == "no_tpu") %>%
    with(., shapiro.test(jmax[curv_meth == "DAT"]))

all_avg_lf_res %>% 
    filter(fit_type == "no_tpu") %>% 
    with(., shapiro.test(jmax[curv_meth == "SS"]))

#All but one deviates from normal




# Looking at histograms for normality, WITH TPU, all data

all_avg_lf_res %>% 
    filter(fit_type == "tpu") %>%
    ggplot(aes(x=vcmax)) + 
    geom_histogram()

all_avg_lf_res %>% 
    filter(fit_type == "tpu") %>%
    ggplot(aes(x=jmax)) + 
    geom_histogram()



# Levene's for homogeneity of variance; WITH TPU, all data
all_avg_lf_res %>% 
    filter(fit_type == "tpu") %>% 
    leveneTest(vcmax ~ curv_meth, data = .)

all_avg_lf_res %>%
    filter(fit_type == "tpu") %>% 
    leveneTest(jmax ~ curv_meth, data = .)

# Both non-significant. The variances are not significantly different from one another. Therefore we have met the assumption of equal variances



# Shapiro-Wilk's test for normality, WITH TPU, all data

### Vcmax
all_avg_lf_res %>% 
    filter(fit_type == "tpu") %>% 
    with(., shapiro.test(vcmax))

all_avg_lf_res %>%
    filter(fit_type == "tpu") %>% 
    with(., shapiro.test(vcmax[curv_meth == "DAT"]))

all_avg_lf_res %>% 
    filter(fit_type == "tpu") %>%
    with(., shapiro.test(vcmax[curv_meth == "SS"]))

### Jmax
all_avg_lf_res %>% 
    filter(fit_type == "tpu") %>%
    with(., shapiro.test(jmax))

all_avg_lf_res %>%
    filter(fit_type == "tpu") %>% 
    with(., shapiro.test(jmax[curv_meth == "DAT"]))

all_avg_lf_res %>%
    filter(fit_type == "tpu") %>%
    with(., shapiro.test(jmax[curv_meth == "SS"]))
#Most deviate from normal




# Testing assumptions for No Overshoot data

# Levene's test; no overshoot

## WITH TPU
nd_complete %>%
    filter(fit_type == "tpu") %>% 
    leveneTest(vcmax ~ curv_meth, data = .)

nd_complete %>%
    filter(fit_type == "tpu") %>%
    leveneTest(jmax ~ curv_meth, data = .)

# WITHOUT TPU
nd_complete %>% 
    filter(fit_type == "no_tpu") %>% 
    leveneTest(vcmax ~ curv_meth, data = .)

nd_complete %>%
    filter(fit_type == "no_tpu") %>% 
    leveneTest(jmax ~ curv_meth, data = .)

# Met assumption of homogeneity of variances



# Shapiro Wilk Test; no overshoot, WITHOUT TPU

### Vcmax
nd_complete %>%
    filter(fit_type == "no_tpu") %>% 
    with(., shapiro.test(vcmax))

nd_complete %>% 
    filter(fit_type == "no_tpu") %>%
    with(., shapiro.test(vcmax[curv_meth == "DAT"]))

nd_complete %>% 
    filter(fit_type == "no_tpu") %>% 
    with(., shapiro.test(vcmax[curv_meth == "SS"]))


### Jmax
nd_complete %>% 
    filter(fit_type == "no_tpu") %>%
    with(., shapiro.test(jmax))

nd_complete %>% 
    filter(fit_type == "no_tpu") %>%
    with(., shapiro.test(jmax[curv_meth == "DAT"]))

nd_complete %>% 
    filter(fit_type == "no_tpu") %>% 
    with(., shapiro.test(jmax[curv_meth == "SS"]))



### Shapiro-Wilk test; no overshoot, WITH TPU

### Vcmax
nd_complete %>% 
    filter(fit_type == "tpu") %>% 
    with(., shapiro.test(vcmax))

nd_complete %>%
    filter(fit_type == "tpu") %>% 
    with(., shapiro.test(vcmax[curv_meth == "DAT"]))

nd_complete %>% 
    filter(fit_type == "tpu") %>%
    with(., shapiro.test(vcmax[curv_meth == "SS"]))

### Jmax
nd_complete %>%
    filter(fit_type == "tpu") %>% 
    with(., shapiro.test(jmax))

nd_complete %>%
    filter(fit_type == "tpu") %>%
    with(., shapiro.test(jmax[curv_meth == "DAT"]))

nd_complete %>% 
    filter(fit_type == "tpu") %>% 
    with(., shapiro.test(jmax[curv_meth == "SS"]))

# They all deviate from normal.


# Checking the symmetric distribution of Wilcoxon -----------------------------------------

# To use wilcoxon paired, we assume the differences between paired samples are distributed symmetrically about the median


### Combine the difference datasets
#diff_lf <- rbind(diff_notpu_lf, diff_tpu_lf)



gghistogram(diff_notpu_lf, x = "vc_diff", bins = 10, add_density = TRUE)
#That should be fine

gghistogram(diff_tpu_lf, x = "vc_diff", bins = 10, add_density = TRUE)
#okay


### Jmax
gghistogram(diff_notpu_lf, x = "j_diff", bins = 10, add_density = TRUE)

gghistogram(diff_tpu_lf, x = "j_diff", bins = 10, add_density = TRUE)
#### not great!; far worse than Vcmax differences



# No Overshoot data

#nd_diff_lf <- rbind(nd_diff_notpu_lf, nd_diff_tpu_lf)

### Vcmax
gghistogram(nd_diff_notpu_lf, x = "vc_diff", bins = 10, add_density = TRUE)
#That should be fine

gghistogram(nd_diff_tpu_lf, x = "vc_diff", bins = 10, add_density = TRUE)
#okay


### Jmax
gghistogram(nd_diff_notpu_lf, x = "j_diff", bins = 10, add_density = TRUE)

gghistogram(nd_diff_tpu_lf, x = "j_diff", bins = 10, add_density = TRUE)
#### not very good; Vcmax is much better


# It's just Jmax where we don't totally meet the Wilcoxon assumptions. 
# We will run both Sign and Wilcoxon for both of these, and make a note of this.



# Run Wilcoxon and Sign tests ------------------------------------------


# All data, Wilcoxon signed rank test on paired samples (and a few sign tests)

### Vcmax Wilcoxon by curve method 
all_avg_lf_res %>%
    group_by(fit_type) %>%
    wilcox_test(data =., vcmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    add_significance()

## Effect size
all_avg_lf_res %>%
    group_by(fit_type) %>% 
    wilcox_effsize(data = ., vcmax ~ curv_meth, paired = TRUE)


### Vcmax Wilcoxon by TPU v. no TPU
all_avg_lf_res %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>% 
    add_significance()

all_avg_lf_res %>%
    group_by(curv_meth) %>% 
    wilcox_effsize(data = ., vcmax ~ fit_type, paired = TRUE)





### Jmax Wilcoxon by curve method
all_avg_lf_res %>%
    group_by(fit_type) %>%
    wilcox_test(data =., jmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    add_significance()

# Note we're running the sign test here in addition!
all_avg_lf_res %>%
    group_by(fit_type) %>%
    sign_test(data =., jmax ~ curv_meth, detailed = TRUE) %>%
    add_significance()

all_avg_lf_res %>%
    group_by(fit_type) %>% 
    wilcox_effsize(data = ., jmax ~ curv_meth, paired = TRUE)



### Jmax Wilcoxon by TPU v. no TPU
all_avg_lf_res %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., jmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    add_significance()

# Sign test here. Note we have a different result.
all_avg_lf_res %>%
    group_by(curv_meth) %>%
    sign_test(data =., jmax ~ fit_type, detailed = TRUE) %>%
    add_significance()

all_avg_lf_res %>%
    group_by(curv_meth) %>% 
    wilcox_effsize(data = ., jmax ~ fit_type, paired = TRUE)




# Wilcoxon Tests, no Overshoot subset


### Vcmax Wilcoxon by curve method
nd_complete %>%
    group_by(fit_type) %>%
    wilcox_test(data =., vcmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    add_significance()

nd_complete %>%
    group_by(fit_type) %>% 
    wilcox_effsize(data = ., vcmax ~ curv_meth, paired = TRUE)


### Vcmax Wilcoxon by TPU v. no TPU

### Note that all the SS curves for which TPU was fit had overshoot, so the datasets are same
nd_complete %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., vcmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    add_significance()




### Jmax Wilcoxon by curve method
nd_complete %>%
    group_by(fit_type) %>%
    wilcox_test(data =., jmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    add_significance()

nd_complete %>%
    group_by(fit_type) %>%
    wilcox_effsize(data = ., jmax ~ curv_meth, paired = TRUE)


### Jmax Wilcoxon by TPU v. no TPU

### Note that all the SS curves for which TPU was fit had overshoot, so the datasets are same
nd_complete %>%
    group_by(curv_meth) %>%
    wilcox_test(data =., jmax ~ fit_type, paired = TRUE, detailed = TRUE) %>%
    add_significance()

#####


# TPU Wilcoxon

tpu_just6_all %>% 
    wilcox_test(data = ., tpu ~ curv_meth, paired = TRUE, detailed = TRUE) %>% 
    add_significance()

tpu_just6_all %>% 
    wilcox_effsize(data = ., tpu ~ curv_meth, paired = TRUE)


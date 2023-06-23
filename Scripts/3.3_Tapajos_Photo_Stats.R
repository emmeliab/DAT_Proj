
# Statistical Analysis of A/Ci curve data
rm(list = ls())

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

wd <- "C://Users/emmel/Desktop/DAT_proj/"
setwd(wd)

##Read in Datasets -------------------------------------------

pho_dat <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_no_TPU.csv"), sep = ",", 
                    header = TRUE, na.strings = 1000) ## TPU values at 1000 are coded as NA
pho_trad <- read.csv(file = paste0(wd, "Results/trad_fits_photo_pars_correct_no_TPU.csv"), sep = ",", 
                     header = TRUE, na.strings = 1000)


pho_dat_tpu <- read.csv(file = paste0(wd, "Results/dat_fits_photo_pars_filt_correct_with_TPU.csv"), sep = ",", 
                        header = TRUE, na.strings = 1000) ## TPU values at 1000 are coded as NA
pho_trad_tpu <- read.csv(file = paste0(wd, "Results/trad_fits_photo_pars_correct_with_TPU.csv"), sep = ",", 
                         header = TRUE, na.strings = 1000)

pho_both <- rbind(pho_dat, pho_trad)
pho_both_tpu <- rbind(pho_dat_tpu, pho_trad_tpu)

#No overshoot data

pho_nd <- pho_dat %>% subset(ID != "K6707L2" & ID != "K6707L2-2" & ID != "K6707L1" & ID != "K6707L1-1" & ID != "K6709L6" & ID != "K6714L2" & ID != "K6714L1" & ID != "K6702L1" & ID != "K6706L2" & ID != "K6706L1" & ID != "K6709L2")
pho_nd_tpu <- pho_dat_tpu %>% subset(ID != "K6707L2" & ID != "K6707L2-2" & ID != "K6707L1" & ID != "K6707L1-1" & ID != "K6709L6" & ID != "K6714L2" & ID != "K6714L1" & ID != "K6702L1" & ID != "K6706L2" & ID != "K6706L1" & ID != "K6709L2")

pho_nd_both <- rbind(pho_nd, pho_trad)
pho_nd_both_tpu <- rbind(pho_trad_tpu, pho_nd_tpu)


#Initial Processing/Grouping ----------------------------------

#Full data, no TPU
pho_leaf <- pho_both %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           curv_meth = Data_point,
           leaf_id = ID)
#CHECK! MAY NOT NEED 'V_TPU'
pho_stat <- select(pho_leaf, 'curv_meth', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU', 'leaf_id', 
                   'leaf_unique', "V_cmax_se", "J_se", "V_TPU_se")
pho_stat$curv_meth[pho_stat$curv_meth == "Before_DAT"] <- "DAT"
pho_stat <- rename(pho_stat,
                   vcmax = Best_Vcmax_25C,
                   jmax = Best_Jmax_25C,
                   tpu = V_TPU,) %>% 
    mutate(fit_type = "no_tpu")
#Describe factor levels: 0 is traditional, 1 is DAT
pho_stat$curv_meth <- factor(pho_stat$curv_meth)


#Full data, with TPU
pho_leaf_tpu <- pho_both_tpu %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           curv_meth = Data_point,
           leaf_id = ID)
pho_stat_tpu <- select(pho_leaf_tpu, 'curv_meth', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU', 
                       'leaf_id', 'leaf_unique', 'V_cmax_se', 'J_se', 'V_TPU_se')
pho_stat_tpu$curv_meth[pho_stat_tpu$curv_meth == "Before_DAT"] <- "DAT"
pho_stat_tpu <- rename(pho_stat_tpu,
                       vcmax = Best_Vcmax_25C,
                       jmax = Best_Jmax_25C,
                       tpu = V_TPU) %>% mutate(fit_type = "tpu")
#Describe factor levels: 0 is traditional, 1 is DAT
pho_stat_tpu$curv_meth <- factor(pho_stat_tpu$curv_meth)


all_results <- rbind(pho_stat, pho_stat_tpu)


#No overshoot, no TPU
pho_nd_leaf <- pho_nd_both %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           curv_meth = Data_point,
           leaf_id = ID)
pho_nd_stat <- select(pho_nd_leaf, 'curv_meth', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU', 'leaf_id', 
                      'leaf_unique', 'V_cmax_se', 'J_se')
pho_nd_stat$curv_meth[pho_nd_stat$curv_meth == "Before_DAT"] <- "DAT"
pho_nd_stat <- rename(pho_nd_stat,
                      vcmax = Best_Vcmax_25C,
                      jmax = Best_Jmax_25C,
                      tpu = V_TPU) %>% mutate(fit_type = "no_tpu")
#Describe factor levels: 0 is traditional, 1 is DAT
pho_nd_stat$curv_meth <- factor(pho_nd_stat$curv_meth)


#No overshoot, with TPU
pho_nd_leaf_tpu <- pho_nd_both_tpu %>%
    mutate(leaf_unique = substring(ID, 1, 7),
           curv_meth = Data_point,
           leaf_id = ID)
pho_nd_stat_tpu <- select(pho_nd_leaf_tpu, 'curv_meth', 'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU',
                          'leaf_id', 'leaf_unique', "V_cmax_se", "J_se")
pho_nd_stat_tpu$curv_meth[pho_nd_stat_tpu$curv_meth == "Before_DAT"] <- "DAT"
pho_nd_stat_tpu <- rename(pho_nd_stat_tpu,
                          vcmax = Best_Vcmax_25C,
                          jmax = Best_Jmax_25C,
                          tpu = V_TPU) %>% mutate(fit_type = "tpu")
#Describe factor levels: 0 is traditional, 1 is DAT
pho_nd_stat_tpu$curv_meth <- factor(pho_nd_stat_tpu$curv_meth)



#group data for further analysis -------------------------------------------

#Merge no-TPU with TPU
tpu_results_grp <- pho_stat_tpu %>%
    group_by(curv_meth, leaf_unique) %>% 
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>% 
    mutate(fit_type = "tpu")

notpu_results_grp <- pho_stat %>%
    group_by(curv_meth, leaf_unique) %>% 
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>%
    mutate(fit_type = "no_tpu")

all_results2 <- rbind(tpu_results_grp, notpu_results_grp)
all_results2$fit_type <- factor(all_results2$fit_type)


#Creating a table with Vcmax and Jmax differences by species ----------------------------
substring(all_results2$leaf_unique, "K67", 1)
treename <- as.numeric(substring(all_results2$leaf_unique, 4, 5))

all_results2$treename <- treename

all_res_dat_tpu <- all_results2 %>% filter(curv_meth == "DAT") %>% filter(fit_type == "tpu")
all_res_trad_tpu <- all_results2 %>% filter(curv_meth == "Traditional") %>% filter(fit_type == "tpu")

all_res_dat_tpu_summ <- all_res_dat_tpu %>% group_by(treename) %>%
    summarize(dat_vcmax = mean(vcmax),
              dat_vc_sd = sd(vcmax),
              dat_jmax = mean(jmax),
              dat_j_sd = sd(jmax))

all_res_trad_tpu_summ <- all_res_trad_tpu %>% group_by(treename) %>%
    summarize(trad_vcmax = mean(vcmax),
              trad_vc_sd = sd(vcmax),
              trad_jmax = mean(jmax),
              trad_j_sd = sd(jmax))


species_summ <- cbind(all_res_dat_tpu_summ, all_res_trad_tpu_summ)
species_summ$vc_diff <- species_summ$trad_vcmax - species_summ$dat_vcmax
species_summ$j_diff <- species_summ$trad_jmax - species_summ$dat_jmax
species_summ$vc_diff_se <- (abs(species_summ$trad_vc_sd) + abs(species_summ$dat_vc_sd))/2
species_summ$j_diff_se <- (abs(species_summ$trad_j_sd) + abs(species_summ$dat_j_sd))/2
species_summ2 <- species_summ[,-6]

species_summ2$treename <- as.character(species_summ2$treename)

rel_can_pos <- c(11, 5.1, 6, 4, 2, 3, 1, 10, 12, 8, 7, 9, 5.2)
species_summ2$rel_can_pos <- rel_can_pos

species_summ2 <- species_summ2[order(species_summ2$rel_can_pos, decreasing = TRUE),]

codebook <- read_csv("Results/id_codebook.csv") %>% arrange(desc(rel_can_pos)) %>% select(-c(overshoot, treeid, rel_can_pos))

species_summ3 <- cbind(species_summ2, codebook) %>% select(-15)

#Note that error bars represent the mean of the absolute error for the differences
vc_diff_hist <- ggplot(data = species_summ3, aes(x = reorder(gen_spec_id, desc(rel_can_pos)), y = vc_diff)) +
    geom_bar(stat="identity", fill = "cadetblue2", color = "grey20") +
    labs(x = NULL,
         y = "TPU-enabled Vcmax Differences (steady-state - DAT)")+
    geom_errorbar(aes(x=gen_spec_id, ymin=vc_diff-vc_diff_se, ymax=vc_diff+vc_diff_se), width=0.3, colour="#CA0068", alpha=0.9, linewidth=0.5)+
    theme_classic(base_family = "serif") +
    geom_hline(yintercept=0, linetype="solid", color="black", linewidth=0.8) +
    ylim(-20, 50)+
    theme(axis.text.y = element_text(face = "italic"))+
    coord_flip()
vc_diff_hist
ggsave(plot = vc_diff_hist, "Figures/vc_diff_hist.png")

j_diff_hist <- ggplot(data = species_summ3, aes(x = reorder(gen_spec_id, desc(rel_can_pos)), y = j_diff)) +
    geom_bar(stat="identity", fill = "cadetblue2", color = "grey20") +
    labs(x = NULL,
         y = "TPU-enabled Jmax Differences (steady-state - DAT)")+
    geom_errorbar(aes(x=gen_spec_id, ymin=j_diff-j_diff_se, ymax=j_diff+j_diff_se), width=0.3, colour="#CA0068", alpha=0.9, size=0.5)+
    theme_classic(base_family = "serif") +
    geom_hline(yintercept=0, linetype="solid", color="black", linewidth=0.8) +
    ylim(-20, 50) +
    theme(axis.text.y = element_text(face = "italic"))+
    coord_flip()
j_diff_hist
ggsave(plot = j_diff_hist, "Figures/j_diff_hist.png")


plot_arranged <- grid.arrange(vc_diff_hist, j_diff_hist)
ggsave(plot = plot_arranged, "Figures/diff_histos.png", width = 4.3, height = 7)

write.csv(species_summ3, "Results/species_diffs_summary_tpu.csv")


#Now for non-TPU fit curves
all_res_dat_notpu <- all_results2 %>% filter(curv_meth == "DAT") %>% filter(fit_type == "no_tpu")
all_res_trad_notpu <- all_results2 %>% filter(curv_meth == "Traditional") %>% filter(fit_type == "no_tpu")

all_res_dat_notpu_summ <- all_res_dat_notpu %>% group_by(treename) %>%
    summarize(dat_vcmax = mean(vcmax),
              dat_vc_sd = sd(vcmax),
              dat_jmax = mean(jmax),
              dat_j_sd = sd(jmax))

all_res_trad_notpu_summ <- all_res_trad_notpu %>% group_by(treename) %>%
    summarize(trad_vcmax = mean(vcmax),
              trad_vc_sd = sd(vcmax),
              trad_jmax = mean(jmax),
              trad_j_sd = sd(jmax))


species_summ_notpu <- cbind(all_res_dat_notpu_summ, all_res_trad_notpu_summ)
species_summ_notpu$vc_diff <- species_summ_notpu$trad_vcmax - species_summ_notpu$dat_vcmax
species_summ_notpu$j_diff <- species_summ_notpu$trad_jmax - species_summ_notpu$dat_jmax
species_summ_notpu$vc_diff_se <- (abs(species_summ_notpu$trad_vc_sd) + abs(species_summ_notpu$dat_vc_sd))/2
species_summ_notpu$j_diff_se <- (abs(species_summ_notpu$trad_j_sd) + abs(species_summ_notpu$dat_j_sd))/2
species_summ2_notpu <- species_summ_notpu[,-6]

species_summ2_notpu$treename <- as.character(species_summ2_notpu$treename)

rel_can_pos <- c(11, 5.1, 6, 4, 2, 3, 1, 10, 12, 8, 7, 9, 5.2)
species_summ2_notpu$rel_can_pos <- rel_can_pos

species_summ2_notpu <- species_summ2_notpu[order(species_summ2_notpu$rel_can_pos, decreasing = TRUE),]

codebook <- read_csv("Results/id_codebook.csv") %>% arrange(desc(rel_can_pos)) %>% select(-c(overshoot, treeid, rel_can_pos))

species_summ3_notpu <- cbind(species_summ2_notpu, codebook) %>% select(-15)

write.csv(species_summ3_notpu, "Results/species_diffs_summary_notpu.csv")

vc_diff_hist_notpu <- ggplot(data = species_summ3_notpu, aes(x = reorder(gen_spec_id, desc(rel_can_pos)), y = vc_diff)) +
    geom_bar(stat="identity", fill = "cadetblue2", color = "grey20") +
    labs(x = NULL,
         y = "No TPU \U0394Vcmax (steady-state - DAT)")+
    geom_errorbar(aes(x=gen_spec_id, ymin=vc_diff-vc_diff_se, ymax=vc_diff+vc_diff_se), width=0.3, colour="#CA0068", alpha=0.9, size=0.5)+
    theme_classic(base_family = "serif") +
    geom_hline(yintercept=0, linetype="solid", color="black", linewidth=0.8) +
    ylim(-20, 30)+
    theme(axis.text.y = element_text(face = "italic", size = rel(1.5)),
          axis.text.x = element_text(size = rel(1.2)),
          axis.title.x = element_text(size = rel(1.5)))+
    coord_flip()
vc_diff_hist_notpu
ggsave(plot = vc_diff_hist_notpu, "Figures/vc_diff_hist.png")

j_diff_hist_notpu <- ggplot(data = species_summ3_notpu, aes(x = reorder(gen_spec_id, desc(rel_can_pos)), y = j_diff)) +
    geom_bar(stat="identity", fill = "cadetblue2", color = "grey20") +
    labs(x = NULL,
         y = "No TPU Jmax Differences (steady-state - DAT)")+
    geom_errorbar(aes(x=gen_spec_id, ymin=j_diff-j_diff_se, ymax=j_diff+j_diff_se), width=0.3, colour="#CA0068", alpha=0.9, size=0.5)+
    theme_classic(base_family = "serif") +
    geom_hline(yintercept=0, linetype="solid", color="black", linewidth=0.8) +
    ylim(-20, 50) +
    theme(axis.text.y = element_text(face = "italic"))+
    coord_flip()
j_diff_hist_notpu
ggsave(plot = j_diff_hist_notpu, "Figures/j_diff_hist.png")

plot_arranged2 <- grid.arrange(vc_diff_hist_notpu, j_diff_hist_notpu)
ggsave(plot = plot_arranged2, "Figures/diff_histos_notpu.png", width = 4.3, height = 7)

gA <- ggplotGrob(vc_diff_hist_notpu)
gB <- ggplotGrob(vc_diff_hist + rremove("y.text"))
gC <- ggplotGrob(j_diff_hist_notpu)
gD <- ggplotGrob(j_diff_hist + rremove("y.text"))

plot_arranged3 <- grid.arrange(arrangeGrob(cbind(gA, gB), arrangeGrob(cbind(gC, gD))))
ggsave(plot = plot_arranged3, "Figures/diff_full_fig.png", width = 8.3, height = 5)


#understanding mean differences --------------------------
all_res_dat_tpu <- all_results2 %>%
    filter(curv_meth == "DAT") %>%
    filter(fit_type == "tpu") %>%
    mutate(dat_vcmax = vcmax,
           dat_jmax = jmax) %>%
    select(-c(vcmax, jmax))
all_res_trad_tpu <- all_results2 %>%
    filter(curv_meth == "Traditional") %>%
    filter(fit_type == "tpu") %>%
    mutate(trad_vcmax = vcmax,
           trad_jmax = jmax) %>%
    select(-c(vcmax, jmax))

res_tpu_summ <- cbind(all_res_dat_tpu, all_res_trad_tpu) %>% select(-c(1, 7, 8, 9, 10)) %>%
    rename(leaf_unique = leaf_unique...2,
           fit_type = fit_type...3,
           treename = treename...4)

res_tpu_summ$vc_diff <- res_tpu_summ$trad_vcmax - res_tpu_summ$dat_vcmax
res_tpu_summ$j_diff <- res_tpu_summ$trad_jmax - res_tpu_summ$dat_jmax

#no tpu
all_res_dat_notpu <- all_results2 %>%
    filter(curv_meth == "DAT") %>%
    filter(fit_type == "no_tpu") %>%
    mutate(dat_vcmax = vcmax,
           dat_jmax = jmax) %>%
    select(-c(vcmax, jmax))
all_res_trad_notpu <- all_results2 %>%
    filter(curv_meth == "Traditional") %>%
    filter(fit_type == "no_tpu") %>%
    mutate(trad_vcmax = vcmax,
           trad_jmax = jmax) %>%
    select(-c(vcmax, jmax))

res_notpu_summ <- cbind(all_res_dat_notpu, all_res_trad_notpu) %>% select(-c(1, 7, 8, 9, 10)) %>%
    rename(leaf_unique = leaf_unique...2,
           fit_type = fit_type...3,
           treename = treename...4)

res_notpu_summ$vc_diff <- res_notpu_summ$trad_vcmax - res_notpu_summ$dat_vcmax
res_notpu_summ$j_diff <- res_notpu_summ$trad_jmax - res_notpu_summ$dat_jmax


#mean differences on leaf basis
mean(res_tpu_summ$vc_diff)
mean(res_tpu_summ$j_diff)
mean(res_notpu_summ$vc_diff)
mean(res_notpu_summ$j_diff)

#mean differences on individual tree basis
mean(species_summ3$vc_diff)
mean(species_summ3$j_diff)
mean(species_summ3_notpu$vc_diff)
mean(species_summ3_notpu$j_diff)

#Difference histograms on a species basis

sp_diff_hist_vc_notpu <- species_summ3_notpu %>% ggplot(aes(x = vc_diff)) +
    stat_function(fun = dnorm, args = list(0, sd(species_summ3_notpu$vc_diff)), color = "red", linetype = "dashed")+
    geom_density(linewidth = 0.8)+
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3)+
    geom_vline(xintercept = mean(species_summ3_notpu$vc_diff), color = "black", alpha = 0.4)+
    xlab("\U0394Vcmax: No TPU") +
    ylab("Density")+
    xlim(-20,50)+
    ylim(0, 0.25)+
    theme_classic()+
    annotate("text", x = 30, y = 0.25, label = paste0("Mean = ", round(mean(species_summ3_notpu$vc_diff), digits = 2)), size = rel(2.7))+
    annotate("text", x = 30, y = 0.23, label = paste0("SD = ", round(sd(species_summ3_notpu$vc_diff), digits = 2)), size = rel(2.7))
sp_diff_hist_vc_notpu

sp_diff_hist_vc_tpu <- species_summ3 %>% ggplot(aes(x = vc_diff)) +
    stat_function(fun = dnorm, args = list(0, sd(species_summ3$vc_diff)), color = "red", linetype = "dashed")+
    geom_density(linewidth = 0.8)+
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3)+
    geom_vline(xintercept = mean(species_summ3$vc_diff), color = "black", alpha = 0.4)+
    xlab("\U0394Vcmax: TPU-enabled") +
    ylab("Density")+
    xlim(-20,50)+
    ylim(0, 0.25)+
    theme_classic()+
    annotate("text", x = 30, y = 0.25, label = paste0("Mean = ", round(mean(species_summ3$vc_diff), digits = 2)), size = rel(2.7))+
    annotate("text", x = 30, y = 0.23, label = paste0("SD = ", round(sd(species_summ3$vc_diff), digits = 2)), size = rel(2.7))
sp_diff_hist_vc_tpu

sp_diff_hist_j_notpu <- species_summ3_notpu %>% ggplot(aes(x = j_diff)) +
    stat_function(fun = dnorm, args = list(0, sd(species_summ3_notpu$j_diff)), linetype = "dashed", color = "red")+
    geom_density(linewidth = 0.8)+
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3)+
    geom_vline(xintercept = mean(species_summ3_notpu$j_diff), color = "black", alpha = 0.4)+
    xlab("\U0394J: No TPU") +
    ylab("Density")+
    xlim(-20,50)+
    ylim(0, 0.25)+
    theme_classic()+
    annotate("text", x = 30, y = 0.25, label = paste0("Mean = ", round(mean(species_summ3_notpu$j_diff), digits = 2)), size = rel(2.7))+
    annotate("text", x = 30, y = 0.23, label = paste0("SD = ", round(sd(species_summ3_notpu$j_diff), digits = 2)), size = rel(2.7))
sp_diff_hist_j_notpu

sp_diff_hist_j_tpu <- species_summ3 %>% ggplot(aes(x = j_diff)) +
    stat_function(fun = dnorm, args = list(0, sd(species_summ3$j_diff)), color = "red", linetype = "dashed")+
    geom_density(linewidth = 0.8)+
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3)+
    geom_vline(xintercept = mean(species_summ3$j_diff), color = "black", alpha = 0.4)+
    xlab("\U0394J: TPU-enabled") +
    ylab("Density")+
    xlim(-20,50)+
    ylim(0, 0.25)+
    theme_classic()+
    annotate("text", x = 30, y = 0.25, label = paste0("Mean = ", round(mean(species_summ3_notpu$j_diff), digits = 2)), size = rel(2.7))+
    annotate("text", x = 30, y = 0.23, label = paste0("SD = ", round(sd(species_summ3$j_diff), digits = 2)), size = rel(2.7))
sp_diff_hist_j_tpu

gW <- ggplotGrob(sp_diff_hist_vc_notpu)
gX <- ggplotGrob(sp_diff_hist_vc_tpu)
gY <- ggplotGrob(sp_diff_hist_j_notpu)
gZ <- ggplotGrob(sp_diff_hist_j_tpu)

diff_arranged <- grid.arrange(arrangeGrob(cbind(gW, gX), arrangeGrob(cbind(gY, gZ))))

ggsave(plot = diff_arranged, "Figures/diff_density_full_fig.png", width = 6.5, height = 5)





## Density distribution of standard errors
#colors blue = "#31688EFF", yellow = "#FDE725FF", yellow text = "#FFBF00"

#### Vcmax no tpu
hist_vc_se_notpu <- ggplot(mapping = aes(x = V_cmax_se)) +
    #stat_function(data = pho_dat, fun = dnorm, args = list(0, sd(pho_dat$V_cmax_se)), 
    #              color = "black", linetype = "dashed")+
    geom_density(data = pho_dat, linewidth = 0.8, color = "#31688EFF")+
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.8)+
    geom_vline(xintercept = mean(pho_dat$V_cmax_se), color = "#31688EFF", alpha = 0.5)+
    geom_density(data = pho_trad, linewidth = 0.8, color = "#FDE725FF")+
    geom_vline(xintercept = mean(pho_trad$V_cmax_se), color = "#FDE725FF", alpha = 0.5)+
    annotate("text", x = 5.5, y = 1.8, label = "DAT", size = rel(2.7), color = "#31688EFF")+
    annotate("text", x = 9, y = 1.8, label = "Steady-State", size = rel(2.7), color = "#FFBF00")+
    annotate("text", x = 3, y = 1.6, label = "Mean", size = rel(2.7))+
    annotate("text", x = 3, y = 1.4, label = "SD", size = rel(2.7)) +
    annotate("text", x = 3, y = 1.2, label = "Range", size = rel(2.7)) +
    annotate("text", x = 5.5, y = 1.6, label = round(mean(pho_dat$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 5.5, y = 1.4, label = round(sd(pho_dat$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 5.5, y = 1.2, 
             label = paste0(round(min(pho_dat$V_cmax_se), digits = 2), ", ", round(max(pho_dat$V_cmax_se), digits = 2)), 
             size = rel(2.7)) +
    annotate("text", x = 9, y = 1.6, label = round(mean(pho_trad$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 9, y = 1.4, label = round(sd(pho_trad$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 9, y = 1.2, 
             label = paste0(round(min(pho_trad$V_cmax_se), digits = 2), ", ", round(max(pho_trad$V_cmax_se), digits = 2)), 
             size = rel(2.7)) +
    xlab("Vcmax SE: No TPU") +
    ylab("Density")+
    xlim(-3,10)+
    ylim(0, 2)+
    theme_classic()
hist_vc_se_notpu


#### Vcmax with TPU
hist_vc_se_tpu <- ggplot(mapping = aes(x = V_cmax_se)) +
    #stat_function(data = pho_dat_tpu, fun = dnorm, args = list(0, sd(pho_dat_tpu$V_cmax_se)), color = "black", linetype = "dashed")+
    geom_density(data = pho_dat_tpu, linewidth = 0.8, color = "#31688EFF")+
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.8)+
    geom_vline(xintercept = mean(pho_dat_tpu$V_cmax_se), color = "#31688EFF", alpha = 0.5)+
    geom_density(data = pho_trad_tpu, linewidth = 0.8, color = "#FDE725FF")+
    geom_vline(xintercept = mean(pho_trad_tpu$V_cmax_se), color = "#FDE725FF", alpha = 0.5)+
    annotate("text", x = 5.5, y = 1.8, label = "DAT", size = rel(2.7), color = "#31688EFF")+
    annotate("text", x = 9, y = 1.8, label = "Steady-State", size = rel(2.7), color = "#FFBF00")+
    annotate("text", x = 3, y = 1.6, label = "Mean", size = rel(2.7))+
    annotate("text", x = 3, y = 1.4, label = "SD", size = rel(2.7)) +
    annotate("text", x = 3, y = 1.2, label = "Range", size = rel(2.7)) +
    annotate("text", x = 5.5, y = 1.6, label = round(mean(pho_dat_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 5.5, y = 1.4, label = round(sd(pho_dat_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 5.5, y = 1.2, label = paste0(round(min(pho_dat_tpu$V_cmax_se), digits = 2), ", ", round(max(pho_dat_tpu$V_cmax_se), digits = 2)), 
             size = rel(2.7)) +
    annotate("text", x = 9, y = 1.6, label = round(mean(pho_trad_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 9, y = 1.4, label = round(sd(pho_trad_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 9, y = 1.2, label = paste0(round(min(pho_trad_tpu$V_cmax_se), digits = 2), ", ", round(max(pho_trad_tpu$V_cmax_se), digits = 2)), 
             size = rel(2.7)) +
    xlab("Vcmax SE: TPU") +
    ylab("Density")+
    xlim(-3,10)+
    ylim(0, 2)+
    theme_classic()
hist_vc_se_tpu



### jmax no tpu
hist_j_se_notpu <- ggplot(mapping = aes(x = J_se)) +
    #stat_function(data = pho_dat, fun = dnorm, args = list(0, sd(pho_dat$J_se)), color = "black", linetype = "dashed")+
    geom_density(data = pho_dat, linewidth = 0.8, color = "#31688EFF")+
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.8)+
    geom_vline(xintercept = mean(pho_dat$J_se), color = "#31688EFF", alpha = 0.5)+
    geom_vline(xintercept = mean(pho_trad$J_se), color = "#FDE725FF", alpha = 0.5)+
    geom_density(data = pho_trad, linewidth = 0.8, color = "#FDE725FF")+
    annotate("text", x = 1.0, y = 18, label = "DAT", size = rel(2.7), color = "#31688EFF")+
    annotate("text", x = 1.7, y = 18, label = "Steady-State", size = rel(2.7), color = "#FFBF00")+
    annotate("text", x = 0.4, y = 16, label = "Mean", size = rel(2.7))+
    annotate("text", x = 0.4, y = 14, label = "SD", size = rel(2.7)) +
    annotate("text", x = 0.4, y = 12, label = "Range", size = rel(2.7)) +
    annotate("text", x = 1.0, y = 16, label = round(mean(pho_dat$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 1.0, y = 14, label = round(sd(pho_dat$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 1.0, y = 12, label = paste0(round(min(pho_dat$J_se), digits = 2), ", ", round(max(pho_dat$J_se), digits = 2)), 
             size = rel(2.7)) +
    annotate("text", x = 1.7, y = 16, label = round(mean(pho_trad$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 1.7, y = 14, label = round(sd(pho_trad$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 1.7, y = 12, label = paste0(round(min(pho_trad$J_se), digits = 2), ", ", round(max(pho_trad$J_se), digits = 2)), 
             size = rel(2.7)) +
    xlab("Jmax SE: No TPU") +
    ylab("Density")+
    xlim(-1,2)+
    ylim(0, 20)+
    theme_classic()
hist_j_se_notpu


#### jmax with TPU
hist_j_se_tpu <- ggplot(mapping = aes(x = J_se)) +
    #stat_function(data = pho_dat_tpu, fun = dnorm, args = list(0, sd(pho_dat_tpu$J_se)), color = "black", linetype = "dashed")+
    geom_density(data = pho_dat_tpu, linewidth = 0.8, color = "#31688EFF")+
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.8)+
    geom_vline(xintercept = mean(pho_dat_tpu$J_se), color = "#31688EFF", alpha = 0.5)+
    geom_density(data = pho_trad_tpu, linewidth = 0.8, color = "#FDE725FF")+
    geom_vline(xintercept = mean(pho_trad_tpu$J_se), color = "#FDE725FF", alpha = 0.5)+
    annotate("text", x = 1.0, y = 18, label = "DAT", size = rel(2.7), color = "#31688EFF")+
    annotate("text", x = 1.7, y = 18, label = "Steady-State", size = rel(2.7), color = "#FFBF00")+
    annotate("text", x = 0.4, y = 16, label = "Mean", size = rel(2.7))+
    annotate("text", x = 0.4, y = 14, label = "SD", size = rel(2.7)) +
    annotate("text", x = 0.4, y = 12, label = "Range", size = rel(2.7)) +
    annotate("text", x = 1.0, y = 16, label = round(mean(pho_dat_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 1.0, y = 14, label = round(sd(pho_dat_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 1.0, y = 12, label = paste0(round(min(pho_dat_tpu$J_se), digits = 2), ", ", round(max(pho_dat_tpu$J_se), digits = 2)), 
             size = rel(2.7)) +
    annotate("text", x = 1.7, y = 16, label = round(mean(pho_trad_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 1.7, y = 14, label = round(sd(pho_trad_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 1.7, y = 12, label = paste0(round(min(pho_trad_tpu$J_se), digits = 2), ", ", round(max(pho_trad_tpu$J_se), digits = 2)), 
             size = rel(2.7)) +
    xlab("Jmax SE: TPU") +
    ylab("Density")+
    xlim(-1,2)+
    ylim(0, 20)+
    theme_classic()
hist_j_se_tpu

gM <- ggplotGrob(hist_vc_se_notpu)
gN <- ggplotGrob(hist_vc_se_tpu)
gO <- ggplotGrob(hist_j_se_notpu)
gP <- ggplotGrob(hist_j_se_tpu)

se_arranged <- grid.arrange(arrangeGrob(cbind(gM, gN), arrangeGrob(cbind(gO, gP))))

ggsave(plot = se_arranged, "Figures/se_density_full_fig.png", width = 6.5, height = 5)




#No-overshoot data, no TPU
grp_pho_nd_dat <- pho_nd_stat %>%
    filter(curv_meth == "DAT") %>%
    group_by(leaf_unique) %>%
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>% 
    mutate(curv_meth = "DAT",
           fit_type = "no_tpu")

grp_pho_nd_trad <- pho_nd_stat %>%
    filter(curv_meth == "Traditional") %>% 
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
    mutate(curv_meth = "Traditional",
           fit_type = "no_tpu")

grp_pho_nd_all <- rbind(grp_pho_nd_dat, grp_pho_nd_trad)


#No-overshoot data, with TPU
grp_pho_nd_dat_tpu <- pho_nd_stat_tpu %>%
    filter(curv_meth == "DAT") %>%
    group_by(leaf_unique) %>%
    summarise(vcmax = mean(vcmax),
              jmax = mean(jmax)) %>% 
    mutate(curv_meth = "DAT",
           fit_type = "tpu")

grp_pho_nd_trad_tpu <- pho_nd_stat_tpu %>%
    filter(curv_meth == "Traditional") %>% 
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
    mutate(curv_meth = "Traditional",
           fit_type = "tpu")

grp_pho_nd_all_tpu <- rbind(grp_pho_nd_dat_tpu, grp_pho_nd_trad_tpu)


nd_complete <- rbind(grp_pho_nd_all, grp_pho_nd_all_tpu)
nd_complete$curv_meth <- factor(nd_complete$curv_meth)
nd_complete$fit_type <- factor(nd_complete$fit_type)

#Just for the TPU analysis

grp_narm_pho_tpu <- pho_stat_tpu %>%
    na.omit(.$tpu) #This leaves 6 trad and 22 DAT

grp_tpu_6trad <- grp_narm_pho_tpu%>%
    filter(curv_meth == "Traditional") %>% 
    group_by(leaf_unique) %>% 
    summarize(vcmax = vcmax,
              jmax = jmax,
              tpu = tpu,
              curv_meth = "Traditional")

grp_tpu_6dat <- grp_narm_pho_tpu %>%
    filter(curv_meth == "DAT") %>%
    filter(leaf_unique %in% c("K6706L1", "K6707L1", "K6707L2", "K6709L6", "K6714L1", "K6714L2")) %>%
    group_by(leaf_unique) %>% 
    summarize(vcmax = mean(vcmax),
              jmax = mean(jmax),
              tpu = mean(tpu),
              curv_meth = "DAT")

tpu_just6_all <- rbind(grp_tpu_6trad, grp_tpu_6dat)
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
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax[curv_meth == "Traditional"]))

all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax))
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
all_results2 %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "Traditional"]))
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
all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax[curv_meth == "Traditional"]))

all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax))
all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
all_results2 %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "Traditional"]))
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

#Levene's test
nd_complete %>% filter(fit_type == "tpu") %>% leveneTest(vcmax ~ curv_meth, data = .)
nd_complete %>% filter(fit_type == "tpu") %>% leveneTest(jmax ~ curv_meth, data = .)
nd_complete %>% filter(fit_type == "no_tpu") %>% leveneTest(vcmax ~ curv_meth, data = .)
nd_complete %>% filter(fit_type == "no_tpu") %>% leveneTest(jmax ~ curv_meth, data = .)
#Met assumption of homogeneity of variances

#Shapiro Wilk Test
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(vcmax[curv_meth == "Traditional"]))

nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "no_tpu") %>% with(., shapiro.test(jmax[curv_meth == "Traditional"]))

nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(vcmax[curv_meth == "Traditional"]))

nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "DAT"]))
nd_complete %>% filter(fit_type == "tpu") %>% with(., shapiro.test(jmax[curv_meth == "Traditional"]))
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
    mutate(differences = Traditional - DAT)

all_res_wide_vcmax %>% filter(fit_type == "no_tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#That should be fine

all_res_wide_vcmax %>% filter(fit_type == "tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#okay

#Jmax
all_res_wide_jmax <- all_results2 %>%
    dcast(., leaf_unique + fit_type ~ curv_meth, value.var="jmax") %>%
    mutate(differences = Traditional - DAT)

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
    mutate(differences = Traditional - DAT)

nd_comp_wide_vcmax %>% filter(fit_type == "no_tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#That should be fine

nd_comp_wide_vcmax %>% filter(fit_type == "tpu") %>% 
    gghistogram(x = "differences", bins = 10, add_density = TRUE)
#okay

#Jmax
nd_comp_wide_jmax <- nd_complete %>%
    dcast(., leaf_unique + fit_type ~ curv_meth, value.var="jmax") %>%
    mutate(differences = Traditional - DAT)

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

#try without Tachi to see what changes
all_results2 %>%
    group_by(fit_type) %>%
    filter(leaf_unique != "K6707L1" & leaf_unique != "K6707L2") %>% 
    wilcox_test(data =., vcmax ~ curv_meth, paired = TRUE, detailed = TRUE) %>%
    add_significance()

all_results2 %>%
    group_by(fit_type) %>% filter(leaf_unique != "K6707L1" & leaf_unique != "K6707L2") %>% wilcox_effsize(data = ., vcmax ~ curv_meth, paired = TRUE)
#Not significant without Tachi!

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


#Boxplots ------------------------------------------------------------------


## TPU v No TPU
#Full data
vcmax_all_TPUvNoTPU <- ggplot(all_results2, aes(x = fit_type, y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), 
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("No TPU", "TPU")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(V[cmax]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=13, family = "serif", colour = "gray10"))
vcmax_all_TPUvNoTPU
ggsave(plot = vcmax_all_TPUvNoTPU, "Figures/box_vcmax_all_TPUvNoTPU.png")


jmax_all_TPUvNoTPU <- ggplot(all_results2, aes(x = fit_type, y = jmax, 
                                               fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("No TPU", "TPU"))+
    theme_classic()+
    labs(x="Fit Type",
         y = expression(J[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=13, family = "serif", color = "grey10"))
jmax_all_TPUvNoTPU
ggsave(plot = jmax_all_TPUvNoTPU, "Figures/box_jmax_all_TPUvNoTPU.png")


#No Overshoot
vcmax_nOS_TPUvNoTPU <- ggplot(nd_complete, aes(x = fit_type, y = vcmax,
                                               fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("No TPU", "TPU")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(V[cmax]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=13, family = "serif", color = "grey10"))
vcmax_nOS_TPUvNoTPU
ggsave(plot = vcmax_nOS_TPUvNoTPU, "Figures/box_vcmax_nOS_TPUvNoTPU.png")

jmax_nOS_TPUvNoTPU <- ggplot(nd_complete, aes(x = fit_type, y = jmax,
                                              fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("No TPU", "TPU")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(J[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=13, family = "serif", color = "grey10"))
jmax_nOS_TPUvNoTPU
ggsave(plot = jmax_nOS_TPUvNoTPU, "Figures/box_jmax_nOS_TPUvNoTPU.png")




### Test TPU v No TPU plot
test_vcmax_all_TPUvNoTPU <- ggplot(all_results2, aes(x = curv_meth, y = vcmax, fill = fit_type)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Fit Method", labels = c("Without TPU", "With TPU"), 
      values = c("skyblue", "red"))+
    scale_x_discrete(labels = c("DAT", "Steady-State")) +
    theme_classic()+
    labs(x="Curve Method",
         y = expression(V[cmax]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=13, family = "serif", colour = "gray10"))
test_vcmax_all_TPUvNoTPU
ggsave(plot = test_vcmax_all_TPUvNoTPU, "Figures/test_box_vcmax_all_TPUvNoTPU.png")




test_jmax_all_TPUvNoTPU <- ggplot(all_results2, aes(x = curv_meth, y = jmax, fill = fit_type)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Fit Method", labels = c("Without TPU", "With TPU"), 
                      values = c("skyblue", "red"))+
    scale_x_discrete(labels = c("DAT", "Steady-State")) +
    theme_classic()+
    labs(x="Curve Method",
         y = expression(J[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=13, family = "serif", colour = "gray10")) 
test_jmax_all_TPUvNoTPU
ggsave(plot = test_jmax_all_TPUvNoTPU, "Figures/test_box_jmax_all_TPUvNoTPU.png")





## Double Boxplots: All v No OS

pho_nd_stat <- mutate(pho_nd_stat, subset = "nOS")
pho_stat <- mutate(pho_stat, subset = "wOS")
all_and_nOS_noTPU <- rbind(pho_stat, pho_nd_stat)


## No TPU
vcmax_AllvNoOS_noTPU <- all_and_nOS_noTPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), 
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(V[cmax]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=13, family = "serif", colour = "gray10")) +
    scale_y_continuous(limits = c(0, 85))
vcmax_AllvNoOS_noTPU
ggsave(plot = vcmax_AllvNoOS_noTPU, "Figures/box_vcmax_AllvNoOS_noTPU.png")


jmax_AllvNoOS_noTPU <- all_and_nOS_noTPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = jmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("All", "No Overshoot"))+
    theme_classic()+
    labs(x="Fit Type",
         y = expression(J[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=13, family = "serif", color = "grey10")) +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_AllvNoOS_noTPU
ggsave(plot = jmax_AllvNoOS_noTPU, "Figures/box_jmax_AllvNoOS_noTPU.png")





pho_nd_stat_tpu <- mutate(pho_nd_stat_tpu, subset = "nOS")
pho_stat_tpu <- mutate(pho_stat_tpu, subset = "wOS")
all_and_nOS_TPU <- rbind(pho_stat_tpu, pho_nd_stat_tpu)


## With TPU
vcmax_AllvNoOS_TPU <- all_and_nOS_TPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), 
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    theme_classic()+
    labs(x="Fit Type",
         y = expression(V[cmax]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=13, family = "serif", colour = "gray10")) +
    scale_y_continuous(limits = c(0, 85))
vcmax_AllvNoOS_TPU
ggsave(plot = vcmax_AllvNoOS_TPU, "Figures/box_vcmax_AllvNoOS_TPU.png")


jmax_AllvNoOS_TPU <- all_and_nOS_TPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = jmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2)+
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF"))+
    scale_x_discrete(labels = c("All", "No Overshoot"))+
    theme_classic()+
    labs(x="Fit Type",
         y = expression(J[max]*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=13, family = "serif", color = "grey10")) + 
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_AllvNoOS_TPU
ggsave(plot = jmax_AllvNoOS_TPU, "Figures/box_jmax_AllvNoOS_TPU.png")







## Single Plots
lab_DATTrad <- c('DAT', 'Steady-State')

#### All data, No TPU
vcmax_all_noTPU <- ggplot(pho_leaf, aes(x=curv_meth, y=Best_Vcmax_25C)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
vcmax_all_noTPU
ggsave(plot = vcmax_all_noTPU, "Figures/pho_box_DvT_vcmax_noTPU.png")


jmax_all_noTPU <- ggplot(pho_leaf, aes(x=curv_meth, y=Best_Jmax_25C, fill = DAT)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(J[max]*" "*(mu*mol~m^{-2}~s^{-1} *"")))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad) +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_all_noTPU
ggsave(plot = jmax_all_noTPU, "Figures/photo_box_DvT_jmax_noTPU.png")




### All data, with TPU
vcmax_all_TPU <- ggplot(pho_leaf_tpu, aes(x=curv_meth, y=Best_Vcmax_25C)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
vcmax_all_TPU
ggsave(plot = vcmax_all_TPU, "Figures/pho_box_DvT_vcmax_TPU.png")


jmax_all_TPU <- ggplot(pho_leaf_tpu, aes(x=curv_meth, y=Best_Jmax_25C, fill = DAT)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(J[max]*" "*(mu*mol~m^{-2}~s^{-1} *"")))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad) +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_all_TPU
ggsave(plot = jmax_all_TPU, "Figures/photo_box_DvT_jmax_TPU.png")





## No Overshoot, without TPU
nd_vcmax_box_noTPU <- ggplot(grp_pho_nd_all, aes(x=curv_meth, y=vcmax, fill = method)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
nd_vcmax_box_noTPU
ggsave(plot = nd_vcmax_box_noTPU, "Figures/pho_box_noOS_DvT_vcmax_noTPU.png")

nd_jmax_box_noTPU <- ggplot(grp_pho_nd_all, aes(x=curv_meth, y=jmax, fill = method)) +
    geom_boxplot(position = position_dodge(1), outlier.shape = 17, outlier.size = 2, fill = "lightgrey")+
    labs(x="Method", y = expression(J[max]*" "*(mu*mol~m^{-2}~s^{-1} *"")))+
    theme_classic()+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=18, family = "serif"),
          axis.title.y=element_text(size=18, family = "serif"),
          axis.text.x=element_text(size=14, family = "serif", colour = "grey10"),
          axis.text.y=element_text(size = 12, family = "serif", colour = "grey10"),
          legend.position="none")+
    scale_x_discrete(labels=lab_DATTrad)
nd_jmax_box_noTPU
ggsave(plot = nd_jmax_box_noTPU, "Figures/photo_box_nOS_DvT_jmax_noTPU.png")




## 1:1 Plots -----------------------------------------------------------------------------------

## 1:1 plots DAT vs Trad, NO TPU

#Vcmax -- work in progress!!
leaf_sub_vcmax <- select(pho_stat, vcmax, curv_meth, leaf_unique, V_cmax_se)
leaf_wide_vcmax <- reshape(leaf_sub_vcmax, idvar = "leaf_unique", timevar = "curv_meth",
                           direction = "wide")
names(leaf_wide_vcmax)[2:5]=c("vcmax_DAT", "vcmax_DAT_se", "vcmax_Trad", "vcmax_Trad_se")
cor1 <- round(cor(leaf_wide_vcmax$vcmax_DAT, leaf_wide_vcmax$vcmax_Trad), 3)
pho_1to1_vcmax_NoTPU <- ggplot(data = leaf_wide_vcmax, mapping = aes(x = vcmax_Trad,
                                                                 y = vcmax_DAT,
                                                                 color = leaf_unique)) +
    geom_point() +
    geom_errorbar(aes(ymin = vcmax_DAT - vcmax_DAT_se, ymax = vcmax_DAT + vcmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = vcmax_Trad - vcmax_Trad_se, xmax = vcmax_Trad + vcmax_Trad_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6)+
    theme_classic()+
    labs(x = expression("Steady-State "*V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})),
         y = expression("DAT "*V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})),
         col = "Unique Leaf") +
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=12, family = "serif"),
          axis.title.y=element_text(size=12, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=8, family = "serif", color = "gray10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 100)) + 
    scale_y_continuous(limits = c(1, 100)) +
    annotate(geom = "text", label = paste0("r = ", cor1), x = 25, y = 75)
pho_1to1_vcmax_NoTPU
ggsave(plot = pho_1to1_vcmax_NoTPU, "Figures/pho_1to1_vcmax_NoTPU.png")

#Jmax
leaf_sub_jmax <- select(pho_stat, jmax, curv_meth, leaf_unique, J_se)
leaf_wide_jmax <- reshape(leaf_sub_jmax, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_jmax)[2:5]=c("jmax_DAT", "jmax_DAT_se", "jmax_Trad", "jmax_Trad_se")
cor2 <- round(cor(leaf_wide_jmax$jmax_DAT, leaf_wide_jmax$jmax_Trad), 3)
pho_1to1_jmax_noTPU <- ggplot(data = leaf_wide_jmax, mapping = aes(x = jmax_Trad,
                                                               y = jmax_DAT,
                                                               color = leaf_unique))+
    geom_point()+
    geom_errorbar(aes(ymin = jmax_DAT - jmax_DAT_se, ymax = jmax_DAT + jmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = jmax_Trad - jmax_Trad_se, xmax = jmax_Trad + jmax_Trad_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6)+
    theme_classic()+
    labs(x = expression("Steady-State "*J[max]*" "*(mu*mol~m^{-2}~s^{-1})),
         y = expression("DAT "*J[max]*" "*(mu*mol~m^{-2}~s^{-1})),
         col = "Unique Leaf")+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=12, family = "serif"),
          axis.title.y=element_text(size=12, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=8, family = "serif", color = "grey10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 130)) + 
    scale_y_continuous(limits = c(1, 130)) +
    annotate(geom = "text", label = paste0("r = ", cor2), x = 100, y = 30)
pho_1to1_jmax_noTPU
ggsave(plot = pho_1to1_jmax_noTPU, "Figures/pho_1to1_jmax_NoTPU.png")



## 1:1 plots DAT vs Trad WITH TPU

#Vcmax
leaf_sub_vcmax_tpu <- select(pho_stat_tpu, vcmax, curv_meth, leaf_unique, V_cmax_se)
leaf_wide_vcmax_tpu <- reshape(leaf_sub_vcmax_tpu, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_vcmax_tpu)[2:5]=c("vcmax_DAT", "vcmax_DAT_se", "vcmax_Trad", "vcmax_Trad_se")
cor3 <- round(cor(leaf_wide_vcmax_tpu$vcmax_DAT, leaf_wide_vcmax_tpu$vcmax_Trad), 3)
#leaf_wide_vcmax_tpu <- subset(leaf_wide_vcmax_tpu, select = -tree_id.Traditional)
pho_1to1_vcmax_tpu <- ggplot(data = leaf_wide_vcmax_tpu, mapping = aes(x = vcmax_Trad,
                                                               y = vcmax_DAT,
                                                               color = leaf_unique))+
    geom_point()+
    geom_errorbar(aes(ymin = vcmax_DAT - vcmax_DAT_se, ymax = vcmax_DAT + vcmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = vcmax_Trad - vcmax_Trad_se, xmax = vcmax_Trad + vcmax_Trad_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, cex = 0.6)+
    theme_classic()+
    labs(x = expression("Steady-State "*V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})),
         y = expression("DAT "*V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})),
         col = "Unique Leaf") +
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=12, family = "serif"),
          axis.title.y=element_text(size=12, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif", color = "gray10"),
          axis.text.y=element_text(size=8, family = "serif", color = "gray10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 100)) + 
    scale_y_continuous(limits = c(1, 100)) +
    annotate(geom = "text", label = paste0("r = ", cor3), x = 25, y = 75)
pho_1to1_vcmax_tpu
ggsave(plot = pho_1to1_vcmax_tpu, "Figures/pho_1to1_vcmax_tpu.png")


#Jmax
leaf_sub_jmax_tpu <- select(pho_stat_tpu, jmax, curv_meth, leaf_unique, J_se)
leaf_wide_jmax_tpu <- reshape(leaf_sub_jmax_tpu, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_jmax_tpu)[2:5]=c("jmax_DAT", "jmax_DAT_se", "jmax_Trad", "jmax_Trad_se")
cor4 <- round(cor(leaf_wide_jmax_tpu$jmax_DAT, leaf_wide_jmax_tpu$jmax_Trad), 3)
pho_1to1_jmax_tpu <- ggplot(data = leaf_wide_jmax_tpu, mapping = aes(x = jmax_Trad,
                                                                       y = jmax_DAT,
                                                                       color = leaf_unique))+
    geom_point()+
    geom_errorbar(aes(ymin = jmax_DAT - jmax_DAT_se, ymax = jmax_DAT + jmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = jmax_Trad - jmax_Trad_se, xmax = jmax_Trad + jmax_Trad_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6)+
    theme_classic()+
    labs(x = expression("Steady-State "*J[max]*" "*(mu*mol~m^{-2}~s^{-1})),
         y = expression("DAT "*J[max]*" "*(mu*mol~m^{-2}~s^{-1})),
         col = "Unique Leaf")+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=12, family = "serif"),
          axis.title.y=element_text(size=12, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=8, family = "serif", color = "grey10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 130)) + 
    scale_y_continuous(limits = c(1, 130)) +
    annotate(geom = "text", label = paste0("r = ", cor4), x = 100, y = 30)
pho_1to1_jmax_tpu
ggsave(plot = pho_1to1_jmax_tpu, "Figures/pho_1to1_datvtrad_jmax_tpu.png")


#TPU
leaf_sub_tpu <- select(pho_stat_tpu, tpu, curv_meth, leaf_unique, V_TPU_se)
leaf_wide_tpu <- reshape(leaf_sub_tpu, idvar = "leaf_unique", timevar = "curv_meth", direction = "wide")
names(leaf_wide_tpu)[2:5]=c("tpu_DAT", "tpu_DAT_se", "tpu_Trad", "tpu_Trad_se")
cor5 <- round(cor(leaf_wide_tpu$tpu_DAT, leaf_wide_tpu$tpu_Trad), 3)
only_tpu_fit <- filter(leaf_wide_tpu, !is.na(tpu_DAT) & !is.na(tpu_Trad))
only_tpu_fit$tpu_DAT_se <- as.numeric(only_tpu_fit$tpu_DAT_se)
only_tpu_fit$tpu_Trad_se <- as.numeric(only_tpu_fit$tpu_Trad_se)
pho_1to1_tpu_tpu <- ggplot(data = only_tpu_fit, mapping = aes(x = tpu_Trad,
                                                             y = tpu_DAT,
                                                             color = leaf_unique))+
    geom_point(cex = 2.5)+
    geom_errorbar(aes(ymin = tpu_DAT - tpu_DAT_se, ymax = tpu_DAT + tpu_DAT_se)) + 
    geom_errorbarh(aes(xmin = tpu_Trad - tpu_Trad_se, xmax = tpu_Trad + tpu_Trad_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6)+
    theme_classic()+
    labs(x = expression("Steady-State TPU "*(mu*mol~m^{-2}~s^{-1})),
              y = expression("DAT TPU "*(mu*mol~m^{-2}~s^{-1})), col = "Unique Leaf")+
    theme(aspect.ratio = 1,
          axis.title.x=element_text(size=12, family = "serif"),
          axis.title.y=element_text(size=12, family = "serif"),
          axis.text.x=element_text(size=8, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=8, family = "serif", color = "grey10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(0, 12), breaks = c(0,3,6,9,12)) + 
    scale_y_continuous(limits = c(0, 12), breaks = c(0,3,6,9,12)) +
    annotate(geom = "text", label = paste0("r = ", cor5), x = 4, y = 8)
pho_1to1_tpu_tpu
ggsave(plot = pho_1to1_tpu_tpu, "Figures/pho_1to1_datvtrad_tpu.png")



# This script creates all figures and tables for the article at DOI: 10.1093/treephys/tpae153
# This work is licensed under CC by 4.0
# Questions can be directed to Loren Albert (corresponding author) at loren.albert@oregonstate.edu, Emmelia Braun (first-author) at emmelia.braun@oregonstate.edu, or Charles Southwick (first author) at charles.southwick@oregonstate.edu


# Load packages
library(tidyverse)
library(gridExtra)
library(here)
library(grid)
library(ggpubr)
library(ggmagnify) # for figure S1; available at https://github.com/hughjonesd/ggmagnify
library(readxl)

theme_set(theme_classic(base_size = 12, base_family = "serif"))


###

# Load in the data --------------------------------------------------------

### ID codebook
ids <- read.csv(here("3_Clean_data/id_codebook.csv")) 
#rename(treeid = ï..treeid) # this line may not be needed, depending on computer encoding


## For the raw data plots, Fig 1, and Fig S2
cmplt.rm_out <- read.csv(here("3_Clean_data/clean_aci_noOutliers.csv"),
                         header = TRUE,
                         fileEncoding = "latin1") %>% 
    # Remove K6709L3 since it was not a paired DAT-SS on the same leaf
    subset(unique_id != "K6709L3")

cmplt.grp <- group_by(cmplt.rm_out, unique_id)


### For 1:1 plots, summary table, and double boxplots
pho_stat <- read.csv(here("5_Results/pho_stat.csv")) %>% 
    left_join(ids, by = "treeid")
pho_stat_tpu <- read.csv(here("5_Results/pho_stat_tpu.csv")) %>% 
    left_join(ids, by = "treeid")


### For difference figs
diff_tpu_lf <- read.csv(here("5_Results/lf_diffs_summ_TPU.csv"))
diff_notpu_lf <- read.csv(here("5_Results/lf_diffs_summ_noTPU.csv"))


### For double boxplots
pho_nd_stat <- read.csv(here("5_Results/pho_nd_stat.csv")) %>% 
    left_join(ids, by = "treeid")
pho_nd_stat_tpu <- read.csv(here("5_Results/pho_nd_stat_tpu.csv")) %>% 
    left_join(ids, by = "treeid")


### For Figure S5: Differences by species
all_diff_tpu_codes <- read.csv(here("5_Results/tree_diffs_summary_TPU.csv"))
all_diff_notpu_codes <- read.csv(here("5_Results/tree_diffs_summary_noTPU.csv"))


### Boostrapped results
boot_res <- read.csv(here("5_Results/boot_res.csv"))


###


# Plot Curves with Clean Data (not fit; exploratory) ---------------------------------

## Plot all ACi curves on one graph by leaf
ggplot(cmplt.grp, mapping = aes(x = Ci, y = A, color = unique_id)) +
    geom_point(mapping = aes(pch = curv_meth)) 


# Make and save plots for each leaf
for (id in unique(cmplt.grp$unique_id)) {
    df1 <- cmplt.grp %>% 
        filter(unique_id == id)
    gg1 <- ggplot(data = df1, 
                  mapping = aes(x = Ci, y = A, 
                                shape = curv_meth,
                                color = curv_meth, fill = curv_meth)) +
        geom_point(size = 3) +
        theme_classic() +
        labs(y = expression(italic("A"[net])*" "*(mu*mol~m^{-2}~s^{-1})), 
             x = expression(italic("C"[i])*" "*(mu*mol~m^{-2}~s^{-1}))) +
        scale_shape_manual(name = "Method", 
                           labels = c("DAT", "Steady-State"), values = c(21, 24)) +
        scale_fill_viridis_d(begin = 0.3, name = "Method", 
                              labels = c("DAT", "Steady-State")) +
        scale_color_manual(values = c(0, "black")) +
        guides(color = "none") +
        ggtitle(id)
    plot(gg1)
    filename1 <- paste0(id, ".png")
    ggsave(filename = filename1, path = paste0(getwd(), "/6_Figures/"), plot = gg1)
}



###

# Figure 1: Overshoot and no-overshoot comparison ---------------------------------

k6717l1 <- ggplot(data = filter(cmplt.grp, unique_id == "K6717L1"), 
                  mapping = aes(x = Ci, y = A,
                                shape = curv_meth,
                                fill = curv_meth, color = curv_meth)) +
    geom_point(size = 3) +
    labs(y = expression("A"[net]*" "*(mu*mol~m^{-2}~s^{-1})), 
         x = expression("C"[i]*" "*(mu*mol~m^{-2}~s^{-1})),
         tag = "(a)",
         title = expression(italic("Aparisthmium cordatum")*", Leaf 1")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none",
          plot.tag = element_text(size = rel(1)),
          title = element_text(size = 10),
          plot.tag.position = c(0.05, 1)) +
    scale_fill_viridis_d(begin = 0.3, name = "Method", 
                          labels = c("DAT", "SS")) +
    scale_color_manual(values = c(0, "black")) +
    scale_shape_manual(name = "Method", 
                       labels = c("DAT", "SS"), values = c(21, 24))

plot(k6717l1)

k6707l2 <- ggplot(data = filter(cmplt.grp, unique_id == "K6707L2"),
                  mapping = aes(x = Ci, y = A, 
                                shape = curv_meth,
                                fill = curv_meth, color = curv_meth)) +
    geom_point(size = 3) +
    labs(y = expression("A"[net]*" "*(mu*mol~m^{-2}~s^{-1})), 
         x = expression("C"[i]*" "*(mu*mol~m^{-2}~s^{-1})),
         tag = "(b)",
         title = expression(italic("Tachigali chrysophylla")*", Leaf 2")) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.tag = element_text(size = rel(1)),
          title = element_text(size = 10),
          plot.tag.position = c(0.05, 1)) +
    scale_color_manual(values = c(0, "black")) +
    scale_fill_viridis_d(begin = 0.3, name = "Method", 
                          labels = c("DAT", "SS")) +
    scale_shape_manual(name = "Method",
                       labels = c("DAT", "SS"), values = c(21, 24)) +
    guides(color = "none")
plot(k6707l2)


g1 <- ggplotGrob(k6717l1)
g2 <- ggplotGrob(k6707l2)

fig1 <- grid.arrange(arrangeGrob(cbind(g1, g2)))


ggsave(plot = fig1, here("6_Figures/figure1.tif"),
       width = 12.5, height = 4.5, units = "cm", dpi = 600, scaling = 0.75)


###

# Figure 2: Vcmax 1:1 and Density of Differences --------------------------------------


### Vcmax 1:1 WITH TPU
leaf_sub_vcmax_tpu <- select(pho_stat_tpu, vcmax, curv_meth, leaf_unique, V_cmax_se,
                             rel_can_pos)
leaf_wide_vcmax_tpu <- reshape(leaf_sub_vcmax_tpu, 
                               idvar = "leaf_unique", 
                               timevar = "curv_meth",
                               direction = "wide") %>% 
    select(-c(rel_can_pos.DAT)) 
names(leaf_wide_vcmax_tpu)[2:6] = c("vcmax_DAT", "vcmax_DAT_se", 
                                    "vcmax_SS", "vcmax_SS_se",
                                    "rel_can_pos")
cor3 <- round(cor(leaf_wide_vcmax_tpu$vcmax_DAT, leaf_wide_vcmax_tpu$vcmax_SS), 3)


pho_1to1_vcmax_tpu <- ggplot(data = leaf_wide_vcmax_tpu, 
                             mapping = aes(x = vcmax_SS,
                                           y = vcmax_DAT,
                                           color = desc(rel_can_pos))) +
    geom_point() +
    geom_errorbar(aes(ymin = vcmax_DAT - vcmax_DAT_se,
                      ymax = vcmax_DAT + vcmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = vcmax_SS - vcmax_SS_se, 
                       xmax = vcmax_SS + vcmax_SS_se)) +
    geom_abline(intercept = 0, slope = 1, 
                linetype = 5, linewidth = 0.6) +
    labs(x = expression("V"["cmax-SS"]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         y = expression("V"["cmax-DAT"]* " " *(mu*mol~m^{-2}~s^{-1})),
         tag = "(a)") +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.93)),
          legend.position = "none",
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(-0.08, 0.97)) +
    scale_x_continuous(limits = c(1, 100)) + 
    scale_y_continuous(limits = c(1, 100)) +
    scale_color_continuous(type = "viridis") +
    annotate(geom = "text", label = paste0("r = ", cor3),
             x = 25, y = 75, size = rel(3.5)) +
    annotate(geom = "text", label = paste0("AMDE = ", round(boot_res[3,8], 2)),
             x = 25, y = 65, size = rel(3.5))
pho_1to1_vcmax_tpu


### Vcmax 1:1 WITHOUT TPU
leaf_sub_vcmax <- select(pho_stat, vcmax, curv_meth, leaf_unique, V_cmax_se, rel_can_pos)
leaf_wide_vcmax <- reshape(leaf_sub_vcmax, 
                           idvar = "leaf_unique", 
                           timevar = "curv_meth",
                           direction = "wide") %>% 
    select(-c(rel_can_pos.DAT))
names(leaf_wide_vcmax)[2:6] = c("vcmax_DAT", "vcmax_DAT_se", 
                                "vcmax_SS", "vcmax_SS_se", 
                                "rel_can_pos")
leaf_wide_vcmax <- arrange(leaf_wide_vcmax, rel_can_pos)
cor1 <- round(cor(leaf_wide_vcmax$vcmax_DAT, leaf_wide_vcmax$vcmax_SS), 3)


pho_1to1_vcmax_NoTPU <- ggplot(data = leaf_wide_vcmax, 
                               mapping = aes(x = vcmax_SS,
                                             y = vcmax_DAT,
                                             color = desc(rel_can_pos))) +
    geom_point() +
    geom_errorbar(aes(ymin = vcmax_DAT - vcmax_DAT_se, 
                      ymax = vcmax_DAT + vcmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = vcmax_SS - vcmax_SS_se, 
                       xmax = vcmax_SS + vcmax_SS_se)) +
    geom_abline(intercept = 0, slope = 1, 
                linetype = 5, linewidth = 0.6) +
    labs(x = expression("V"["cmax-SS"]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         y = expression("V"["cmax-DAT"]* " " *(mu*mol~m^{-2}~s^{-1})),
        # col = "Unique Leaf",
         tag = "(c)") +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.93)),
          legend.position = "none",
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(-0.08, 0.97)) +
    scale_x_continuous(limits = c(1, 100)) + 
    scale_y_continuous(limits = c(1, 100)) +
    scale_color_continuous(type = "viridis") +
    annotate(geom = "text", label = paste0("r = ", cor1), 
             x = 25, y = 75, size = rel(3.5)) +
    annotate(geom = "text", label = paste0("AMDE = ", round(boot_res[1,8], 2)),
             x = 25, y = 65, size = rel(3.5))
pho_1to1_vcmax_NoTPU


### Vcmax Density of differences WITH TPU
sp_diff_hist_vc_tpu <- diff_tpu_lf %>%
    ggplot(aes(x = vc_diff)) +
    stat_function(fun = dnorm, 
                  args = list(0, sd(diff_tpu_lf$vc_diff)), 
                  color = "red", linetype = "dashed") +
    geom_density(linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = mean(diff_tpu_lf$vc_diff), color = "black", alpha = 0.4) +
    labs(x = expression("SS - DAT "*"V"["cmax"]* " " *(mu*mol~m^{-2}~s^{-1})),
         y = "Density",
         tag = "(b)") +
    xlim(-30, 65) +
    ylim(0, 0.20) +
    theme(axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.93)),
          plot.tag = element_text(size = rel(.9))) +
    annotate("text", x = 37, y = 0.18, 
             label = paste0("Mean = ", round(mean(diff_tpu_lf$vc_diff), digits = 2)), 
             size = rel(3.5)) +
    annotate("text", x = 37, y = 0.15, 
             label = paste0("Median = ", round(median(diff_tpu_lf$vc_diff),
                                               digits = 2)),
             size = rel(3.5)) +
    annotate("text", x = 37, y = 0.12, 
             label = paste0("SD = ", round(sd(diff_tpu_lf$vc_diff), digits = 2)), 
             size = rel(3.5))
sp_diff_hist_vc_tpu


### Vcmax density of differences WITHOUT TPU
sp_diff_hist_vc_notpu <- diff_notpu_lf %>% 
    ggplot(aes(x = vc_diff)) +
    stat_function(fun = dnorm, 
                  args = list(0, sd(diff_notpu_lf$vc_diff)), 
                  color = "red", linetype = "dashed")+
    geom_density(linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = mean(diff_notpu_lf$vc_diff), color = "black", alpha = 0.4) +
    labs(x = expression("SS - DAT "*"V"["cmax"]* " " *(mu*mol~m^{-2}~s^{-1})),
         y = "Density",
         tag = "(d)") +
    xlim(-30, 65) +
    ylim(0, 0.20) +
    theme(axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.93)),
          plot.tag = element_text(size = rel(0.9))) +
    annotate("text", x = 37, y = 0.18, 
             label = paste0("Mean = ", round(mean(diff_notpu_lf$vc_diff), digits = 2)), 
             size = rel(3.5)) +
    annotate("text", x = 37, y = 0.15, 
             label = paste0("Median = ", round(median(diff_notpu_lf$vc_diff), digits = 2)),
             size = rel(3.5)) +
    annotate("text", x = 37, y = 0.12, 
             label = paste0("SD = ", round(sd(diff_notpu_lf$vc_diff), digits = 2)), 
             size = rel(3.5))
sp_diff_hist_vc_notpu


gQ <- ggplotGrob(pho_1to1_vcmax_tpu)
gR <- ggplotGrob(sp_diff_hist_vc_tpu)
gS <- ggplotGrob(pho_1to1_vcmax_NoTPU)
gT <- ggplotGrob(sp_diff_hist_vc_notpu)
TPU <- textGrob("With TPU", rot = 90,
                 gp = gpar(fontfamily = "serif", fontface = "bold", cex = 0.9))
noTPU <- textGrob("Without TPU", rot = 90, 
                  gp = gpar(fontfamily = "serif", fontface = "bold", cex = 0.9))

fig2 <- grid.arrange(arrangeGrob(gQ, left = TPU), gR, arrangeGrob(gS, left = noTPU), gT, nrow = 2) 


ggsave(plot = fig2, here("6_Figures/figure2.tif"),
       width = 17.16, height = 14, units = "cm", dpi = 600,
       scaling = 0.98)

###

# Figure 3: Jmax 1:1 and Density of Differences ----------------------------------------------------------------


### Jmax 1:1 WITHOUT TPU
leaf_sub_jmax <- select(pho_stat, jmax, curv_meth, leaf_unique, J_se, rel_can_pos)
leaf_wide_jmax <- reshape(leaf_sub_jmax,
                          idvar = "leaf_unique",
                          timevar = "curv_meth", 
                          direction = "wide") %>% 
    select(-c(rel_can_pos.DAT))
names(leaf_wide_jmax)[2:6]=c("jmax_DAT", "jmax_DAT_se", "jmax_SS", "jmax_SS_se",
                             "rel_can_pos")
cor2 <- round(cor(leaf_wide_jmax$jmax_DAT, leaf_wide_jmax$jmax_SS), 3)

pho_1to1_jmax_noTPU <- ggplot(data = leaf_wide_jmax, 
                              mapping = aes(x = jmax_SS,
                                            y = jmax_DAT,
                                            color = desc(rel_can_pos))) +
    geom_point() +
    geom_errorbar(aes(ymin = jmax_DAT - jmax_DAT_se, 
                      ymax = jmax_DAT + jmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = jmax_SS - jmax_SS_se, 
                       xmax = jmax_SS + jmax_SS_se)) +
    geom_abline(intercept = 0, slope = 1, 
                linetype = 5, linewidth = 0.6) +
    labs(x = expression("J"["max-SS"]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         y = expression("J"["max-DAT"]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         tag = "(c)") +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size = 8, color = "grey10"),
          axis.text.y = element_text(size = 8, color = "grey10"),
          axis.title = element_text(size = rel(0.93)),
          legend.position = "none",
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(-0.08, 0.97)) +
    scale_x_continuous(limits = c(1, 130)) + 
    scale_y_continuous(limits = c(1, 130)) +
    scale_color_continuous(type = "viridis") +
    annotate(geom = "text", label = paste0("r = ", cor2), 
             x = 35, y = 100, size = rel(3.5)) +
    annotate(geom = "text", label = paste0("AMDE = ", round(boot_res[2,8], 2)),
             x = 35, y = 90, size = rel(3.5))
pho_1to1_jmax_noTPU


### Jmax 1:1 WITH TPU
leaf_sub_jmax_tpu <- select(pho_stat_tpu, jmax, curv_meth, leaf_unique, J_se, rel_can_pos)
leaf_wide_jmax_tpu <- reshape(leaf_sub_jmax_tpu, 
                              idvar = "leaf_unique",
                              timevar = "curv_meth", 
                              direction = "wide") %>% 
    select(-c(rel_can_pos.DAT))
names(leaf_wide_jmax_tpu)[2:6] = c("jmax_DAT", "jmax_DAT_se", "jmax_SS", "jmax_SS_se",
                                   "rel_can_pos")
cor4 <- round(cor(leaf_wide_jmax_tpu$jmax_DAT, leaf_wide_jmax_tpu$jmax_SS), 3)

pho_1to1_jmax_tpu <- ggplot(data = leaf_wide_jmax_tpu,
                            mapping = aes(x = jmax_SS,
                                          y = jmax_DAT,
                                          color = desc(rel_can_pos))) +
    geom_point() +
    geom_errorbar(aes(ymin = jmax_DAT - jmax_DAT_se, ymax = jmax_DAT + jmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = jmax_SS - jmax_SS_se, xmax = jmax_SS + jmax_SS_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6) +
    labs(x = expression("J"["max-SS"]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         y = expression("J"["max-DAT"]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         tag = "(a)") +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size = 8, color = "grey10"),
          axis.text.y = element_text(size = 8, color = "grey10"),
          axis.title = element_text(size = rel(0.93)),
          legend.position = "none",
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(-0.08, 0.97)) +
    scale_x_continuous(limits = c(1, 130)) + 
    scale_y_continuous(limits = c(1, 130)) +
    scale_color_continuous(type = "viridis") +
    annotate(geom = "text", label = paste0("r = ", cor4), 
             x = 35, y = 100, size = rel(3.5)) +
    annotate(geom = "text", label = paste0("AMDE = ", round(boot_res[4,8], 2)),
             x = 35, y = 90, size = rel(3.5))
pho_1to1_jmax_tpu


### Jmax WITHOUT TPU
sp_diff_hist_j_notpu <- diff_notpu_lf %>% 
    ggplot(aes(x = j_diff)) +
    stat_function(fun = dnorm, 
                  args = list(0, sd(diff_notpu_lf$j_diff)), 
                  linetype = "dashed", color = "red") +
    geom_density(linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = mean(diff_notpu_lf$j_diff), color = "black", alpha = 0.4) +
    labs(x = expression("SS - DAT "*"J"["max"]* " " * (mu*mol~m^{-2}~s^{-1})),
         y = "Density",
         tag = "(d)") +
    xlim(-30, 65) +
    ylim(0, 0.20) +
    theme(axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.93)),
          plot.tag = element_text(size = rel(0.9))) +
    annotate("text", x = 37, y = 0.18, 
             label = paste0("Mean = ", round(mean(diff_notpu_lf$j_diff), digits = 2)),
             size = rel(3.5)) +
    annotate("text", x = 37, y = 0.15, 
             label = paste0("Median = ", round(median(diff_notpu_lf$j_diff), 
                                               digits = 2)),
             size = rel(3.5)) +
    annotate("text", x = 37, y = 0.12, 
             label = paste0("SD = ", round(sd(diff_notpu_lf$j_diff), digits = 2)), 
             size = rel(3.5))
sp_diff_hist_j_notpu



### Jmax WITH TPU
sp_diff_hist_j_tpu <- diff_tpu_lf %>% 
    ggplot(aes(x = j_diff)) +
    stat_function(fun = dnorm, 
                  args = list(0, sd(diff_tpu_lf$j_diff)), 
                  color = "red", linetype = "dashed") +
    geom_density(linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = mean(diff_tpu_lf$j_diff), color = "black", alpha = 0.4) +
    labs(x = expression("SS - DAT "*"J"["max"]* " " * (mu*mol~m^{-2}~s^{-1})),
         y = "Density",
         tag = "(b)") +
    xlim(-30, 65) +
    ylim(0, 0.20) +
    theme(axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.93)),
          plot.tag = element_text(size = rel(0.9))) +
    annotate("text", x = 37, y = 0.18, 
             label = paste0("Mean = ", round(mean(diff_tpu_lf$j_diff), digits = 2)),
             size = rel(3.5)) +
    annotate("text", x = 37, y = 0.12, 
             label = paste0("SD = ", round(sd(diff_tpu_lf$j_diff), digits = 2)),
             size = rel(3.5)) +
    annotate("text", x = 37, y = 0.15, 
             label = paste0("Median = ", round(median(diff_tpu_lf$j_diff), digits = 2)),
             size = rel(3.5))
sp_diff_hist_j_tpu



gQ2 <- ggplotGrob(pho_1to1_jmax_tpu)
gR2 <- ggplotGrob(sp_diff_hist_j_tpu)
gS2 <- ggplotGrob(pho_1to1_jmax_noTPU)
gT2 <- ggplotGrob(sp_diff_hist_j_notpu)
TPU <- textGrob("With TPU", rot = 90,
                gp = gpar(fontfamily = "serif", fontface = "bold", cex = 0.9))
noTPU <- textGrob("Without TPU", rot = 90, 
                  gp = gpar(fontfamily = "serif", fontface = "bold", cex = 0.9))


fig3 <- grid.arrange(arrangeGrob(gQ2, left = TPU), 
                     gR2, arrangeGrob(gS2, left = noTPU), gT2, nrow = 2) 

ggsave(plot = fig3, "6_Figures/figure3.tif", 
       width = 17.16, height = 14, units = "cm", dpi = 600,
       scaling = 0.98)

###

# Figure 4: Density distribution of standard errors ---------------------------------

# colors: blue = "#31688EFF", yellow = "#FDE725FF", yellow text = "#FFBF00"

## Read in the data

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


#### Vcmax WITHOUT TPU
hist_vc_se_notpu <- ggplot(mapping = aes(x = V_cmax_se)) +
    # Density and mean of DAT SE
    geom_density(data = pho_dat, linewidth = 0.8, color = "#31688EFF") +
    geom_vline(xintercept = mean(pho_dat$V_cmax_se), color = "#31688EFF", alpha = 0.5) +
    # Density and mean of SS SE
    geom_density(data = pho_SS, linewidth = 0.8, color = "#FDE725FF") +
    geom_vline(xintercept = mean(pho_SS$V_cmax_se), color = "#FDE725FF", alpha = 0.5) +
    # Add in a zero line
    geom_vline(xintercept = 0, 
               color = "black", linetype = "dashed", alpha = 0.8) +
    # Add in the statistics
    annotate("text", x = 1.8, y = 1.8, label = "DAT", 
             size = rel(3), color = "#31688EFF") +
    annotate("text", x = 3.8, y = 1.8, label = "SS",
             size = rel(3), color = "#FFBF00") +
    annotate("text", x = -0.95, y = 1.6, label = "Mean SE", size = rel(3)) +
    annotate("text", x = -0.95, y = 1.4, label = "SD SE", size = rel(3)) +
    annotate("text", x = -0.95, y = 1.2, label = "Range SE", size = rel(3)) +
    annotate("text", x = 1.8, y = 1.6,
             label = round(mean(pho_dat$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 1.8, y = 1.4, 
             label = round(sd(pho_dat$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 1.8, y = 1.2, 
             label = paste0(round(min(pho_dat$V_cmax_se), digits = 2), ", ", 
                            round(max(pho_dat$V_cmax_se), digits = 2)), 
             size = rel(3)) +
    annotate("text", x = 3.8, y = 1.6, 
             label = round(mean(pho_SS$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 3.8, y = 1.4, 
             label = round(sd(pho_SS$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 3.8, y = 1.2, 
             label = paste0(round(min(pho_SS$V_cmax_se), digits = 2), ", ", 
                            round(max(pho_SS$V_cmax_se), digits = 2)), 
             size = rel(3)) +
    labs(x = expression(V[cmax]* " SE " *(mu*mol~m^{-2}~s^{-1})),
         y = "Density",
         title = "Without TPU",
         tag = "(b)") +
    theme(axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.tag = element_text(size = rel(1)),
          plot.tag.position = c(0, 0.95)) +
    xlim(-1.5,5) +
    ylim(0, 2)
hist_vc_se_notpu


#### Vcmax WITH TPU
hist_vc_se_tpu <- ggplot(mapping = aes(x = V_cmax_se)) +
    # Density and mean of DAT SE
    geom_density(data = pho_dat_tpu, linewidth = 0.8, color = "#31688EFF") +
    geom_vline(xintercept = mean(pho_dat_tpu$V_cmax_se),
               color = "#31688EFF", alpha = 0.5) +
    # Density and mean of SS SE
    geom_density(data = pho_SS_tpu, linewidth = 0.8, color = "#FDE725FF") +
    geom_vline(xintercept = mean(pho_SS_tpu$V_cmax_se), 
               color = "#FDE725FF", alpha = 0.5) +
    # Add in a zero line
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.8) +
    # Add in the statistics
    annotate("text", x = 1.8, y = 1.8, 
             label = "DAT", size = rel(3), color = "#31688EFF") +
    annotate("text", x = 3.8, y = 1.8,
             label = "SS", size = rel(3), color = "#FFBF00") +
    annotate("text", x = -0.95, y = 1.6, 
             label = "Mean SE", size = rel(3)) +
    annotate("text", x = -0.95, y = 1.4, 
             label = "SD SE", size = rel(3)) +
    annotate("text", x = -0.95, y = 1.2, 
             label = "Range SE", size = rel(3)) +
    annotate("text", x = 1.8, y = 1.6, 
             label = round(mean(pho_dat_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 1.8, y = 1.4, 
             label = round(sd(pho_dat_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(3))+
    annotate("text", x = 1.8, y = 1.2, 
             label = paste0(round(min(pho_dat_tpu$V_cmax_se), digits = 2), ", ", 
                            round(max(pho_dat_tpu$V_cmax_se), digits = 2)), 
             size = rel(3)) +
    annotate("text", x = 3.8, y = 1.6, 
             label = round(mean(pho_SS_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 3.8, y = 1.4, 
             label = round(sd(pho_SS_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 3.8, y = 1.2, 
             label = paste0(round(min(pho_SS_tpu$V_cmax_se), digits = 2), ", ", 
                            round(max(pho_SS_tpu$V_cmax_se), digits = 2)), 
             size = rel(3)) +
    labs(x = expression(V[cmax]* " SE " *(mu*mol~m^{-2}~s^{-1})),
         y = "Density",
         title = "With TPU",
         tag = "(a)") +
    theme(axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          plot.tag = element_text(size = rel(1)),
          plot.tag.position = c(0.05, 0.95)) +
    xlim(-1.5,5) +
    ylim(0, 2)
hist_vc_se_tpu



### Jmax WITHOUT tpu
hist_j_se_notpu <- ggplot(mapping = aes(x = J_se)) +
    # Density and mean of DAT SE
    geom_density(data = pho_dat, linewidth = 0.8, color = "#31688EFF") +
    geom_vline(xintercept = mean(pho_dat$J_se), color = "#31688EFF", alpha = 0.5) +
    # Density and mean of SS SE
    geom_vline(xintercept = mean(pho_SS$J_se), color = "#FDE725FF", alpha = 0.5) +
    geom_density(data = pho_SS, linewidth = 0.8, color = "#FDE725FF") +
    # Add in a zero line
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.8) +
    # Add in the statistics
    annotate("text", x = 0.6, y = 18, label = "DAT", 
             size = rel(3), color = "#31688EFF") + 
    annotate("text", x = 1.5, y = 18, label = "SS", size = rel(3), color = "#FFBF00") +
    annotate("text", x = -0.4, y = 16, label = "Mean SE", size = rel(3)) +
    annotate("text", x = -0.4, y = 14, label = "SD SE", size = rel(3)) +
    annotate("text", x = -0.4, y = 12, label = "Range SE", size = rel(3)) +
    annotate("text", x = 0.6, y = 16, 
             label = round(mean(pho_dat$J_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 0.6, y = 14, 
             label = round(sd(pho_dat$J_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 0.6, y = 12, 
             label = paste0(round(min(pho_dat$J_se), digits = 2), ", ", 
                            round(max(pho_dat$J_se), digits = 2)), 
             size = rel(3)) +
    annotate("text", x = 1.5, y = 16, 
             label = round(mean(pho_SS$J_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 1.5, y = 14, 
             label = round(sd(pho_SS$J_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 1.5, y = 12, 
             label = paste0(round(min(pho_SS$J_se), digits = 2), ", ", 
                            round(max(pho_SS$J_se), digits = 2)), 
             size = rel(3)) +
    labs(x = expression(J[max]* " SE " *(mu*mol~m^{-2}~s^{-1})), 
         y = "Density",
         tag = "(d)") +
    theme(axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          plot.tag = element_text(size = rel(1))) +
    xlim(-.75,2) +
    ylim(0, 20)
hist_j_se_notpu


#### Jmax WITH TPU
hist_j_se_tpu <- ggplot(mapping = aes(x = J_se)) +
    # Density and mean of DAT SE
    geom_density(data = pho_dat_tpu, linewidth = 0.8, color = "#31688EFF") +
    geom_vline(xintercept = mean(pho_dat_tpu$J_se), color = "#31688EFF", alpha = 0.5) +
    # Density and mean of SS SE
    geom_density(data = pho_SS_tpu, linewidth = 0.8, color = "#FDE725FF") +
    geom_vline(xintercept = mean(pho_SS_tpu$J_se), color = "#FDE725FF", alpha = 0.5) +
    # Add in a zero line
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.8) +
    # Add in the statistics
    annotate("text", x = 0.6, y = 18, label = "DAT", 
             size = rel(3), color = "#31688EFF") +
    annotate("text", x = 1.5, y = 18, label = "SS", size = rel(3), color = "#FFBF00") +
    annotate("text", x = -0.4, y = 16, label = "Mean SE", size = rel(3)) +
    annotate("text", x = -0.4, y = 14, label = "SD SE", size = rel(3)) +
    annotate("text", x = -0.4, y = 12, label = "Range SE", size = rel(3)) +
    annotate("text", x = 0.6, y = 16, 
             label = round(mean(pho_dat_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 0.6, y = 14, 
             label = round(sd(pho_dat_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 0.6, y = 12, 
             label = paste0(round(min(pho_dat_tpu$J_se), digits = 2), ", ", 
                            round(max(pho_dat_tpu$J_se), digits = 2)), 
             size = rel(3)) +
    annotate("text", x = 1.5, y = 16, 
             label = round(mean(pho_SS_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 1.5, y = 14, 
             label = round(sd(pho_SS_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(3)) +
    annotate("text", x = 1.5, y = 12, 
             label = paste0(round(min(pho_SS_tpu$J_se), digits = 2), ", ", 
                            round(max(pho_SS_tpu$J_se), digits = 2)), 
             size = rel(3)) +
    labs(x = expression(J[max]* " SE " *(mu*mol~m^{-2}~s^{-1})), 
         y = "Density",
         tag = "(c)") +
    theme(axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          plot.tag = element_text(size = rel(1)),
          plot.tag.position = c(0.05, 1)) +
    xlim(-0.75,2) +
    ylim(0, 20)
hist_j_se_tpu


gM <- ggplotGrob(hist_vc_se_tpu)
gN <- ggplotGrob(hist_vc_se_notpu + rremove("ylab"))
gO <- ggplotGrob(hist_j_se_tpu)
gP <- ggplotGrob(hist_j_se_notpu + rremove("ylab"))

se_arranged <- grid.arrange(arrangeGrob(cbind(gM, gN)), arrangeGrob(cbind(gO, gP)),
                            heights = c(4.25,4))

ggsave(plot = se_arranged, here("6_Figures/figure4.tif"),
       width = 17.16, height = 14, units = "cm", dpi = 600,
       scaling = 0.9)

###

# Table S1: Summary stats for SS, TPU-enabled parameters------------------------------------------

pho_stat_tpu_ss <- rename(pho_stat_tpu,
                          vcmax_se = V_cmax_se,
                          jmax_se = J_se) %>% 
    filter(curv_meth == "SS")


ss_keyparams_tpu <- pho_stat_tpu_ss %>%
    group_by(treeid) %>% 
    summarize(mean_vcmax = mean(vcmax),
              sd_vcmax = sd(vcmax),
              mean_jmax = mean(jmax),
              sd_jmax = sd(jmax))

# join with codes
ss_keyparams_tpu_codes <- left_join(ss_keyparams_tpu, ids, by = "treeid")

# Sort the dataframe by rel_can_pos in ascending order
sorted_df <- ss_keyparams_tpu_codes[order(ss_keyparams_tpu_codes$rel_can_pos), ]

###

# Figure S1: DAT split curve ---------------------------------------------------------

## Note have to copy and paste excel values as values, not formulas, for them to be read into R

dat_test <- read_excel(path = here("1_Raw_data/2022-07-07-1518_DAT_Insitu.xlsx"), skip = 14, sheet = "Measurements")
dat_test <- dat_test[-1,]

dat_test$A <- as.numeric(dat_test$A)
dat_test$Ci <- as.numeric(dat_test$Ci)
dat_test$date <- as_datetime(dat_test$date)


# Make and save plots for each leaf
split_plot <- ggplot(data = dat_test, 
                     mapping = aes(x = Ci, y = A, color = date)) +
    geom_point(size = 3) +
    theme(axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.93))) +
    labs(y = expression("A"[net]*" "*(mu*mol~m^{-2}~s^{-1})), 
         x = expression("C"[i]*" "*(mu*mol~m^{-2}~s^{-1}))) +
    guides(color = "none")


ggsave(split_plot, file = here("6_Figures/figureS1.tif"),
       width = 8.23, height = 6, units = "cm", scale = 2, dpi = 600)

###
# Figure S2: Backwards points ---------------------------------------------

k6715l1.1 <- ggplot(data = filter(cmplt.grp, unique_id == "K6715L1-1"), 
                  mapping = aes(x = Ci, y = A)) +
    geom_point(size = 3) +
    labs(y = expression("A"[net]*" "*(mu*mol~m^{-2}~s^{-1})), 
         x = expression("C"[i]*" "*(mu*mol~m^{-2}~s^{-1}))) +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.93))) +
    geom_magnify(from = c(-5, 100, -1, 3), to = c(500, 1000, 0, 7), axes = "xy",
                 shadow = TRUE)
k6715l1.1

ggsave(plot = k6715l1.1, file = here("6_figures/figureS2.tif"), dpi = 600,
       height = 7, width = 8.23, units = "cm", scale = 1.3)

###



# Figure S3: TPU 1:1 -------------------------------------------------------

# TPU 1:1
leaf_sub_tpu <- select(pho_stat_tpu, tpu, curv_meth, leaf_unique, V_TPU_se, rel_can_pos)
leaf_wide_tpu <- reshape(leaf_sub_tpu,
                         idvar = "leaf_unique",
                         timevar = "curv_meth", 
                         direction = "wide") %>% 
    na.omit() %>% 
    select(-c("rel_can_pos.DAT"))
names(leaf_wide_tpu)[2:6] = c("tpu_DAT", "tpu_DAT_se", "tpu_SS", "tpu_SS_se", "rel_can_pos")
cor5 <- round(cor(leaf_wide_tpu$tpu_DAT, leaf_wide_tpu$tpu_SS), 3)
only_tpu_fit <- filter(leaf_wide_tpu, !is.na(tpu_DAT) & !is.na(tpu_SS))
only_tpu_fit$tpu_DAT_se <- as.numeric(only_tpu_fit$tpu_DAT_se)
only_tpu_fit$tpu_SS_se <- as.numeric(only_tpu_fit$tpu_SS_se)
only_tpu_fit <- arrange(only_tpu_fit, rel_can_pos)


pho_1to1_tpu_tpu <- ggplot(data = only_tpu_fit,
                           mapping = aes(x = tpu_SS,
                                         y = tpu_DAT,
                                         color = desc(rel_can_pos))) +
    geom_point(cex = 2.5) +
    geom_errorbar(aes(ymin = tpu_DAT - tpu_DAT_se, ymax = tpu_DAT + tpu_DAT_se)) + 
    geom_errorbarh(aes(xmin = tpu_SS - tpu_SS_se, xmax = tpu_SS + tpu_SS_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6) +
    labs(x = expression("TPU"[SS] * " " * (mu*mol~m^{-2}~s^{-1})),
         y = expression("TPU"[DAT] * " " * (mu*mol~m^{-2}~s^{-1})),
         col = "Unique Leaf") +
    theme(aspect.ratio = 1,
          axis.text.x = element_text(size = 8, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.93)),
          legend.position = "none") +
    scale_x_continuous(limits = c(0, 12), breaks = c(0,3,6,9,12)) + 
    scale_y_continuous(limits = c(0, 12), breaks = c(0,3,6,9,12)) +
    scale_color_continuous(type = "viridis") +
    #geom_text(aes(label = rel_can_pos), hjust = 2) +
    annotate(geom = "text", label = paste0("r = ", cor5), x = 3, y = 9, size = rel(3.5)) +
    annotate(geom = "text", label = paste0("AMDE = ", round(boot_res[5,8], 2)),
             x = 3, y = 8, size = rel(3.5))
pho_1to1_tpu_tpu

ggsave(plot = pho_1to1_tpu_tpu, here("6_Figures/figureS3.tif"), dpi = 600,
       width = 8.23, height = 8.23, units = "cm", scaling = 0.99)


###


# Figure S4: Double Boxplots: All v No OS --------------------------------------------

#### TPU-omitted Boxplots
pho_nd_stat <- mutate(pho_nd_stat, subset = "nOS")
pho_stat <- mutate(pho_stat, subset = "wOS")
all_and_nOS_noTPU <- rbind(pho_stat, pho_nd_stat)


## Vcmax WITHOUT TPU
vcmax_AllvNoOS_noTPU <- all_and_nOS_noTPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 1.6, linewidth = 0.3) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "SS"), 
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    labs(x = "Dataset",
         y = expression("V"[cmax]*" "*(mu*mol~m^{-2}~s^{-1})),
         tag = "(b)",
         title = "Without TPU") +
    annotate("text", label = "**", x = 1, y = 120, color = "red") +
    annotate("text", label = "*", x = 2, y = 120, color = "red") +
    theme(axis.text.x = element_text(size = 9, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.85)),
          plot.title = element_text(size = rel(0.93),face = "bold", hjust = 0.5),
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(0, 0.95),
          legend.title = element_text(size = rel(0.9)),
          legend.text = element_text(size = rel(0.88))) +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
vcmax_AllvNoOS_noTPU

### Jmax WITHOUT TPU
jmax_AllvNoOS_noTPU <- all_and_nOS_noTPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = jmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 1.6, linewidth = 0.3) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "SS"),
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    labs(x = "Dataset",
         y = expression(J[max]*" "*(mu*mol~m^{-2}~s^{-1})),
         tag = "(d)") +
    annotate("text", label = "****", x = 1, y = 120, color = "red") +
    annotate("text", label = "***", x = 2, y = 120, color = "red") +
    theme(axis.text.x = element_text(size = 9, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.85)),
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(0, 1.05),
          legend.title = element_text(size = rel(0.9)),
          legend.text = element_text(size = rel(0.88))) +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_AllvNoOS_noTPU




#### TPU-enabled boxplots

pho_nd_stat_tpu <- mutate(pho_nd_stat_tpu, subset = "nOS")
pho_stat_tpu <- mutate(pho_stat_tpu, subset = "wOS")
all_and_nOS_TPU <- rbind(pho_stat_tpu, pho_nd_stat_tpu)


### Vcmax WITH TPU
vcmax_AllvNoOS_TPU <- all_and_nOS_TPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 1.6, linewidth = 0.3) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "SS"), 
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    labs(x = "Dataset",
         y = expression(V[cmax]*" "*(mu*mol~m^{-2}~s^{-1})),
         title = "With TPU",
         tag = "(a)") +
    annotate("text", label = "*", x = 2, y = 120, color = "red") +
    theme(axis.text.x = element_text(size = 9, color = "gray10"),
          axis.text.y = element_text(size = 8, colour = "gray10"),
          axis.title = element_text(size = rel(0.85)),
          plot.title = element_text(size = rel(0.93), face = "bold", hjust = 0.5),
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(0.15, 0.95),
          legend.title = element_text(size = rel(0.9)),
          legend.text = element_text(size = rel(0.88))) +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
vcmax_AllvNoOS_TPU


### Jmax WITH TPU
jmax_AllvNoOS_TPU <- all_and_nOS_TPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = jmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 1.6, linewidth = 0.3) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "SS"),
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    labs(x = "Dataset",
         y = expression(J[max]*" "*(mu*mol~m^{-2}~s^{-1})),
         tag = "(c)") +
    annotate("text", label = "****", x = 1, y = 120, color = "red") +
    annotate("text", label = "***", x = 2, y = 120, color = "red") +
    theme(axis.text.x = element_text(size = 9, color = "gray10"),
          axis.text.y = element_text(size = 8, color = "gray10"),
          axis.title = element_text(size = rel(0.85)),
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(0.15, 1.05),
          legend.title = element_text(size = rel(0.9)),
          legend.text = element_text(size = rel(0.88))) + 
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_AllvNoOS_TPU


gW <- ggplotGrob(vcmax_AllvNoOS_TPU + rremove("legend") + rremove("xlab"))
gX <- ggplotGrob(vcmax_AllvNoOS_noTPU + rremove("ylab")  + rremove("legend") + rremove("xlab")) 
gY <- ggplotGrob(jmax_AllvNoOS_TPU + rremove("legend"))
gZ <- ggplotGrob(jmax_AllvNoOS_noTPU + rremove("ylab") + rremove("legend") )

figS3_grid <- grid.arrange(arrangeGrob(cbind(gW, gX)), arrangeGrob(cbind(gY, gZ)),
                           heights = c(1.05, 0.95))

legend <- cowplot::get_legend(vcmax_AllvNoOS_TPU)
figS3 <- grid.arrange(
    figS3_grid,
    legend,
   widths = c(4, 1)
)



ggsave(plot = figS3, here("6_Figures/figureS4.tif"),
       width = 12.5, height = 10, units = "cm", dpi = 600, scale = 0.9)




###

# Figure S5: Species difference histograms, sorted by relative canopy height ----------

## Note that error bars represent the mean of the absolute error for the differences

## Vcmax WITH TPU
vc_diff_hist <- ggplot(data = all_diff_tpu_codes, 
                       aes(x = reorder(gen_spec_id, desc(rel_can_pos.x)),
                           y = vc_diff,
                           fill = desc(rel_can_pos.x))) +
    geom_bar(stat = "identity", color = "grey20") +
    labs(x = NULL,
         y = expression("SS - DAT "*V[cmax]* " " *(mu*mol~m^{-2}~s^{-1})),
         title = "With TPU",
         tag = "(a)") +
    geom_errorbar(aes(x = gen_spec_id, 
                      ymin = vc_diff - vc_diff_se,
                      ymax = vc_diff + vc_diff_se),
                  width = 0.3, colour = "#CA0068", alpha = 0.9, linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    ylim(-5, 45) +
    theme(axis.text.y = element_text(size = 10, face = "italic"),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = rel(0.9)),
          plot.title = element_text(size = rel(0.93), face = "bold", hjust = 0.5),
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(0.4, 0.9),
          legend.position = "none") +
    scale_fill_continuous(type = "viridis") +
    coord_flip()
vc_diff_hist


# Jmax WITH TPU
j_diff_hist <- ggplot(data = all_diff_tpu_codes, 
                      aes(x = reorder(gen_spec_id, desc(rel_can_pos.x)),
                          y = j_diff,
                          fill = desc(rel_can_pos.x))) +
    geom_bar(stat = "identity", color = "grey20") +
    labs(x = NULL,
         y = expression("SS - DAT "*J[max]*" "*(mu*mol~m^{-2}~s^{-1})),
         tag = "(c)") +
    geom_errorbar(aes(x = gen_spec_id, 
                      ymin = j_diff - j_diff_se,
                      ymax = j_diff + j_diff_se),
                  width = 0.3, colour = "#CA0068", alpha = 0.9, linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    ylim(-5, 45) +
    theme(axis.text.y = element_text(size = 10, face = "italic"),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = rel(0.9)),
          plot.title = element_text(size = rel(0.93), face = "bold", hjust = 0.5),
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(0.4, 0.95),
          legend.position = "none") +
    scale_fill_continuous(type = "viridis") +
    coord_flip()
j_diff_hist



## Vcmax WITHOUT TPU
vc_diff_hist_notpu <- ggplot(data = all_diff_notpu_codes, 
                             aes(x = reorder(gen_spec_id, desc(rel_can_pos.x)),
                                 y = vc_diff,
                                 fill = desc(rel_can_pos.x))) +
    geom_bar(stat = "identity", color = "grey20") +
    labs(x = NULL,
         y = expression("SS - DAT "*V[cmax]* " " *(mu*mol~m^{-2}~s^{-1})),
         title = "Without TPU",
         tag = "(b)") +
    geom_errorbar(aes(x = gen_spec_id, 
                      ymin = vc_diff - vc_diff_se, 
                      ymax = vc_diff + vc_diff_se),
                  width = 0.3, colour = "#CA0068", alpha = 0.9, linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    ylim(-5, 45) +
    theme(axis.text.y = element_text(size = 10, face = "italic"),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = rel(0.9)),
          plot.title = element_text(size = rel(0.93), face = "bold", hjust = 0.5),
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(0.08, 0.9),
          legend.position = "none") +
    scale_fill_continuous(type = "viridis") +
    coord_flip()
vc_diff_hist_notpu

## Jmax WITHOUT TPU
j_diff_hist_notpu <- ggplot(data = all_diff_notpu_codes,
                            aes(x = reorder(gen_spec_id, desc(rel_can_pos.x)),
                                y = j_diff,
                                fill = desc(rel_can_pos.x))) +
    geom_bar(stat = "identity", color = "grey20") +
    labs(x = NULL,
         y = expression("SS - DAT "*J[max]* " " *(mu*mol~m^{-2}~s^{-1})),
         tag = "(d)") +
    geom_errorbar(aes(x = gen_spec_id, 
                      ymin = j_diff - j_diff_se,
                      ymax = j_diff + j_diff_se), 
                  width = 0.3, colour = "#CA0068", alpha = 0.9, linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    ylim(-5, 45) +
    theme(axis.text.y = element_text(size = 10, face = "italic"),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = rel(0.9)),
          plot.title = element_text(size = rel(0.93), face = "bold", hjust = 0.5),
          plot.tag = element_text(size = rel(0.9)),
          plot.tag.position = c(0.08, 0.95),
          legend.position = "none") +
    scale_fill_continuous(type = "viridis") +
    coord_flip()
j_diff_hist_notpu



gA <- ggplotGrob(vc_diff_hist)
gB <- ggplotGrob(vc_diff_hist_notpu + rremove("y.text"))
gC <- ggplotGrob(j_diff_hist)
gD <- ggplotGrob(j_diff_hist_notpu + rremove("y.text"))

figS4 <- grid.arrange(arrangeGrob(cbind(gA, gB)), arrangeGrob(cbind(gC, gD)), 
                      heights = c(4.5,4))

ggsave(plot = figS4, here("6_Figures/figureS5.tif"), 
       width = 17.16, height = 12, units = "cm", 
       dpi = 600, scale = 1.1)


###


# TPU v. No TPU Boxplots --------------------------------------------------------------

### Not in manuscript - just visualization
# 
# # Full data
# 
# ### Vcmax ALL data
# vcmax_all_TPUvNoTPU <- ggplot(all_avg_lf_res, 
#                               aes(x = factor(fit_type, levels = c("tpu", "no_tpu")),
#                                   y = vcmax, fill = curv_meth)) +
#     geom_boxplot(outlier.shape = 17, outlier.size = 2) +
#     scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), 
#                       values = c("#31688EFF", "#FDE725FF")) +
#     scale_x_discrete(labels = c("TPU-enabled", "TPU-omitted")) +
#     theme_classic() +
#     labs(x = "Fit Type",
#          y = expression(italic(V[cmax])* " " *(mu*mol~m^{-2}~s^{-1}))) +
#     theme(axis.title.x = element_text(size = 16, family = "serif"),
#           axis.title.y = element_text(size = 16, family = "serif"),
#           axis.text.x = element_text(size = 13, family = "serif", color = "gray10"),
#           axis.text.y = element_text(size = 13, family = "serif", colour = "gray10"))
# vcmax_all_TPUvNoTPU
# # ggsave(plot = vcmax_all_TPUvNoTPU, here("6_Figures/box_vcmax_all_TPUvNoTPU.png"))
# 
# ### Jmax ALL data
# jmax_all_TPUvNoTPU <- ggplot(all_avg_lf_res, 
#                              aes(x = factor(fit_type, levels = c("tpu", "no_tpu")),
#                                  y = jmax, fill = curv_meth)) +
#     geom_boxplot(outlier.shape = 17, outlier.size = 2) +
#     scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
#                       values = c("#31688EFF", "#FDE725FF")) +
#     scale_x_discrete(labels = c("TPU-enabled", "TPU-omitted")) +
#     theme_classic() +
#     labs(x="Fit Type",
#          y = expression(italic(J[max])* " " *(mu*mol~m^{-2}~s^{-1}))) +
#     theme(axis.title.x=element_text(size=16, family = "serif"),
#           axis.title.y=element_text(size=16, family = "serif"),
#           axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
#           axis.text.y=element_text(size=13, family = "serif", color = "grey10"))
# jmax_all_TPUvNoTPU
# # ggsave(plot = jmax_all_TPUvNoTPU, here("6_Figures/box_jmax_all_TPUvNoTPU.png"))
# 
# 
# 
# # No Overshoot curves
# 
# ### Vcmax no overhoot data
# vcmax_nOS_TPUvNoTPU <- ggplot(nd_complete, 
#                               aes(x = factor(fit_type, levels = c("tpu", "no_tpu")), 
#                                   y = vcmax, fill = curv_meth)) +
#     geom_boxplot(outlier.shape = 17, outlier.size = 2) +
#     scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
#                       values = c("#31688EFF", "#FDE725FF")) +
#     scale_x_discrete(labels = c("TPU-enabled", "TPU-omitted")) +
#     theme_classic() +
#     labs(x = "Fit Type",
#          y = expression(italic(V[cmax])*" "*(mu*mol~m^{-2}~s^{-1}))) +
#     theme(axis.title.x = element_text(size = 16, family = "serif"),
#           axis.title.y = element_text(size = 16, family = "serif"),
#           axis.text.x = element_text(size = 13, family = "serif", color = "grey10"),
#           axis.text.y = element_text(size = 13, family = "serif", color = "grey10"))
# vcmax_nOS_TPUvNoTPU
# # ggsave(plot = vcmax_nOS_TPUvNoTPU, here("6_Figures/box_vcmax_nOS_TPUvNoTPU.png"))
# 
# 
# ### Jmax no overshoot data
# jmax_nOS_TPUvNoTPU <- ggplot(nd_complete, 
#                              aes(x = factor(fit_type, levels = c("tpu", "no_tpu")),
#                                  y = jmax, fill = curv_meth)) +
#     geom_boxplot(outlier.shape = 17, outlier.size = 2) +
#     scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
#                       values = c("#31688EFF", "#FDE725FF")) +
#     scale_x_discrete(labels = c("TPU-enabled", "TPU-omitted")) +
#     theme_classic() +
#     labs(x = "Fit Type",
#          y = expression(italic(J[max])*" "*(mu*mol~m^{-2}~s^{-1}))) +
#     theme(axis.title.x = element_text(size = 16, family = "serif"),
#           axis.title.y = element_text(size = 16, family = "serif"),
#           axis.text.x = element_text(size = 13, family = "serif", color = "grey10"),
#           axis.text.y = element_text(size = 13, family = "serif", color = "grey10"))
# jmax_nOS_TPUvNoTPU
# # ggsave(plot = jmax_nOS_TPUvNoTPU, here("6_Figures/box_jmax_nOS_TPUvNoTPU.png"))
# 




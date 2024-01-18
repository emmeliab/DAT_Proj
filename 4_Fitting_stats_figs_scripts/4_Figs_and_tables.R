### All the in-text and supplemental figures and tables

# Load packages
library(tidyverse)
library(gridExtra)

# Load in the data and filter
cmplt.rm_out <- read.csv(here("3_Clean_data/clean_aci_noOutliers.csv"),
                         header = TRUE)

cmplt.grp <- group_by(cmplt.rm_out, unique_id)




# Plot Curves with Clean Data (not fit) -----------------------------------------------------

## Plot all ACi curves on one graph by leaf
ggplot(cmplt.grp, mapping = aes(x = Ci, y = A, color = unique_id)) +
    geom_point(mapping = aes(pch = curv_meth)) +
    theme_classic()


# Make and save plots for each leaf
for (id in unique(cmplt.grp$unique_id)) {
    df1 <- cmplt.grp %>% 
        filter(unique_id == id)
    gg1 <- ggplot(data = df1, 
                  mapping = aes(x = Ci, y = A, 
                                shape = curv_meth, color = curv_meth)) +
        geom_point(size = 3) +
        theme_classic() +
        labs(y = expression(italic("A"[net])*" "*(mu*mol~m^{-2}~s^{-1})), 
             x = expression(italic("C"[i])*" "*(mu*mol~m^{-2}~s^{-1}))) +
        scale_shape_manual(name = "Method", 
                           labels = c("DAT", "Steady-State"), values = c(19, 17)) +
        scale_color_viridis_d(begin = 0.3, name = "Method", 
                              labels = c("DAT", "Steady-State")) +
        ggtitle(id)
    plot(gg1)
    filename1 <- paste0(id, ".png")
    ggsave(filename = filename1, path = paste0(getwd(), "/6_Figures/"), plot = gg1)
}





# Figure 1: Overshoot and no-overshoot comparison ---------------------------------


k6717l1 <- ggplot(data = filter(cmplt.grp, unique_id == "K6717L1"), 
                  mapping = aes(x = Ci, y = A,
                                shape = curv_meth, color = curv_meth)) +
    geom_point(size = 3) +
    theme_classic() +
    labs(y = expression(italic("A"[net])*" "*(mu*mol~m^{-2}~s^{-1})), 
         x = expression(italic("C"[i])*" "*(mu*mol~m^{-2}~s^{-1})),
         tag = "a",
         title = expression(italic("Aparisthmium cordatum")*", Leaf 1")) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    scale_shape_manual(name = "Method", 
                       labels = c("DAT", "Steady-State"), values = c(19, 17)) +
    scale_color_viridis_d(begin = 0.3, name = "Method", 
                          labels = c("DAT", "Steady-State"))
plot(k6717l1)

k6707l2 <- ggplot(data = filter(cmplt.grp, unique_id == "K6707L2"), mapping = aes(x = Ci, y = A, shape = curv_meth, color = curv_meth)) +
    geom_point(size = 3) +
    theme_classic() +
    labs(y = expression(italic("A"[net])*" "*(mu*mol~m^{-2}~s^{-1})), 
         x = expression(italic("C"[i])*" "*(mu*mol~m^{-2}~s^{-1})),
         tag = "b",
         title = expression(italic("Tachigali chrysophylla")*", Leaf 2")) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_shape_manual(name = "Method",
                       labels = c("DAT", "Steady-State"), values = c(19, 17)) +
    scale_color_viridis_d(begin = 0.3, name = "Method", 
                          labels = c("DAT", "Steady-State"))
plot(k6707l2)

g1 <- ggplotGrob(k6717l1)
g2 <- ggplotGrob(k6707l2)


fig1 <- grid.arrange(arrangeGrob(cbind(g1, g2)), heights = c(3,2))
ggsave(plot = fig1, here("6_Figures/figure1.png"), width = 8, height = 5)




# Figure 2: Vcmax 1:1 and Density of Differences ----------------------------------------------------------------

### Vcmax 1:1 WITH TPU
leaf_sub_vcmax_tpu <- select(pho_stat_tpu, vcmax, curv_meth, leaf_unique, V_cmax_se)
leaf_wide_vcmax_tpu <- reshape(leaf_sub_vcmax_tpu, 
                               idvar = "leaf_unique", 
                               timevar = "curv_meth",
                               direction = "wide")
names(leaf_wide_vcmax_tpu)[2:5]=c("vcmax_DAT", "vcmax_DAT_se", "vcmax_SS", "vcmax_SS_se")
cor3 <- round(cor(leaf_wide_vcmax_tpu$vcmax_DAT, leaf_wide_vcmax_tpu$vcmax_SS), 3)
#leaf_wide_vcmax_tpu <- subset(leaf_wide_vcmax_tpu, select = -tree_id.SS)

pho_1to1_vcmax_tpu <- ggplot(data = leaf_wide_vcmax_tpu, 
                             mapping = aes(x = vcmax_SS,
                                           y = vcmax_DAT,
                                           color = leaf_unique))+
    geom_point()+
    geom_errorbar(aes(ymin = vcmax_DAT - vcmax_DAT_se, ymax = vcmax_DAT + vcmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = vcmax_SS - vcmax_SS_se, xmax = vcmax_SS + vcmax_SS_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6) +
    theme_classic() +
    labs(x = expression(italic("V")[italic("cmax-SS")]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         y = expression(italic("V")[italic("cmax-DAT")]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         col = "Unique Leaf") +
    theme(aspect.ratio = 1,
          axis.title.x = element_text(size = 14, family = "serif"),
          axis.title.y = element_text(size = 14, family = "serif"),
          axis.text.x = element_text(size = 8, family = "serif", color = "gray10"),
          axis.text.y = element_text(size = 8, family = "serif", color = "gray10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 100)) + 
    scale_y_continuous(limits = c(1, 100)) +
    annotate(geom = "text", label = paste0("r = ", cor3), x = 25, y = 75, size = rel(3))
pho_1to1_vcmax_tpu
# ggsave(plot = pho_1to1_vcmax_tpu, here("6_Figures/pho_1to1_vcmax_tpu.png"))


### Vcmax 1:1 WITHOUT TPU
leaf_sub_vcmax <- select(pho_stat, vcmax, curv_meth, leaf_unique, V_cmax_se)
leaf_wide_vcmax <- reshape(leaf_sub_vcmax, 
                           idvar = "leaf_unique", 
                           timevar = "curv_meth",
                           direction = "wide")
names(leaf_wide_vcmax)[2:5] = c("vcmax_DAT", "vcmax_DAT_se", "vcmax_SS", "vcmax_SS_se")
cor1 <- round(cor(leaf_wide_vcmax$vcmax_DAT, leaf_wide_vcmax$vcmax_SS), 3)

pho_1to1_vcmax_NoTPU <- ggplot(data = leaf_wide_vcmax, 
                               mapping = aes(x = vcmax_SS,
                                             y = vcmax_DAT,
                                             color = leaf_unique)) +
    geom_point() +
    geom_errorbar(aes(ymin = vcmax_DAT - vcmax_DAT_se, ymax = vcmax_DAT + vcmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = vcmax_SS - vcmax_SS_se, xmax = vcmax_SS + vcmax_SS_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6) +
    theme_classic() +
    labs(x = expression(italic("V")[italic("cmax-SS")]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         y = expression(italic("V")[italic("cmax-DAT")]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         col = "Unique Leaf") +
    theme(aspect.ratio = 1,
          axis.title.x = element_text(size = 14, family = "serif"),
          axis.title.y = element_text(size = 14, family = "serif"),
          axis.text.x = element_text(size = 8, family = "serif", color = "gray10"),
          axis.text.y = element_text(size = 8, family = "serif", color = "gray10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 100)) + 
    scale_y_continuous(limits = c(1, 100)) +
    annotate(geom = "text", label = paste0("r = ", cor1), x = 25, y = 75, size = rel(3))
pho_1to1_vcmax_NoTPU
# ggsave(plot = pho_1to1_vcmax_NoTPU, here("6_Figures/pho_1to1_vcmax_NoTPU.png"))


### Vcmax Density of differences WITH TPU
sp_diff_hist_vc_tpu <- all_diff_tpu2 %>%
    ggplot(aes(x = vc_diff)) +
    stat_function(fun = dnorm, 
                  args = list(0, sd(all_diff_tpu2$vc_diff)), 
                  color = "red", linetype = "dashed") +
    geom_density(linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = mean(all_diff_tpu2$vc_diff), color = "black", alpha = 0.4) +
    xlab(expression("SS - DAT "*italic("V")[italic("cmax")]* " " *(mu*mol~m^{-2}~s^{-1}))) +
    ylab("Density") +
    xlim(-30, 50) +
    ylim(0, 0.25) +
    theme_classic() +
    theme(axis.title.x = element_text(size = 14, family = "serif"),
          axis.title.y = element_text(size = 14, family = "serif"),
          axis.text.x = element_text(size = 8, family = "serif", color = "gray10"),
          axis.text.y = element_text(size = 8, family = "serif", color = "gray10")) +
    annotate("text", x = 30, y = 0.18, 
             label = paste0("Mean = ", round(mean(all_diff_tpu2$vc_diff), digits = 2)), 
             size = rel(3)) +
    annotate("text", x = 30, y = 0.15, 
             label = paste0("SD = ", round(sd(all_diff_tpu2$vc_diff), digits = 2)), 
             size = rel(3))
sp_diff_hist_vc_tpu


### Vcmax density of differences WITHOUT TPU
sp_diff_hist_vc_notpu <- all_diff_notpu2 %>% 
    ggplot(aes(x = vc_diff)) +
    stat_function(fun = dnorm, 
                  args = list(0, sd(all_diff_notpu2$vc_diff)), 
                  color = "red", linetype = "dashed")+
    geom_density(linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = mean(all_diff_notpu2$vc_diff), color = "black", alpha = 0.4) +
    xlab(expression("SS - DAT "*italic("V")[italic("cmax")]*  " " *(mu*mol~m^{-2}~s^{-1}))) +
    ylab("Density") +
    xlim(-30, 50) +
    ylim(0, 0.25) +
    theme_classic() +
    theme(axis.title.x = element_text(size = 14, family = "serif"),
          axis.title.y = element_text(size = 14, family = "serif"),
          axis.text.x = element_text(size = 8, family = "serif", color = "gray10"),
          axis.text.y = element_text(size = 8, family = "serif", color = "gray10")) +
    annotate("text", x = 30, y = 0.18, 
             label = paste0("Mean = ", round(mean(all_diff_notpu2$vc_diff), digits = 2)), 
             size = rel(3)) +
    annotate("text", x = 30, y = 0.15, 
             label = paste0("SD = ", round(sd(all_diff_notpu2$vc_diff), digits = 2)), 
             size = rel(3))
sp_diff_hist_vc_notpu


gQ <- ggplotGrob(pho_1to1_vcmax_tpu)
gR <- ggplotGrob(sp_diff_hist_vc_tpu)
gS <- ggplotGrob(pho_1to1_vcmax_NoTPU)
gT <- ggplotGrob(sp_diff_hist_vc_notpu)

fig2 <- grid.arrange(arrangeGrob(cbind(gQ, gR), arrangeGrob(cbind(gS, gT))))
ggsave(plot = fig2, here("6_Figures/figure2.png"), width = 6.5, height = 5)



# Figure 3: Jmax 1:1 and Density of Differences ----------------------------------------------------------------


### Jmax 1:1 WITHOUT TPU
leaf_sub_jmax <- select(pho_stat, jmax, curv_meth, leaf_unique, J_se)
leaf_wide_jmax <- reshape(leaf_sub_jmax,
                          idvar = "leaf_unique",
                          timevar = "curv_meth", 
                          direction = "wide")
names(leaf_wide_jmax)[2:5]=c("jmax_DAT", "jmax_DAT_se", "jmax_SS", "jmax_SS_se")
cor2 <- round(cor(leaf_wide_jmax$jmax_DAT, leaf_wide_jmax$jmax_SS), 3)

pho_1to1_jmax_noTPU <- ggplot(data = leaf_wide_jmax, 
                              mapping = aes(x = jmax_SS,
                                            y = jmax_DAT,
                                            color = leaf_unique))+
    geom_point() +
    geom_errorbar(aes(ymin = jmax_DAT - jmax_DAT_se, ymax = jmax_DAT + jmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = jmax_SS - jmax_SS_se, xmax = jmax_SS + jmax_SS_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6) +
    theme_classic() +
    labs(x = expression(italic("J")[italic("max-SS")]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         y = expression(italic("J")[italic("max-DAT")]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         col = "Unique Leaf") +
    theme(aspect.ratio = 1,
          axis.title.x = element_text(size = 14, family = "serif"),
          axis.title.y = element_text(size = 14, family = "serif"),
          axis.text.x = element_text(size = 8, family = "serif", color = "grey10"),
          axis.text.y = element_text(size = 8, family = "serif", color = "grey10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 130)) + 
    scale_y_continuous(limits = c(1, 130)) +
    annotate(geom = "text", label = paste0("r = ", cor2), x = 100, y = 30, size = rel(3))
pho_1to1_jmax_noTPU
# ggsave(plot = pho_1to1_jmax_noTPU, here("6_Figures/pho_1to1_jmax_NoTPU.png"))


### Jmax 1:1 WITH TPU
leaf_sub_jmax_tpu <- select(pho_stat_tpu, jmax, curv_meth, leaf_unique, J_se)
leaf_wide_jmax_tpu <- reshape(leaf_sub_jmax_tpu, 
                              idvar = "leaf_unique",
                              timevar = "curv_meth", 
                              direction = "wide")
names(leaf_wide_jmax_tpu)[2:5] = c("jmax_DAT", "jmax_DAT_se", "jmax_SS", "jmax_SS_se")
cor4 <- round(cor(leaf_wide_jmax_tpu$jmax_DAT, leaf_wide_jmax_tpu$jmax_SS), 3)

pho_1to1_jmax_tpu <- ggplot(data = leaf_wide_jmax_tpu,
                            mapping = aes(x = jmax_SS,
                                          y = jmax_DAT,
                                          color = leaf_unique)) +
    geom_point() +
    geom_errorbar(aes(ymin = jmax_DAT - jmax_DAT_se, ymax = jmax_DAT + jmax_DAT_se)) + 
    geom_errorbarh(aes(xmin = jmax_SS - jmax_SS_se, xmax = jmax_SS + jmax_SS_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6) +
    theme_classic() +
    labs(x = expression(italic("J")[italic("max-SS")]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         y = expression(italic("J")[italic("max-DAT")]* " " *(mu*mol~m^{-2}~s^{-1} *"")),
         col = "Unique Leaf") +
    theme(aspect.ratio = 1,
          axis.title.x = element_text(size = 12, family = "serif"),
          axis.title.y = element_text(size = 12, family = "serif"),
          axis.text.x = element_text(size = 8, family = "serif", color = "grey10"),
          axis.text.y = element_text(size = 8, family = "serif", color = "grey10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(1, 130)) + 
    scale_y_continuous(limits = c(1, 130)) +
    annotate(geom = "text", label = paste0("r = ", cor4), x = 100, y = 30, size = rel(3))
pho_1to1_jmax_tpu
# ggsave(plot = pho_1to1_jmax_tpu, here("6_Figures/pho_1to1_datvSS_jmax_tpu.png"))


### Jmax WITHOUT TPU
sp_diff_hist_j_notpu <- all_diff_notpu2 %>% 
    ggplot(aes(x = j_diff)) +
    stat_function(fun = dnorm, 
                  args = list(0, sd(all_diff_notpu2$j_diff)), 
                  linetype = "dashed", color = "red")+
    geom_density(linewidth = 0.8)+
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3)+
    geom_vline(xintercept = mean(all_diff_notpu2$j_diff), color = "black", alpha = 0.4)+
    xlab(expression("SS - DAT "*italic("J")[italic("max")]* " " * (mu*mol~m^{-2}~s^{-1}))) +
    ylab("Density")+
    xlim(-30, 50)+
    ylim(0, 0.25)+
    theme_classic()+
    theme(axis.title.x = element_text(size = 14, family = "serif"),
          axis.title.y = element_text(size = 14, family = "serif"),
          axis.text.x = element_text(size = 8, family = "serif", color = "gray10"),
          axis.text.y = element_text(size = 8, family = "serif", color = "gray10"))+
    annotate("text", x = 30, y = 0.18, 
             label = paste0("Mean = ", round(mean(all_diff_notpu2$j_diff), digits = 2)), 
             size = rel(3))+
    annotate("text", x = 30, y = 0.15, 
             label = paste0("SD = ", round(sd(all_diff_notpu2$j_diff), digits = 2)), 
             size = rel(3))
sp_diff_hist_j_notpu



### Jmax WITH TPU
sp_diff_hist_j_tpu <- all_diff_tpu2 %>% 
    ggplot(aes(x = j_diff)) +
    stat_function(fun = dnorm, 
                  args = list(0, sd(all_diff_tpu2$j_diff)), 
                  color = "red", linetype = "dashed") +
    geom_density(linewidth = 0.8) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed", alpha = 0.3) +
    geom_vline(xintercept = mean(all_diff_tpu2$j_diff), color = "black", alpha = 0.4) +
    xlab(expression("SS - DAT "*italic("J")[italic("max")]* " " * (mu*mol~m^{-2}~s^{-1}))) +
    ylab("Density") +
    xlim(-30, 50) +
    ylim(0, 0.25) +
    theme_classic() +
    theme(axis.title.x = element_text(size = 14, family = "serif"),
          axis.title.y = element_text(size = 14, family = "serif"),
          axis.text.x = element_text(size = 8, family = "serif", color = "gray10"),
          axis.text.y = element_text(size = 8, family = "serif", color = "gray10")) +
    annotate("text", x = 30, y = 0.18, 
             label = paste0("Mean = ", round(mean(all_diff_tpu2$j_diff), digits = 2)),
             size = rel(3)) +
    annotate("text", x = 30, y = 0.15, 
             label = paste0("SD = ", round(sd(all_diff_tpu2$j_diff), digits = 2)),
             size = rel(3))
sp_diff_hist_j_tpu





gQ2 <- ggplotGrob(pho_1to1_jmax_tpu)
gR2 <- ggplotGrob(sp_diff_hist_j_tpu)
gS2 <- ggplotGrob(pho_1to1_jmax_noTPU)
gT2 <- ggplotGrob(sp_diff_hist_j_notpu)

fig3 <- grid.arrange(arrangeGrob(cbind(gQ2, gR2), arrangeGrob(cbind(gS2, gT2))))
ggsave(plot = fig3, "Figures/figure3.png", width = 6.5, height = 5)


# Figure 4 Histograms of differences by species, sorted by relative canopy height ----------

## Note that error bars represent the mean of the absolute error for the differences

### Read in the data
all_diff_tpu_codes <- read.csv(here("5_Results/species_diffs_summary_TPU.csv"))
all_diff_notpu_codes <- read.csv(here("5_Results/species_diffs_summary_noTPU.csv"))


## Vcmax WITH TPU
vc_diff_hist <- ggplot(data = all_diff_tpu_codes, 
                       aes(x = reorder(gen_spec_id, desc(rel_can_pos)),
                           y = vc_diff)) +
    geom_bar(stat = "identity", fill = "cadetblue2", color = "grey20") +
    labs(x = NULL,
         y = expression("SS - DAT "*italic("V")[italic("cmax")]* " " *(mu*mol~m^{-2}~s^{-1}))) +
    geom_errorbar(aes(x = gen_spec_id, 
                      ymin = vc_diff - vc_diff_se,
                      ymax = vc_diff + vc_diff_se),
                  width = 0.3, colour = "#CA0068", alpha = 0.9, linewidth = 0.5) +
    theme_classic(base_family = "serif") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    ylim(-20, 50)+
    theme(axis.text.y = element_text(face = "italic"))+
    coord_flip()
vc_diff_hist
#ggsave(plot = vc_diff_hist, here("6_Figures/vc_diff_hist.png"))


# Jmax WITH TPU
j_diff_hist <- ggplot(data = all_diff_tpu_codes, 
                      aes(x = reorder(gen_spec_id, desc(rel_can_pos)),
                          y = j_diff)) +
    geom_bar(stat = "identity", fill = "cadetblue2", color = "grey20") +
    labs(x = NULL,
         y = expression("SS - DAT "*italic("J")[italic("max")]* " " *(mu*mol~m^{-2}~s^{-1}))) +
    geom_errorbar(aes(x = gen_spec_id, 
                      ymin = j_diff - j_diff_se,
                      ymax = j_diff + j_diff_se),
                  width = 0.3, colour = "#CA0068", alpha = 0.9, size = 0.5)+
    theme_classic(base_family = "serif") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    ylim(-20, 50) +
    theme(axis.text.y = element_text(face = "italic"))+
    coord_flip()
j_diff_hist
#ggsave(plot = j_diff_hist, here("6_Figures/j_diff_hist.png"))


# plot_arranged <- grid.arrange(vc_diff_hist, j_diff_hist)
# ggsave(plot = plot_arranged, here("6_Figures/diff_histos.png"), width = 4.3, height = 7)



## Vcmax WITHOUT TPU
vc_diff_hist_notpu <- ggplot(data = all_diff_notpu_codes, 
                             aes(x = reorder(gen_spec_id, desc(rel_can_pos)),
                                 y = vc_diff)) +
    geom_bar(stat = "identity", fill = "cadetblue2", color = "grey20") +
    labs(x = NULL,
         y = expression("SS - DAT "*italic("V")[italic("cmax")]* " " *(mu*mol~m^{-2}~s^{-1}))) +
    geom_errorbar(aes(x = gen_spec_id, 
                      ymin = vc_diff - vc_diff_se, 
                      ymax = vc_diff + vc_diff_se),
                  width = 0.3, colour = "#CA0068", alpha = 0.9, linewidth = 0.5)+
    theme_classic(base_family = "serif") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    ylim(-20, 50) +
    theme(axis.text.y = element_text(face = "italic")) +
    coord_flip()
vc_diff_hist_notpu
#ggsave(plot = vc_diff_hist, here("6_Figures/vc_diff_hist_notpu.png"))


## Jmax WITHOUT TPU
j_diff_hist_notpu <- ggplot(data = all_diff_notpu_codes,
                            aes(x = reorder(gen_spec_id, desc(rel_can_pos)),
                                y = j_diff)) +
    geom_bar(stat = "identity", fill = "cadetblue2", color = "grey20") +
    labs(x = NULL,
         y = expression("SS - DAT "*italic("J")[italic("max")]* " " *(mu*mol~m^{-2}~s^{-1}))) +
    geom_errorbar(aes(x = gen_spec_id, 
                      ymin = j_diff - j_diff_se,
                      ymax = j_diff + j_diff_se), 
                  width = 0.3, colour = "#CA0068", alpha = 0.9, size = 0.5)+
    theme_classic(base_family = "serif") +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.8) +
    ylim(-20, 50) +
    theme(axis.text.y = element_text(face = "italic"))+
    coord_flip()
j_diff_hist_notpu
#ggsave(plot = j_diff_hist_notpu, here("6_Figures/j_diff_hist_notpu.png"))



gA <- ggplotGrob(vc_diff_hist)
gB <- ggplotGrob(vc_diff_hist_notpu + rremove("y.text"))
gC <- ggplotGrob(j_diff_hist)
gD <- ggplotGrob(j_diff_hist_notpu + rremove("y.text"))

plot_arranged3 <- grid.arrange(arrangeGrob(cbind(gA, gB), arrangeGrob(cbind(gC, gD))))
ggsave(plot = plot_arranged3, here("6_Figures/figure4.png"), width = 8.3, height = 5)




# Table S1: Summary stats ------------------------------------------
## NEED TO REARRANGE! CHARLIE MADE THIS ON 1/11/24

pho_stat_tpu <- select(pho_leaf_tpu, 'curv_meth', 
                       'Best_Vcmax_25C', 'Best_Jmax_25C', 'V_TPU', 
                       'ID', 'leaf_unique', 'V_cmax_se', 'J_se', 'V_TPU_se',
                       treeid # see if we can just use this column
                       )


pho_stat_tpu_ss <- rename(pho_stat_tpu,
                          vcmax = Best_Vcmax_25C,
                          jmax = Best_Jmax_25C,
                          tpu = V_TPU,
                          vcmax_se = V_cmax_se,
                          jmax_se = J_se,
                          leaf_id = ID) %>% 
    mutate(fit_type = "tpu") %>% 
    filter(curv_meth == "SS")

treename <- as.numeric(substring(pho_stat_tpu_ss$leaf_unique, 4, 5))  ##### Pretty sure we can use treeid for this

pho_stat_tpu_ss$treename <- treename

ss_keyparams_tpu <- pho_stat_tpu_ss %>%
    group_by(treename) %>% 
    summarize(mean_vcmax = mean(vcmax),
              sd_vcmax = sd(vcmax),
              mean_jmax = mean(jmax),
              sd_jmax = sd(jmax))

ss_keyparams_tpu_codes <- left_join(ss_keyparams_tpu, codebook, by = "treename")

######### END AWKWARD SECTION


########## Don't we need this for the other ones, not just SS TPU?




# Figure S3: TPU 1:1 -------------------------------------------------------

# TPU 1:1
leaf_sub_tpu <- select(pho_stat_tpu, tpu, curv_meth, leaf_unique, V_TPU_se)
leaf_wide_tpu <- reshape(leaf_sub_tpu,
                         idvar = "leaf_unique",
                         timevar = "curv_meth", 
                         direction = "wide")
names(leaf_wide_tpu)[2:5] = c("tpu_DAT", "tpu_DAT_se", "tpu_SS", "tpu_SS_se")
cor5 <- round(cor(leaf_wide_tpu$tpu_DAT, leaf_wide_tpu$tpu_SS), 3)
only_tpu_fit <- filter(leaf_wide_tpu, !is.na(tpu_DAT) & !is.na(tpu_SS))
only_tpu_fit$tpu_DAT_se <- as.numeric(only_tpu_fit$tpu_DAT_se)
only_tpu_fit$tpu_SS_se <- as.numeric(only_tpu_fit$tpu_SS_se)


pho_1to1_tpu_tpu <- ggplot(data = only_tpu_fit,
                           mapping = aes(x = tpu_SS,
                                         y = tpu_DAT,
                                         color = leaf_unique)) +
    geom_point(cex = 2.5) +
    geom_errorbar(aes(ymin = tpu_DAT - tpu_DAT_se, ymax = tpu_DAT + tpu_DAT_se)) + 
    geom_errorbarh(aes(xmin = tpu_SS - tpu_SS_se, xmax = tpu_SS + tpu_SS_se)) +
    geom_abline(intercept = 0, slope = 1, linetype = 5, linewidth = 0.6) +
    theme_classic() +
    labs(x = expression("TPU"[SS] * " " * (mu*mol~m^{-2}~s^{-1})),
         y = expression("TPU"[DAT] * " " * (mu*mol~m^{-2}~s^{-1})), col = "Unique Leaf") +
    theme(aspect.ratio = 1,
          axis.title.x = element_text(size = 12, family = "serif"),
          axis.title.y = element_text(size = 12, family = "serif"),
          axis.text.x = element_text(size = 8, family = "serif", color = "grey10"),
          axis.text.y = element_text(size = 8, family = "serif", color = "grey10"),
          legend.position = "none") +
    scale_x_continuous(limits = c(0, 12), breaks = c(0,3,6,9,12)) + 
    scale_y_continuous(limits = c(0, 12), breaks = c(0,3,6,9,12)) +
    annotate(geom = "text", label = paste0("r = ", cor5), x = 4, y = 8, size = rel(3))
pho_1to1_tpu_tpu
# ggsave(plot = pho_1to1_tpu_tpu, here("6_Figures/figureS3.png"))







# Figure S4: Double Boxplots: All v No OS --------------------------------------------

#### Again, get the data that needs to be read into these


#### TPU-omitted Boxplots
pho_nd_stat <- mutate(pho_nd_stat, subset = "nOS")
pho_stat <- mutate(pho_stat, subset = "wOS")
all_and_nOS_noTPU <- rbind(pho_stat, pho_nd_stat)


## Vcmax WITHOUT TPU
vcmax_AllvNoOS_noTPU <- all_and_nOS_noTPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "SS"), 
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    theme_classic() +
    labs(x="Fit Type",
         y = expression(italic(V[cmax])*" "*(mu*mol~m^{-2}~s^{-1}))) +
    theme(axis.title.x = element_text(size = 16, family = "serif"),
          axis.title.y = element_text(size = 16, family = "serif"),
          axis.text.x = element_text(size = 13, family = "serif", color = "gray10"),
          axis.text.y = element_text(size = 13, family = "serif", colour = "gray10")) +
    scale_y_continuous(limits = c(0, 85))
vcmax_AllvNoOS_noTPU
#gsave(plot = vcmax_AllvNoOS_noTPU, here("6_Figures/box_vcmax_AllvNoOS_noTPU.png"))


### Jmax WITHOUT TPU
jmax_AllvNoOS_noTPU <- all_and_nOS_noTPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = jmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "SS"),
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    theme_classic() +
    labs(x = "Fit Type",
         y = expression(italic(J[max])*" "*(mu*mol~m^{-2}~s^{-1})))+
    theme(axis.title.x = element_text(size = 16, family = "serif"),
          axis.title.y = element_text(size = 16, family = "serif"),
          axis.text.x = element_text(size = 13, family = "serif", color = "grey10"),
          axis.text.y = element_text(size = 13, family = "serif", color = "grey10")) +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_AllvNoOS_noTPU
# ggsave(plot = jmax_AllvNoOS_noTPU, here("6_Figures/box_jmax_AllvNoOS_noTPU.png"))




#### TPU-enabled boxplots

pho_nd_stat_tpu <- mutate(pho_nd_stat_tpu, subset = "nOS")
pho_stat_tpu <- mutate(pho_stat_tpu, subset = "wOS")
all_and_nOS_TPU <- rbind(pho_stat_tpu, pho_nd_stat_tpu)


### Vcmax WITH TPU
vcmax_AllvNoOS_TPU <- all_and_nOS_TPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "SS"), 
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    theme_classic() +
    labs(x = "Fit Type",
         y = expression(italic(V[cmax])*" "*(mu*mol~m^{-2}~s^{-1}))) +
    theme(axis.title.x = element_text(size = 16, family = "serif"),
          axis.title.y = element_text(size = 16, family = "serif"),
          axis.text.x = element_text(size = 13, family = "serif", color = "gray10"),
          axis.text.y = element_text(size = 13, family = "serif", colour = "gray10")) +
    scale_y_continuous(limits = c(0, 85))
vcmax_AllvNoOS_TPU
# ggsave(plot = vcmax_AllvNoOS_TPU, here("6_Figures/box_vcmax_AllvNoOS_TPU.png"))


### Jmax WITH TPU
jmax_AllvNoOS_TPU <- all_and_nOS_TPU %>% 
    mutate(subset = fct_reorder(subset, subset , .fun = "length", .desc = TRUE)) %>% 
    ggplot(aes(x = subset, y = jmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "SS"),
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("All", "No Overshoot")) +
    theme_classic()+
    labs(x = "Fit Type",
         y = expression(italic(J[max])*" "*(mu*mol~m^{-2}~s^{-1}))) +
    theme(axis.title.x = element_text(size = 16, family = "serif"),
          axis.title.y = element_text(size = 16, family = "serif"),
          axis.text.x = element_text(size = 13, family = "serif", color = "grey10"),
          axis.text.y = element_text(size = 13, family = "serif", color = "grey10")) + 
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_AllvNoOS_TPU
# ggsave(plot = jmax_AllvNoOS_TPU, here("6_Figures/box_jmax_AllvNoOS_TPU.png"))


gW <- ggplotGrob(vcmax_AllvNoOS_TPU + rremove("xlab") + rremove("legend"))
gX <- ggplotGrob(vcmax_AllvNoOS_noTPU + rremove("ylab") + rremove("xlab"))
gY <- ggplotGrob(jmax_AllvNoOS_TPU + rremove("legend"))
gZ <- ggplotGrob(jmax_AllvNoOS_noTPU + rremove("ylab"))

figS4 <- grid.arrange(arrangeGrob(cbind(gW, gX), arrangeGrob(cbind(gY, gZ))))

ggsave(plot = figS4, here("6_Figures/figureS4.png"), width = 6.5, height = 5)



# Figure S5: Density distribution of standard errors ---------------------------------

# colors: blue = "#31688EFF", yellow = "#FDE725FF", yellow text = "#FFBF00"

### Read in the data
pho_dat <- read.csv(file = here("5_Results/DAT_photo_pars_crct_noTPU.csv"), 
                    sep = ",", 
                    header = TRUE, na.strings = 1000) %>% 
    ## TPU values at 1000 are coded as NA
    subset(ID != "K6709L3")

pho_SS <- read.csv(file = here("5_Results/SS_photo_pars_crct_noTPU.csv"),
                   sep = ",", 
                   header = TRUE, na.strings = 1000) %>% 
    subset(ID != "K6709L3")



#### Vcmax WITHOUT TPU
hist_vc_se_notpu <- ggplot(mapping = aes(x = V_cmax_se)) +
    #stat_function(data = pho_dat, fun = dnorm, args = list(0, sd(pho_dat$V_cmax_se)), 
    #              color = "black", linetype = "dashed") +
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
    annotate("text", x = 5.5, y = 1.8,
             label = "DAT", size = rel(2.7), color = "#31688EFF") +
    annotate("text", x = 9, y = 1.8, 
             label = "SS", size = rel(2.7), color = "#FFBF00") +
    annotate("text", x = 2.5, y = 1.6, 
             label = "Mean SE", size = rel(2.7)) +
    annotate("text", x = 2.5, y = 1.4, 
             label = "SD SE", size = rel(2.7)) +
    annotate("text", x = 2.5, y = 1.2, 
             label = "Range SE", size = rel(2.7)) +
    annotate("text", x = 5.5, y = 1.6,
             label = round(mean(pho_dat$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 5.5, y = 1.4, 
             label = round(sd(pho_dat$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 5.5, y = 1.2, 
             label = paste0(round(min(pho_dat$V_cmax_se), digits = 2), ", ", 
                            round(max(pho_dat$V_cmax_se), digits = 2)), 
             size = rel(2.7)) +
    annotate("text", x = 9, y = 1.6, 
             label = round(mean(pho_SS$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 9, y = 1.4, 
             label = round(sd(pho_SS$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 9, y = 1.2, 
             label = paste0(round(min(pho_SS$V_cmax_se), digits = 2), ", ", 
                            round(max(pho_SS$V_cmax_se), digits = 2)), 
             size = rel(2.7)) +
    xlab(expression(italic("V")[italic("cmax")]* " SE " *(mu*mol~m^{-2}~s^{-1} *""))) +
    ylab("Density") +
    xlim(-3,10) +
    ylim(0, 2) +
    theme_classic()
hist_vc_se_notpu


#### Vcmax WITH TPU
hist_vc_se_tpu <- ggplot(mapping = aes(x = V_cmax_se)) +
    # stat_function(data = pho_dat_tpu, fun = dnorm, 
    #               args = list(0, sd(pho_dat_tpu$V_cmax_se)), color = "black", 
    #               linetype = "dashed") +
    # Density and mean of DAT SE
    geom_density(data = pho_dat_tpu, linewidth = 0.8, color = "#31688EFF") +
    geom_vline(xintercept = mean(pho_dat_tpu$V_cmax_se), color = "#31688EFF", alpha = 0.5) +
    # Density and mean of SS SE
    geom_density(data = pho_SS_tpu, linewidth = 0.8, color = "#FDE725FF") +
    geom_vline(xintercept = mean(pho_SS_tpu$V_cmax_se), color = "#FDE725FF", alpha = 0.5) +
    # Add in a zero line
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.8) +
    # Add in the statistics
    annotate("text", x = 5.5, y = 1.8, 
             label = "DAT", size = rel(2.7), color = "#31688EFF") +
    annotate("text", x = 9, y = 1.8,
             label = "SS", size = rel(2.7), color = "#FFBF00") +
    annotate("text", x = 2.5, y = 1.6, 
             label = "Mean SE", size = rel(2.7)) +
    annotate("text", x = 2.5, y = 1.4, 
             label = "SD SE", size = rel(2.7)) +
    annotate("text", x = 2.5, y = 1.2, 
             label = "Range SE", size = rel(2.7)) +
    annotate("text", x = 5.5, y = 1.6, 
             label = round(mean(pho_dat_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7)) +
    annotate("text", x = 5.5, y = 1.4, 
             label = round(sd(pho_dat_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7))+
    annotate("text", x = 5.5, y = 1.2, 
             label = paste0(round(min(pho_dat_tpu$V_cmax_se), digits = 2), ", ", 
                            round(max(pho_dat_tpu$V_cmax_se), digits = 2)), 
             size = rel(2.7)) +
    annotate("text", x = 9, y = 1.6, 
             label = round(mean(pho_SS_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7)) +
    annotate("text", x = 9, y = 1.4, 
             label = round(sd(pho_SS_tpu$V_cmax_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7)) +
    annotate("text", x = 9, y = 1.2, 
             label = paste0(round(min(pho_SS_tpu$V_cmax_se), digits = 2), ", ", 
                            round(max(pho_SS_tpu$V_cmax_se), digits = 2)), 
             size = rel(2.7)) +
    xlab(expression(italic("V")[italic("cmax")]* " SE " *(mu*mol~m^{-2}~s^{-1} *""))) +
    ylab("Density")+
    xlim(-3,10)+
    ylim(0, 2)+
    theme_classic()
hist_vc_se_tpu



### Jmax WITHOUT tpu
hist_j_se_notpu <- ggplot(mapping = aes(x = J_se)) +
    # stat_function(data = pho_dat, fun = dnorm, 
    #               args = list(0, sd(pho_dat$J_se)), color = "black", 
    #               linetype = "dashed") +
    # Density and mean of DAT SE
    geom_density(data = pho_dat, linewidth = 0.8, color = "#31688EFF") +
    geom_vline(xintercept = mean(pho_dat$J_se), color = "#31688EFF", alpha = 0.5) +
    # Density and mean of SS SE
    geom_vline(xintercept = mean(pho_SS$J_se), color = "#FDE725FF", alpha = 0.5) +
    geom_density(data = pho_SS, linewidth = 0.8, color = "#FDE725FF") +
    # Add in a zero line
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.8) +
    # Add in the statistics
    annotate("text", x = 1.0, y = 18, label = "DAT", size = rel(2.7), color = "#31688EFF")+ 
    annotate("text", x = 1.7, y = 18, label = "SS", size = rel(2.7), color = "#FFBF00") +
    annotate("text", x = 0.4, y = 16, label = "Mean SE", size = rel(2.7)) +
    annotate("text", x = 0.4, y = 14, label = "SD SE", size = rel(2.7)) +
    annotate("text", x = 0.4, y = 12, label = "Range SE", size = rel(2.7)) +
    annotate("text", x = 1.0, y = 16, 
             label = round(mean(pho_dat$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7)) +
    annotate("text", x = 1.0, y = 14, 
             label = round(sd(pho_dat$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7)) +
    annotate("text", x = 1.0, y = 12, 
             label = paste0(round(min(pho_dat$J_se), digits = 2), ", ", 
                            round(max(pho_dat$J_se), digits = 2)), 
             size = rel(2.7)) +
    annotate("text", x = 1.7, y = 16, 
             label = round(mean(pho_SS$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7)) +
    annotate("text", x = 1.7, y = 14, 
             label = round(sd(pho_SS$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7)) +
    annotate("text", x = 1.7, y = 12, 
             label = paste0(round(min(pho_SS$J_se), digits = 2), ", ", 
                            round(max(pho_SS$J_se), digits = 2)), 
             size = rel(2.7)) +
    xlab(expression(italic("J")[italic("max")]* " SE " *(mu*mol~m^{-2}~s^{-1} *""))) +
    ylab("Density") +
    xlim(-1,2) +
    ylim(0, 20) +
    theme_classic()
hist_j_se_notpu


#### Jmax WITH TPU
hist_j_se_tpu <- ggplot(mapping = aes(x = J_se)) +
    # stat_function(data = pho_dat_tpu, fun = dnorm, 
    #               args = list(0, sd(pho_dat_tpu$J_se)), color = "black", 
    #               linetype = "dashed") +
    # Density and mean of DAT SE
    geom_density(data = pho_dat_tpu, linewidth = 0.8, color = "#31688EFF") +
    geom_vline(xintercept = mean(pho_dat_tpu$J_se), color = "#31688EFF", alpha = 0.5) +
    # Density and mean of SS SE
    geom_density(data = pho_SS_tpu, linewidth = 0.8, color = "#FDE725FF") +
    geom_vline(xintercept = mean(pho_SS_tpu$J_se), color = "#FDE725FF", alpha = 0.5) +
    # Add in a zero line
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.8) +
    # Add in the statistics
    annotate("text", x = 1.0, y = 18, label = "DAT", size = rel(2.7), color = "#31688EFF") +
    annotate("text", x = 1.7, y = 18, label = "SS", size = rel(2.7), color = "#FFBF00") +
    annotate("text", x = 0.4, y = 16, label = "Mean SE", size = rel(2.7)) +
    annotate("text", x = 0.4, y = 14, label = "SD SE", size = rel(2.7)) +
    annotate("text", x = 0.4, y = 12, label = "Range SE", size = rel(2.7)) +
    annotate("text", x = 1.0, y = 16, 
             label = round(mean(pho_dat_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7)) +
    annotate("text", x = 1.0, y = 14, 
             label = round(sd(pho_dat_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7)) +
    annotate("text", x = 1.0, y = 12, 
             label = paste0(round(min(pho_dat_tpu$J_se), digits = 2), ", ", 
                            round(max(pho_dat_tpu$J_se), digits = 2)), 
             size = rel(2.7)) +
    annotate("text", x = 1.7, y = 16, 
             label = round(mean(pho_SS_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7)) +
    annotate("text", x = 1.7, y = 14, 
             label = round(sd(pho_SS_tpu$J_se, na.rm = TRUE), digits = 2), 
             size = rel(2.7)) +
    annotate("text", x = 1.7, y = 12, 
             label = paste0(round(min(pho_SS_tpu$J_se), digits = 2), ", ", 
                            round(max(pho_SS_tpu$J_se), digits = 2)), 
             size = rel(2.7)) +
    xlab(expression(italic("J")[italic("max")]* " SE " *(mu*mol~m^{-2}~s^{-1} *""))) +
    ylab("Density") +
    xlim(-1,2) +
    ylim(0, 20) +
    theme_classic()
hist_j_se_tpu


gM <- ggplotGrob(hist_vc_se_tpu)
gN <- ggplotGrob(hist_vc_se_notpu + rremove("ylab"))
gO <- ggplotGrob(hist_j_se_tpu)
gP <- ggplotGrob(hist_j_se_notpu + rremove("ylab"))

se_arranged <- grid.arrange(arrangeGrob(cbind(gM, gN), arrangeGrob(cbind(gO, gP))))

ggsave(plot = se_arranged, "Figures/figureS5.png", width = 6.5, height = 5)





# TPU v. No TPU Boxplots --------------------------------------------------------------

# Full data

### Vcmax ALL data
vcmax_all_TPUvNoTPU <- ggplot(all_results2, 
                              aes(x = factor(fit_type, levels = c("tpu", "no_tpu")),
                                  y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"), 
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("TPU-enabled", "TPU-omitted")) +
    theme_classic() +
    labs(x = "Fit Type",
         y = expression(italic(V[cmax])* " " *(mu*mol~m^{-2}~s^{-1}))) +
    theme(axis.title.x = element_text(size = 16, family = "serif"),
          axis.title.y = element_text(size = 16, family = "serif"),
          axis.text.x = element_text(size = 13, family = "serif", color = "gray10"),
          axis.text.y = element_text(size = 13, family = "serif", colour = "gray10"))
vcmax_all_TPUvNoTPU
# ggsave(plot = vcmax_all_TPUvNoTPU, here("6_Figures/box_vcmax_all_TPUvNoTPU.png"))

### Jmax ALL data
jmax_all_TPUvNoTPU <- ggplot(all_results2, 
                             aes(x = factor(fit_type, levels = c("tpu", "no_tpu")),
                                 y = jmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("TPU-enabled", "TPU-omitted")) +
    theme_classic() +
    labs(x="Fit Type",
         y = expression(italic(J[max])* " " *(mu*mol~m^{-2}~s^{-1}))) +
    theme(axis.title.x=element_text(size=16, family = "serif"),
          axis.title.y=element_text(size=16, family = "serif"),
          axis.text.x=element_text(size=13, family = "serif", color = "grey10"),
          axis.text.y=element_text(size=13, family = "serif", color = "grey10"))
jmax_all_TPUvNoTPU
# ggsave(plot = jmax_all_TPUvNoTPU, here("6_Figures/box_jmax_all_TPUvNoTPU.png"))



# No Overshoot curves

### Vcmax no overhoot data
vcmax_nOS_TPUvNoTPU <- ggplot(nd_complete, 
                              aes(x = factor(fit_type, levels = c("tpu", "no_tpu")), 
                                  y = vcmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("TPU-enabled", "TPU-omitted")) +
    theme_classic() +
    labs(x = "Fit Type",
         y = expression(italic(V[cmax])*" "*(mu*mol~m^{-2}~s^{-1}))) +
    theme(axis.title.x = element_text(size = 16, family = "serif"),
          axis.title.y = element_text(size = 16, family = "serif"),
          axis.text.x = element_text(size = 13, family = "serif", color = "grey10"),
          axis.text.y = element_text(size = 13, family = "serif", color = "grey10"))
vcmax_nOS_TPUvNoTPU
# ggsave(plot = vcmax_nOS_TPUvNoTPU, here("6_Figures/box_vcmax_nOS_TPUvNoTPU.png"))


### Jmax no overshoot data
jmax_nOS_TPUvNoTPU <- ggplot(nd_complete, 
                             aes(x = factor(fit_type, levels = c("tpu", "no_tpu")),
                                 y = jmax, fill = curv_meth)) +
    geom_boxplot(outlier.shape = 17, outlier.size = 2) +
    scale_fill_manual(name = "Curve Method", labels = c("DAT", "Steady-state"),
                      values = c("#31688EFF", "#FDE725FF")) +
    scale_x_discrete(labels = c("TPU-enabled", "TPU-omitted")) +
    theme_classic() +
    labs(x = "Fit Type",
         y = expression(italic(J[max])*" "*(mu*mol~m^{-2}~s^{-1}))) +
    theme(axis.title.x = element_text(size = 16, family = "serif"),
          axis.title.y = element_text(size = 16, family = "serif"),
          axis.text.x = element_text(size = 13, family = "serif", color = "grey10"),
          axis.text.y = element_text(size = 13, family = "serif", color = "grey10"))
jmax_nOS_TPUvNoTPU
# ggsave(plot = jmax_nOS_TPUvNoTPU, here("6_Figures/box_jmax_nOS_TPUvNoTPU.png"))






# Single boxplots ---------------------------------------------------------

####### Can we just take this out?????, think it's a relic



lab_DATSS <- c('DAT', 'Steady-State')

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
    scale_x_discrete(labels=lab_DATSS)
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
    scale_x_discrete(labels=lab_DATSS) +
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
    scale_x_discrete(labels=lab_DATSS)
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
    scale_x_discrete(labels=lab_DATSS) +
    scale_y_continuous(limits = c(0, 130), breaks = c(0, 40,  80,  120))
jmax_all_TPU
ggsave(plot = jmax_all_TPU, "Figures/photo_box_DvT_jmax_TPU.png")



################## One instance of 'method'

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
    scale_x_discrete(labels=lab_DATSS)
nd_vcmax_box_noTPU
ggsave(plot = nd_vcmax_box_noTPU, "Figures/pho_box_noOS_DvT_vcmax_noTPU.png")

################ another instance of method. Can we use curv_meth?

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
    scale_x_discrete(labels=lab_DATSS)
nd_jmax_box_noTPU
ggsave(plot = nd_jmax_box_noTPU, "Figures/photo_box_nOS_DvT_jmax_noTPU.png")









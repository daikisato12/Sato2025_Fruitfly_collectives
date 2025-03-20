#### load packages ####
targetPackages <- c('tidyverse','arrow','patchwork','ggpubr','gg.layers')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)
# install.packages("remotes")
# remotes::install_github("rpkgs/gg.layers")

#### Figure 3 ####
##### load dataset #####
dfm_s5min_2995_5990_freezing_duration <- read_parquet("../data/3_mutants/dfm_s5min_2995_5990_freezing_duration.parquet") %>%
  ungroup()

dfm_motion_cue_exit_coeff_trial2 <- read_parquet("../data/3_mutants/dfm_motion_cue_exit_coeff_trial2.parquet") %>%
  ungroup()

##### stat #####
stats::wilcox.test(freezing_duration ~ strain, 
                   data = dfm_s5min_2995_5990_freezing_duration %>%
                     filter(strain %in% c("IMPTNT_R9B08", "TNT_R9B08"),
                            sex == "Female", n_inds == "Group"))
stats::wilcox.test(motion_cue_exit_intercept ~ strain, 
                   data = dfm_motion_cue_exit_coeff_trial2 %>%
                     filter(strain %in% c("IMPTNT_R9B08", "TNT_R9B08"),
                            sex == "Female"))

# stats::wilcox.test(freezing_duration ~ strain, 
#                    data = dfm_s5min_2995_5990_freezing_duration %>%
#                      filter(strain %in% c("Ptp99ARNAi_NA", "Ptp99ARNAi_R9B08"),
#                             sex == "Female", n_inds == "Group"))
stats::wilcox.test(motion_cue_exit_intercept ~ strain, 
                   data = dfm_motion_cue_exit_coeff_trial2 %>%
                     filter(strain %in% c("NA_R9B08", "Ptp99ARNAi_R9B08"),
                            sex == "Female"))
stats::wilcox.test(motion_cue_exit_intercept ~ strain, 
                   data = dfm_motion_cue_exit_coeff_trial2 %>%
                     filter(strain %in% c("Ptp99ARNAi_NA", "Ptp99ARNAi_R9B08"),
                            sex == "Female"))

# stats::wilcox.test(freezing_duration ~ strain, 
#                    data = dfm_s5min_2995_5990_freezing_duration %>%
#                      filter(strain %in% c("IMPTNT_Ptp99A", "TNT_Ptp99A"),
#                             sex == "Female", n_inds == "Group"))
stats::wilcox.test(motion_cue_exit_intercept ~ strain, 
                   data = dfm_motion_cue_exit_coeff_trial2 %>%
                     filter(strain %in% c("IMPTNT_Ptp99A", "TNT_Ptp99A"),
                            sex == "Female"))


##### Figure 3a #####
###### freezing duration Ptp99ARNAi_R9B08 ######
my_comparisons_Ptp99ARNAi_R9B08 <- list( c("NA_R9B08", "Ptp99ARNAi_R9B08"), c("Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA") )
g_freezing_Ptp99ARNAi_R9B08 <- 
  ggplot(dfm_s5min_2995_5990_freezing_duration %>%
           filter(strain %in% c("NA_R9B08", "Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA"),
                  sex == "Female") %>%
           transform(strain = factor(strain, levels = c("NA_R9B08", "Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA"))),
         aes(x = strain, y = freezing_duration, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_Ptp99ARNAi_R9B08) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("R9B08-GAL4"), italic("R9B08>Ptp99A-RNAi"), italic("Ptp99A-RNAi"))) +
  scale_y_continuous(breaks = seq(0, 15, 5), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 18)) +
  ylab("Freezing duration (s)") +
  facet_grid(n_inds ~ .) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

g_freezing_Ptp99ARNAi_R9B08

###### freezing duration TNT_Ptp99A ######
my_comparisons_TNT_Ptp99A <- list( c("IMPTNT_Ptp99A", "TNT_Ptp99A") )
g_freezing_TNT_Ptp99A <- 
  ggplot(dfm_s5min_2995_5990_freezing_duration %>%
           filter(strain %in% c("IMPTNT_Ptp99A", "TNT_Ptp99A"),
                  sex == "Female") %>% #, n_inds == "Single"
           transform(strain = factor(strain, levels = c("IMPTNT_Ptp99A", "TNT_Ptp99A"))),
         aes(x = strain, y = freezing_duration, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_TNT_Ptp99A) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("Ptp99A>IMPTNT"), italic("Ptp99A>TNT"))) +
  scale_y_continuous(breaks = seq(0, 15, 5), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 16)) +
  ylab("Freezing duration (s)") +
  facet_grid(n_inds ~ .) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

g_freezing_TNT_Ptp99A

###### freezing duration TNT_R9B08 ######
my_comparisons_TNT_R9B08 <- list( c("IMPTNT_R9B08", "TNT_R9B08") )
g_freezing_TNT_R9B08 <- ggplot(dfm_s5min_2995_5990_freezing_duration %>%
                                  filter(strain %in% c("IMPTNT_R9B08", "TNT_R9B08"),
                                         sex == "Female") %>% #, n_inds == "Single"
                                  transform(strain = factor(strain, levels = c("IMPTNT_R9B08", "TNT_R9B08"))),
                               aes(x = strain, y = freezing_duration, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_TNT_R9B08) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("R9B08>IMPTNT"), italic("R9B08>TNT"))) +
  scale_y_continuous(breaks = seq(0, 15, 5), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 16)) +
  ylab("Freezing duration (s)") +
  facet_grid(n_inds ~ .) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

g_freezing_TNT_R9B08

##### Figure 3b #####
dfm_plot <- 
  dfm_motion_cue_exit_coeff_trial2 %>%
  pivot_longer(cols = contains("motion_cue"), names_to = "var", values_to = "value")

###### visual responsiveness Ptp99ARNAi_R9B08 ######
my_comparisons_Ptp99ARNAi_R9B08 <- list( c("NA_R9B08", "Ptp99ARNAi_R9B08"), c("Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA") )
g_visual_responsiveness_Ptp99ARNAi_R9B08 <- 
  ggplot(dfm_plot %>%
           filter(strain %in% c("NA_R9B08", "Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA"),
                  sex == "Female", var == "motion_cue_exit_intercept") %>%
           transform(strain = factor(strain, levels = c("NA_R9B08", "Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA"))),
         aes(x = strain, y = value, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_Ptp99ARNAi_R9B08) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("R9B08-GAL4"), italic("R9B08>Ptp99A-RNAi"), italic("Ptp99A-RNAi"))) +
  coord_cartesian(ylim = c(-4.5, -2.2)) +
  ylab("Visual responsiveness") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none")

g_visual_responsiveness_Ptp99ARNAi_R9B08


###### visual responsiveness TNT_Ptp99A ######
my_comparisons_TNT_Ptp99A <- list( c("IMPTNT_Ptp99A", "TNT_Ptp99A") )
g_visual_responsiveness_TNT_Ptp99A <- 
  ggplot(dfm_plot %>%
           filter(strain %in% c("IMPTNT_Ptp99A", "TNT_Ptp99A"),
                  sex == "Female", var == "motion_cue_exit_intercept") %>%
           transform(strain = factor(strain, levels = c("IMPTNT_Ptp99A", "TNT_Ptp99A"))),
         aes(x = strain, y = value, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_TNT_Ptp99A) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("Ptp99A>IMPTNT"), italic("Ptp99A>TNT"))) +
  coord_cartesian(ylim = c(-4.5, -2.2)) +
  ylab("Visual responsiveness") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none")

g_visual_responsiveness_TNT_Ptp99A


###### visual responsiveness TNT_R9B08 ######
my_comparisons_TNT_R9B08 <- list( c("IMPTNT_R9B08", "TNT_R9B08") )
g_visual_responsiveness_TNT_R9B08 <- 
  ggplot(dfm_plot %>%
           filter(strain %in% c("IMPTNT_R9B08", "TNT_R9B08"),
                  sex == "Female", var == "motion_cue_exit_intercept") %>%
           transform(strain = factor(strain, levels = c("IMPTNT_R9B08", "TNT_R9B08"))),
         aes(x = strain, y = value, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_TNT_R9B08) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("R9B08>IMPTNT"), italic("R9B08>TNT"))) +
  coord_cartesian(ylim = c(-4.5, -2.2)) +
  ylab("Visual responsiveness") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none")

g_visual_responsiveness_TNT_R9B08


g3ab <- 
  g_freezing_TNT_R9B08 + g_freezing_TNT_Ptp99A + g_freezing_Ptp99ARNAi_R9B08 +
  g_visual_responsiveness_TNT_R9B08 + g_visual_responsiveness_TNT_Ptp99A + g_visual_responsiveness_Ptp99ARNAi_R9B08 +
  plot_layout(nrow = 2, ncol = 3, widths = c(0.66, 0.66, 1), heights = c(1.5, 1))
g3ab
ggsave("../figures/Figure3ab.pdf", g3ab, w = 6, h = 6)



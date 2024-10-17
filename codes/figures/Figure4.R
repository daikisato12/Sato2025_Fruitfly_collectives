#### load packages ####
library(tidyverse)
library(gtools)
library(ggpubr)
library(arrow)


#### Figure 4a ####
##### load dataset #####
dfd_f10min_speed_sgd_strain <- 
  read_parquet("../data/2_mixed_strain/dfd_f10min_speed_sgd_strain.parquet") %>%
  ungroup()

dfd_f5min_speed_sgd_ave_strain <- 
  dfd_f10min_speed_sgd_strain %>%
  filter(seconds_total < 300) %>%
  group_by(sex, strain, n_inds, var, type, mixed, alpha) %>%
  dplyr::summarize(speed_f5min_ave = mean(speed, na.rm = T))

dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain <- read_parquet("../data/2_mixed_strain/dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain.parquet") %>%
  ungroup()

dfd_s5min_2995_5990_speed_sgd_normbyf5minave_allmean_pairedwilcox_res <- data.frame()
comp_list <- c("Single_Mean-Group_Mean", "Group_Mean-Group_Mixed")
for (v in comp_list){
  var1 <- str_split(v, "-")[[1]][1]
  var2 <- str_split(v, "-")[[1]][2]
  
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_allmean_pairedwilcox_res <- 
    bind_rows(dfd_s5min_2995_5990_speed_sgd_normbyf5minave_allmean_pairedwilcox_res,
              dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain %>%
                filter(var %in% c(var1, var2)) %>%
                group_by(stim_time) %>%
                rstatix::wilcox_test(speed_norm_mean ~ var, paired = TRUE) %>%
                rstatix::adjust_pvalue(method = "bonferroni") %>%
                rstatix::add_significance("p.adj")
    )
}
dfd_s5min_2995_5990_speed_sgd_normbyf5minave_allmean_pairedwilcox_res <- replace(dfd_s5min_2995_5990_speed_sgd_normbyf5minave_allmean_pairedwilcox_res,
                                                                                 dfd_s5min_2995_5990_speed_sgd_normbyf5minave_allmean_pairedwilcox_res=="ns",
                                                                                 NA)

##### make plot #####
g_d_f10min_speed_gd_onlymix_normbyf5minave_allmean <- 
  ggplot(dfd_f10min_speed_sgd_strain %>%
           filter(n_inds == "Group",
                  type == "Mixed",
                  !str_detect(strain, "norpA")) %>%
           group_by(alpha, seconds_total) %>%
           dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
           left_join(dfd_f5min_speed_sgd_ave_strain %>%
                       filter(n_inds == "Group",
                              type == "Mixed",
                              !str_detect(strain, "norpA")) %>%
                       group_by(alpha) %>%
                       dplyr::summarize(speed_f5min_ave = mean(speed_f5min_ave, na.rm = TRUE)), 
                     by = "alpha") %>%
           mutate(speed_normbyf5minave = speed / speed_f5min_ave),
         aes(x = seconds_total, 
             y = speed_normbyf5minave, 
             alpha = alpha,
             group = alpha)) +
  geom_line(color = "#874c70") +
  geom_vline(xintercept=300, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=315, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=330, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=345, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=360, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=375, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=390, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=405, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=420, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=435, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=450, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=465, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=480, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=495, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=510, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=525, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=540, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=555, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=570, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=585, linetype="dotted", alpha=.6, linewidth=.4) +  
  xlab("Time (min)") +
  ylab("Moving speed relative to \n the average during the initial 5 min") +
  scale_x_continuous(breaks = seq(0, 600, 60), labels = seq(0, 10, 1)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.22, 0.2),
        legend.key = element_blank(),
        legend.background = element_blank())
g_d_f10min_speed_gd_onlymix_normbyf5minave_allmean


g_d_s5min_stim_speed_gd_2strains_normbyf5minave_allmean <- 
  ggplot(dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain %>%
           filter(n_inds == "Group", !str_detect(strain, "norpA")) %>%
           group_by(stim_time, var, type, mixed, alpha) %>%
           dplyr::summarize(speed_norm_mean = mean(speed_norm_mean, na.rm = TRUE)) %>%
           filter(type == "Mixed"),
         aes(x = stim_time, 
             y = speed_norm_mean, 
             shape = mixed,
             alpha = alpha,
             group = var)) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = -Inf, ymax = Inf, alpha = 0.4) +
  geom_point(size = 3, color = "#874c70") +
  geom_line(color = "#874c70") +
  annotate(geom = "text",
           x = dfd_s5min_2995_5990_speed_sgd_normbyf5minave_allmean_pairedwilcox_res %>%
             filter(group2 == "Group_Mixed") %>%
             pull(stim_time) + 0.25,
           y = 1.01, 
           label = dfd_s5min_2995_5990_speed_sgd_normbyf5minave_allmean_pairedwilcox_res %>%
             filter(group2 == "Group_Mixed") %>%
             pull(p.adj.signif),
           angle = 90) +
  scale_x_continuous(breaks = seq(0, 15, 3)) +
  scale_shape_manual(values = c(17, 16), guide = "none") +
  scale_alpha_manual(values = c(0.4, 1)) +
  coord_cartesian(xlim = c(0, 15)) +
  xlab("Time after stimulus (s)") +
  ylab("Moving speed relative to \n the average during the initial 5 min") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.7, 0.2),
        legend.key = element_blank(),
        legend.background = element_blank())
g_d_s5min_stim_speed_gd_2strains_normbyf5minave_allmean


g_fig_gd_speed <- g_d_f10min_speed_gd_onlymix_normbyf5minave_allmean +
  g_d_s5min_stim_speed_gd_2strains_normbyf5minave_allmean +
  patchwork::plot_layout(ncol = 2,
                         width = c(0.5, 0.3))
g_fig_gd_speed
ggsave("../figures/Figure4a.pdf", g_fig_gd_speed, w=8, h=3)



#### Figure 4c ####
##### load dataset #####
dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial <- 
  read_parquet("../data/2_mixed_strain/dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial.parquet") %>%
  ungroup()

for_sec <- 3
a2 <- for_sec / 0.5
a <- (29 - a2) / a2

dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging <- 
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial %>%
  filter(n_inds == "Group", type == "Mixed", stim_time >= 0) %>%
  mutate(foraging_success = if_else(stim_time < for_sec, speed * (-a), speed)) %>%
  group_by(strain, sex, n_inds, trial, var, type, mixed, alpha) %>%
  dplyr::summarize(foraging_success = sum(foraging_success))

dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging2 <- 
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial %>%
  filter(n_inds == "Group", stim_time >= 0) %>%
  mutate(foraging_success = if_else(stim_time < for_sec, speed * (-a), speed)) %>%
  group_by(strain, sex, n_inds, trial, var, type, mixed, alpha) %>%
  dplyr::summarize(foraging_success = sum(foraging_success))

dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging4 <- 
  dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging2 %>%
  filter(n_inds == "Group", type == "Mixed") %>%
  group_by(strain, sex, n_inds, var, type, mixed, alpha) %>%
  dplyr::summarize(foraging_success = mean(foraging_success, na.rm = TRUE)) %>%
  ungroup()

##### make plot #####
g_d_s5min_stim_speed_raw_foraging <-  
  dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging4 %>%
  filter(!str_detect(strain, "norpA")) %>%
  ggplot(aes(x = alpha, y = foraging_success, 
             col = alpha)) +
  geom_boxplot() + 
  geom_path(aes(group = strain), color = "gray", linewidth = 0.4) +
  geom_point(aes(shape = mixed), show.legend = F, size = 2) + 
  ggpubr::stat_compare_means(paired = TRUE) +
  scale_color_manual(values = c("#9e8896", "#874c70")) +
  scale_shape_manual(values = c(17, 16)) +
  ylab("Virtual fitness") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        legend.position = "none")
g_d_s5min_stim_speed_raw_foraging
ggsave("../figures/Figure4c.pdf", g_d_s5min_stim_speed_raw_foraging, w=2, h=3)


#### Figure 4d ####
##### load dataset #####
dfd_s5min_2995_5990_freezing_stopduration_group <-
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial %>% #dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial #dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial
  filter(stim_time != -0.5, n_inds == "Group", mixed == "Single strain", type != "Mixed", speed_norm_mean > 1) %>% #speed_norm_mean #speed_normbystimminus05
  group_by(strain, trial, sex, n_inds, var, type) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(strain, trial, sex, n_inds),
              names_from = type, values_from = stim_time,
              values_fill = 14) %>%
  group_by(strain) %>%
  dplyr::summarize(Strain1 = mean(Strain1, na.rm = TRUE),
                   Strain2 = mean(Strain2, na.rm = TRUE)) %>%
  mutate(diff_stop_duration = abs(Strain1 - Strain2))

##### make plot #####
g_d_s5min_stim_speed_normbyf5minave_foraging_freezing_stopduration_diff <- 
  ggplot(
    dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging_tr %>%
      inner_join(dfd_s5min_2995_5990_freezing_stopduration_group, by = "strain"), 
    aes(x = diff_stop_duration, y = overyield)) +
  geom_point() +
  stat_smooth(linewidth = 2, color= "grey", method = "lm") +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(
                          after_stat(rr.label),
                          stat(p.value.label),
                          sep = "~~~")),
                        label.x = "left",
                        label.y = "top",
                        parse = TRUE, size=4) +
  xlab("Difference in freezing duration") +
  ylab("Diversity effect on fitness") +
  theme_bw()
g_d_s5min_stim_speed_normbyf5minave_foraging_freezing_stopduration_diff
ggsave("../figures/Figure4d.pdf", g_d_s5min_stim_speed_normbyf5minave_foraging_freezing_stopduration_diff, w=2.5, h=3)


#### Figure 4e ####
##### load dataset #####
df_go <- list.files(paste0("../data/2_mixed_strain/gwas/enrichment/"),  #clusterProfiler_1kbp_foraging_success5_2000_regress_pi_mean_plus
                     pattern = "*.tsv", full.names = TRUE, recursive=T) %>% 
  map_df(~ data.table::fread(.) %>% mutate(Description = as.character(Description),
                                           ID = as.character(ID),
                                           GeneRatio = as.character(GeneRatio),
                                           BgRatio = as.character(BgRatio),
                                           geneID = as.character(geneID),
                                           Type = as.character(Type))) %>%
  mutate(Description = str_to_sentence(Description),
         GeneRatio = str_split(GeneRatio, "/") %>% 
           map_chr(1) %>% 
           as.numeric() / str_split(GeneRatio, "/") %>% 
           map_chr(2) %>% 
           as.numeric()) %>%
  group_by(Type) %>%
  arrange(qvalue) %>%
  slice_head(n = 10) %>%
  ungroup()

##### make plot #####
g_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- 
  ggplot(df_go, 
         aes(y = reorder(Description, GeneRatio), 
             x = GeneRatio, col = qvalue)) +
  geom_point(aes(shape = Type, size = Count)) +
  tidytext::scale_y_reordered() +
  scale_color_viridis_c(name = expression(paste(italic(q), "-value"))) +
  scale_shape_manual(values = c(17, 16, 15), labels = c("Biological process", "Cellular component", "Molecular Function")) +
  xlab("Gene ratio") +
  ylab("GO term") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank())

g_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus
ggsave("../figures/Figure4e.pdf", g_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, w=6, h=3)


#### Figure 4f ####
##### load dataset #####
df_gwas_test <- read.delim("../data/1_single_strain/gwas/selection/df_gwas_test.tsv")
df_DE_test <- read.delim("../data/2_mixed_strain/gwas/selection/df_DE_test.tsv")

##### make plot #####
g_gwas <- 
  ggplot(df_gwas_test %>%
           transform(data = factor(data, levels = c("Not-associated", "Associated"))),
         aes(x = data, y = TajimaD, color = data)) +
  geom_boxplot(outlier.color = NA) + 
  geom_jitter(size = 2, height = 0, width = 0.2, alpha = 0.2, shape = 16) + #show.legend = F) + 
  ggpubr::stat_compare_means() +
  scale_color_manual(values = c("#5a5359", "#a25768")) +
  coord_cartesian(ylim = c(-0.25, 0.1)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
g_gwas


g_ghas <- 
  ggplot(df_DE_test %>%
           filter(cor == "plus") %>%
           transform(data = factor(data, levels = c("others", "DE involved"))),
         aes(x = data, y = TajimaD, col = data)) +
  geom_boxplot(outlier.color = NA) + 
  geom_jitter(size = 2, height = 0, width = 0.2, alpha = 0.2, shape = 16) + #show.legend = F) + 
  ggpubr::stat_compare_means() +
  scale_color_manual(values = c("#5a5359", "#874c70")) +
  coord_cartesian(ylim = c(-0.25, 0.1)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
g_ghas

g_selection <- g_gwas + g_ghas 

ggsave("../figures/Figure4f.pdf", g_selection, w = 3, h = 3)

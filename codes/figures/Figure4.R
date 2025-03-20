#### load packages ####
targetPackages <- c('tidyverse','gtools','ggpubr','arrow','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### Figure 4a ####
##### load dataset #####
dfd_f10min_speed_sgd_strain <- 
  read_parquet("../data/2_mixed_strain/parquet/dfd_f10min_speed_sgd_strain.parquet") %>%
  ungroup() %>%
  filter(!str_detect(strain, "norpA"))

dfd_f5min_speed_sgd_ave_strain <- 
  dfd_f10min_speed_sgd_strain %>%
  filter(seconds_total < 300) %>%
  group_by(sex, strain, n_inds, var, type, mixed, alpha) %>%
  dplyr::summarize(speed_f5min_ave = mean(speed, na.rm = T))

dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain <- read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain.parquet") %>%
  ungroup() %>%
  filter(!str_detect(strain, "norpA"))

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
        legend.background = element_blank(),
        axis.text = element_text(color = "black"))
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
        legend.background = element_blank(),
        axis.text = element_text(color = "black"))
g_d_s5min_stim_speed_gd_2strains_normbyf5minave_allmean


g_fig_gd_speed <- g_d_f10min_speed_gd_onlymix_normbyf5minave_allmean +
  g_d_s5min_stim_speed_gd_2strains_normbyf5minave_allmean +
  patchwork::plot_layout(ncol = 2,
                         width = c(0.5, 0.3))
g_fig_gd_speed
ggsave("../figures/Figure4a.pdf", g_fig_gd_speed, w=8, h=3)


#### Figure 4b ####
##### load dataset #####
dfd_fig4 <- read_delim("../data/2_mixed_strain/dfd_fig4.tsv")

##### stat #####
df_for_stat <- dfd_fig4 %>% 
  filter(thr_sec == 3) %>%
  pivot_wider(id_cols = strain, names_from = "alpha", values_from = performance)
stats::wilcox.test(Pair(df_for_stat$Expected, df_for_stat$Observed) ~ 1,
                   data = df_for_stat)

##### make plot #####
g_d_s5min_stim_speed_normbyf5minave_performance_sec3 <-
  dfd_fig4 %>%
  filter(thr_sec == 3) %>%
  ggplot(aes(x = alpha, y = performance,
             col = alpha)) +#col = t.test.p.value)) +
  geom_boxplot() +
  geom_path(aes(group = strain), color = "gray", linewidth = 0.4) +
  geom_point(aes(shape = alpha), show.legend = F, size = 2) +
  ggpubr::stat_compare_means(paired = TRUE, size = 1.8) +
  scale_color_manual(values = c("#9e8896", "#874c70")) +#values = viridis(3)[1:2]) +
  scale_shape_manual(values = c("Expected" = 17, "Observed" = 16)) +#values = viridis(3)[1:2]) +
  ylab("Behavioral performance") +
  # facet_wrap( ~ thr_sec, nrow = 2) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none")
g_d_s5min_stim_speed_normbyf5minave_performance_sec3
ggsave("../figures/Figure4b.pdf", g_d_s5min_stim_speed_normbyf5minave_performance_sec3, w=2, h=3)


#### Figure 4c ####
##### load dataset #####
dfd_s5min_2995_5990_gd_freezing_duration_dist <- 
  read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_gd_freezing_duration_dist.parquet") %>%
  ungroup()

##### make plot #####
g_d_s5min_stim_speed_normbyf5minave_foraging_freezing_stopduration_diff <- 
  ggplot(
    dfd_s5min_2995_5990_gd_freezing_duration_dist %>%
      filter(thr_sec == "Threshold: 3 s"),
    aes(x = freezing_duration_dist, y = performance)) +
  geom_point() +
  stat_smooth(linewidth = 2, color= "grey", method = "lm") +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(
                          after_stat(rr.label),
                          stat(p.value.label),
                          sep = "~~~")),
                        label.x = "left",
                        label.y = "top",
                        parse = TRUE, size = 4) +
  xlab("Difference in freezing duration (s)") +
  ylab("Behavioral performance") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"))
g_d_s5min_stim_speed_normbyf5minave_foraging_freezing_stopduration_diff
ggsave("../figures/Figure4c.pdf", g_d_s5min_stim_speed_normbyf5minave_foraging_freezing_stopduration_diff, w=2, h=3)



#### Figure 4d ####
##### load dataset #####
df_spider_group_attack_dist <- read_parquet("../data/6_spiderACI/Experiment2_group/parquet/df_spider_group_attack_dist.parquet") %>%
  group_by(id, trial, id_trial, order, stop_duration, social_influence) %>%
  dplyr::summarize(num_attack = mean(num_attack, na.rm = TRUE),
                   Fly_distance = mean(Fly_distance, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(stop_duration = case_when(stop_duration == "1-1-1-1-1-1s" ~ "1 s",
                                          stop_duration == "1-1-1-3-3-3s" ~ "1 / 3 s",
                                          stop_duration == "1-1-1-6-6-6s" ~ "1 / 6 s",
                                          stop_duration == "1-1-1-12-12-12s" ~ "1 / 12 s",
                                          stop_duration == "3-3-3-3-3-3s" ~ "3 s",
                                          stop_duration == "3-3-3-6-6-6s" ~ "3 / 6 s",
                                          stop_duration == "3-3-3-12-12-12s" ~ "3 / 12 s",
                                          stop_duration == "6-6-6-6-6-6s" ~ "6 s",
                                          stop_duration == "6-6-6-12-12-12s" ~ "6 / 12 s",
                                          stop_duration == "12-12-12-12-12-12s" ~ "12 s")) %>%
  transform(stop_duration = factor(stop_duration, levels = c("1 s", "3 s", "6 s", "12 s",
                                                             "1 / 3 s", "1 / 6 s", "1 / 12 s",
                                                             "3 / 6 s", "3 / 12 s", "6 / 12 s")))

# PCA
spider_group_pca_result <- prcomp(df_spider_group_attack_dist %>% 
                                    dplyr::select(Fly_distance, num_attack), scale = TRUE)
# 各PCの説明力（分散説明率）を表示
summary(spider_group_pca_result)
# 負荷量（loadings_spider_group）を取得してデータフレームに変換
loadings_spider_group <- as.data.frame(spider_group_pca_result$rotation[, 1:2])  # PC1, PC2の負荷量
loadings_spider_group$variables <- rownames(loadings_spider_group)  # 変数名を列として追加

spider_group_pca_data <- as.data.frame(spider_group_pca_result$x) %>%
  bind_cols(df_spider_group_attack_dist)

N = length(unique(spider_group_pca_data$id_trial))
n = length(unique(spider_group_pca_data$id))


spider_group_pca_data_long <-
  spider_group_pca_data %>%
  pivot_longer(cols = c(PC1, PC2), names_to = "var", values_to = "value")

spider_group_pca_data_long_corrected <- 
  spider_group_pca_data_long %>%
  group_by(id_trial, var) %>%
  transmute(order, id_trial, var, corrected_value = value - lm(value ~ order)$coefficients["order"] * order) %>%
  right_join(spider_group_pca_data_long) %>%
  dplyr::select(colnames(spider_group_pca_data_long), everything())

s_list <- c(1, 3, 6, 12)
spider_group_pca_data_long_corrected_mixed <- data.frame()
for (i in 1:3){
  # i <- 1
  for (j in (i+1):4){
    # if (i == 1 & j == 2){ next }
    # if (i == 1 & j == 5){ next }
    # j <- 4
    print(paste(i, j, sep = "_"))
    
    stop_duration_tmp1 <- paste0(s_list[i], " s")
    stop_duration_tmp2 <- paste0(s_list[j], " s")
    stop_duration_tmp <- paste0(s_list[i], " / ", s_list[j], " s")
    
    spider_group_pca_data_long_corrected_tmp <-
      spider_group_pca_data_long_corrected %>%
      filter(stop_duration %in% c(stop_duration_tmp1, stop_duration_tmp2)) %>%
      pivot_wider(id_cols = c(social_influence, id, trial, id_trial, var), 
                  names_from = stop_duration, values_from = corrected_value) %>%
      na.omit() %>%
      dplyr::rename(value_strain1 = all_of(stop_duration_tmp1),
                    value_strain2 = all_of(stop_duration_tmp2)) %>%
      dplyr::mutate(value_strainmean = (value_strain1 + value_strain2) / 2) %>%
      
      left_join(spider_group_pca_data_long_corrected %>%
                  filter(stop_duration == stop_duration_tmp)) %>%
      dplyr::mutate(value_diff_mean = corrected_value - value_strainmean)
    
    spider_group_pca_data_long_corrected_mixed <-
      bind_rows(spider_group_pca_data_long_corrected_mixed,
                spider_group_pca_data_long_corrected_tmp)
    
  }
}

spider_group_pca_data_long_corrected_mixed_plot <-
  spider_group_pca_data_long_corrected_mixed %>%
  filter(var == "PC1") %>%
  pivot_longer(cols = c("value_strainmean", "corrected_value"), 
               names_to = "type", values_to = "score") %>%
  dplyr::mutate(type = if_else(type == "value_strainmean", "Expected", "Observed")) %>%
  arrange(id, trial, order)

##### stat #####
df_for_stat_no <- spider_group_pca_data_long_corrected_mixed_plot %>% 
  filter(social_influence == "No interaction") %>%
  dplyr::select(c(id_trial, stop_duration, type, score)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(id_trial, stop_duration), 
              names_from = type, 
              values_from = score)
stats::wilcox.test(Pair(df_for_stat_no$Expected, df_for_stat_no$Observed) ~ 1,
                   data = df_for_stat_no)

df_for_stat_yes <- spider_group_pca_data_long_corrected_mixed_plot %>% 
  filter(social_influence == "Social interaction") %>%
  dplyr::select(c(id_trial, stop_duration, type, score)) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(id_trial, stop_duration), 
              names_from = type, 
              values_from = score)
stats::wilcox.test(Pair(df_for_stat_yes$Expected, df_for_stat_yes$Observed) ~ 1,
                   data = df_for_stat_yes)

##### make plot #####
col_list_all <-
  c(pals::tol.rainbow(5)[2:5], 
    pals::tol.rainbow(length(unique(spider_group_pca_data$stop_duration))+3)[6:10],
    pals::tol.rainbow(length(unique(spider_group_pca_data$stop_duration))+3)[12])

g_spider_group_pca <-
  ggplot(spider_group_pca_data, aes(x = PC1, y = PC2, color = stop_duration)) +
  stat_ellipse(type = "norm", level = 0.68, linetype = "solid", linewidth = 0.6, alpha = 0.6) +
  geom_point(aes(shape = social_influence), alpha = 0.6) + #, alpha = 0.6
  # 変数ベクトルの描画
  geom_segment(data = loadings_spider_group, aes(x = 0, y = 0, xend = PC1, yend = PC2), # スケーリングを適宜調整
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  ggrepel::geom_text_repel(data = loadings_spider_group, aes(x = PC1, y = PC2, label = variables), # スケーリングを適宜調整
                           color = "black", size = 4) +
  annotate("text", x = Inf, y = -2, , check_overlap = TRUE,
           label = substitute(expr = paste(italic("N"), " = ", N_val, " (", italic("n"), " = ", n_val, ")"),
                              env = list(N_val = N, n_val = n)), 
           hjust = 1.1, vjust = 1.5, size = 3, fontface = "plain")+
  # ggrepel::geom_text_repel(aes(label = id_trial), color = "black") +
  scale_color_manual(values = col_list_all) +
  scale_shape_manual(values = c(16, 15)) +
  xlab(paste0("PC1 (", round(summary(spider_group_pca_result)$importance[2,1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(spider_group_pca_result)$importance[2,2] * 100, 1), "%)")) +
  # facet_wrap(~ social_influence) +
  theme_classic() +
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) +
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_blank(),
    legend.position = c(0.8, 0.8),
    legend.spacing.y = unit(0.1, "cm"),   # レジェンド項目間の横スペース
    legend.box.margin = margin(c(0,0,0,0)),
    legend.box.spacing = unit(0, "cm"),
    legend.box = "vertical",
    legend.direction = "vertical",
    legend.key = element_blank(),
    legend.key.size = unit(0.3, "cm"),
    legend.background = element_blank()) +#,
  # legend.margin = margin(c(1,1,0.1,1)),
  # legend.spacing = unit(0, "cm")) +
  guides(color = guide_legend(ncol = 2))  # レジェンドの列数を2列に設定
g_spider_group_pca


col_list <- col_list_all[5:10]

g_spider_group_pca_DE_pair <-
  ggplot(spider_group_pca_data_long_corrected_mixed_plot,
         aes(x = type, y = score)) +
  geom_boxplot(aes(color = type), outlier.shape = NA) +
  geom_path(aes(group = interaction(id_trial, social_influence, stop_duration), 
                color = stop_duration), 
            linewidth = 0.4, alpha = 0.6) +
  geom_point(aes(shape = type, color = type), show.legend = F, size = 2) +
  ggpubr::stat_compare_means(paired = TRUE, size = 1.8) +
  scale_color_manual(values = c(col_list, "#9e8896", "#874c70")) +#values = viridis(3)[1:2]) + "#e6b422", "#478384", "#4d4c61",
  scale_shape_manual(values = c("Expected" = 17, "Observed" = 16)) +#values = viridis(3)[1:2]) +
  ylab("Corrected PC1 score") +
  facet_grid( ~ social_influence) + #sex+
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        strip.background = element_blank(),
        legend.position = "none")
g_spider_group_pca_DE_pair

final_plot <- g_spider_group_pca + g_spider_group_pca_DE_pair + 
  plot_layout(ncol = 2, widths = c(1, 1.2))
final_plot
ggsave("../figures/Figure4d.pdf", final_plot, w = 5, h = 3.0)



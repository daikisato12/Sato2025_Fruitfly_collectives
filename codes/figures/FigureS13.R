#### load packages ####
targetPackages <- c('tidyverse','arrow','patchwork','emmeans')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### Figure S13a-b ####
##### load dataset #####
df_spider_single <- read_parquet("../data/6_spiderACI/Experiment1_single/parquet/df_spider_single.parquet")

df_spider_single_prefix <- df_spider_single %>%
  group_by(prefix, order, stop_duration, stop_threshold, speed_ave_set, id, trial, id_trial) %>% 
  dplyr::summarize()

df_spider_single_freezing <- df_spider_single %>%
  right_join(df_spider_single_prefix) %>%
  dplyr::mutate(Fly_stopped2 = if_else(Fly_stopped == 0, "Moving", "Freezing")) %>%
  group_by(prefix, order, stop_duration, stop_threshold, speed_ave_set, id, trial, id_trial, Fly_stopped2) %>%
  dplyr::summarize(n = n()) %>%
  pivot_wider(id_cols = c(prefix, order, stop_duration, stop_threshold, speed_ave_set, id, trial, id_trial), 
              names_from = Fly_stopped2, values_from = n, values_fill = 0) %>%
  dplyr::mutate(Total = Moving + Freezing,
                Freezing_ratio = Freezing / Total) %>%
  transform(speed_ave_set = factor(speed_ave_set, levels = c("6.0mmps", "12.0mmps", "24.0mmps")),
            stop_duration = factor(stop_duration, levels = c("0.0s", "1.0s", "3.0s", "6.0s", "12.0s")),
            stop_threshold = factor(stop_threshold, levels = c("50px", "100px")),
            id_trial = factor(id_trial, levels = unique(.$id_trial) %>% gtools::mixedsort())) %>%
  arrange(id, trial, order)

df_spider_single_speed <- df_spider_single %>%
  group_by(prefix, order, stop_duration, stop_threshold, speed_ave_set, id, trial, id_trial) %>% 
  dplyr::summarize(Fly_speed = mean(Fly_speed, na.rm = TRUE),
                   Fly_distance = sum(Fly_travelled_dist_diff, na.rm = TRUE),
                   Spider_speed = mean(Spider_speed, na.rm = TRUE)) %>%
  transform(speed_ave_set = factor(speed_ave_set, levels = c("6.0mmps", "12.0mmps", "24.0mmps")),
            stop_duration = factor(stop_duration, levels = c("0.0s", "1.0s", "3.0s", "6.0s", "12.0s")),
            stop_threshold = factor(stop_threshold, levels = c("50px", "100px")),
            id_trial = factor(id_trial, levels = unique(.$id_trial) %>% gtools::mixedsort())) %>%
  arrange(id, trial, order)


df_spider_single_aggression_tmp <- df_spider_single %>%
  right_join(df_spider_single_prefix) %>%
  # dplyr::filter(prefix == "20241030194711_360s_3.0s_50px_24.0mmps_m4_male_1") %>%
  dplyr::mutate(
    match_condition = if_else(
      seconds_total != 0 & Spider_speed > 30 & Distance > lead(Distance, default = Inf), TRUE, FALSE
    ),
    Fly_stopped_1_in_next_3 = if_else(
      match_condition,
      map_lgl(seq_along(Fly_stopped), ~ any(Fly_stopped[(.x + 1):min(.x + 3, n())] == 1, na.rm = TRUE)),
      NA
    ),
    Fly_catched_1_in_next_3 = if_else(
      match_condition,
      map_lgl(seq_along(Distance), ~ any(Distance[(.x + 1):min(.x + 3, n())] < 10, na.rm = TRUE)),
      NA
    ),
    row_id = row_number()
  ) %>%
  dplyr::select(-row_id)

df_spider_single_aggression_num_attack <- df_spider_single_aggression_tmp %>%
  dplyr::filter(match_condition == TRUE) %>%
  group_by(prefix, order, stop_duration, stop_threshold, speed_ave_set, id, trial, id_trial, Fly_catched_1_in_next_3) %>%
  dplyr::summarize(n = n()) %>%
  right_join(df_spider_single_prefix) %>%
  dplyr::mutate(Fly_catched_1_in_next_3 = case_when(Fly_catched_1_in_next_3 == TRUE ~ "Success",
                                                    TRUE ~ "Failed"),
                n = if_else(is.na(n), 0, n)) %>%
  pivot_wider(id_cols = c(prefix, order, stop_duration, stop_threshold, speed_ave_set, id, trial, id_trial), 
              names_from = Fly_catched_1_in_next_3, values_from = n, values_fill = 0) %>%
  dplyr::mutate(Total = Success + Failed) %>%
  pivot_longer(cols = c("Success", "Failed", "Total"), names_to = "type", values_to = "num_attack") %>%
  transform(speed_ave_set = factor(speed_ave_set, levels = c("6.0mmps", "12.0mmps", "24.0mmps")),
            stop_duration = factor(stop_duration, levels = c("0.0s", "1.0s", "3.0s", "6.0s", "12.0s")),
            stop_threshold = factor(stop_threshold, levels = c("50px", "100px")),
            id_trial = factor(id_trial, levels = unique(.$id_trial) %>% gtools::mixedsort())) %>%
  filter(type == "Success")

df_spider_single_all <- left_join(df_spider_single_speed, df_spider_single_freezing) %>%
  left_join(df_spider_single_aggression_num_attack) %>%
  pivot_longer(cols = c("Fly_speed", "Fly_distance", "Freezing_ratio", "num_attack"), 
               names_to = "var", values_to = "value")

##### Figure S13a #####
g_spider_single_var_all <-
  ggplot(df_spider_single_all %>%
           filter(var != "Fly_distance") %>%
           dplyr::mutate(stop_duration = parse_number(as.character(stop_duration)) %>% as.character()) %>%
           transform(stop_duration = factor(stop_duration, levels = unique(.$stop_duration) %>% gtools::mixedsort())), 
         aes(x = stop_duration, y = value, color = stop_duration)) +
  stat_summary(fun.data = "mean_se", size = 0.3, linewidth = 0.3) +
  scale_color_manual(values = pals::tol.rainbow(5)) +
  facet_grid(var ~ stop_threshold+speed_ave_set, scales = "free_y") +
  xlab("Freezing duration (s)") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        strip.background = element_blank())

g_spider_single_var_all
ggsave("../figures/FigureS13a.pdf", g_spider_single_var_all, w = 6, h = 4)

##### Figure S13b #####
num_id <- length(unique(paste0(df_spider_single_prefix$id, df_spider_single_prefix$trial)))
g_spider_single_num_attack_order <-
  ggplot(df_spider_single_aggression_num_attack %>%
           filter(type == "Success"),
         aes(x = order, y = num_attack)) +
  stat_smooth(color = "grey", se = FALSE, method = "lm") +
  geom_point(aes(color = stop_duration, shape = speed_ave_set, alpha = stop_threshold), size = 2) +
  geom_text(
    data = df_spider_single_aggression_num_attack %>%
      group_by(id_trial) %>%
      dplyr::summarize(x = 30, y = 9, .groups = "drop"),
    aes(x = x, y = y, label = id_trial),
    inherit.aes = FALSE,
    hjust = 1, vjust = 0, size = 3.5
  ) +
  scale_color_manual(values = pals::tol.rainbow(5)) +
  scale_alpha_manual(values = c(0.4, 0.8)) +
  scale_shape_manual(values = c(16, 17, 15)) +
  facet_wrap(~ id_trial, ncol = 5, strip.position = "top") +
  xlab("Order of experiment") +
  ylab("Number of attacks") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        strip.text = element_blank(),
        strip.background = element_blank())

g_spider_single_num_attack_order
ggsave("../figures/FigureS13b.pdf", g_spider_single_num_attack_order, w = 8, h = 4)


#### Figure S13c ####
##### load dataset #####
df_spider_single_attack_dist <- read_parquet("../data/6_spiderACI/Experiment1_single/parquet/df_spider_single_attack_dist.parquet") %>%
  dplyr::mutate(stop_duration = paste0(parse_number(as.character(stop_duration)), " s")) %>%
  transform(stop_duration = factor(stop_duration, levels = c("0 s", "1 s", "3 s", "6 s", "12 s")))

# PCA
spider_single_pca_result <- prcomp(df_spider_single_attack_dist %>% 
                                     dplyr::select(Fly_distance, num_attack), scale = TRUE)
# 各PCの説明力（分散説明率）を表示
summary(spider_single_pca_result)
# 負荷量（loadings_spider_single）を取得してデータフレームに変換
loadings_spider_single <- as.data.frame(spider_single_pca_result$rotation[, 1:2])  # PC1, PC2の負荷量
loadings_spider_single$variables <- rownames(loadings_spider_single)  # 変数名を列として追加

spider_single_pca_data <- as.data.frame(spider_single_pca_result$x) %>%
  bind_cols(df_spider_single_attack_dist)

N = length(unique(spider_single_pca_data$id_trial))
n = length(unique(spider_single_pca_data$id))

##### make plot #####
g_spider_single_pca <-
  ggplot(spider_single_pca_data, aes(x = PC1, y = PC2, color = stop_duration)) +
  stat_ellipse(type = "norm", level = 0.68, linetype = "solid", linewidth = 0.6) +
  geom_point(shape = 16) + #, alpha = 0.6
  # 変数ベクトルの描画
  geom_segment(data = loadings_spider_single, aes(x = 0, y = 0, xend = PC1, yend = PC2), # スケーリングを適宜調整
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  ggrepel::geom_text_repel(data = loadings_spider_single, aes(x = PC1, y = PC2, label = variables), # スケーリングを適宜調整
                           color = "black", size = 4) +
  annotate("text", x = Inf, y = Inf, check_overlap = TRUE,
           label = substitute(expr = paste(italic("N"), " = ", N_val, " (", italic("n"), " = ", n_val, ")"),
                              env = list(N_val = N, n_val = n)), 
           hjust = 1.1, vjust = 1.5, size = 3, fontface = "plain")+
  # ggrepel::geom_text_repel(aes(label = id_trial), color = "black") +
  scale_color_manual(values = pals::tol.rainbow(5)) +
  xlab(paste0("PC1 (", round(summary(spider_single_pca_result)$importance[2,1] * 100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(spider_single_pca_result)$importance[2,2] * 100, 1), "%)")) +
  theme_classic() +
  theme(#legend.position = "bottom",
    axis.text = element_text(color = "black"),
    axis.ticks = element_line(color = "black"),
    legend.title = element_blank(),
    legend.position = c(0.85, 0.2),
    legend.key = element_blank(),
    legend.key.size = unit(0.4, "cm"),
    legend.background = element_blank()) +#,
  # legend.margin = margin(c(1,1,0.1,1)),
  # legend.spacing = unit(0, "cm")) +
  guides(color = guide_legend(ncol = 1))  # レジェンドの列数を2列に設定
g_spider_single_pca

g_spider_single_pc2 <- ggplot(spider_single_pca_data %>%
                                dplyr::mutate(stop_duration = as.character(parse_number(as.character(stop_duration)))) %>%
                                transform(stop_duration = factor(stop_duration, levels = c("0", "1", "3", "6", "12"))), 
                              aes(x = stop_duration, 
                                  y = PC2, 
                                  color = stop_duration)) +
  # geom_boxplot() +
  # geom_smooth((aes(group = 1)), color = "grey", method = "lm", formula=y ~ poly(x, 2, raw=TRUE)) +
  stat_summary(aes(y = PC2, group = 1), 
               fun = mean, geom = "line", linewidth = 1, color = "grey") +
  geom_point(shape = 16, alpha = .4) +
  stat_summary(fun.data = "mean_se", size = 0.6) +
  scale_color_manual(values = pals::tol.rainbow(5)) +
  xlab("Freezing duration (s)") +
  ylab("PC2 score") +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.ticks = element_line(color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.spacing = unit(1, "lines"))
g_spider_single_pc2

final_plot <- g_spider_single_pca + g_spider_single_pc2 + 
  plot_layout(ncol = 2, widths = c(1, 0.6))
final_plot
ggsave("../figures/FigureS13c.pdf", final_plot, w = 4.5, h = 2.8)


##### stat #####
car::Anova(lme4::lmer(PC2 ~ stop_duration + (1|id) + (1|trial) + (1|order), 
                      data = spider_single_pca_data))
model_lmer <- lme4::lmer(PC2 ~ stop_duration + (1|id) + (1|trial) + (1|order), 
                         data = spider_single_pca_data)
# emmeansオブジェクトを作成
emmeans_results <- emmeans::emmeans(model_lmer, ~ stop_duration)
# plot(emmeans_results, comparison=T)
# ペアワイズ比較を行い、多重比較の調整をTukey法で実施
pairwise_results <- pairs(emmeans_results, adjust = "tukey")
emmeans::pwpp(emmeans_results, sort = FALSE) + theme_bw()
print(pairwise_results)



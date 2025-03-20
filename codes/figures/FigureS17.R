#### load packages ####
targetPackages <- c('tidyverse','arrow','patchwork','gtools')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### load dataset ####
df_spider_group <- read_parquet("../data/6_spiderACI/Experiment2_group/parquet/df_spider_group.parquet")
order_id_trial <- unique(df_spider_group$id_trial) %>% gtools::mixedsort()

df_spider_group_prefix <- df_spider_group %>%
  group_by(prefix, order, stop_duration, social_influence, id, trial, id_trial) %>% 
  dplyr::summarize()

df_spider_group_freezing <- df_spider_group %>%
  right_join(df_spider_group_prefix) %>%
  dplyr::mutate(Fly_stopped2 = if_else(Fly_stopped == 0, "Moving", "Freezing")) %>%
  group_by(prefix, order, stop_duration, social_influence, id, trial, id_trial, Fly_stopped2) %>%
  dplyr::summarize(n = n()) %>%
  pivot_wider(id_cols = c(prefix, order, stop_duration, social_influence, id, trial, id_trial), 
              names_from = Fly_stopped2, values_from = n, values_fill = 0) %>%
  dplyr::mutate(Total = Moving + Freezing,
                Freezing_ratio = Freezing / Total) %>%
  arrange(id, trial, order)

df_spider_group_speed <- df_spider_group %>%
  group_by(prefix, order, stop_duration, social_influence, id, trial, id_trial) %>% 
  dplyr::summarize(Fly_speed = mean(Fly_speed, na.rm = TRUE),
                   Fly_distance = sum(Fly_travelled_dist_diff, na.rm = TRUE),
                   Spider_speed = mean(Spider_speed, na.rm = TRUE)) %>%
  arrange(id, trial, order)

df_spider_group_aggression_tmp <- df_spider_group %>%
  right_join(df_spider_group_prefix) %>%
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
  dplyr::select(-row_id) # 不要な列を削除

df_spider_group_aggression_num_attack <- df_spider_group_aggression_tmp %>%
  dplyr::filter(match_condition == TRUE) %>%
  group_by(prefix, order, stop_duration, social_influence, id, trial, id_trial, Fly_catched_1_in_next_3) %>%
  dplyr::summarize(n = n()) %>%
  right_join(df_spider_group_prefix) %>%
  dplyr::mutate(Fly_catched_1_in_next_3 = case_when(Fly_catched_1_in_next_3 == TRUE ~ "Success",
                                                    TRUE ~ "Failed"),
                n = if_else(is.na(n), 0, n)) %>%
  pivot_wider(id_cols = c(prefix, order, stop_duration, social_influence, id, trial, id_trial), 
              names_from = Fly_catched_1_in_next_3, values_from = n, values_fill = 0) %>%
  dplyr::mutate(Total = Success + Failed) %>%
  pivot_longer(cols = c("Success", "Failed", "Total"), names_to = "type", values_to = "num_attack") %>%
  filter(type == "Success") %>%
  left_join(df_spider_group_speed) %>%
  left_join(df_spider_group_freezing) %>%
  dplyr::mutate(social_influence = case_when(social_influence == "no" ~ "Unaligned speed", 
                                             social_influence == "social" ~ "Aligned speed"),
                stop_duration = case_when(stop_duration == "1-1-1-1-1-1s" ~ "1",
                                          stop_duration == "1-1-1-3-3-3s" ~ "1 / 3",
                                          stop_duration == "1-1-1-6-6-6s" ~ "1 / 6",
                                          stop_duration == "1-1-1-12-12-12s" ~ "1 / 12",
                                          stop_duration == "3-3-3-3-3-3s" ~ "3",
                                          stop_duration == "3-3-3-6-6-6s" ~ "3 / 6",
                                          stop_duration == "3-3-3-12-12-12s" ~ "3 / 12",
                                          stop_duration == "6-6-6-6-6-6s" ~ "6",
                                          stop_duration == "6-6-6-12-12-12s" ~ "6 / 12",
                                          stop_duration == "12-12-12-12-12-12s" ~ "12"),
                fill = case_when(str_detect(stop_duration, "/") & stop_duration == "1 / 3" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[6],
                                 str_detect(stop_duration, "/") & stop_duration == "1 / 6" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[7],
                                 str_detect(stop_duration, "/") & stop_duration == "1 / 12" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[8],
                                 str_detect(stop_duration, "/") & stop_duration == "3 / 6" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[9],
                                 str_detect(stop_duration, "/") & stop_duration == "3 / 12" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[10],
                                 str_detect(stop_duration, "/") & stop_duration == "6 / 12" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[12],
                                 TRUE ~ "white"
                )) %>%
  transform(stop_duration = factor(stop_duration, 
                                   levels = c("1", "1 / 3", "1 / 6", "1 / 12",
                                              "3", "3 / 6", "3 / 12", 
                                              "6", "6 / 12", "12")),
            social_influence = factor(social_influence, levels = c("Unaligned speed", 
                                                                   "Aligned speed")),
            id_trial = factor(id_trial, levels = order_id_trial)) %>%
  arrange(id_trial)

#### Figure S17a ####
num_id <- length(unique(paste0(df_spider_group_prefix$id, df_spider_group_prefix$trial)))
g_spider_group_num_attack_order <-
  ggplot(df_spider_group_aggression_num_attack %>%
           filter(type == "Success"),
         aes(x = order, y = num_attack)) +
  stat_smooth(color = "grey", se = FALSE, method = "lm") +
  geom_point(aes(color = stop_duration, shape = social_influence), alpha = 0.8, size = 2) +
  geom_text(
    data = df_spider_group_aggression_num_attack %>%
      group_by(id_trial) %>%
      dplyr::summarize(x = 20, y = 11, .groups = "drop"),
    aes(x = x, y = y, label = id_trial),
    inherit.aes = FALSE,
    hjust = 1, vjust = 0, size = 2.5
  ) +
  scale_color_manual(values = col_list_group_all) +
  scale_shape_manual(values = c("Unaligned speed" = 16,
                                "Aligned speed" = 17)) +
  facet_wrap(~ id_trial, ncol = 5, strip.position = "top") +
  xlab("Order of experiment") +
  ylab("Number of attacks") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        strip.text = element_blank(),
        strip.background = element_blank())

g_spider_group_num_attack_order
ggsave("../figures/FigureS17a.pdf", g_spider_group_num_attack_order, w = 6, h = 3)


#### Figure S17bc ####
##### load dataset #####
df_spider_group_attack_dist <- read_parquet("../data/6_spiderACI/Experiment2_group/parquet/df_spider_group_attack_dist.parquet") %>%
  group_by(id, trial, id_trial, order, stop_duration, social_influence) %>%
  dplyr::summarize(num_attack = mean(num_attack, na.rm = TRUE),
                   Fly_distance = mean(Fly_distance, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(social_influence = case_when(social_influence == "No interaction" ~ "Unaligned speed", 
                                             social_influence == "Social interaction" ~ "Aligned speed"),
                stop_duration = case_when(stop_duration == "1-1-1-1-1-1s" ~ "1",
                                          stop_duration == "1-1-1-3-3-3s" ~ "1 / 3",
                                          stop_duration == "1-1-1-6-6-6s" ~ "1 / 6",
                                          stop_duration == "1-1-1-12-12-12s" ~ "1 / 12",
                                          stop_duration == "3-3-3-3-3-3s" ~ "3",
                                          stop_duration == "3-3-3-6-6-6s" ~ "3 / 6",
                                          stop_duration == "3-3-3-12-12-12s" ~ "3 / 12",
                                          stop_duration == "6-6-6-6-6-6s" ~ "6",
                                          stop_duration == "6-6-6-12-12-12s" ~ "6 / 12",
                                          stop_duration == "12-12-12-12-12-12s" ~ "12"),
                fill = case_when(str_detect(stop_duration, "/") & stop_duration == "1 / 3" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[6],
                                 str_detect(stop_duration, "/") & stop_duration == "1 / 6" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[7],
                                 str_detect(stop_duration, "/") & stop_duration == "1 / 12" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[8],
                                 str_detect(stop_duration, "/") & stop_duration == "3 / 6" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[9],
                                 str_detect(stop_duration, "/") & stop_duration == "3 / 12" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[10],
                                 str_detect(stop_duration, "/") & stop_duration == "6 / 12" ~ pals::tol.rainbow(length(unique(.$stop_duration))+3)[12],
                                 TRUE ~ "white"
                )) %>%
  transform(stop_duration = factor(stop_duration, 
                                   levels = c("1", "1 / 3", "1 / 6", "1 / 12",
                                              "3", "3 / 6", "3 / 12", 
                                              "6", "6 / 12", "12")),
            social_influence = factor(social_influence, levels = c("Unaligned speed", 
                                                                   "Aligned speed")),
            id_trial = factor(id_trial, levels = order_id_trial))

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

N_spider_group = length(unique(spider_group_pca_data$id_trial))
n_spider_group = length(unique(spider_group_pca_data$id))


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
    
    stop_duration_tmp1 <- s_list[i]
    stop_duration_tmp2 <- s_list[j]
    stop_duration_tmp <- paste0(s_list[i], " / ", s_list[j])
    
    spider_group_pca_data_long_corrected_tmp <-
      spider_group_pca_data_long_corrected %>%
      filter(stop_duration %in% c(stop_duration_tmp1, stop_duration_tmp2)) %>%
      pivot_wider(id_cols = c(social_influence, id, trial, id_trial, var), 
                  names_from = stop_duration, values_from = corrected_value) %>%
      na.omit() %>%
      dplyr::rename(value_strain1 = as.character(stop_duration_tmp1),
                    value_strain2 = as.character(stop_duration_tmp2)) %>%
      dplyr::mutate(value_strainmean = (value_strain1 + value_strain2) / 2) %>%
      
      left_join(spider_group_pca_data_long_corrected %>%
                  filter(stop_duration == stop_duration_tmp)) %>%
      dplyr::mutate(value_diff_mean = corrected_value - value_strainmean)
    
    spider_group_pca_data_long_corrected_mixed <-
      bind_rows(spider_group_pca_data_long_corrected_mixed,
                spider_group_pca_data_long_corrected_tmp)
    
  }
}


##### make plot #####
g_spider_group_pc1 <- ggplot(spider_group_pca_data, 
                              aes(x = stop_duration, 
                                  y = PC1, 
                                  color = stop_duration,
                                  fill = fill)) +
  # geom_boxplot() +
  # geom_smooth((aes(group = 1)), color = "grey", method = "lm", formula=y ~ poly(x, 2, raw=TRUE)) +
  stat_summary(aes(y = PC1, group = 1), 
               fun = mean, geom = "line", linewidth = 1, color = "grey") +
  geom_point(shape = 16, alpha = .4) +
  stat_summary(fun.data = "mean_se", shape = 21, size = 0.6) +
  scale_color_manual(values = col_list_group_all) +
  scale_fill_identity() +
  xlab("Freezing duration (s)") +
  ylab("PC1 Score") +
  facet_wrap( ~ social_influence) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks = element_line(color = "black"),
        # plot.margin = unit(c(0, 0, 0, 0), "cm"),
        # panel.spacing = unit(1, "lines"),
        strip.background = element_blank())

g_spider_group_pc1

g_spider_group_pc1_DE <-
  ggplot(spider_group_pca_data_long_corrected_mixed %>%
           filter(var == "PC1"),
         aes(x = stop_duration, y = value_diff_mean, color = fill, fill = fill)) +
  # geom_histogram(aes(x = Spider_speed), binwidth = 0.2) +
  # coord_cartesian(xlim = c(0, 10))
  # geom_boxplot(outlier.shape = NA) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  stat_summary(fun.data = "mean_se", shape = 21, size = 0.6) +
  # geom_point(aes(color = id_trial), size = 2, alpha = 0.6, shape = 16) +
  # geom_line(aes(group = id_trial, color = id_trial), linewidth = 0.4) +
  # ggpubr::stat_compare_means(paired = TRUE, size = 1.8) +
  scale_color_identity() +
  scale_fill_identity() +
  coord_cartesian(ylim = c(-1, 1)) +
  facet_wrap( ~ social_influence) + 
  xlab("Freezing duration (s)") +
  ylab("Diversity effect in PC1 score") +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1))

g_spider_group_pc1_DE

g_spider_group_pc1_all <- 
  g_spider_group_pc1 + 
  g_spider_group_pc1_DE +
  plot_layout(nrow = 2)
g_spider_group_pc1_all
ggsave("../figures/FigureS17bc.pdf", g_spider_group_pc1_all, w = 6, h = 4.5)

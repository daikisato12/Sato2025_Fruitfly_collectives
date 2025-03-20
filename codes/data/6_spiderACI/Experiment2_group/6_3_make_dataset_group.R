#### load packages ####
targetPackages <- c('tidyverse','data.table','arrow','car','emmeans','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### load dataset and convert to parquet ####
datadir <- "../data/6_spiderACI/Experiment2_group/moddata/tracks" # data directory for tsv files
# ファイル一覧を取得
file_list <- list.files(datadir, pattern = "*.tsv", full.names = TRUE, recursive=T)
# 各ファイルを読み込み、ファイル名の接頭辞を列に追加
df_spider_group_tmp <- file_list %>%
  map_df(~ {
    # ファイル名から接頭辞を抽出（拡張子前までの部分）
    prefix <- tools::file_path_sans_ext(basename(.x)) %>%
      str_remove("_track_modified")
    
    # データを読み込み、新しい列 'filename_prefix' に接頭辞を追加
    fread(.x) %>%
      mutate(prefix = prefix)
  }) %>%
  dplyr::mutate(Spider_X = if_else(Spider_X != 500, Spider_X, NA),
                Spider_Y = if_else(Spider_Y != 0, Spider_Y, NA))

arena_diameter <- 84 # mm
camera_height <- 480 # px

df_spider_group <- df_spider_group_tmp %>%
  dplyr::mutate(Fly_Y = camera_height - Fly_Y,
                Fly_Direction_Radians = -Fly_Direction_Radians,
                Spider_Y = camera_height - Spider_Y) %>%
  group_by(prefix, Fly_ID) %>%
  dplyr::mutate(pos_x_max = max(.$Fly_X, na.rm = T), # 411 px = 84 mm -> 50 px = 10.22 mm, 100 px = 20.44 mm
                pos_x_min = min(.$Fly_X, na.rm = T),
                pos_y_max = max(.$Fly_Y, na.rm = T),
                pos_y_min = min(.$Fly_Y, na.rm = T)) %>%
  dplyr::mutate(Fly_X = ((Fly_X - pos_x_min) / (pos_x_max - pos_x_min)) * arena_diameter, #convert to mm
                Fly_Y = ((Fly_Y - pos_y_min) / (pos_y_max - pos_y_min)) * arena_diameter, #convert to mm
                Fly_X_1framebefore = lag(Fly_X),
                Fly_Y_1framebefore = lag(Fly_Y),
                Fly_travelled_dist_diff = ( (Fly_X - Fly_X_1framebefore)**2 + (Fly_Y - Fly_Y_1framebefore)**2 )**(1/2),
                Spider_X = ((Spider_X - pos_x_min) / (pos_x_max - pos_x_min)) * arena_diameter, #convert to mm
                Spider_Y = ((Spider_Y - pos_y_min) / (pos_y_max - pos_y_min)) * arena_diameter, #convert to mm
                Spider_X_1framebefore = lag(Spider_X),
                Spider_Y_1framebefore = lag(Spider_Y),
                Spider_travelled_dist_diff = ( (Spider_X - Spider_X_1framebefore)**2 + (Spider_Y - Spider_Y_1framebefore)**2)**(1/2),
                Distance = ( (Spider_X - Fly_X)**2 + (Spider_Y - Fly_Y)**2 )**(1/2)) %>%
  ungroup() %>%
  dplyr::mutate(id = str_split(prefix, "_") %>% map_chr(1),
                trial = str_split(prefix, "_") %>% map_chr(2) %>% as.integer(),
                order = str_split(prefix, "_") %>% map_chr(3) %>% as.integer(),
                stop_duration = str_split(prefix, "_") %>% map_chr(4),
                stop_threshold = str_split(prefix, "_") %>% map_chr(5),
                speed_ave_set = str_split(prefix, "_") %>% map_chr(6),
                social_influence = str_split(prefix, "_") %>% map_chr(7),
                id_trial = paste0(id, "_", trial),
                id_num = parse_number(id)
  ) %>%
  dplyr::rename(seconds_total = Time) %>%
  group_by(prefix, stop_duration, stop_threshold, speed_ave_set, social_influence, id_num, id, trial, id_trial, order, Fly_ID, seconds_total = as.integer(seconds_total*4)/4) %>%
  dplyr::summarize(Fly_pos_x = mean(Fly_X, na.rm = TRUE),
                   Fly_pos_y = mean(Fly_Y, na.rm = TRUE),
                   Fly_travelled_dist_diff = sum(Fly_travelled_dist_diff, na.rm = TRUE),
                   Fly_direction = mean(Fly_Direction_Radians, na.rm = TRUE),
                   Fly_speed = sum(Fly_travelled_dist_diff, na.rm = TRUE)/sum(Time_interval, na.rm = TRUE),
                   Fly_stopped = mean(Stopped, na.rm = TRUE),
                   Spider_pos_x = mean(Spider_X, na.rm = TRUE),
                   Spider_pos_y = mean(Spider_Y, na.rm = TRUE),
                   Spider_travelled_dist_diff = sum(Spider_travelled_dist_diff, na.rm = TRUE),
                   Spider_speed = sum(Spider_travelled_dist_diff, na.rm = TRUE)/sum(Time_interval, na.rm = TRUE),
                   Distance = mean(Distance, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(id_num, trial, order, Fly_ID, seconds_total) %>%
  dplyr::select(!id_num)

write_parquet(df_spider_group, "../data/6_spiderACI/Experiment2_group/parquet/df_spider_group.parquet")
# df_spider_group <- read_parquet("../data/6_spiderACI/Experiment2_group/parquet/df_spider_group.parquet")

df_spider_group_prefix <- df_spider_group %>%
  group_by(id, trial, id_trial, order, stop_duration, stop_threshold, speed_ave_set, social_influence, prefix) %>% 
  dplyr::summarize()

num_id <- length(unique(paste0(df_spider_group_prefix$id, df_spider_group_prefix$trial)))

#### travelled distance ####
df_spider_group_speed <- df_spider_group %>%
  group_by(id, trial, id_trial, order, Fly_ID, stop_duration, stop_threshold, speed_ave_set, social_influence, prefix) %>% 
  dplyr::summarize(Fly_speed = mean(Fly_speed, na.rm = TRUE),
                   Fly_distance = sum(Fly_travelled_dist_diff, na.rm = TRUE),
                   Spider_speed = mean(Spider_speed, na.rm = TRUE)) %>%
  dplyr::mutate(Fly_ID = as.character(Fly_ID),
                id_num = parse_number(id),
                social_influence = case_when(social_influence == "no" ~ "No interaction",
                                             social_influence == "weak" ~ "Weak interaction",
                                             social_influence == "strong" ~ "Strong interaction")) %>%
  arrange(id_num, trial, order) %>%
  dplyr::select(!id_num) %>%
  transform(stop_duration = factor(stop_duration,
                                   levels = c("1-1-1-1-1-1s", "1-1-1-3-3-3s", "1-1-1-6-6-6s", "1-1-1-12-12-12s",
                                              "3-3-3-3-3-3s", "3-3-3-6-6-6s", "3-3-3-12-12-12s",
                                              "6-6-6-6-6-6s", "6-6-6-12-12-12s",
                                              "12-12-12-12-12-12s")),
            social_influence = factor(social_influence, levels = c("No interaction", 
                                                                   "Weak interaction",
                                                                   "Strong interaction")),
            id_trial = factor(id_trial, levels = unique(.$id_trial) %>% gtools::mixedsort()))


#### number of attacks ####
df_spider_group_aggression_tmp <- df_spider_group %>%
  group_by(prefix, stop_duration, stop_threshold, speed_ave_set, social_influence, id, trial, id_trial, Fly_ID) %>% 
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
  ungroup() %>%
  dplyr::select(-row_id) # 不要な列を削除

df_spider_group_aggression_num_attack <- df_spider_group_aggression_tmp %>%
  dplyr::filter(match_condition == TRUE) %>%
  group_by(id, trial, id_trial, order, Fly_ID, Fly_catched_1_in_next_3, stop_duration, stop_threshold, speed_ave_set, social_influence, prefix) %>%
  dplyr::summarize(n = n()) %>%
  right_join(df_spider_group_prefix) %>%
  dplyr::mutate(Fly_catched_1_in_next_3 = case_when(Fly_catched_1_in_next_3 == TRUE ~ "Success",
                                                    TRUE ~ "Failed"),
                n = if_else(is.na(n), 0, n)) %>%
  arrange(prefix, stop_duration, stop_threshold, speed_ave_set, social_influence, id, trial, id_trial, order, Fly_ID) %>%
  pivot_wider(id_cols = c(id, trial, id_trial, order, stop_duration, stop_threshold, speed_ave_set, social_influence, prefix), 
              names_from = c(Fly_ID, Fly_catched_1_in_next_3), values_from = n, values_fill = 0) %>%
  dplyr::select(!NA_Failed) %>%
  pivot_longer(cols = matches("^[0-5]+_"), names_sep = "_", names_to = c("Fly_ID", "Fly_catched_1_in_next_3")) %>%
  pivot_wider(id_cols = c(id, trial, id_trial, order, Fly_ID, stop_duration, stop_threshold, speed_ave_set, social_influence, prefix), 
              names_from = Fly_catched_1_in_next_3, values_from = value, values_fill = 0) %>%
  dplyr::mutate(Total = Success + Failed) %>%
  pivot_longer(cols = c("Success", "Failed", "Total"), names_to = "type", values_to = "num_attack") %>%
  dplyr::mutate(id_num = parse_number(id),
                social_influence = case_when(social_influence == "no" ~ "No interaction",
                                             social_influence == "weak" ~ "Weak interaction",
                                             social_influence == "strong" ~ "Strong interaction")) %>%
  arrange(id_num, trial, order) %>%
  dplyr::select(!id_num) %>%
  transform(stop_duration = factor(stop_duration,
                                   levels = c("1-1-1-1-1-1s", "1-1-1-3-3-3s", "1-1-1-6-6-6s", "1-1-1-12-12-12s",
                                              "3-3-3-3-3-3s", "3-3-3-6-6-6s", "3-3-3-12-12-12s",
                                              "6-6-6-6-6-6s", "6-6-6-12-12-12s",
                                              "12-12-12-12-12-12s")),
            social_influence = factor(social_influence, levels = c("No interaction", 
                                                                   "Weak interaction",
                                                                   "Strong interaction")),
            id_trial = factor(id_trial, levels = unique(.$id_trial) %>% gtools::mixedsort()))


df_spider_group_attack_dist <-
  left_join(df_spider_group_speed, 
            df_spider_group_aggression_num_attack %>%
              filter(type == "Success"))

write_parquet(df_spider_group_attack_dist, "../data/6_spiderACI/Experiment2_group/parquet/df_spider_group_attack_dist.parquet")


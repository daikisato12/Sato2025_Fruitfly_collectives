#### load packages ####
targetPackages <- c('tidyverse','arrow','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)


#### Figure S8 ####
##### load simulation dataset #####
###### groupsinglemove ######
files_groupsinglemove <- list.files("../data/5_simulation/1_grid_parameters/groupsinglemove/", pattern = "df*", full.names = TRUE, recursive=T) %>% 
  gtools::mixedsort()
param_t_groupsinglemove <- files_groupsinglemove %>%
  basename() %>%
  str_split("_") %>%
  map_chr(4)
names(param_t_groupsinglemove) <- seq(1:length(param_t_groupsinglemove))
param_a_groupsinglemove <- files_groupsinglemove %>%
  basename() %>%
  str_split("_") %>%
  map_chr(5)
names(param_a_groupsinglemove) <- seq(1:length(param_a_groupsinglemove))

df_groupsinglemove_pos <- map_dfr(files_groupsinglemove, arrow::read_parquet, .id = 'id') %>%
  dplyr::mutate(param_t = param_t_groupsinglemove[id] %>% as.numeric(),
                param_a = param_a_groupsinglemove[id] %>% as.numeric(),
                social_factor = 0,
                id = as.numeric(id)) %>%
  ungroup() %>%
  arrange(stim_time, id, rep)

###### groupsocialcontagion ######
files_groupsocialcontagion <- list.files("../data/5_simulation/1_grid_parameters/groupsocialcontagion", pattern = "df*", full.names = TRUE, recursive=T) %>% 
  gtools::mixedsort()
param_t_groupsocialcontagion <- files_groupsocialcontagion %>%
  basename() %>%
  str_split("_") %>%
  map_chr(4)
names(param_t_groupsocialcontagion) <- seq(1:length(param_t_groupsocialcontagion))
param_a_groupsocialcontagion <- files_groupsocialcontagion %>%
  basename() %>%
  str_split("_") %>%
  map_chr(5)
names(param_a_groupsocialcontagion) <- seq(1:length(param_a_groupsocialcontagion))
social_factor_groupsocialcontagion <- files_groupsocialcontagion %>%
  basename() %>%
  str_split("_") %>%
  map_chr(6) %>%
  str_remove(".parquet")
names(social_factor_groupsocialcontagion) <- seq(1:length(social_factor_groupsocialcontagion))

df_groupsocialcontagion_pos <- map_dfr(files_groupsocialcontagion, arrow::read_parquet, .id = 'id') %>%
  dplyr::mutate(param_t = param_t_groupsocialcontagion[id] %>% as.numeric(),
                param_a = param_a_groupsocialcontagion[id] %>% as.numeric(),
                social_factor = social_factor_groupsocialcontagion[id] %>% as.numeric(),
                id = as.numeric(id)) %>%
  ungroup() %>%
  arrange(stim_time, id, rep)

###### merge ######
df_groupmerge_pos <-
  df_groupsocialcontagion_pos %>%
  group_by(param_t, param_a, social_factor, stim_time) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  dplyr::mutate(type = "Group w interaction") %>%
  bind_rows(df_groupsinglemove_pos %>%
              group_by(param_t, param_a, social_factor, stim_time) %>%
              dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
              dplyr::mutate(type = "Group w/o interaction")) %>%
  dplyr::mutate(alpha = if_else(str_detect(type, "w/o"), 0.5, 1)) %>%
  transform(type = factor(type, levels = c("Group w/o interaction",
                                           "Group w interaction")))



##### load original fly dataset #####
df_f5min_speed_ave <- read_parquet("../data/1_single_strain/parquet/df_f5min_speed_ave.parquet") %>%
  ungroup()

df_s5min_2995_5990_speed <- read_parquet("../data/1_single_strain/parquet/df_s5min_2995_5990_speed.parquet") %>%
  ungroup()

df_s5min_2995_5990_speed_sg_normbyf5minave_strain <-
  df_s5min_2995_5990_speed %>%
  group_by(strain, stim_time, n_inds, sex) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  left_join(df_f5min_speed_ave %>%
              group_by(strain, n_inds, sex) %>%
              dplyr::summarize(speed_f5min_ave = mean(speed_f5min_ave, na.rm = TRUE)), 
            by = c("strain", "n_inds", "sex")) %>%
  group_by(strain, n_inds, sex, stim_time) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE),
                   speed_f5min_ave = mean(speed_f5min_ave, na.rm = TRUE)) %>%
  mutate(speed_normbyf5minave = speed / speed_f5min_ave)

##### calculate distance from real movement pattern of fly groups #####
df_groupmerge_realfly <-
  df_groupmerge_pos %>%
  left_join(df_s5min_2995_5990_speed_sg_normbyf5minave_strain %>%
              filter(sex == "Female", n_inds == "Group") %>%
              filter(strain %in% c("DGRP88", "DGRP101", "DGRP136", "DGRP161",
                                   "DGRP189", "DGRP309", "DGRP324", "DGRP357",
                                   "DGRP358", "DGRP360", "DGRP399", "DGRP786",
                                   "DGRP852", "DGRP855")) %>%
              group_by(stim_time) %>%
              dplyr::summarize(speed_fly = speed_normbyf5minave))

df_groupmerge_realfly_dist <-
  df_groupmerge_realfly %>%
  dplyr::mutate(dist = abs(speed_fly - speed)) %>%
  group_by(param_t, param_a, social_factor) %>%
  dplyr::summarize(dist = sum(dist, na.rm = TRUE)) %>%
  ungroup() %>%
  group_by(social_factor) %>%
  arrange(dist) %>%
  dplyr::mutate(top = if_else(row_number() <= 5, "yes", "no")) %>%
  ungroup() %>%
  transform(top = factor(top, levels = c("no", "yes")))

g_groupmerge_realfly_dist <-
  ggplot(df_groupmerge_realfly_dist %>%
           dplyr::mutate(param_a = as.factor(param_a),
                         social_factor = paste0("Social factor: ", social_factor)), 
         aes(x = param_t, y = param_a, fill = dist)) +
  geom_tile(aes(color = top), linewidth = 2) +
  # geom_text(aes(label = round(dist, 3)), size = 2, color = "white") +
  scale_color_manual(values = c("no" = "transparent", "yes" = "white")) +
  scale_fill_viridis_c(option = "F") +
  # scale_fill_gradientn("value", colours = rev(brewer.pal(9, "Spectral")), na.value = "white") +
  scale_x_continuous(breaks = seq(0, 4.5, by = 0.5)) +
  xlab("Parameter t") +
  ylab("Parameter a") +
  coord_cartesian(xlim = c(0, 4.5)) +
  facet_wrap(~ social_factor, ncol = 4) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 11, color = "black"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_text(color = "black"))
g_groupmerge_realfly_dist
ggsave("../figures/FigureS8.pdf", g_groupmerge_realfly_dist, w = 10, h = 8)


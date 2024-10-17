#### load libraries ####
library(tidyverse)
library(data.table)
library(arrow)

#### single ####
dt = 0.02
nrep_single = 1000

files_single <- list.files("../data/5_simulation/single/", pattern = "*pos*", full.names = TRUE, recursive=T)
names(files_single) <- str_remove(basename(files_single), "\\.tsv") %>% parse_number()

df_single_pos <- map_dfr(files_single, fread, .id = 'rep') %>%
  mutate(rep = as.numeric(rep),
         id = as.factor(id),
         seconds_diff = dt) %>%
  filter(rep < nrep_single)

df_single_pos2 <- df_single_pos %>% 
  dplyr::arrange(rep, id, frame) %>%
  mutate(
    seconds_total = frame * dt, 
    pos_x_1framebefore = lag(pos_x),
    pos_y_1framebefore = lag(pos_y),
    travelled_dist_diff = ( (pos_x - pos_x_1framebefore)**2 + (pos_y - pos_y_1framebefore)**2 )**(1/2) 
  ) %>%
  filter(frame >= 60) %>%
  group_by(rep, seconds_total = as.integer(seconds_total*2)/2, id) %>%
  dplyr::summarize(pos_x = mean(pos_x),
                   pos_y = mean(pos_y),
                   speed = sum(travelled_dist_diff, na.rm = TRUE) / sum(seconds_diff), 
                   travelled_dist_diff = sum(travelled_dist_diff, na.rm = TRUE)) %>%
  mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5)
write_parquet(df_single_pos2, "../data/5_simulation/df_single_pos2.parquet")


#### group_socialcontagion ####
dt = 0.02
nrep_group_socialcontagion = 200

files_group_socialcontagion <- list.files("../data/5_simulation/group_socialcontagion/", pattern = "*pos*", full.names = TRUE, recursive=T)
names(files_group_socialcontagion) <- str_remove(basename(files_group_socialcontagion), "\\.tsv") %>% parse_number()

df_group_socialcontagion_pos <- map_dfr(files_group_socialcontagion, fread, .id = 'rep') %>%
  mutate(rep = as.numeric(rep),
         id = as.factor(id),
         seconds_diff = dt) %>%
  filter(rep < nrep_group_socialcontagion)

df_group_socialcontagion_pos2 <- df_group_socialcontagion_pos %>% 
  dplyr::arrange(rep, id, frame) %>%
  mutate(
    seconds_total = frame * dt, 
    pos_x_1framebefore = lag(pos_x),
    pos_y_1framebefore = lag(pos_y),
    travelled_dist_diff = ( (pos_x - pos_x_1framebefore)**2 + (pos_y - pos_y_1framebefore)**2 )**(1/2) 
  ) %>%
  filter(frame >= 60) %>%
  group_by(rep, seconds_total = as.integer(seconds_total*2)/2, id) %>%
  dplyr::summarize(pos_x = mean(pos_x),
                   pos_y = mean(pos_y),
                   speed = sum(travelled_dist_diff, na.rm = TRUE) / sum(seconds_diff), #mean(travelled_distance_per_seconds), 
                   travelled_dist_diff = sum(travelled_dist_diff, na.rm = TRUE)) %>%
  mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5)
write_parquet(df_group_socialcontagion_pos2, "../data/5_simulation/df_group_socialcontagion_pos2.parquet")


#### group_socialcontagion_diversity ####
dt = 0.02
list_strains <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
nrep_group_socialcontagion_diversity = 10

files_group_socialcontagion_diversity <- list.files("../data/5_simulation/group_socialcontagion_diversity", pattern = "*pos*", full.names = TRUE, recursive=T) %>% 
  gtools::mixedsort()
index_strain <- 
  data.frame(strain1 = files_group_socialcontagion_diversity %>%
               basename() %>%
               str_split("_") %>%
               map_chr(3),
             strain2 = files_group_socialcontagion_diversity %>%
               basename() %>%
               str_split("_") %>%
               map_chr(4),
             rep = files_group_socialcontagion_diversity %>%
               basename() %>%
               str_split("_") %>%
               map_chr(5) %>%
               parse_number()) %>%
  tibble::rowid_to_column() %>%
  filter(strain1 %in% list_strains,
         strain2 %in% list_strains,
         rep < nrep_group_socialcontagion_diversity) %>%
  pull(rowid)

df_group_socialcontagion_diversity_pos <- map_dfr(files_group_socialcontagion_diversity[index_strain], fread) %>%
  mutate(strain = paste0(strain1, "_", strain2),
         id = as.factor(id),
         seconds_diff = dt)

df_group_socialcontagion_diversity_pos2 <- df_group_socialcontagion_diversity_pos %>% 
  dplyr::arrange(strain, strain1, strain2, rep, id, frame) %>%
  mutate(
    seconds_total = frame * dt, 
    pos_x_1framebefore = lag(pos_x),
    pos_y_1framebefore = lag(pos_y),
    travelled_dist_diff = ( (pos_x - pos_x_1framebefore)**2 + (pos_y - pos_y_1framebefore)**2 )**(1/2) 
  ) %>%
  filter(frame >= 60) %>%
  group_by(strain, strain1, strain2, rep, seconds_total = as.integer(seconds_total*2)/2, id) %>%
  dplyr::summarize(pos_x = mean(pos_x),
                   pos_y = mean(pos_y),
                   speed = sum(travelled_dist_diff, na.rm = TRUE) / sum(seconds_diff), #mean(travelled_distance_per_seconds), 
                   travelled_dist_diff = sum(travelled_dist_diff, na.rm = TRUE)) %>%
  mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5)
write_parquet(df_group_socialcontagion_diversity_pos2, "../data/5_simulation/df_group_socialcontagion_diversity_pos2.parquet")


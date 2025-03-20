#### load packages ####
targetPackages <- c('tidyverse','data.table','slider','gtools','sf','arrow')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### load functions ####
# function to add nearest neighbor distance for each row (each individual per frame)
add_nnd <- function(dat){
  dat <- dat %>%
    st_as_sf(., coords = c("pos_x","pos_y")) %>%
    st_distance() %>%
    data.frame(dat,
               dist_ave = apply(., 1, FUN = function(x) {if(length(x)==1 && x==0){NA}else{mean(x[x > 0])}}),
               nn_id = apply(., 1, FUN = function(x) {if(length(x)==1 && x==0){NA}else{paste0("Fly",which(x == min(x[x > 0])))}}),
               nnd = apply(., 1, FUN = function(x) {if(length(x)==1 && x==0){NA}else{min(x[x > 0])}})
    )
  return(dat)
}

calc_angle <- function(M,N){
  atan2(N[2],N[1]) - atan2(M[2],M[1]) 
}

correct_N <- function(dat){
  dat <- dat %>%
    dplyr::mutate(N2 = n()) %>%
    data.frame(., trial = seq(.$N2))
  return (dat)
}

add_samplenum_trial <- function(dat){
  
  dat_n <- dat %>%
    dplyr::select(filename, strain, sex, age, n_inds) %>%
    dplyr::distinct() %>%
    dplyr::group_by(strain, sex, age, n_inds) %>%
    dplyr::summarize(strain_n = n())
  
  dat_ig_n_1  <-  dat_n %>%
    dplyr::mutate(strain_i_n1 = lag(strain_n),
                  strain_ig_n1 = paste0(strain, " ", sex, " (n = ", strain_i_n1, "/", strain_n, ")")) #%>%
  
  dat_ig_n_2  <-  dat_n %>%
    dplyr::mutate(strain_i_n2 = lead(strain_n),
                  strain_ig_n2 = paste0(strain, " ", sex, " (n = ", strain_n, "/", strain_i_n2, ")")) %>%
    dplyr::select(-strain_n) %>%
    right_join(., dat_ig_n_1, by = c("strain", "sex", "age", "n_inds")) %>%
    dplyr::mutate(strain_ig_n = if_else(n_inds == 1, strain_ig_n2, strain_ig_n1)) %>%
    dplyr::group_by(strain, sex, age) %>%
    dplyr::summarize(strain_ig_n = strain_ig_n) %>%
    dplyr::distinct()
  
  dat_group_n <- dat %>%
    dplyr::select(filename, strain, sex, age, n_inds, id) %>%
    dplyr::distinct() %>%
    dplyr::group_by(strain, sex, age, n_inds) %>%
    dplyr::summarize(strain_ind_n = n())
  
  dat <- dat %>% left_join(dat_n, by = c("strain", "sex", "age", "n_inds")) %>%
    left_join(dat_group_n, by = c("strain", "sex", "age", "n_inds")) %>%
    left_join(dat_ig_n_2, by = c("strain", "sex", "age")) %>%
    dplyr::mutate(strain_n = paste0(strain, " ", sex, " (n = ", strain_n, ")"),
                  strain_ind_n = paste0(strain, " ", sex, " (n = ", strain_ind_n, ")"))
  
  dat <- dat %>% 
    group_nest(prefix, strain, sex, age, n_inds, N) %>% 
    dplyr::select(prefix, strain, sex, age, n_inds, N) %>%
    dplyr::distinct() %>%
    group_nest(strain, sex, age, n_inds) %>% 
    dplyr::mutate(data = map(data, correct_N)) %>% 
    unnest(data) %>%
    dplyr::select(strain, sex, age, n_inds, N, trial) %>%
    left_join(dat, ., by=c("strain", "sex", "age", "n_inds", "N"))
  
  return(dat)
}

count_visual_num_ind <- function(dat){
  length_arena_data = 100 #diameter of arena in data
  length_arena_real = 30 #mm
  a <- length_arena_data * 2/2 / length_arena_real #assuming 2mm as body length -> a = 2/2 mm
  b <- length_arena_data * 1/2 / length_arena_real #assuming 1mm as body width -> b = 1/2 mm
  set_angle <- 180
  count_li = c()
  motion_cue_li = c()
  for (i in 1:nrow(dat)){
    count <- 0
    motion_cue <- 0
    for (j in 1:nrow(dat)){
      if (i != j){
        vec1 <- c(dat[i,"angle_x_end"] %>% as.numeric() - dat[i,"pos_x"] %>% as.numeric(),
                  dat[i,"angle_y_end"] %>% as.numeric() - dat[i,"pos_y"] %>% as.numeric())
        vec2 <- c(dat[j,"pos_x"] %>% as.numeric() - dat[i,"pos_x"] %>% as.numeric(),
                  dat[j,"pos_y"] %>% as.numeric() - dat[i,"pos_y"] %>% as.numeric())
        angle_tmp <- abs(calc_angle(vec1, vec2)*180/pi)
        if (angle_tmp < set_angle / 2){ #count number of flies behind the focused individual
          count <- count + 1
        }
        
        psi <- atan2(dat[j,"pos_y"] %>% as.numeric() - dat[i,"pos_y"] %>% as.numeric(),
                     dat[j,"pos_x"] %>% as.numeric() - dat[i,"pos_x"] %>% as.numeric()) + pi/2
        psi <- if_else(psi > pi, psi - pi, psi)
        phi <- atan2(dat[j,"angle_y_end"] %>% as.numeric() - dat[j,"pos_y"] %>% as.numeric(),
                     dat[j,"angle_x_end"] %>% as.numeric() - dat[j,"pos_x"] %>% as.numeric())
        x0 <- a**2/((a**2+b**2*tan(psi-phi)^2)**(1/2))
        y0 <- b**2*tan(psi-phi)/((a**2+b**2*tan(psi-phi)^2)**(1/2))
        size = 2*a**2*b**2/((b**4*x0**2+a**4*y0**2)**(1/2))
        
        dist <- norm(vec2, type="2")
        visual_angle <- 2*atan(size/(2*dist))
        motion_cue <- motion_cue + dat[j,"speed"] %>% as.numeric() * visual_angle %>% as.double()
      }
    }
    count_li <- c(count_li, count) %>% as.integer()
    motion_cue_li <- c(motion_cue_li, motion_cue) %>% as.double()
  }
  dat <- data.frame(dat, 
                    count_visual_num_ind = count_li,
                    motion_cue = motion_cue_li)
}

calc_glm_coeff <- function(dat){
  glm_res <- dat %>%
    glm(formula = change_posture ~ motion_cue_diff3, family = "binomial")
  return(c(glm_res$coefficients[1],glm_res$coefficients[2]))
}


# a list to centralize the coordinate
cent_mask_x <- c(-5, -3, -2, 2, -5, -3, -2, 2, -5, -3, -2, 2)
cent_mask_y <- c(-5, -5, -5, -5, -1, -1, -1, -1, 4, 4, 4, 4)
names(cent_mask_x) <- paste0("no",seq(1,12))
names(cent_mask_y) <- paste0("no",seq(1,12))


#### load dataset and convert to parquet ####
datadir <- "$DATADIR" # data directory for tsv files

tbl_fread <- 
  list.files(paste0(datadir, "/1_single_strain/tsv/"), 
             pattern = "*.tsv", full.names = TRUE, recursive=T) %>% 
  map_df(~fread(.))

# set angle
num_roll_pos <- 10
num_roll_angle <- 50
slice_lines <- num_roll_pos + num_roll_angle
unit = 1 #making vector length = 1

df_tmp <- tbl_fread %>%
  dplyr::mutate(filename = paste0(prefix, "_900s_", place, "-", strain, "-", sex, "-", age, "-", n_inds, "-", N),
                prefix = as.character(prefix),
                date = str_sub(prefix, start=1, end=8),
                pos_x = pos_x_wma - cent_mask_x[place],
                pos_y = 160 - (pos_y_wma - cent_mask_y[place])) %>%
  filter(filename %in% include_samples) %>%

  # convert position to (0-100, 0-100)
  group_by(place) %>%
  mutate(pos_x_max = max(pos_x, na.rm = T),
         pos_x_min = min(pos_x, na.rm = T),
         pos_y_max = max(pos_y, na.rm = T),
         pos_y_min = min(pos_y, na.rm = T)) %>%
  ungroup() %>%
  arrange(filename, date, prefix, place, strain, sex, age, n_inds, N, id, seconds_total) %>%
  mutate(pos_x = (pos_x - pos_x_min) * 100 / (pos_x_max - pos_x_min),
         pos_y = (pos_y - pos_y_min) * 100 / (pos_y_max - pos_y_min),
         pos_x_1framebefore = lag(pos_x),
         pos_y_1framebefore = lag(pos_y),
         travelled_dist_diff = ( (pos_x - pos_x_1framebefore)**2 + (pos_y - pos_y_1framebefore)**2 )**(1/2) ) %>%
  
  # additional variables
  mutate(
    vec_norm = unit/(1+tan(angle_diff_based)**2)**(1/2), # when making vector same size
    angle_x_end = if_else(cos(angle_diff_based) > 0, pos_x + vec_norm, pos_x - vec_norm),
    angle_y_end = if_else(angle_diff_based > 0, 
                          (pos_y + abs(tan(angle_diff_based))*vec_norm),
                          (pos_y - abs(tan(angle_diff_based))*vec_norm)),
    angle = 180 * atan2(angle_y_end - pos_y, angle_x_end - pos_x) / pi,
    angle_stim = if_else(angle + 270 > 360, angle - 90, angle + 270),
    vec_x = angle_x_end - pos_x,
    vec_y = angle_y_end - pos_y) %>%
  
  # remove unnecessary columns
  group_by(filename, date, prefix, strain, sex, age, n_inds, N, place, id) %>%
  slice(slice_lines:n()) %>%
  ungroup() %>%
  dplyr::arrange(filename, date, prefix, strain, sex, age, n_inds, N, place, frame, id) %>%
  dplyr::select(c(filename, date, prefix, strain, sex, age, n_inds, N, place,
                  frame, timestamp, seconds_diff, seconds_total, stimuli,
                  id, pos_x, pos_y, travelled_dist_diff, angle_diff_based,
                  vec_norm, angle_x_end, angle_y_end, angle, angle_stim))

df <- df_tmp %>%
  group_by(filename, date, prefix, place, strain, sex, age, n_inds, N, id, seconds_total = as.integer(seconds_total*2)/2) %>%
  dplyr::summarize(pos_x = mean(pos_x),
                   pos_y = mean(pos_y),
                   angle_x_end = mean(angle_x_end),
                   angle_y_end = mean(angle_y_end),
                   speed = sum(travelled_dist_diff)/sum(seconds_diff),
                   travelled_dist_diff = sum(travelled_dist_diff),
  ) %>%
  dplyr::mutate(speed_ma = slide_vec(.x = speed, .f = mean, .before = 5),
                angle = 180 * atan2(angle_y_end - pos_y, angle_x_end - pos_x) / pi,
                angle_stim = if_else(angle + 270 > 360, angle - 90, angle + 270),
                posture = if_else(speed*30/100 < 4, "stop", "walk"), #slower than 4mm/s -> stopping (Zacaris et al. 2018 Nat Commmun)
                strain = as.factor(strain),
  ) %>%
  transform(strain = factor(strain, levels = mixedsort(levels(.$strain)))) %>%
  group_nest(filename, place, seconds_total) %>%
  dplyr::mutate(data = map(data, add_nnd)) %>%
  unnest(data) %>%
  dplyr::select(filename, prefix, place, strain, sex, age, n_inds, N, seconds_total, id, posture, everything()) %>%
  dplyr::select(-c(matches("^X\\d|\\.")))
write_parquet(df, paste0(datadir, "/1_single_strain/parquet/df_0.parquet"))


#### make df ####
df <- add_samplenum_trial(df) %>%
  group_by(strain, sex, age, n_inds, trial, id) %>%
  dplyr::summarize() %>%
  dplyr::filter(!(prefix == "20211011091236" & place == "no10" & trial == 2),
                !(prefix == "20211013084910" & place == "no5" & trial == 1)) %>%
  dplyr::mutate(vec_x = angle_x_end - pos_x,
                vec_y = angle_y_end - pos_y,
                angle_diff_based = atan2(vec_y, vec_x)) %>%
  group_by(filename, seconds_total) %>%
  dplyr::mutate(sum_vec_x = ifelse(n_inds != 1, sum(vec_x), NA),
                sum_vec_y = ifelse(n_inds != 1, sum(vec_y), NA),
                sum_r = ifelse(n_inds != 1, abs(sum(exp(1i*angle_diff_based))), NA), 
                other_ave_vec_x = (sum_vec_x - vec_x) / 5,
                other_ave_vec_y = (sum_vec_y - vec_y) / 5,
                orientation_order = sum_r / 6,
                sex = str_to_title(sex)) %>%
  ungroup()
write_parquet(df, paste0(datadir, "/data/1_single_strain/parquet/df.parquet"))
# df <- read_parquet(paste0(datadir, "/data/1_single_strain/parquet/df.parquet")) %>%
#   ungroup()

## data first 5 min
df_f5min <- df %>%
  filter(seconds_total < 300)
write_parquet(df_f5min, paste0(datadir, "/data/1_single_strain/parquet/df_f5min.parquet"))
# df_f5min <- read_parquet(paste0(datadir, "/data/1_single_strain/parquet/df_f5min.parquet")) %>%
#   ungroup()

## data first 10 min
df_f10min <- df %>%
  filter(seconds_total < 600)
write_parquet(df_f10min, paste0(datadir, "/data/1_single_strain/parquet/df_f10min.parquet"))
# df_f10min <- read_parquet(paste0(datadir, "/data/1_single_strain/parquet/df_f10min.parquet")) %>%
#   ungroup()

## data middle 5 min
df_s5min <- df %>%
  filter(298.5 < seconds_total, seconds_total < 600)
write_parquet(df_s5min, paste0(datadir, "/data/1_single_strain/parquet/df_s5min.parquet"))
# df_s5min <- read_parquet(paste0(datadir, "/data/1_single_strain/parquet/df_s5min.parquet")) %>%
#   ungroup()

## add stim info
df_s5min_stim <- df_s5min %>%
  mutate(stimuli = case_when(
    as.integer(seconds_total) %% 15 == 14 & seconds_total - as.integer(seconds_total) == 0 ~ "-1.0", 
    as.integer(seconds_total) %% 15 == 0 & seconds_total - as.integer(seconds_total) == 0.5 ~ "+0.5",
    as.integer(seconds_total) %% 15 == 10 & seconds_total - as.integer(seconds_total) == 0 ~ "+10.0",)) %>%
  dplyr::filter(!is.na(stimuli)) %>%
  mutate(n_inds = if_else(n_inds == 1, "Single", "Group"),
         time = as.POSIXct(prefix, format = format('%Y%m%d%H%M%S')) %>% format("%H:%M")) %>%
  transform(stimuli = factor(stimuli, levels=c("-1.0", "+0.5", "+10.0")),
            n_inds = factor(n_inds, levels=c("Single", "Group")),
            strain = factor(strain, levels = unique(.$strain) %>% as.character() %>% mixedsort())) %>%
  mutate(stim_no = as.integer((seconds_total-299) %/% 15 + 1),
         posture_1s = lead(posture),
         change_posture = case_when(posture == "walk" & posture_1s == "stop" ~ "ws",
                                    posture == "walk" & posture_1s == "walk" ~ "ww",
                                    posture == "stop" & posture_1s == "stop" ~ "ss",
                                    posture == "stop" & posture_1s == "walk" ~ "sw")) %>%
  filter(seconds_total != 599)

write_parquet(df_s5min_stim, paste0(datadir, "/data/1_single_strain/parquet/df_s5min_stim.parquet"))
# df_s5min_stim <- read_parquet(paste0(datadir, "/data/1_single_strain/parquet/df_s5min_stim.parquet")) %>%
#   ungroup()

## s5min (from 299.5-599)
df_s5min_2995_5990 <- df %>%
  filter(299 < seconds_total, seconds_total < 599.5)
write_parquet(df_s5min_2995_5990, paste0(datadir, "/data/1_single_strain/parquet/df_s5min_2995_5990.parquet"))
# df_s5min_2995_5990 <- read_parquet(paste0(datadir, "/data/1_single_strain/parquet/df_s5min_2995_5990.parquet")) %>%
#   ungroup()


#### additional dataset ####
##### moving speed #####
df_f5min_speed <- df_f5min %>%
  mutate(n_inds = if_else(n_inds == 1, "Single", "Group"),
         speed = speed * 30 /100) %>%
  group_by(prefix, date, place, strain, n_inds, sex, trial, id) %>%
  dplyr::summarize(speed = mean(speed, na.rm = T)) %>%
  mutate(time = as.POSIXct(prefix, format = format('%Y%m%d%H%M%S')) %>% format("%H:%M")) %>%
  arrange(strain, n_inds, sex, trial, id) %>%
  transform(n_inds = factor(n_inds, levels = c("Single", "Group")),
            strain = factor(strain, levels = unique(.$strain) %>% as.character() %>% gtools::mixedsort()))
write_parquet(df_f5min_speed, "../data/1_single_strain/parquet/df_f5min_speed.parquet")
# df_f5min_speed <- read_parquet("../data/1_single_strain/parquet/df_f5min_speed.parquet") %>%
#   ungroup()

df_f5min_speed_ave <- df_f5min_speed %>%
  group_by(strain, n_inds, sex, trial, id) %>%
  dplyr::summarize(speed_f5min_ave = mean(speed, na.rm = T))
write_parquet(df_f5min_speed_ave, "../data/1_single_strain/parquet/df_f5min_speed_ave.parquet")
# df_f5min_speed_ave <- read_parquet("../data/1_single_strain/parquet/df_f5min_speed_ave.parquet") %>%
#   ungroup()

df_f10min_speed <- df_f10min %>%
  mutate(n_inds = if_else(n_inds == 1, "Single", "Group"),
         speed = speed * 30 /100) %>%
  group_by(strain, n_inds, sex, seconds_total) %>%
  dplyr::summarize(speed = mean(speed, na.rm = T)) %>%
  transform(n_inds = factor(n_inds, levels = c("Single", "Group")),
            strain = factor(strain, levels = unique(.$strain) %>% as.character() %>% gtools::mixedsort()))
write_parquet(df_f10min_speed, "../data/1_single_strain/parquet/df_f10min_speed.parquet")
# df_f10min_speed <- read_parquet("../data/1_single_strain/parquet/df_f10min_speed.parquet") %>%
#   ungroup()

df_s5min_2995_5990_speed <- df_s5min_2995_5990 %>%
  mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5,
         speed = speed * 30 / 100,
         n_inds = if_else(n_inds == 1, "Single", "Group"),
         time = as.POSIXct(prefix, format = format('%Y%m%d%H%M%S')) %>% format("%H:%M")) %>%
  group_by(strain, sex, n_inds, trial, id, stim_time, date, time, filename, prefix, place) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  transform(n_inds = factor(n_inds, levels = c("Single", "Group")),
            strain = factor(strain, levels = unique(.$strain) %>% as.character() %>% mixedsort()))
write_parquet(df_s5min_2995_5990_speed, "../data/1_single_strain/parquet/df_s5min_2995_5990_speed.parquet")
# df_s5min_2995_5990_speed <- read_parquet("../data/1_single_strain/parquet/df_s5min_2995_5990_speed.parquet") %>%
#   ungroup()


##### nnd with pseudo-group #####
include_samples <- read.csv("../data/1_single_strain/included_samples.txt", header = FALSE) %>%
  pull(V1)
sex <- c("male", "female")
df_f5min_rand <- data.frame()
for (sex in sex){
  file_random_list <- list()
  for (i in 1:20){
    set.seed(i)
    place <- sample(1:12, 1)
    pattern <- paste0("no",place,"-.*-",sex,"-3to5-1-")
    file_random <- sample_n(as_tibble(include_samples) %>%
                              filter(str_detect(value, "DGRP"), str_detect(value, pattern)), 6)
    file_random_list[[i]] <- c(file_random$value)
    df_f5min_test <- filter(df_f5min, filename %in% file_random_list[[i]]) %>% ungroup()
    df_f5min_test2 <- df_f5min_test %>%
      group_by(seconds_total) %>%
      dplyr::mutate(speed_mean = mean(speed)) %>%
      ungroup() %>%
      filter(speed_mean == min(.$speed_mean)) %>%
      group_nest(seconds_total) %>%
      dplyr::mutate(data = map(data, add_nnd)) %>%
      unnest(data)
    df_f5min_rand <- bind_rows(df_f5min_rand, 
                               data.frame(sex = sex, 
                                          strain = "random", 
                                          prefix = sample(df_f5min_test2$prefix, 1),
                                          place = sample(df_f5min_test2$place, 1),
                                          nnd = mean(df_f5min_test2$nnd.1) * 30 / 100))
  }
}
df_f5min_rand <- df_f5min_rand %>%
  dplyr::mutate(sex = str_to_title(sex)) %>%
  bind_cols(data.frame(trial = rep(seq(1:20), 2)))
write_parquet(df_f5min_rand, "../data/1_single_strain/parquet/df_f5min_rand.parquet")
# df_f5min_rand <- read_parquet("../data/1_single_strain/parquet/df_f5min_rand.parquet") %>%
#   ungroup()

df_f5min_nnd <- df_f5min %>%
  filter(n_inds == 6) %>%
  group_by(prefix, date, place, strain, sex, trial, seconds_total) %>%
  dplyr::mutate(speed_mean = mean(speed * 30 / 100, na.rm = T),
                nnd = mean(nnd * 30 / 100, na.rm = T)) %>%
  ungroup() %>%
  group_by(prefix, date, place, strain, sex, trial) %>%
  slice(which.min(speed_mean)) %>%
  arrange(strain, sex, trial)

df_f5min_nnd_rand <- 
  bind_rows(df_f5min_nnd, 
            df_f5min_rand) %>%
  mutate(date = str_sub(prefix, start=1, end=8),
         time = str_sub(prefix, start=9, end=12) %>% as.POSIXct(format = format('%H%M'))) %>%
  dplyr::select(date, time, prefix, place, strain, sex, trial, nnd) %>%
  arrange(strain, sex, trial)
write_parquet(df_f5min_nnd_rand, "../data/1_single_strain/parquet/df_f5min_nnd_rand.parquet")
# df_f5min_nnd_rand <- read_parquet("../data/1_single_strain/parquet/df_f5min_nnd_rand.parquet") %>%
#   ungroup()


###### nnd female scaled ######
df_out_nnd_female_scaled_GREML <- df_f5min_nnd %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Female") %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  group_by(strain) %>%
  dplyr::summarize(nnd = mean(nnd, na.rm=T)) %>%
  mutate(row = strain) %>%
  dplyr::select(row, strain, nnd) %>%
  mutate_at(vars(starts_with("nnd")), ~(predict(bestNormalize(.)) %>% as.vector))
write_delim(df_out_nnd_female_scaled_GREML, "../data/1_single_strain/gwas/pheno/df_out_nnd_female_scaled_GREML.txt", delim = " ", col_names = F)

###### nnd male scaled ######
df_out_nnd_male_scaled_GREML <- df_f5min_nnd %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Male") %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  group_by(strain) %>%
  dplyr::summarize(nnd = mean(nnd, na.rm=T)) %>%
  mutate(row = strain) %>%
  dplyr::select(row, strain, nnd) %>%
  mutate_at(vars(starts_with("nnd")), ~(predict(bestNormalize(.)) %>% as.vector))
write_delim(df_out_nnd_male_scaled_GREML, "../data/1_single_strain/gwas/pheno/df_out_nnd_male_scaled_GREML.txt", delim = " ", col_names = F)


##### motion cue exit #####
## calculate motion cue
df_s5min_stim_vis <- df_s5min_stim %>%
  group_nest(strain, sex, n_inds, seconds_total, stimuli, stim_no, place, prefix, date, time) %>%
  dplyr::mutate(data = map(data, count_visual_num_ind)) %>%
  unnest(data) %>%
  arrange(strain, n_inds, sex, trial, id, stim_no, stimuli)
write_parquet(df_s5min_stim_vis, "../data/1_single_strain/parquet/df_s5min_stim_vis_0.5.parquet")
# df_s5min_stim_vis <- read_parquet("../data/1_single_strain/parquet/df_s5min_stim_vis_0.5.parquet") %>%
#   ungroup()

df_s5min_stim_freez_vis <- df_s5min_stim_vis %>%
  dplyr::filter(n_inds == "Group") %>%
  transform(change_posture = factor(change_posture, levels=c("ww", "ws", "sw", "ss"))) %>%
  mutate(motion_cue_diff = lead(motion_cue) - motion_cue,
         motion_cue_diff2 = motion_cue - lead(motion_cue),
         motion_cue_diff3 = motion_cue - lag(motion_cue),
         motion_cue_next = lead(motion_cue))

df_motion_cue_exit <- df_s5min_stim_freez_vis %>%
  filter(stimuli == "+10.0", change_posture %in% c("ss", "sw")) %>%
  mutate(change_posture = if_else(change_posture == "sw", 1, 0),
         posture_num = if_else(posture == "walk", 1, 0))

df_motion_cue_exit_coeff <- df_motion_cue_exit %>%
  group_nest(strain, sex) %>%
  dplyr::mutate(data = map(data, calc_glm_coeff)) %>%
  unnest(data) %>%
  group_by(strain, sex) %>%
  dplyr::mutate(var = if_else(row_number() == 1, "motion_cue_exit_intercept", "motion_cue_exit_coeff")) %>%
  pivot_wider(names_from = var, values_from = data)

write_parquet(df_motion_cue_exit_coeff, "../data/1_single_strain/parquet/df_motion_cue_exit_coeff.parquet")
# df_motion_cue_exit_coeff <- read_parquet("../data/1_single_strain/parquet/df_motion_cue_exit_coeff.parquet") %>%
#   ungroup()

df_motion_cue_exit_coeff_trial <- df_motion_cue_exit %>%
  group_nest(strain, sex, trial, date, time, prefix, place) %>%
  dplyr::mutate(data = map(data, calc_glm_coeff)) %>%
  unnest(data) %>%
  group_by(strain, sex, trial, date, time, prefix, place) %>%
  dplyr::mutate(var = if_else(row_number() == 1, "motion_cue_exit_intercept", "motion_cue_exit_coeff")) %>%
  pivot_wider(names_from = var, values_from = data)

write_parquet(df_motion_cue_exit_coeff_trial, "../data/1_single_strain/parquet/df_motion_cue_exit_coeff_trial.parquet")
# df_motion_cue_exit_coeff_trial <- read_parquet("../data/1_single_strain/parquet/df_motion_cue_exit_coeff_trial.parquet") %>%
#   ungroup()

###### motion_cue_exit_intecept female ######
df_out_motion_cue_exit_intercept_female_GREML <- df_motion_cue_exit_coeff %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Female") %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  dplyr::select(row, strain, motion_cue_exit_intercept)
write_delim(df_out_motion_cue_exit_intercept_female_GREML, "../data/1_single_strain/gwas/pheno/df_out_motion_cue_exit_intercept_female_GREML.txt", delim = " ", col_names = F)

###### motion_cue_exit_intecept male ######
df_out_motion_cue_exit_intercept_male_GREML <- df_motion_cue_exit_coeff %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Male") %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  dplyr::select(row, strain, motion_cue_exit_intercept)
write_delim(df_out_motion_cue_exit_intercept_male_GREML, "../data/1_single_strain/gwas/pheno/df_out_motion_cue_exit_intercept_male_GREML.txt", delim = " ", col_names = F)


###### motion_cue_exit_intecept female scaled ######
set.seed(12345)
df_out_motion_cue_exit_intercept_female_scaled_GREML <- df_motion_cue_exit_coeff %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Female") %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  dplyr::select(row, strain, motion_cue_exit_intercept) %>% 
  as.data.frame() %>%
  mutate_at(vars(starts_with("motion_cue_exit_intercept")), ~(predict(bestNormalize(.)) %>% as.vector))
write_delim(df_out_motion_cue_exit_intercept_female_scaled_GREML, "../data/1_single_strain/gwas/pheno/df_out_motion_cue_exit_intercept_female_scaled_GREML.txt", delim = " ", col_names = F)

###### motion_cue_exit_intecept male scaled ######
set.seed(12345)
df_out_motion_cue_exit_intercept_male_scaled_GREML <- df_motion_cue_exit_coeff %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Male") %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  dplyr::select(row, strain, motion_cue_exit_intercept) %>% 
  as.data.frame() %>%
  mutate_at(vars(starts_with("motion_cue_exit_intercept")), ~(predict(bestNormalize(.)) %>% as.vector))
write_delim(df_out_motion_cue_exit_intercept_male_scaled_GREML, "../data/1_single_strain/gwas/pheno/df_out_motion_cue_exit_intercept_male_scaled_GREML.txt", delim = " ", col_names = F)



##### motion cue alternative metric #####
df_s5min_2995_5990_vis <- df_s5min_2995_5990 %>%
  mutate(stim_no = as.integer((seconds_total - 299) %/% 15 + 1),
         stim_time = (seconds_total + 0.5) %% 15 - 0.5,
         n_inds = if_else(n_inds == 1, "Single", "Group"),
         time = as.POSIXct(prefix, format = format('%Y%m%d%H%M%S')) %>% format("%H:%M"),
         speed = speed * 30 / 100) %>%
  group_nest(filename, stim_no, stim_time) %>%
  dplyr::mutate(data = map(data, count_visual_num_ind)) %>%
  unnest(data) %>%
  arrange(strain, n_inds, sex, trial, id, stim_no, stim_time)

write_parquet(df_s5min_2995_5990_vis, "../data/1_single_strain/parquet/df_s5min_2995_5990_vis_0.5.parquet")
# df_s5min_2995_5990_vis <- read_parquet("../data/1_single_strain/parquet/df_s5min_2995_5990_vis_0.5.parquet") %>%
#   ungroup()

df_s5min_2995_5990_vis_speed <-
  df_s5min_2995_5990_vis %>%
  group_by(strain, sex, n_inds, trial, id, stim_no, stim_time, date, time, filename, prefix, place) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE),
                   motion_cue  = mean(motion_cue, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(df_f5min_speed_ave, 
            by = c("strain", "n_inds", "sex", "trial", "id"))
write_parquet(df_s5min_2995_5990_vis_speed, "../data/1_single_strain/parquet/df_s5min_2995_5990_vis_speed.parquet")
# df_s5min_2995_5990_vis_speed <- read_parquet("../data/1_single_strain/parquet/df_s5min_2995_5990_vis_speed.parquet") %>%
#   ungroup()

df_s5min_2995_5990_vis_speed_filt0.5_no <- 
  df_s5min_2995_5990_vis_speed %>%
  dplyr::filter(n_inds == "Group") %>%
  ungroup() %>%
  dplyr::filter(stim_no <= 20, stim_time > 0) %>%
  dplyr::group_by(strain, sex, n_inds, trial, id, stim_no) %>%
  dplyr::mutate(motion_cue_cumsum = cumsum(motion_cue)) %>%
  ungroup()

df_s5min_2995_5990_vis_speed_filt0.5_no_exit_abs4mmsec <-
  df_s5min_2995_5990_vis_speed_filt0.5_no %>%
  dplyr::filter(speed > 4) %>%
  group_by(strain, sex, n_inds, trial, id, stim_no) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(strain, sex, n_inds, id, stim_no),
              names_from = trial, values_from = stim_time,
              values_fill = 14) %>%
  pivot_longer(cols = !c(strain, sex, n_inds, id, stim_no),
               names_to = "trial",
               values_to = "stim_time") %>%
  dplyr::mutate(trial = as.numeric(trial)) %>%
  left_join(df_s5min_2995_5990_vis_speed_filt0.5_no,
            by = c("strain", "sex", "n_inds", "trial", "id", "stim_no", "stim_time")) %>%
  dplyr::mutate(motion_cue_cumsum_ave = motion_cue_cumsum / stim_time) %>%
  arrange(strain, sex, n_inds, trial, stim_no, id) %>%
  dplyr::select(strain, sex, n_inds, trial, id, stim_no, everything())

write_parquet(df_s5min_2995_5990_vis_speed_filt0.5_no_exit_abs4mmsec, "../data/1_single_strain/parquet/df_s5min_2995_5990_vis_speed_filt0.5_no_exit_abs4mmsec.parquet")
# df_s5min_2995_5990_vis_speed_filt0.5_no_exit_abs4mmsec <- read_parquet("../data/1_single_strain/parquet/df_s5min_2995_5990_vis_speed_filt0.5_no_exit_abs4mmsec.parquet") %>%
#   ungroup()

###### motion_cue_stop_freezing_abs female scaled ######
set.seed(12345)
df_out_motion_cue_stop_freezing_abs_female_scaled_GREML <- df_s5min_2995_5990_vis_speed_filt0.5_no_exit_abs4mmsec %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Female") %>%
  group_by(strain) %>%
  dplyr::summarize(motion_cue_stop_freezing_abs = mean(motion_cue_cumsum, na.rm=T)) %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  dplyr::select(row, strain, motion_cue_stop_freezing_abs) %>%
  as.data.frame() %>%
  # mutate_at(vars(starts_with("motion_cue_stop_freezing_abs")), ~(scale(.) %>% as.vector))
  mutate_at(vars(starts_with("motion_cue_stop_freezing_abs")), 
            ~(bestNormalize::bestNormalize(.) %>% predict() %>% as.vector))
write_delim(df_out_motion_cue_stop_freezing_abs_female_scaled_GREML, "../data/1_single_strain/gwas/pheno/df_out_motion_cue_stop_freezing_abs_female_scaled_GREML.txt", delim = " ", col_names = F)

###### motion_cue_stop_freezing_abs male scaled ######
set.seed(12345)
df_out_motion_cue_stop_freezing_abs_male_scaled_GREML <- df_s5min_2995_5990_vis_speed_filt0.5_no_exit_abs4mmsec %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Male") %>%
  group_by(strain) %>%
  dplyr::summarize(motion_cue_stop_freezing_abs = mean(motion_cue_cumsum, na.rm=T)) %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  dplyr::select(row, strain, motion_cue_stop_freezing_abs) %>%
  as.data.frame() %>%
  # mutate_at(vars(starts_with("motion_cue_stop_freezing_abs")), ~(scale(.) %>% as.vector))
  mutate_at(vars(starts_with("motion_cue_stop_freezing_abs")), 
            ~(bestNormalize::bestNormalize(.) %>% predict() %>% as.vector))
write_delim(df_out_motion_cue_stop_freezing_abs_male_scaled_GREML, "../data/1_single_strain/gwas/pheno/df_out_motion_cue_stop_freezing_abs_male_scaled_GREML.txt", delim = " ", col_names = F)


##### freezing duration #####
df_s5min_2995_5990_speed_sg_normbyf5minave_strain_trial <-
  df_s5min_2995_5990_speed %>%
  group_by(strain, sex, n_inds, trial, id, stim_time, date, time, prefix, place) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  left_join(df_f5min_speed_ave, 
            by = c("strain", "n_inds", "sex", "trial", "id")) %>%
  group_by(strain, n_inds, sex, trial, id, stim_time, date, time, prefix, place) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE),
                   speed_f5min_ave = mean(speed_f5min_ave, na.rm = TRUE)) %>%
  mutate(speed_normbyf5minave = speed / speed_f5min_ave)

df_s5min_2995_5990_freezing_duration <-
  df_s5min_2995_5990_speed_sg_normbyf5minave_strain_trial %>% 
  filter(stim_time != -0.5, speed_normbyf5minave > 1) %>% 
  group_by(strain, sex, n_inds, trial, id) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(strain, sex, n_inds, id),
              names_from = trial, values_from = stim_time,
              values_fill = 14) %>%
  pivot_longer(cols = !c(strain, sex, n_inds, id), names_to = "trial", values_to = "freezing_duration") %>%
  mutate(trial = as.numeric(trial)) %>%
  group_by(strain, n_inds, sex, trial) %>%
  dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm = T)) %>%
  ungroup() %>%
  arrange(strain, sex, n_inds, trial) %>%
  left_join(df_s5min_2995_5990_speed_sg_normbyf5minave_strain_trial %>%
              group_by(strain, sex, n_inds, trial) %>%
              dplyr::slice(1) %>%
              dplyr::select(strain, sex, n_inds, trial, date, time, prefix, place))

write_parquet(df_s5min_2995_5990_freezing_duration, "../data/1_single_strain/parquet/df_s5min_2995_5990_freezing_duration.parquet")
# df_s5min_2995_5990_freezing_duration <- read_parquet("../data/1_single_strain/parquet/df_s5min_2995_5990_freezing_duration.parquet") %>%
#   ungroup()


###### freezing_duration_single female scaled ######
set.seed(12345)
df_out_freezing_duration_single_female_scaled_GREML <- df_s5min_2995_5990_freezing_duration %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Female", n_inds == "Single") %>%
  group_by(strain) %>%
  dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm=T)) %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  dplyr::select(row, strain, freezing_duration) %>%
  as.data.frame() %>%
  # mutate_at(vars(starts_with("freezing_duration")), ~(scale(.) %>% as.vector))
  mutate_at(vars(starts_with("freezing_duration")), 
            ~(bestNormalize::bestNormalize(.) %>% predict() %>% as.vector))
write_delim(df_out_freezing_duration_single_female_scaled_GREML, "../data/1_single_strain/gwas/pheno/df_out_freezing_duration_single_female_scaled_GREML.txt", delim = " ", col_names = F)

###### freezing_duration_single male scaled ######
set.seed(12345)
df_out_freezing_duration_single_male_scaled_GREML <- df_s5min_2995_5990_freezing_duration %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Male", n_inds == "Single") %>%
  group_by(strain) %>%
  dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm=T)) %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  dplyr::select(row, strain, freezing_duration) %>%
  as.data.frame() %>%
  # mutate_at(vars(starts_with("freezing_duration")), ~(scale(.) %>% as.vector))
  mutate_at(vars(starts_with("freezing_duration")), 
            ~(bestNormalize::bestNormalize(.) %>% predict() %>% as.vector))
write_delim(df_out_freezing_duration_single_male_scaled_GREML, "../data/1_single_strain/gwas/pheno/df_out_freezing_duration_single_male_scaled_GREML.txt", delim = " ", col_names = F)

###### freezing_duration_group female scaled ######
df_out_freezing_duration_group_female_scaled_GREML <- df_s5min_2995_5990_freezing_duration %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Female", n_inds == "Group") %>%
  group_by(strain) %>%
  dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm=T)) %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  dplyr::select(row, strain, freezing_duration) %>%
  as.data.frame() %>%
  # mutate_at(vars(starts_with("freezing_duration")), ~(scale(.) %>% as.vector))
  mutate_at(vars(starts_with("freezing_duration")), 
            ~(bestNormalize::bestNormalize(.) %>% predict() %>% as.vector))
write_delim(df_out_freezing_duration_group_female_scaled_GREML, "../data/1_single_strain/gwas/pheno/df_out_freezing_duration_group_female_scaled_GREML.txt", delim = " ", col_names = F)

###### freezing_duration_group male scaled ######
set.seed(12345)
df_out_freezing_duration_group_male_scaled_GREML <- df_s5min_2995_5990_freezing_duration %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","DGRP208_norpA"), sex == "Male", n_inds == "Group") %>%
  group_by(strain) %>%
  dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm=T)) %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  dplyr::select(row, strain, freezing_duration) %>%
  as.data.frame() %>%
  # mutate_at(vars(starts_with("freezing_duration")), ~(scale(.) %>% as.vector))
  mutate_at(vars(starts_with("freezing_duration")), 
            ~(bestNormalize::bestNormalize(.) %>% predict() %>% as.vector))
write_delim(df_out_freezing_duration_group_male_scaled_GREML, "../data/1_single_strain/gwas/pheno/df_out_freezing_duration_group_male_scaled_GREML.txt", delim = " ", col_names = F)


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
  # ref: https://oshiete.goo.ne.jp/qa/6059674.html
  length_arena_data = 100 #diameter of arena in data
  length_arena_real = 30 #mm
  a <- length_arena_data * 2/2 / length_arena_real #assuming 2mm as body length -> a = 2/2 mm
  b <- length_arena_data * 1/2 / length_arena_real #assuming 1mm as body width -> b = 1/2 mm
  #  dat <- dfm_test
  set_angle <- 180
  count_li = c()
  motion_cue_li = c()
  for (i in 1:nrow(dat)){
    # i <- 3
    count <- 0
    motion_cue <- 0
    for (j in 1:nrow(dat)){
      # j <- 4
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


calc_glm_coeff_500trial <- function(dat, trial){
  # dat <- dfm_motion_cue_exit %>%
  #   filter(strain %in% c("Dop1R1_KO"), sex == "Female")
  res <- data.frame()
  for(i in 1:trial){
    set.seed(i)
    glm_res <- dat %>%
      sample_n(500, replace = TRUE) %>%
      glm(formula = change_posture ~ motion_cue_diff3, family = "binomial")
    res <- bind_rows(res,
                     data.frame(trial = i, 
                                motion_cue_exit_intercept = glm_res$coefficients[1],
                                motion_cue_exit_coeff = glm_res$coefficients[2]))
  }
  return(res)
}


# a list to centralize the coordinate
cent_mask_x <- c(-5, -3, -2, 2, -5, -3, -2, 2, -5, -3, -2, 2)
cent_mask_y <- c(-5, -5, -5, -5, -1, -1, -1, -1, 4, 4, 4, 4)
names(cent_mask_x) <- paste0("no",seq(1,12))
names(cent_mask_y) <- paste0("no",seq(1,12))


#### load dataset and convert to parquet ####
datadir <- "$DATADIR" # data directory for tsv files

tbl_fread <- 
  list.files(paste0(datadir, "/3_mutants/tsv/"), 
             pattern = "*.tsv", full.names = TRUE, recursive=T) %>% 
  map_df(~fread(.))

# set angle
num_roll_pos <- 10
num_roll_angle <- 50
slice_lines <- num_roll_pos + num_roll_angle
unit = 1 #making vector length = 1

dfm_tmp <- tbl_fread %>%
  dplyr::mutate(filename = paste0(prefix, "_900s_", place, "-", strain, "-", sex, "-", age, "-", n_inds, "-", N),
                prefix = as.character(prefix),
                date = str_sub(prefix, start=1, end=8),
                pos_x = pos_x_wma - cent_mask_x[place],
                pos_y = 160 - (pos_y_wma - cent_mask_y[place])) %>%
  # filter(!filename %in% exclude_samples) %>%
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
  mutate(# vec_norm = (1+tan(angle_diff_based)**2)**(1/2),
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

dfm <- dfm_tmp %>%
  group_by(filename, date, prefix, place, strain, sex, age, n_inds, N, id, seconds_total = as.integer(seconds_total*2)/2) %>%
  dplyr::summarize(pos_x = mean(pos_x),
                   pos_y = mean(pos_y),
                   angle_x_end = mean(angle_x_end),
                   angle_y_end = mean(angle_y_end),
                   speed = sum(travelled_dist_diff)/sum(seconds_diff), #mean(travelled_distance_per_seconds), 
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

write_parquet(dfm, paste0(datadir, "/3_mutants/parquet/dfm_0.parquet"))


dfm <- add_samplenum_trial(dfm) %>%
  group_by(strain, sex, age, n_inds, trial, id) %>%
  dplyr::summarize() %>%
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

## data first 5 min
dfm_f5min <- dfm %>%
  filter(seconds_total < 300)
write_parquet(dfm_f5min, "../data/3_mutants/parquet/dfm_f5min.parquet")
# dfm_f5min <- read_parquet("../data/3_mutants/parquet/dfm_f5min.parquet") %>%
#   ungroup()

## data first 10 min
dfm_f10min <- dfm %>%
  filter(seconds_total < 600)
write_parquet(dfm_f10min, "../data/3_mutants/parquet/dfm_f10min.parquet")
# dfm_f10min <- read_parquet("../data/3_mutants/parquet/dfm_f10min.parquet") %>%
#   ungroup()

## data middle 5 min
dfm_s5min <- dfm %>%
  filter(298.5 < seconds_total, seconds_total < 600)
write_parquet(dfm_s5min, "../data/3_mutants/parquet/dfm_s5min.parquet")
# dfm_s5min <- read_parquet("../data/3_mutants/parquet/dfm_s5min.parquet") %>%
#   ungroup()

## add stim info
dfm_s5min_stim <- dfm_s5min %>%
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

write_parquet(dfm_s5min_stim, "../data/3_mutants/parquet/dfm_s5min_stim.parquet")
# dfm_s5min_stim <- read_parquet("../data/3_mutants/parquet/dfm_s5min_stim.parquet") %>%
#   ungroup()

## s5min (from 299.5-599) for foraging success
dfm_s5min_2995_5990 <- dfm %>%
  filter(299 < seconds_total, seconds_total < 599.5)

write_parquet(dfm_s5min_2995_5990, "../data/3_mutants/parquet/dfm_s5min_2995_5990.parquet")
# dfm_s5min_2995_5990 <- read_parquet("../data/3_mutants/parquet/dfm_s5min_2995_5990.parquet") %>%
#   ungroup()




#### additional dataset ####
##### dfm_s5min_2995_5990_freezing_duration #####
dfm_s5min_2995_5990_speed <- 
  dfm %>%
  filter(299 < seconds_total, seconds_total < 599.5) %>%
  mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5,
         speed = speed * 30 / 100,
         n_inds = if_else(n_inds == 1, "Single", "Group"),
         time = as.POSIXct(prefix, format = format('%Y%m%d%H%M%S')) %>% format("%H:%M")) %>%
  group_by(strain, sex, n_inds, trial, id, stim_time, date, time, filename, prefix, place) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  transform(n_inds = factor(n_inds, levels = c("Single", "Group")),
            strain = factor(strain, levels = unique(.$strain) %>% as.character() %>% mixedsort()))


dfm_s5min_2995_5990_speed_sg_normbyf5minave_strain_trial <-
  dfm_s5min_2995_5990_speed %>%
  group_by(strain, stim_time, n_inds, sex, trial, id, date, time, prefix, place) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  left_join(dfm_f5min_speed_ave, 
            by = c("strain", "n_inds", "sex", "trial", "id")) %>%
  group_by(strain, n_inds, sex, trial, id, stim_time, date, time, prefix, place) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE),
                   speed_f5min_ave = mean(speed_f5min_ave, na.rm = TRUE)) %>%
  mutate(speed_normbyf5minave = speed / speed_f5min_ave)

dfm_s5min_2995_5990_freezing_duration <-
  dfm_s5min_2995_5990_speed_sg_normbyf5minave_strain_trial %>% #dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial #dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial
  filter(stim_time != -0.5, speed_normbyf5minave > 1) %>% #speed_norm_mean #speed_normbystimminus05
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
  left_join(dfm_s5min_2995_5990_speed_sg_normbyf5minave_strain_trial %>%
              group_by(strain, sex, n_inds, trial) %>%
              dplyr::slice(1) %>%
              dplyr::select(strain, sex, n_inds, trial, date, time, prefix, place))
write_parquet(dfm_s5min_2995_5990_freezing_duration, "../data/3_mutants/dfm_s5min_2995_5990_freezing_duration.parquet")


##### dfm_motion_cue_exit_coeff_trial2 #####
dfm_s5min_stim_vis <- dfm_s5min_stim %>%
  group_nest(strain, sex, n_inds, seconds_total, stimuli, stim_no, place, prefix, date, time) %>%
  dplyr::mutate(data = map(data, count_visual_num_ind)) %>%
  unnest(data) %>%
  arrange(strain, n_inds, sex, trial, id, stim_no, stimuli)

dfm_s5min_stim_freez_vis <- dfm_s5min_stim_vis %>%
  dplyr::filter(n_inds == "Group") %>%
  transform(change_posture = factor(change_posture, levels=c("ww", "ws", "sw", "ss"))) %>%
  mutate(motion_cue_diff = lead(motion_cue) - motion_cue,
         motion_cue_diff2 = motion_cue - lead(motion_cue),
         motion_cue_diff3 = motion_cue - lag(motion_cue),
         motion_cue_next = lead(motion_cue))

dfm_motion_cue_exit <- dfm_s5min_stim_freez_vis %>%
  filter(stimuli == "+10.0", change_posture %in% c("ss", "sw")) %>%
  mutate(change_posture = if_else(change_posture == "sw", 1, 0),
         posture_num = if_else(posture == "walk", 1, 0))

dfm_motion_cue_exit_coeff_trial2 <- 
  dfm_motion_cue_exit %>%
  group_nest(strain, sex) %>%
  dplyr::mutate(data = map2(data, 20, calc_glm_coeff_500trial)) %>%
  unnest(data)
write_parquet(dfm_motion_cue_exit_coeff_trial2, "../data/3_mutants/tracked_data/dfm_motion_cue_exit_coeff_trial2.parquet")
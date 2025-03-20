#### load packages ####
targetPackages <- c('tidyverse','arrow','data.table','slider','gtools','sf')
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


# a list to centralize the coordinate
cent_mask_x <- c(-5, -3, -2, 2, -5, -3, -2, 2, -5, -3, -2, 2)
cent_mask_y <- c(-5, -5, -5, -5, -1, -1, -1, -1, 4, 4, 4, 4)
names(cent_mask_x) <- paste0("no",seq(1,12))
names(cent_mask_y) <- paste0("no",seq(1,12))


#### load dataset and convert to parquet ####
datadir <- "$DATADIR" # data directory for tsv files

tbl_fread <- 
  list.files(paste0(datadir, "/2_mixed_strain/tsv/"), 
             pattern = "*.tsv", full.names = TRUE, recursive=T) %>% 
  map_df(~fread(.))

# set angle
num_roll_pos <- 10
num_roll_angle <- 50
slice_lines <- num_roll_pos + num_roll_angle
unit = 1 #making vector length = 1

dfd_tmp <- tbl_fread %>%
  dplyr::mutate(filename = paste0(prefix, "_900s_", place, "-", strain, "-", sex, "-", age, "-", n_inds, "-", N),
                prefix = as.character(prefix),
                date = str_sub(prefix, start=1, end=8),
                pos_x = pos_x_wma - cent_mask_x[place],
                pos_y = 160 - (pos_y_wma - cent_mask_y[place])) %>%
  # filter(!filename %in% exclude_samples) %>%
  filter(filename %in% include_samples) %>%
  
  # convert position to (0-100, 0-100)
  group_by(place) %>%
  dplyr::mutate(pos_x_max = max(pos_x, na.rm = T),
         pos_x_min = min(pos_x, na.rm = T),
         pos_y_max = max(pos_y, na.rm = T),
         pos_y_min = min(pos_y, na.rm = T)) %>%
  ungroup() %>%
  arrange(filename, date, prefix, place, strain, sex, age, n_inds, N, id, seconds_total) %>%
  dplyr::mutate(pos_x = (pos_x - pos_x_min) * 100 / (pos_x_max - pos_x_min),
         pos_y = (pos_y - pos_y_min) * 100 / (pos_y_max - pos_y_min),
         pos_x_1framebefore = lag(pos_x),
         pos_y_1framebefore = lag(pos_y),
         travelled_dist_diff = ( (pos_x - pos_x_1framebefore)**2 + (pos_y - pos_y_1framebefore)**2 )**(1/2) ) %>%
  
  # additional variables
  dplyr::mutate(# vec_norm = (1+tan(angle_diff_based)**2)**(1/2),
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

dfd <- dfd_tmp %>%
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
write_parquet(dfd, paste0(datadir, "/2_mixed_strain/parquet/dfd_0.parquet"))


dfd <- add_samplenum_trial(dfd) %>%
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
write_parquet(dfd, paste0(datadir, "/2_mixed_strain/parquet/dfd.parquet"))
# dfd <- read_parquet(paste0(datadir, "2_mixed_strain/parquet/dfd.parquet")) %>%
#   ungroup()

## data first 5 min
dfd_f5min <- dfd %>%
  filter(seconds_total < 300)
write_parquet(dfd_f5min, "../data/2_mixed_strain/parquet/dfd_f5min.parquet")
# dfd_f5min <- read_parquet("../data/2_mixed_strain/parquet/dfd_f5min.parquet") %>%
#   ungroup()

## data first 10 min
dfd_f10min <- dfd %>%
  filter(seconds_total < 600)
write_parquet(dfd_f10min, "../data/2_mixed_strain/parquet/dfd_f10min.parquet")
# dfd_f10min <- read_parquet("../data/2_mixed_strain/parquet/dfd_f10min.parquet") %>%
#   ungroup()

## data middle 5 min
dfd_s5min <- dfd %>%
  filter(298.5 < seconds_total, seconds_total < 600)
write_parquet(dfd_s5min, "../data/2_mixed_strain/parquet/dfd_s5min.parquet")
# dfd_s5min <- read_parquet("../data/2_mixed_strain/parquet/dfd_s5min.parquet") %>%
#   ungroup()

## add stim info
dfd_s5min_stim <- dfd_s5min %>%
  dplyr::mutate(stimuli = case_when(
    as.integer(seconds_total) %% 15 == 14 & seconds_total - as.integer(seconds_total) == 0 ~ "-1.0", 
    as.integer(seconds_total) %% 15 == 0 & seconds_total - as.integer(seconds_total) == 0.5 ~ "+0.5",
    as.integer(seconds_total) %% 15 == 10 & seconds_total - as.integer(seconds_total) == 0 ~ "+10.0",)) %>%
  dplyr::filter(!is.na(stimuli)) %>%
  dplyr::mutate(n_inds = if_else(n_inds == 1, "Single", "Group"),
         time = as.POSIXct(prefix, format = format('%Y%m%d%H%M%S')) %>% format("%H:%M")) %>%
  transform(stimuli = factor(stimuli, levels=c("-1.0", "+0.5", "+10.0")),
            n_inds = factor(n_inds, levels=c("Single", "Group")),
            strain = factor(strain, levels = unique(.$strain) %>% as.character() %>% mixedsort())) %>%
  dplyr::mutate(stim_no = as.integer((seconds_total-299) %/% 15 + 1),
         posture_1s = lead(posture),
         change_posture = case_when(posture == "walk" & posture_1s == "stop" ~ "ws",
                                    posture == "walk" & posture_1s == "walk" ~ "ww",
                                    posture == "stop" & posture_1s == "stop" ~ "ss",
                                    posture == "stop" & posture_1s == "walk" ~ "sw")) %>%
  filter(seconds_total != 599)

write_parquet(dfd_s5min_stim, "../data/2_mixed_strain/parquet/dfd_s5min_stim.parquet")
# dfd_s5min_stim <- read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_stim.parquet") %>%
#   ungroup()

## s5min (from 299.5-599) for foraging success
dfd_s5min_2995_5990 <- dfd %>%
  filter(299 < seconds_total, seconds_total < 599.5)

write_parquet(dfd_s5min_2995_5990, "../data/2_mixed_strain/parquet/dfd_s5min_2995_5990.parquet")
# dfd_s5min_2995_5990 <- read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990.parquet") %>%
#   ungroup()

list_strain_d <- unique(dfd_s5min_2995_5990$strain)
write.table(list_strain_d, "../data/2_mixed_strain/list_strain_d.txt")
list_strain_d <- read.csv("../data/2_mixed_strain/list_strain_d.txt") %>%
  pull(x)

#### additional dataset ####
##### moving speed #####
###### dfd_f10min_speed_trial ######
dfd_f10min_speed_trial <- dfd_f10min %>%
  dplyr::mutate(n_inds = if_else(n_inds == 1, "Single", "Group"),
         speed = speed * 30 /100) %>%
  group_by(strain, n_inds, sex, trial, seconds_total) %>%
  dplyr::summarize(speed = mean(speed, na.rm = T)) %>%
  transform(n_inds = factor(n_inds, levels = c("Single", "Group")),
            strain = factor(strain, levels = unique(.$strain) %>% as.character() %>% mixedsort()))

write_parquet(dfd_f10min_speed_trial, "../data/2_mixed_strain/parquet/dfd_f10min_speed_trial.parquet")
dfd_f10min_speed_trial <- read_parquet("../data/2_mixed_strain/parquet/dfd_f10min_speed_trial.parquet") %>%
  ungroup()

###### dfd_f10min_speed_trial_id ######
dfd_f10min_speed_trial_id <- dfd_f10min %>%
  dplyr::mutate(n_inds = if_else(n_inds == 1, "Single", "Group"),
         speed = speed * 30 /100) %>%
  group_by(strain, n_inds, sex, trial, id, seconds_total) %>%
  dplyr::summarize(speed = mean(speed, na.rm = T)) %>%
  transform(n_inds = factor(n_inds, levels = c("Single", "Group")),
            strain = factor(strain, levels = unique(.$strain) %>% as.character() %>% mixedsort()))

write_parquet(dfd_f10min_speed_trial_id, "../data/2_mixed_strain/parquet/dfd_f10min_speed_trial_id.parquet")
dfd_f10min_speed_trial_id <- read_parquet("../data/2_mixed_strain/parquet/dfd_f10min_speed_trial_id.parquet") %>%
  ungroup()

###### df_f10min_speed_sgd_strainmean_trial ######
## calculate mean of two strains ##
df_f10min_speed_sgd_strainmean_trial <- data.frame()
for (i in list_strain_d){
  # i <- list_strain_d[1]
  strain1 <- str_split(i, "_")[[1]][1]
  strain2 <- str_split(i, "_")[[1]][2]
  
  df_f10min_speed_sgd_strainmean_trial <- bind_rows(
    df_f10min_speed_sgd_strainmean_trial,
    bind_rows(
      df_f10min_speed_trial %>%
        filter(sex == "Female", strain == eval(strain1)),
      df_f10min_speed_trial %>%
        filter(sex == "Female", strain == eval(strain2))
    ) %>%
      pivot_wider(id_cols = c(sex, n_inds, trial, seconds_total), names_from = strain, values_from = speed) %>%
      dplyr::mutate(Mean = (get(strain1) + get(strain2)) / 2,
             strain = i) %>%
      dplyr::rename(strain1_speed = strain1, 
                    strain2_speed = strain2) %>%
      dplyr::mutate(strain1 = strain1,
             strain2 = strain2)
  )
}

write_parquet(df_f10min_speed_sgd_strainmean_trial, "../data/2_mixed_strain/parquet/df_f10min_speed_sgd_strainmean_trial.parquet")
df_f10min_speed_sgd_strainmean_trial <- read_parquet("../data/2_mixed_strain/parquet/df_f10min_speed_sgd_strainmean_trial.parquet") %>%
  ungroup()

###### dfd_f10min_speed_sgd_strain_trial ######
dfd_f10min_speed_sgd_strain_trial <- 
  dfd_f10min_speed_trial %>%
  dplyr::rename(Mixed = speed) %>%
  right_join(df_f10min_speed_sgd_strainmean_trial) %>%
  pivot_longer(cols = c(Mixed, strain1_speed, strain2_speed, Mean),
               names_to = "var",
               values_to = "speed") %>%
  dplyr::mutate(var = paste0(n_inds, "_", var)) %>%
  # dplyr::mutate(n_inds2 = n_inds) %>%
  pivot_wider(id_cols = c(strain, sex, trial, seconds_total),
              names_from = var,
              values_from = speed) %>%
  dplyr::select(!Single_Mixed) %>%
  pivot_longer(cols = !c(strain, sex, trial, seconds_total),
               names_to = "var",
               values_to = "speed") %>%
  dplyr::mutate(var2 = var) %>%
  separate(var2, into = c("n_inds", "type")) %>%
  dplyr::mutate(type = if_else(str_detect(var, "strain"),
                        type,
                        "Mixed") %>%
           str_to_title(),
         mixed = if_else(str_detect(var, "Mixed"), 
                         "Mixed strain", 
                         "Single strain"),
         alpha = if_else(str_detect(var, "Mean"), 
                         "Expected", 
                         "Observed"),
         var = case_when(var == "Single_strain1_speed" ~ "Single_Strain1",
                         var == "Single_strain2_speed" ~ "Single_Strain2",
                         var == "Group_strain1_speed" ~ "Group_Strain1",
                         var == "Group_strain2_speed" ~ "Group_Strain2",
                         TRUE ~ var)) %>%
  transform(type = factor(type, levels = c("Strain1", "Mixed", "Strain2")),
            mixed = factor(mixed, levels = c("Single strain", "Mixed strain")),
            alpha = factor(alpha, levels = c("Expected", "Observed")),
            var = factor(var, levels = c("Single_Strain1", "Single_Strain2", "Single_Mean",
                                         "Group_Strain1", "Group_Strain2", "Group_Mean",
                                         "Group_Mixed")))


write_parquet(dfd_f10min_speed_sgd_strain_trial, "../data/2_mixed_strain/parquet/dfd_f10min_speed_sgd_strain_trial.parquet")
dfd_f10min_speed_sgd_strain_trial <- read_parquet("../data/2_mixed_strain/parquet/dfd_f10min_speed_sgd_strain_trial.parquet") %>%
  ungroup()


###### dfd_f10min_speed_sgd_strain ######
dfd_f10min_speed_sgd_strain <-
  dfd_f10min_speed_sgd_strain_trial %>%
  group_by(strain, sex, seconds_total, n_inds, var, type, mixed, alpha) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  transform(strain = factor(strain, levels = unique(.$strain) %>% as.character() %>% mixedsort()))
write_parquet(dfd_f10min_speed_sgd_strain, "../data/2_mixed_strain/parquet/dfd_f10min_speed_sgd_strain.parquet")
dfd_f10min_speed_sgd_strain <- read_parquet("../data/2_mixed_strain/parquet/dfd_f10min_speed_sgd_strain.parquet") %>%
  ungroup()


###### dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain ######
dfd_f5min_speed_sgd_ave_strain <- 
  dfd_f10min_speed_sgd_strain %>%
  filter(seconds_total < 300) %>%
  group_by(sex, strain, n_inds, var, type, mixed, alpha) %>%
  dplyr::summarize(speed_f5min_ave = mean(speed, na.rm = T))

dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain <- 
  dfd_f10min_speed_sgd_strain %>%
  filter(299 < seconds_total, seconds_total < 599.5) %>%
  dplyr::mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5) %>%
  group_by(strain, sex, trial, stim_time, n_inds, var, type, mixed, alpha) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  inner_join(dfd_f5min_speed_sgd_ave_strain) %>%
  dplyr::mutate(speed_norm_mean = speed / speed_f5min_ave) %>%
  group_by(strain, stim_time, sex, var, n_inds, type, mixed, alpha) %>%
  dplyr::summarize(speed_norm_mean = mean(speed_norm_mean, na.rm = TRUE),
                   speed = mean(speed, na.rm = TRUE),
                   speed_f5min_ave = mean(speed_f5min_ave, na.rm = TRUE)) %>%
  ungroup() %>%
  transform(var = factor(var, levels = c("Single_Strain1", "Single_Strain2", "Single_Mean",
                                         "Group_Strain1", "Group_Strain2", "Group_Mean",
                                         "Group_Mixed"))) %>%
  # re-calculate the speed_normbyf5min by getting mean of 2 strains
  dplyr::mutate(speed_norm_mean2 = (lag(speed_norm_mean, 1) + dplyr::lag(speed_norm_mean, 2)) / 2,
         speed_norm_mean = case_when(str_detect(var, "Mean") ~ speed_norm_mean2,
                                     TRUE ~ speed_norm_mean)) %>%
  dplyr::select(!speed_norm_mean2)

write_parquet(dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain, "../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain.parquet")
dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain <- read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain.parquet") %>%
  ungroup()


###### dfd_f5min_speed_sgd_ave_strain_trial ######
dfd_f5min_speed_sgd_ave_strain_trial <- 
  dfd_f10min_speed_sgd_strain_trial %>%
  filter(seconds_total < 300) %>%
  group_by(strain, sex, n_inds, trial, var, type, mixed, alpha) %>%
  dplyr::summarize(speed_f5min_ave = mean(speed, na.rm = T))

###### dfd_s5min_2995_5990_speed_sgd_strain_trial ######
dfd_s5min_2995_5990_speed_sgd_strain_trial <- 
  dfd_f10min_speed_sgd_strain_trial %>%
  filter(299 < seconds_total, seconds_total < 599.5) %>%
  dplyr::mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5) %>%
  group_by(strain, sex, trial, stim_time, n_inds, var, type, mixed, alpha) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE))
write_parquet(dfd_s5min_2995_5990_speed_sgd_strain_trial, "../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_sgd_strain_trial.parquet")
dfd_s5min_2995_5990_speed_sgd_strain_trial <- read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_sgd_strain_trial.parquet") %>%
  ungroup()


###### dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial ######
dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial <- 
  dfd_s5min_2995_5990_speed_sgd_strain_trial %>%
  inner_join(dfd_f5min_speed_sgd_ave_strain_trial) %>%
  dplyr::mutate(speed_norm_mean = speed / speed_f5min_ave) %>%
  group_by(strain, sex, trial, stim_time, var, n_inds, type, mixed, alpha) %>%
  dplyr::summarize(speed_norm_mean = mean(speed_norm_mean, na.rm = TRUE),
                   speed = mean(speed, na.rm = TRUE),
                   speed_f5min_ave = mean(speed_f5min_ave, na.rm = TRUE)) %>%
  ungroup() %>%
  transform(var = factor(var, levels = c("Single_Strain1", "Single_Strain2", "Single_Mean",
                                         "Group_Strain1", "Group_Strain2", "Group_Mean",
                                         "Group_Mixed"))) %>%
  # re-calculate the speed_normbyf5min by getting mean of 2 strains
  dplyr::mutate(speed_norm_mean2 = (lag(speed_norm_mean, 1) + dplyr::lag(speed_norm_mean, 2)) / 2,
         speed_norm_mean = case_when(str_detect(var, "Mean") ~ speed_norm_mean2,
                                     TRUE ~ speed_norm_mean)) %>%
  dplyr::select(!speed_norm_mean2)
write_parquet(dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial, "../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial.parquet")
dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial <- read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial.parquet") %>%
  ungroup()

###### speed_normbystimminus05 ######
dfd_s5min_2995_5990_speed_sgd_stimminus05_ave_strain_trial <-
  dfd_f10min_speed_sgd_strain_trial %>%
  filter(299 < seconds_total, seconds_total < 599.5,
         (seconds_total + 0.5) %% 15 == 0) %>%
  group_by(strain, trial, sex, n_inds, var, type, mixed, alpha) %>%
  dplyr::summarize(speed_stimminus05_ave = mean(speed, na.rm = TRUE))

dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial <-
  dfd_f10min_speed_sgd_strain_trial %>%
  filter(299 < seconds_total, seconds_total < 599.5) %>%
  dplyr::mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5) %>%
  inner_join(dfd_s5min_2995_5990_speed_sgd_stimminus05_ave_strain_trial) %>%
  dplyr::mutate(speed_normbystimminus05 = speed / speed_stimminus05_ave) %>%
  dplyr::mutate(speed_normbystimminus05_2 = (lag(speed_normbystimminus05, 1) + lag(speed_normbystimminus05, 2)) / 2,
         speed_normbystimminus05 = case_when(str_detect(var, "Mean") ~ speed_normbystimminus05_2,
                                             TRUE ~ speed_normbystimminus05)) %>%
  dplyr::select(!speed_normbystimminus05_2)
write_parquet(dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial, "../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial.parquet")
dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial <- read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial.parquet") %>%
  ungroup()


#### behavioral performance (balance between vigilance and exploration) ####
##### not separate freezing / exploration #####
t_test_1sample_performance <- function(dat){
  # dat <- dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance %>%
  #   filter(strain == "DGRP88_DGRP101")
  pval <-
    t.test(dat %>%
             filter(alpha == "Observed") %>%
             pull(performance),
           mu = dat %>%
             filter(alpha == "Expected") %>%
             pull(performance) %>%
             mean()
    )$p.value
  return(pval)
}

dfd_fig4 <- data.frame()
dfd_figS14 <- data.frame()
for (for_sec in c(0.5, 1, 2, 3, 4, 5, 6, 12)){
  # for_sec <- 3
  a2 <- for_sec / 0.5
  a <- (28 - a2) / a2
  
  ###### s5min_2995_5990_speed_normbyf5minave ######
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance <- 
    dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial %>%
    filter(n_inds == "Group", type == "Mixed", stim_time > 0, !str_detect(strain, "norpA")) %>%
    # dplyr::mutate(performance = if_else(stim_time <= for_sec, speed * (-a), speed)) %>% #raw speed
    dplyr::mutate(performance = if_else(stim_time <= for_sec, speed_norm_mean * (-a), speed_norm_mean)) %>% #normbyf5min
    group_by(strain, sex, n_inds, trial, var, type, mixed, alpha) %>%
    dplyr::summarize(performance = sum(performance))
  
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance2 <- 
    dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial %>%
    filter(n_inds == "Group", stim_time > 0, !str_detect(strain, "norpA")) %>%
    # dplyr::mutate(performance = if_else(stim_time <= for_sec, speed * (-a), speed)) %>% #raw speed
    dplyr::mutate(performance = if_else(stim_time <= for_sec, speed_norm_mean * (-a), speed_norm_mean)) %>% #normbyf5min
    group_by(strain, sex, n_inds, trial, var, type, mixed, alpha) %>%
    dplyr::summarize(performance = sum(performance))
  
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance_tr <-
    dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance2 %>%
    filter(n_inds == "Group") %>%
    group_by(strain, sex, n_inds, var, type, mixed, alpha) %>%
    dplyr::summarize(performance = mean(performance, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::mutate(strain1 = dplyr::lag(performance, 3),
           strain2 = dplyr::lag(performance, 2),
           mean_2strains =  dplyr::lag(performance, 1)) %>%
    filter(var == "Group_Mixed") %>%
    dplyr::mutate(overyield = performance - mean_2strains,
           tr_overyield = if_else(strain1 > strain2, performance - strain1, performance - strain2),
           tr_underyield = if_else(strain1 > strain2, performance - strain2, performance - strain1),
           col = case_when(tr_overyield > 0 ~ "Transgressive overyielding",
                           tr_overyield <= 0 & overyield > 0 ~ "Overyielding",
                           tr_underyield >= 0 & overyield < 0 ~ "Underyielding",
                           tr_underyield < 0 ~ "Transgressive underyielding")) %>%
    transform(col = factor(col, levels = c("Transgressive overyielding", "Overyielding", "Underyielding", "Transgressive underyielding")))
  
  # write.table(dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance_tr,
  #             paste0("/Users/dsato/Dropbox/研究室/投稿論文/2021/Drosophila/2023/data/2_mixed_strain/gwas/pheno/dfd_s5min_2995_5990_speed_normbyf5minave_performance_sec", for_sec, ".tsv"), sep = "\t", row.names = F, quote = F)
  
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance_t_test <-
    dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance %>%
    group_by(strain) %>%
    nest() %>%
    dplyr::mutate(t.test.p.value = map(data, t_test_1sample_performance)) %>%
    unnest(data) %>%
    dplyr::mutate(t.test.p.value = as.numeric(t.test.p.value)) %>%
    dplyr::select(strain, t.test.p.value) %>%
    distinct() %>%
    dplyr::mutate(col2 = if_else(t.test.p.value < 0.05, "sig", "insig"))
  
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance3 <-
    dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance2 %>%
    group_by(strain, sex, n_inds, var, type, mixed, alpha) %>%
    dplyr::summarize(mean = mean(performance, na.rm = T),
                     error = sd(performance, na.rm = T)) %>%
    left_join(dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance_tr %>%
                dplyr::select(!c(var, type, mixed, alpha))) %>%
    left_join(dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance_t_test %>%
                dplyr::select(strain, t.test.p.value, col2)) %>%
    transform(strain = factor(strain, levels = unique(.$strain) %>% gtools::mixedsort())) %>%
    dplyr::mutate(thr_sec = for_sec)
  
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance4 <- 
    dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance2 %>%
    filter(n_inds == "Group", type == "Mixed") %>%
    group_by(strain, sex, n_inds, var, type, mixed, alpha) %>%
    dplyr::summarize(performance = mean(performance, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::mutate(thr_sec = for_sec)
  
  # write_parquet(dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance4,
  #               "../../data/2_mixed_strain/tracked_data/dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance4.parquet")
  
  dfd_fig4 <- bind_rows(dfd_fig4,
                        dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance4)
  
  dfd_figS14 <- bind_rows(dfd_figS14,
                         dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance3)
}

###### save dataset ######
write.table(dfd_fig4, "../data/2_mixed_strain/dfd_fig4.tsv", 
            sep = "\t", row.names = F, quote = F)
# dfd_fig4 <- read.delim("../data/2_mixed_strain/dfd_fig4.tsv")
write.table(dfd_figS14, "../data/2_mixed_strain/dfd_figS14.tsv", 
            sep = "\t", row.names = F, quote = F)
# dfd_figS14 <- read.delim("../data/2_mixed_strain/dfd_figS14.tsv")

for_sec <- 3
dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance_sec <-
  dfd_fig4 %>%
  filter(thr_sec == for_sec) %>%
  dplyr::mutate(DE = performance - lag(performance)) %>%
  filter(alpha == "Observed")

write.table(dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance_sec, 
            paste0("../data/2_mixed_strain/gwas/pheno/dfd_s5min_2995_5990_speed_normbyf5minave_performance_sec", for_sec, ".tsv"), sep = "\t", row.names = F, quote = F)



##### separate freezing / exploration #####
t_test_1sample_freezing_performance <- function(dat){
  # dat <- dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance %>%
  #   filter(strain == "DGRP88_DGRP101")
  pval <-
    t.test(dat %>%
             filter(alpha == "Observed") %>%
             pull(speed),
           mu = dat %>%
             filter(alpha == "Expected") %>%
             pull(speed) %>%
             mean()
    )$p.value
  return(pval)
}

dfd_fig4_2 <- data.frame()
for (for_sec in c(0.5, 1, 2, 3, 4, 5, 6, 12)){
  # for_sec <- 3
  # a2 <- for_sec / 0.5
  # a <- (29 - a2) / a2
  
  ###### s5min_2995_5990_speed_normbyf5minave ######
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance <- 
    dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial %>%
    filter(n_inds == "Group", type == "Mixed", stim_time > 0, !str_detect(strain, "norpA")) %>%
    dplyr::mutate(time = if_else(stim_time <= for_sec, "Before", "After")) %>%
    group_by(strain, sex, n_inds, trial, var, type, mixed, alpha, time) %>%
    dplyr::summarize(speed = mean(speed_norm_mean, na.rm = TRUE)) #normbyf5min
    # dplyr::summarize(speed = mean(speed, na.rm = TRUE)) #raw
  
  
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance2 <- 
    dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial %>%
    filter(n_inds == "Group", stim_time > 0, !str_detect(strain, "norpA")) %>%
    dplyr::mutate(time = if_else(stim_time <= for_sec, "Before", "After")) %>%
    group_by(strain, sex, n_inds, trial, var, type, mixed, alpha, time) %>%
    dplyr::summarize(speed = mean(speed_norm_mean, na.rm = TRUE)) #normbyf5min
    # dplyr::summarize(speed = mean(speed, na.rm = TRUE)) #raw
  
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance_tr <-
    dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance2 %>%
    # filter(n_inds == "Group") %>%
    group_by(strain, sex, n_inds, trial, var, type, mixed, alpha, time) %>%
    dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::mutate(strain1 = dplyr::lag(speed, 6),
           strain2 = dplyr::lag(speed, 4),
           mean_2strains =  dplyr::lag(speed, 2)) %>%
    filter(var == "Group_Mixed") %>%
    dplyr::mutate(overyield = speed - mean_2strains,
           tr_overyield = if_else(strain1 > strain2, speed - strain1, speed - strain2),
           tr_underyield = if_else(strain1 > strain2, speed - strain2, speed - strain1),
           col = case_when(tr_overyield > 0 ~ "Transgressive overyielding",
                           tr_overyield <= 0 & overyield > 0 ~ "Overyielding",
                           tr_underyield >= 0 & overyield < 0 ~ "Underyielding",
                           tr_underyield < 0 ~ "Transgressive underyielding")) %>%
    transform(col = factor(col, levels = c("Transgressive overyielding", "Overyielding", "Underyielding", "Transgressive underyielding")))
  
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance_t_test <-
    dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance %>%
    group_by(strain, time) %>%
    nest() %>%
    dplyr::mutate(t.test.p.value = map(data, t_test_1sample_freezing_performance)) %>%
    unnest(data) %>%
    dplyr::mutate(t.test.p.value = as.numeric(t.test.p.value)) %>%
    dplyr::select(strain, t.test.p.value) %>%
    distinct() %>%
    dplyr::mutate(col2 = if_else(t.test.p.value < 0.05, "sig", "insig"))
  
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance3 <- 
    dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance2 %>%
    filter(n_inds == "Group", type == "Mixed") %>%
    group_by(strain, sex, n_inds, var, type, mixed, alpha, time) %>%
    dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
    ungroup() %>%
    dplyr::mutate(thr_sec = for_sec)
  
  dfd_fig4_2 <- bind_rows(dfd_fig4_2,
                          dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial_performance3)
}
###### save dataset ######
write.table(dfd_fig4_2, "../data/2_mixed_strain/dfd_fig4_2.tsv", 
            sep = "\t", row.names = F, quote = F)
# dfd_fig4_2 <- read.delim("../data/2_mixed_strain/dfd_fig4_2.tsv")


#### correlation with distance of freezing duration ####
dfd_s5min_2995_5990_speed <- dfd_s5min_2995_5990 %>%
  dplyr::mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5,
         speed = speed * 30 / 100,
         n_inds = if_else(n_inds == 1, "Single", "Group"),
         time = as.POSIXct(prefix, format = format('%Y%m%d%H%M%S')) %>% format("%H:%M")) %>%
  group_by(strain, sex, n_inds, trial, id, stim_time, date, time, filename, prefix, place) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  transform(n_inds = factor(n_inds, levels = c("Single", "Group")),
            strain = factor(strain, levels = unique(.$strain) %>% as.character() %>% gtools::mixedsort()))

write_parquet(dfd_s5min_2995_5990_speed, "../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed.parquet")
# dfd_s5min_2995_5990_speed <- read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed.parquet") %>%
#   ungroup()


dfd_s5min_2995_5990_speed_d_normbyf5minave_strain_trial <-
  dfd_s5min_2995_5990_speed %>%
  group_by(strain, sex, n_inds, trial, id, stim_time, date, time, prefix, place) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  left_join(dfd_f5min_speed %>%
              dplyr::rename(speed_f5min_ave = speed), 
            by = c("strain", "n_inds", "sex", "trial", "id", "date", "time", "prefix", "place")) %>%
  group_by(strain, n_inds, sex, trial, id, stim_time, date, time, prefix, place) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE),
                   speed_f5min_ave = mean(speed_f5min_ave, na.rm = TRUE)) %>%
  dplyr::mutate(speed_normbyf5minave = speed / speed_f5min_ave)

write_parquet(dfd_s5min_2995_5990_speed_d_normbyf5minave_strain_trial, "../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_d_normbyf5minave_strain_trial.parquet")
# dfd_s5min_2995_5990_speed_d_normbyf5minave_strain_trial <- 
#   read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_d_normbyf5minave_strain_trial.parquet") %>%
#   ungroup()


dfd_s5min_2995_5990_d_freezing_duration <-
  dfd_s5min_2995_5990_speed_d_normbyf5minave_strain_trial %>% #dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial #dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial
  filter(stim_time != -0.5, speed_normbyf5minave > 1) %>% #speed_norm_mean #speed_normbystimminus05
  group_by(strain, sex, n_inds, trial, id) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(strain, sex, n_inds, id),
              names_from = trial, values_from = stim_time,
              values_fill = 14) %>%
  pivot_longer(cols = !c(strain, sex, n_inds, id), names_to = "trial", values_to = "freezing_duration") %>%
  dplyr::mutate(trial = as.numeric(trial)) %>%
  group_by(strain, n_inds, sex, trial) %>%
  dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm = T)) %>%
  ungroup() %>%
  arrange(strain, sex, n_inds, trial) %>%
  left_join(dfd_s5min_2995_5990_speed_d_normbyf5minave_strain_trial %>%
              group_by(strain, sex, n_inds, trial) %>%
              dplyr::slice(1) %>%
              dplyr::select(strain, sex, n_inds, trial, date, time, prefix, place))
write_parquet(df_s5min_2995_5990_freezing_duration, "../data/2_mixed_strain/parquet/df_s5min_2995_5990_freezing_duration.parquet")
# df_s5min_2995_5990_freezing_duration <- read_parquet("../data/1_single_strain/parquet/df_s5min_2995_5990_freezing_duration.parquet") %>%
#   ungroup()

## calculate mean of two strains ##
dfd_s5min_2995_5990_g_freezing_duration <- data.frame()
for (i in list_strain_d){
  # i <- list_strain_d[1]
  strain1 <- str_split(i, "_")[[1]][1]
  strain2 <- str_split(i, "_")[[1]][2]
  
  dfd_s5min_2995_5990_g_freezing_duration <- bind_rows(
    dfd_s5min_2995_5990_g_freezing_duration,
    bind_rows(
      df_s5min_2995_5990_freezing_duration %>%
        filter(sex == "Female", n_inds == "Group", strain == eval(strain1)),
      df_s5min_2995_5990_freezing_duration %>%
        filter(sex == "Female", n_inds == "Group", strain == eval(strain2))
    ) %>%
      pivot_wider(id_cols = c(sex, n_inds, trial), 
                  names_from = strain, values_from = freezing_duration) %>%
      dplyr::mutate(Mean = (get(strain1) + get(strain2)) / 2,
             strain = i) %>%
      dplyr::rename(strain1_freezing_duration = strain1, 
                    strain2_freezing_duration = strain2) %>%
      dplyr::mutate(strain1 = strain1,
             strain2 = strain2)
  )
}
dfd_s5min_2995_5990_gd_freezing_duration <-
  dfd_s5min_2995_5990_d_freezing_duration %>%
  group_by(strain, sex, n_inds) %>%
  dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm = TRUE)) %>%
  dplyr::mutate(alpha = "Observed") %>%
  bind_rows(dfd_s5min_2995_5990_g_freezing_duration %>%
              group_by(strain, sex, n_inds) %>%
              dplyr::summarize(freezing_duration = mean(Mean, na.rm = TRUE)) %>%
              dplyr::mutate(alpha = "Expected")) %>%
  arrange(strain, sex, n_inds, alpha)


dfd_s5min_2995_5990_gd_freezing_duration_dist <-
  dfd_figS14 %>% 
  filter(type == "Mixed", alpha == "Observed") %>%
  left_join(dfd_s5min_2995_5990_g_freezing_duration %>%
              dplyr::mutate(freezing_duration_dist = abs(strain1_freezing_duration - strain2_freezing_duration)) %>%
              group_by(strain) %>%
              summarize(freezing_duration_dist = mean(freezing_duration_dist, na.rm = TRUE))
  ) %>%
  dplyr::mutate(thr_sec = paste0("Threshold: ", thr_sec, " s")) %>%
  transform(thr_sec = factor(thr_sec, levels = unique(.$thr_sec) %>% gtools::mixedsort()))
write_parquet(dfd_s5min_2995_5990_gd_freezing_duration_dist, "../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_gd_freezing_duration_dist.parquet")



#### correlation with distance of visual responsiveness ####
df_motion_cue_exit_coeff <- read_parquet("../data/1_single_strain/parquet/df_motion_cue_exit_coeff.parquet") %>%
  ungroup()

dfd_s5min_2995_5990_g_motion_cue_exit_intercept <- data.frame()
for (i in list_strain_d){
  # i <- list_strain_d[1]
  strain1 <- str_split(i, "_")[[1]][1]
  strain2 <- str_split(i, "_")[[1]][2]
  
  dfd_s5min_2995_5990_g_motion_cue_exit_intercept <- bind_rows(
    dfd_s5min_2995_5990_g_motion_cue_exit_intercept,
    bind_rows(
      df_motion_cue_exit_coeff %>%
        filter(sex == "Female", strain == eval(strain1)),
      df_motion_cue_exit_coeff %>%
        filter(sex == "Female", strain == eval(strain2))
    ) %>%
      pivot_wider(id_cols = sex, 
                  names_from = strain, values_from = motion_cue_exit_intercept) %>%
      dplyr::mutate(strain = i) %>%
      dplyr::rename(strain1_motion_cue_exit_intercept = strain1, 
                    strain2_motion_cue_exit_intercept = strain2) %>%
      dplyr::mutate(strain1 = strain1,
             strain2 = strain2)
  )
}

dfd_s5min_2995_5990_gd_motion_cue_exit_intercept_dist <-
  dfd_figS14 %>% 
  filter(type == "Mixed", alpha == "Observed") %>%
  left_join(dfd_s5min_2995_5990_g_motion_cue_exit_intercept %>%
              dplyr::mutate(motion_cue_exit_intercept_dist = abs(strain1_motion_cue_exit_intercept - strain2_motion_cue_exit_intercept)) %>%
              group_by(strain) %>%
              summarize(motion_cue_exit_intercept_dist = mean(motion_cue_exit_intercept_dist, na.rm = TRUE))
  ) %>%
  dplyr::mutate(thr_sec = paste0("Threshold: ", thr_sec, " s")) %>%
  transform(thr_sec = factor(thr_sec, levels = unique(.$thr_sec) %>% gtools::mixedsort()))
write_parquet(dfd_s5min_2995_5990_gd_motion_cue_exit_intercept_dist, "../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_gd_motion_cue_exit_intercept_dist.parquet")


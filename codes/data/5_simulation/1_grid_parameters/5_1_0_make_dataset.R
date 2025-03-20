#### load packages ####
targetPackages <- c('tidyverse','data.table','arrow')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### parameters ####
args <- commandArgs(trailingOnly=TRUE)
#print(args)
type <- args[1]
#print(type)
cond <- setdiff(args, type) %>%
  paste(collapse = "_")
#print(cond)
dt <- 0.02
nrep <- as.integer(args[2])
#print(nrep)

#### load dataset ####
resdir <- paste0("../data/5_simulation/1_grid_parameters/raw", type, "_", cond)
files <- list.files(resdir, pattern = "*pos*", full.names = TRUE, recursive=T)
names(files) <- str_remove(basename(files), "\\.tsv") %>% parse_number()
df_pos <- map_dfr(files, fread, .id = 'rep') %>%
  dplyr::mutate(rep = as.numeric(rep),
                id = as.factor(id),
                seconds_diff = dt)

df_pos2 <- df_pos %>%
  dplyr::arrange(rep, id, frame) %>%
  dplyr::mutate(
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
  dplyr::mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5)

df_pos3_speed_beforestimave <- df_pos2 %>%
  filter(seconds_total < 300) %>%
  group_by(rep, id) %>%
  dplyr::summarize(speed_beforestimave = mean(speed, na.rm = TRUE)) %>%
  dplyr::select(rep, id, speed_beforestimave)

df_pos3_speed_normbybeforestimave <- df_pos2 %>%
  inner_join(df_pos3_speed_beforestimave, by = c("rep","id")) %>%
  dplyr::mutate(speed_normbybeforestimave = speed / speed_beforestimave) %>%
  filter(seconds_total >= 299.5 & seconds_total < 599.5) %>%
  group_by(stim_time, rep) %>%
  dplyr::summarize(speed = mean(speed_normbybeforestimave, na.rm = TRUE))

newdir <- paste0("../data/5_simulation/1_grid_parameters/parquet", type)
dir.create(newdir, showWarnings = FALSE)
filename <- paste0(newdir, "/df_", type, "_", cond, ".parquet")
write_parquet(df_pos3_speed_normbybeforestimave, filename)
unlink(resdir, recursive=TRUE)

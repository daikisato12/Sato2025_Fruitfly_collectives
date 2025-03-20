#### load packages ####
targetPackages <- c('tidyverse','arrow','ggtext','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### load functions ####
log_func <- function(t, a, b){
  p <- if_else(1 - (a*log(t) + b) < 0, 0, 1 - (a*log(t) + b))
  return(p)
}

#### Figure S6a ####
##### load dataset #####
log_param_a = c(0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
log_param_b = 1.02 - 2.45 * log_param_a
log_param <- data.frame(
  a = log_param_a,
  b = log_param_b
) %>%
  dplyr::mutate(name = paste0("a = ", a, ", b = ", b))

##### make plot #####
g_parameter <- 
  ggplot(data.frame(X = c(0, 15)), 
         aes(x = X)) +
  mapply(
    function(a, b, co) stat_function(fun=log_func, args=list(a, b), 
                                     aes(color = co)),
    log_param_a, log_param_b, as.factor(log_param_a)
  ) + 
  # annotate(geom = "text", x = 9.5, y = 0.38,
  #          label = label) +
  # scale_color_discrete(name = paste("Freezing parameter a")) +
  ggtext::geom_richtext(x = 9.2, y = 0.39, 
                family = "Helvetica",
                size = 3,
                fill = NA,
                label.color = NA, 
                label = "Freezing parameter *a*") +
  xlab("Time after stimulus (s)") +
  ylab("Probability to freeze") +
  theme_bw() +
  theme(legend.key.size = unit(0.4, "cm"),
        # legend.margin = margin(unit(1, "cm")),
        legend.position = c(0.6, 0.7),
        legend.title = element_blank(),
        legend.background = element_blank()) +
  guides(col = guide_legend(ncol = 3))
g_parameter
ggsave("../figures/FigureS6a.pdf", g_parameter, w = 3, h = 3)


#### Figure S6b ####
##### load dataset #####
###### single ######
df_single_pos2 <- read_parquet("../data/5_simulation/2_approx_parameters/parquet/df_single_pos.parquet")

df_single_pos3_speed_beforestimave <- df_single_pos2 %>%
  filter(seconds_total < 300) %>%
  group_by(rep, id) %>%
  dplyr::summarize(speed_beforestimave = mean(speed, na.rm = TRUE)) %>%
  dplyr::select(rep, id, speed_beforestimave)

df_single_pos3_speed_normbybeforestimave <- df_single_pos2 %>%
  inner_join(df_single_pos3_speed_beforestimave, by = c("rep","id")) %>%
  dplyr::mutate(speed_normbybeforestimave = speed / speed_beforestimave)

df_single_pos3_speed_stimminus05ave <- df_single_pos2 %>%
  filter(seconds_total >= 299.5, seconds_total < 599.5, stim_time == -0.5) %>%
  group_by(rep, id) %>%
  dplyr::summarize(speed_stimminus05ave = mean(speed, na.rm = TRUE)) %>%
  dplyr::select(rep, id, speed_stimminus05ave)

df_single_pos3_speed_normbystimminus05ave <- df_single_pos2 %>%
  filter(seconds_total >= 299.5, seconds_total < 599.5) %>%
  inner_join(df_single_pos3_speed_stimminus05ave, by = c("rep","id")) %>%
  dplyr::mutate(speed_normbystimminus05ave = speed / speed_stimminus05ave) %>%
  group_by(rep, stim_time, id) %>%
  dplyr::summarize(speed_normbystimminus05ave = mean(speed_normbystimminus05ave, na.rm = TRUE))


###### group w/o interaction ######
df_group_singlemove_pos2 <- read_parquet("../data/5_simulation/2_approx_parameters/parquet/df_group_singlemove_pos.parquet")

df_group_singlemove_pos3_speed_beforestimave <- df_group_singlemove_pos2 %>%
  filter(seconds_total < 300) %>%
  group_by(rep, id) %>%
  dplyr::summarize(speed_beforestimave = mean(speed, na.rm = TRUE)) %>%
  dplyr::select(rep, id, speed_beforestimave)

df_group_singlemove_pos3_speed_normbybeforestimave <- df_group_singlemove_pos2 %>%
  inner_join(df_group_singlemove_pos3_speed_beforestimave, by = c("rep","id")) %>%
  dplyr::mutate(speed_normbybeforestimave = speed / speed_beforestimave)

df_group_singlemove_pos3_speed_stimminus05ave <- df_group_singlemove_pos2 %>%
  filter(seconds_total >= 299.5, seconds_total < 599.5, stim_time == -0.5) %>%
  group_by(rep, id) %>%
  dplyr::summarize(speed_stimminus05ave = mean(speed, na.rm = TRUE)) %>%
  dplyr::select(rep, id, speed_stimminus05ave)

df_group_singlemove_pos3_speed_normbystimminus05ave <- df_group_singlemove_pos2 %>%
  filter(seconds_total >= 299.5, seconds_total < 599.5) %>%
  inner_join(df_group_singlemove_pos3_speed_stimminus05ave, by = c("rep","id")) %>%
  dplyr::mutate(speed_normbystimminus05ave = speed / speed_stimminus05ave) %>%
  group_by(rep, stim_time, id) %>%
  dplyr::summarize(speed_normbystimminus05ave = mean(speed_normbystimminus05ave, na.rm = TRUE))


###### group w interaction ######
df_group_socialcontagion_pos2 <- read_parquet("../data/5_simulation/2_approx_parameters/parquet/df_group_socialcontagion_pos.parquet")

df_group_socialcontagion_pos3_speed_beforestimave <- df_group_socialcontagion_pos2 %>%
  filter(seconds_total < 300) %>%
  group_by(rep, id) %>%
  dplyr::summarize(speed_beforestimave = mean(speed, na.rm = TRUE)) %>%
  dplyr::select(rep, id, speed_beforestimave)

df_group_socialcontagion_pos3_speed_normbybeforestimave <- df_group_socialcontagion_pos2 %>%
  inner_join(df_group_socialcontagion_pos3_speed_beforestimave, by = c("rep","id")) %>%
  dplyr::mutate(speed_normbybeforestimave = speed / speed_beforestimave)

df_group_socialcontagion_pos3_speed_stimminus05ave <- df_group_socialcontagion_pos2 %>%
  filter(seconds_total >= 299.5, seconds_total < 599.5, stim_time == -0.5) %>%
  group_by(rep, id) %>%
  dplyr::summarize(speed_stimminus05ave = mean(speed, na.rm = TRUE)) %>%
  dplyr::select(rep, id, speed_stimminus05ave)

df_group_socialcontagion_pos3_speed_normbystimminus05ave <- df_group_socialcontagion_pos2 %>%
  filter(seconds_total >= 299.5, seconds_total < 599.5) %>%
  inner_join(df_group_socialcontagion_pos3_speed_stimminus05ave, by = c("rep","id")) %>%
  dplyr::mutate(speed_normbystimminus05ave = speed / speed_stimminus05ave) %>%
  group_by(rep, stim_time, id) %>%
  dplyr::summarize(speed_normbystimminus05ave = mean(speed_normbystimminus05ave, na.rm = TRUE))


###### merge single and group ######
df_group_merge_speed_normbybeforestimave <- 
  inner_join(
    df_group_socialcontagion_pos3_speed_normbybeforestimave %>%
      group_by(seconds_total, stim_time) %>%
      dplyr::summarize(`Group w interaction` = mean(speed_normbybeforestimave, na.rm = TRUE)),
    df_single_pos3_speed_normbybeforestimave %>%
      group_by(seconds_total, stim_time) %>%
      dplyr::summarize(`Single` = mean(speed_normbybeforestimave, na.rm = TRUE))) %>%
  inner_join(
    df_group_singlemove_pos3_speed_normbybeforestimave %>%
      group_by(seconds_total, stim_time) %>%
      dplyr::summarize(`Group w/o interaction` = mean(speed_normbybeforestimave, na.rm = TRUE))) %>%
  pivot_longer(cols = !c(seconds_total, stim_time), names_to = "var", values_to = "speed") %>%
  dplyr::mutate(alpha = if_else(str_detect(var, "w/o"), 0.5, 1),
                color = if_else(str_detect(var, "Single"), "Single", "Group")) %>%
  transform(var = factor(var, levels = c("Single", "Group w/o interaction", "Group w interaction")),
            color = factor(color, levels = c("Single", "Group")))

##### make plot #####
g_group_merge_speed_normbybeforestimave <- 
  ggplot(df_group_merge_speed_normbybeforestimave,
         aes(x = seconds_total, y = speed)) +
  geom_vline(xintercept = 300, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 315, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 330, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 345, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 360, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 375, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 390, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 405, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 420, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 435, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 450, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 465, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 480, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 495, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 510, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 525, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 540, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 555, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 570, linetype = "dotted", col = "grey") +
  geom_vline(xintercept = 585, linetype = "dotted", col = "grey") +
  geom_line(aes(col = color, alpha = alpha, group = var)) +
  xlab("Time (min)") +
  ylab("Mean moving speed relative to\n the average during the initial 5 min") +
  scale_color_manual(values = c("#A2A0A2", "#A93439")) +#values = viridis(3)[1:2]) +
  scale_alpha_identity() +
  scale_x_continuous(breaks = seq(0, 600, by = 60), labels = seq(0, 10, by = 1)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.16, 0.2),
        legend.box = "horizontal",
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank())
g_group_merge_speed_normbybeforestimave

df_sim_singlemove_s5min_2995_5990_speed_sg_normbyf5minave_allmean_pairedwilcox_res <- 
  bind_rows(
    df_single_pos3_speed_normbybeforestimave %>%
      filter(seconds_total >= 299.5 & seconds_total < 599.5) %>%
      group_by(stim_time, rep) %>%
      dplyr::summarize(speed = mean(speed_normbybeforestimave, na.rm = TRUE)) %>%
      dplyr::mutate(type = "Single"),
    df_group_socialcontagion_pos3_speed_normbybeforestimave %>%
      filter(seconds_total >= 299.5 & seconds_total < 599.5) %>%
      group_by(stim_time, rep) %>%
      dplyr::summarize(speed = mean(speed_normbybeforestimave, na.rm = TRUE)) %>%
      dplyr::mutate(type = "Group w interaction")) %>%
  bind_rows(
    df_group_singlemove_pos3_speed_normbybeforestimave %>%
      filter(seconds_total >= 299.5 & seconds_total < 599.5) %>%
      group_by(stim_time, rep) %>%
      dplyr::summarize(speed = mean(speed_normbybeforestimave, na.rm = TRUE)) %>%
      dplyr::mutate(type = "Group w/o interaction")) %>%
  group_by(stim_time) %>%
  rstatix::wilcox_test(speed ~ type, paired = FALSE) %>%
  rstatix::adjust_pvalue(method = "bonferroni") %>%
  rstatix::add_significance("p.adj")

df_sim_singlemove_s5min_2995_5990_speed_sg_normbyf5minave_allmean_pairedwilcox_res <- 
  replace(df_sim_singlemove_s5min_2995_5990_speed_sg_normbyf5minave_allmean_pairedwilcox_res,
          df_sim_singlemove_s5min_2995_5990_speed_sg_normbyf5minave_allmean_pairedwilcox_res=="ns",
          NA)

g_group_merge_speed_normbybeforestimave_s5min <- 
  ggplot(df_group_merge_speed_normbybeforestimave %>%
           filter(seconds_total >= 299.5, seconds_total < 599.5) %>%
           group_by(stim_time, var, alpha, color) %>%
           dplyr::summarize(speed = mean(speed, na.rm = TRUE)),
         aes(x = stim_time, y = speed)) +
  geom_line(aes(col = color, alpha = alpha, group = var), linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = -Inf, ymax = Inf, alpha = 0.4) +
  annotate(geom = "text",
           x = df_sim_singlemove_s5min_2995_5990_speed_sg_normbyf5minave_allmean_pairedwilcox_res %>%
             filter(group2 == "Group w/o interaction") %>%
             pull(stim_time) + 0.25,
           y = 1.01,
           label = df_sim_singlemove_s5min_2995_5990_speed_sg_normbyf5minave_allmean_pairedwilcox_res %>%
             filter(group2 == "Group w/o interaction") %>%
             pull(p.adj.signif),
           angle = 90) +
  scale_color_manual(values = c("#A2A0A2", "#A93439")) +#values = viridis(3)[1:2]) +
  scale_alpha_identity() +
  scale_x_continuous(breaks = seq(0, 15, by = 3)) +
  xlab("Time after stimulus (s)") +
  ylab("Mean moving speed relative to\n the average during the initial 5 min") +
  coord_cartesian(xlim = c(-0.5, 15)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.7, 0.2),
        legend.box = "horizontal",
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank())
g_group_merge_speed_normbybeforestimave_s5min

g_merge_single_group_simulation <- 
  g_group_merge_speed_normbybeforestimave +
  g_group_merge_speed_normbybeforestimave_s5min +
  patchwork::plot_layout(ncol = 2,
                         width = c(0.65, 0.3))
g_merge_single_group_simulation
ggsave("../figures/FigureS6b.pdf", g_merge_single_group_simulation, w = 8, h = 3)


#### Figure S6c1 ####
##### load dataset #####
###### group_singlemove_mix ######
df_group_singlemove_diversity_pos2 <- read_parquet("../data/5_simulation/2_approx_parameters/parquet/df_group_singlemove_diversity_pos.parquet")

df_group_singlemove_diversity_pos3_speed_beforestimave <- 
  df_group_singlemove_diversity_pos2 %>%
  filter(seconds_total < 300) %>%
  group_by(strain, strain1, strain2, rep, id) %>%
  dplyr::summarize(speed_beforestimave = mean(speed, na.rm = TRUE)) %>%
  dplyr::select(strain, strain1, strain2, rep, id, speed_beforestimave)

df_group_singlemove_diversity_pos3_speed_normbybeforestimave <- 
  df_group_singlemove_diversity_pos2 %>%
  inner_join(df_group_singlemove_diversity_pos3_speed_beforestimave, by = c("strain", "strain1", "strain2", "rep", "id")) %>%
  dplyr::mutate(speed_normbybeforestimave = speed / speed_beforestimave)

df_group_singlemove_diversity_pos3_speed_stimminus05ave <- 
  df_group_singlemove_diversity_pos2 %>%
  filter(seconds_total >= 299.5, stim_time == -0.5) %>%
  group_by(strain, strain1, strain2, rep, id) %>%
  dplyr::summarize(speed_stimminus05ave = mean(speed, na.rm = TRUE)) %>%
  dplyr::select(strain, strain1, strain2, rep, id, speed_stimminus05ave)

df_group_singlemove_diversity_pos3_speed_normbystimminus05ave <- 
  df_group_singlemove_diversity_pos2 %>%
  filter(seconds_total >= 299.5, seconds_total < 599.5) %>%
  inner_join(df_group_singlemove_diversity_pos3_speed_stimminus05ave, by = c("strain", "strain1", "strain2", "rep", "id")) %>%
  dplyr::mutate(speed_normbystimminus05ave = speed / speed_stimminus05ave) %>%
  group_by(strain, strain1, strain2, rep, stim_time, id) %>%
  dplyr::summarize(speed_normbystimminus05ave = mean(speed_normbystimminus05ave, na.rm = TRUE))

###### calculate diversity effect ######
df_group_singlemove_diversity_mean <-
  df_group_singlemove_diversity_pos2 %>%
  group_by(strain, strain1, strain2, seconds_total) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE))

list_strain_d <- unique(df_group_singlemove_diversity_pos3_speed_beforestimave$strain)
## calculate mean of two strains ##
dfd_group_singlemove_diversity_mean <- data.frame()
for (i in list_strain_d){
  strain1 <- str_split(i, "_")[[1]][1]
  strain2 <- str_split(i, "_")[[1]][2]
  if(strain1 != strain2){
    dfd_group_singlemove_diversity_mean <- bind_rows(
      dfd_group_singlemove_diversity_mean,
      bind_rows(
        df_group_singlemove_diversity_mean[df_group_singlemove_diversity_mean$strain==paste0(strain1, "_", strain1),] %>%
          dplyr::mutate(strain = "strain1"),
        df_group_singlemove_diversity_mean[df_group_singlemove_diversity_mean$strain==paste0(strain2, "_", strain2),] %>%
          dplyr::mutate(strain = "strain2")
      ) %>%
        pivot_wider(id_cols = seconds_total, names_from = strain, values_from = speed) %>%
        dplyr::mutate(Mean = (strain1 + strain2) / 2,
               strain = i) %>%
        dplyr::rename(strain1_speed = strain1, 
                      strain2_speed = strain2) %>%
        dplyr::mutate(strain1 = strain1,
               strain2 = strain2)
    )
  }
}
dfd_group_singlemove_diversity_mean2 <-
  inner_join(df_group_singlemove_diversity_mean %>%
               ungroup() %>%
               dplyr::select(strain, seconds_total, speed) %>%
               dplyr::rename(Mixed_Observed = speed),
             dfd_group_singlemove_diversity_mean %>%
               dplyr::rename(Mixed_Expected = Mean,
                             Strain1_Observed = strain1_speed,
                             Strain2_Observed = strain2_speed), 
             by = c("strain", "seconds_total")) %>%
  dplyr::select(!c(strain1, strain2)) %>%
  pivot_longer(cols = contains("ed"), names_sep = "_", names_to = c("mixed", "type"), values_to = "speed")

dfd_group_singlemove_diversity_mean2_beforestimave <- 
  dfd_group_singlemove_diversity_mean2 %>%
  filter(seconds_total < 300) %>%
  group_by(strain, mixed, type) %>%
  dplyr::summarize(speed_beforestimave = mean(speed, na.rm = TRUE))

dfd_group_singlemove_diversity_mean2_normbybeforestimave <- 
  dfd_group_singlemove_diversity_mean2 %>%
  inner_join(dfd_group_singlemove_diversity_mean2_beforestimave, by = c("strain", "mixed", "type")) %>%
  dplyr::mutate(speed_normbybeforestimave = speed / speed_beforestimave,
         speed_normbybeforestimave = case_when(type == "Expected" ~ (lag(speed_normbybeforestimave, 2) + lag(speed_normbybeforestimave, 1)) / 2,
                                               TRUE ~ speed_normbybeforestimave),
         stim_time = (seconds_total + 0.5) %% 15 - 0.5)

dfd_group_singlemove_diversity_freezing_stopdurationdiff <-
  dfd_group_singlemove_diversity_mean2_normbybeforestimave %>% # 
  filter(seconds_total > 299) %>%
  group_by(strain, stim_time, mixed, type) %>%
  dplyr::summarize(speed_normbybeforestimave = mean(speed_normbybeforestimave, na.rm = TRUE)) %>%
  filter(stim_time != -0.5, mixed != "Mixed", speed_normbybeforestimave > 1) %>% #speed_norm_mean #speed_normbystimminus05
  group_by(strain, mixed, type) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(strain, type),
              names_from = mixed, values_from = stim_time,
              values_fill = 14) %>%
  group_by(strain) %>%
  dplyr::summarize(Strain1 = mean(Strain1, na.rm = TRUE),
                   Strain2 = mean(Strain2, na.rm = TRUE)) %>%
  dplyr::mutate(diff_stop_duration = abs(Strain1 - Strain2))

###### Behavioral performance ######
for_sec <- 3
a2 <- for_sec / 0.5
a <- (29 - a2) / a2

dfd_group_singlemove_diversity_mean2_performance <- 
  dfd_group_singlemove_diversity_mean2 %>%
  dplyr::mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5) %>%
  filter(mixed == "Mixed", stim_time >= 0, seconds_total >= 300) %>%
  group_by(strain, type, mixed, stim_time) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(performance = if_else(stim_time < for_sec, speed * (-a), speed)) %>%
  group_by(strain, type, mixed) %>%
  dplyr::summarize(performance = sum(performance)) %>%
  ungroup() %>%
  dplyr::mutate(Overyield = performance - lag(performance),
         strain1 = str_split(strain, "_") %>% map_chr(1),
         strain2 = str_split(strain, "_") %>% map_chr(2))

dfd_group_singlemove_diversity_mean2_performance_rel <- 
  dfd_group_singlemove_diversity_mean2_performance %>%
  filter(type == "Expected") %>%
  dplyr::mutate(tmp = "a") %>%
  group_by(tmp) %>%
  dplyr::summarize(mean_perf = mean(performance, na.rm = TRUE)) %>%
  left_join(dfd_group_singlemove_diversity_mean2_performance %>%
              dplyr::mutate(tmp = "a")) %>%
  dplyr::mutate(rel_performance = performance / mean_perf)

##### stat #####
df_for_stat <- dfd_group_singlemove_diversity_mean2_performance_rel %>% 
  rstatix::wilcox_test(rel_performance ~ type, 
                       paired = TRUE, detailed = TRUE)

##### make plot #####
g_mix_group_singlemove_diversity_performance <- 
  ggplot(dfd_group_singlemove_diversity_mean2_performance_rel,
         aes(x = type, y = rel_performance, color = type)) +
  geom_boxplot() + 
  geom_path(aes(group = strain), color = "gray", linewidth = 0.4) +
  geom_point(aes(shape = type), show.legend = F, size = 2) + 
  ggpubr::stat_compare_means(paired = TRUE, size = 2.5) +
  scale_color_manual(values = c("#9e8896", "#874c70"), guide = "none") +#values = viridis(3)[1:2]) +
  scale_shape_manual(values = c(17, 16)) +
  scale_alpha_manual(values = c(0.4, 1)) +
  coord_cartesian(ylim = c(0, 2.7)) +
  ylab("Relative behavioral performance") +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())
g_mix_group_singlemove_diversity_performance



#### Figure S6c2 ####
##### load dataset #####
###### group_socialcontagion_mix ######
df_group_socialcontagion_diversity_pos2 <- read_parquet("../data/5_simulation/2_approx_parameters/parquet/df_group_socialcontagion_diversity_pos.parquet")

df_group_socialcontagion_diversity_pos3_speed_beforestimave <- 
  df_group_socialcontagion_diversity_pos2 %>%
  filter(seconds_total < 300) %>%
  group_by(strain, strain1, strain2, rep, id) %>%
  dplyr::summarize(speed_beforestimave = mean(speed, na.rm = TRUE)) %>%
  dplyr::select(strain, strain1, strain2, rep, id, speed_beforestimave)

df_group_socialcontagion_diversity_pos3_speed_normbybeforestimave <- 
  df_group_socialcontagion_diversity_pos2 %>%
  inner_join(df_group_socialcontagion_diversity_pos3_speed_beforestimave, by = c("strain", "strain1", "strain2", "rep", "id")) %>%
  dplyr::mutate(speed_normbybeforestimave = speed / speed_beforestimave)

df_group_socialcontagion_diversity_pos3_speed_stimminus05ave <- 
  df_group_socialcontagion_diversity_pos2 %>%
  filter(seconds_total >= 299.5, stim_time == -0.5) %>%
  group_by(strain, strain1, strain2, rep, id) %>%
  dplyr::summarize(speed_stimminus05ave = mean(speed, na.rm = TRUE)) %>%
  dplyr::select(strain, strain1, strain2, rep, id, speed_stimminus05ave)

df_group_socialcontagion_diversity_pos3_speed_normbystimminus05ave <- 
  df_group_socialcontagion_diversity_pos2 %>%
  filter(seconds_total >= 299.5, seconds_total < 599.5) %>%
  inner_join(df_group_socialcontagion_diversity_pos3_speed_stimminus05ave, by = c("strain", "strain1", "strain2", "rep", "id")) %>%
  dplyr::mutate(speed_normbystimminus05ave = speed / speed_stimminus05ave) %>%
  group_by(strain, strain1, strain2, rep, stim_time, id) %>%
  dplyr::summarize(speed_normbystimminus05ave = mean(speed_normbystimminus05ave, na.rm = TRUE))

###### calculate diversity effect ######
df_group_socialcontagion_diversity_mean <-
  df_group_socialcontagion_diversity_pos2 %>%
  group_by(strain, strain1, strain2, seconds_total) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE))

list_strain_d <- unique(df_group_socialcontagion_diversity_pos3_speed_beforestimave$strain)
## calculate mean of two strains ##
dfd_group_socialcontagion_diversity_mean <- data.frame()
for (i in list_strain_d){
  strain1 <- str_split(i, "_")[[1]][1]
  strain2 <- str_split(i, "_")[[1]][2]
  if(strain1 != strain2){
    dfd_group_socialcontagion_diversity_mean <- bind_rows(
      dfd_group_socialcontagion_diversity_mean,
      bind_rows(
        df_group_socialcontagion_diversity_mean[df_group_socialcontagion_diversity_mean$strain==paste0(strain1, "_", strain1),] %>%
          dplyr::mutate(strain = "strain1"),
        df_group_socialcontagion_diversity_mean[df_group_socialcontagion_diversity_mean$strain==paste0(strain2, "_", strain2),] %>%
          dplyr::mutate(strain = "strain2")
      ) %>%
        pivot_wider(id_cols = seconds_total, names_from = strain, values_from = speed) %>%
        dplyr::mutate(Mean = (strain1 + strain2) / 2,
                      strain = i) %>%
        dplyr::rename(strain1_speed = strain1, 
                      strain2_speed = strain2) %>%
        dplyr::mutate(strain1 = strain1,
                      strain2 = strain2)
    )
  }
}
dfd_group_socialcontagion_diversity_mean2 <-
  inner_join(df_group_socialcontagion_diversity_mean %>%
               ungroup() %>%
               dplyr::select(strain, seconds_total, speed) %>%
               dplyr::rename(Mixed_Observed = speed),
             dfd_group_socialcontagion_diversity_mean %>%
               dplyr::rename(Mixed_Expected = Mean,
                             Strain1_Observed = strain1_speed,
                             Strain2_Observed = strain2_speed), 
             by = c("strain", "seconds_total")) %>%
  dplyr::select(!c(strain1, strain2)) %>%
  pivot_longer(cols = contains("ed"), names_sep = "_", names_to = c("mixed", "type"), values_to = "speed")

dfd_group_socialcontagion_diversity_mean2_beforestimave <- 
  dfd_group_socialcontagion_diversity_mean2 %>%
  filter(seconds_total < 300) %>%
  group_by(strain, mixed, type) %>%
  dplyr::summarize(speed_beforestimave = mean(speed, na.rm = TRUE))

dfd_group_socialcontagion_diversity_mean2_normbybeforestimave <- 
  dfd_group_socialcontagion_diversity_mean2 %>%
  inner_join(dfd_group_socialcontagion_diversity_mean2_beforestimave, by = c("strain", "mixed", "type")) %>%
  dplyr::mutate(speed_normbybeforestimave = speed / speed_beforestimave,
                speed_normbybeforestimave = case_when(type == "Expected" ~ (lag(speed_normbybeforestimave, 2) + lag(speed_normbybeforestimave, 1)) / 2,
                                                      TRUE ~ speed_normbybeforestimave),
                stim_time = (seconds_total + 0.5) %% 15 - 0.5)

dfd_group_socialcontagion_diversity_freezing_stopdurationdiff <-
  dfd_group_socialcontagion_diversity_mean2_normbybeforestimave %>% # 
  filter(seconds_total > 299) %>%
  group_by(strain, stim_time, mixed, type) %>%
  dplyr::summarize(speed_normbybeforestimave = mean(speed_normbybeforestimave, na.rm = TRUE)) %>%
  filter(stim_time != -0.5, mixed != "Mixed", speed_normbybeforestimave > 1) %>% #speed_norm_mean #speed_normbystimminus05
  group_by(strain, mixed, type) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(strain, type),
              names_from = mixed, values_from = stim_time,
              values_fill = 14) %>%
  group_by(strain) %>%
  dplyr::summarize(Strain1 = mean(Strain1, na.rm = TRUE),
                   Strain2 = mean(Strain2, na.rm = TRUE)) %>%
  dplyr::mutate(diff_stop_duration = abs(Strain1 - Strain2))

###### Behavioral performance ######
for_sec <- 3
a2 <- for_sec / 0.5
a <- (29 - a2) / a2

dfd_group_socialcontagion_diversity_mean2_performance <- 
  dfd_group_socialcontagion_diversity_mean2 %>%
  dplyr::mutate(stim_time = (seconds_total + 0.5) %% 15 - 0.5) %>%
  filter(mixed == "Mixed", stim_time >= 0, seconds_total >= 300) %>%
  group_by(strain, type, mixed, stim_time) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(performance = if_else(stim_time < for_sec, speed * (-a), speed)) %>%
  group_by(strain, type, mixed) %>%
  dplyr::summarize(performance = sum(performance)) %>%
  ungroup() %>%
  dplyr::mutate(Overyield = performance - lag(performance),
                strain1 = str_split(strain, "_") %>% map_chr(1),
                strain2 = str_split(strain, "_") %>% map_chr(2))


dfd_group_socialcontagion_diversity_mean2_performance_rel <- 
  dfd_group_socialcontagion_diversity_mean2_performance %>%
  filter(type == "Expected") %>%
  dplyr::mutate(tmp = "a") %>%
  group_by(tmp) %>%
  dplyr::summarize(mean_perf = mean(performance, na.rm = TRUE)) %>%
  left_join(dfd_group_socialcontagion_diversity_mean2_performance %>%
              dplyr::mutate(tmp = "a")) %>%
  dplyr::mutate(rel_performance = performance / mean_perf)

##### stat #####
df_for_stat <- dfd_group_socialcontagion_diversity_mean2_performance_rel %>% 
  rstatix::wilcox_test(rel_performance ~ type, 
                       paired = TRUE, detailed = TRUE)

##### make plot #####
g_mix_group_socialcontagion_diversity_performance <- 
  ggplot(dfd_group_socialcontagion_diversity_mean2_performance_rel,
         aes(x = type, y = rel_performance, color = type)) +
  geom_boxplot() + 
  geom_path(aes(group = strain), color = "gray", linewidth = 0.4) +
  geom_point(aes(shape = type), show.legend = F, size = 2) + 
  ggpubr::stat_compare_means(paired = TRUE, size = 2.5) +
  scale_color_manual(values = c("#9e8896", "#874c70"), guide = "none") +#values = viridis(3)[1:2]) +
  scale_shape_manual(values = c(17, 16)) +
  scale_alpha_manual(values = c(0.4, 1)) +
  coord_cartesian(ylim = c(0, 2.7)) +
  ylab("Relative behavioral performance") +
  theme_bw() +
  theme(legend.title = element_blank(),
        axis.title.x = element_blank())
g_mix_group_socialcontagion_diversity_performance


g_merge_simulation <- 
  g_group_merge_speed_normbybeforestimave +
  g_group_merge_speed_normbybeforestimave_s5min +
  g_mix_group_singlemove_diversity_performance +
  g_mix_group_socialcontagion_diversity_performance +
  patchwork::plot_layout(ncol = 4,
                         width = c(0.5, 0.3, 0.15, 0.15))
g_merge_simulation
ggsave("../figures/FigureS6bc.pdf", g_merge_simulation, w=9, h=2.5)


##### merge #####
g_merge_simulation <- 
  g_group_merge_speed_normbybeforestimave +
  g_group_merge_speed_normbybeforestimave_s5min +
  g_mix_group_singlemove_diversity_performance +
  patchwork::plot_layout(ncol = 3,
                         width = c(0.65, 0.3, 0.25))
g_merge_simulation
# ggsave("../figures/FigureS6bc.pdf", g_merge_simulation, w=8.5, h=2.5)

#### load packages ####
library(tidyverse)
library(arrow)
# library(car)
# library(lmerTest)
# library(ggsci)
# library(ggpmisc)
library(patchwork)

#### Figure 1c ####
##### load dataset #####
df_f5min_speed <- read_parquet("../data/1_single_strain/df_f5min_speed.parquet") %>%
  ungroup()

df_f10min_speed <- read_parquet("../data/1_single_strain/df_f10min_speed.parquet") %>%
  ungroup()

df_f5min_speed_ave <- read_parquet("../data/1_single_strain/df_f5min_speed_ave.parquet") %>%
  ungroup()

df_s5min_2995_5990_speed <- read_parquet("../data/1_single_strain/df_s5min_2995_5990_speed.parquet") %>%
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

##### make plot #####
g_f10min_sg_speed_normbyf5minave_allmean <- 
  ggplot(df_f10min_speed %>%
           filter(strain != "norpA") %>%
           left_join(df_f5min_speed_ave %>%
                       filter(strain != "norpA") %>%
                       group_by(n_inds, sex) %>%
                       dplyr::summarize(speed_f5min_ave = mean(speed_f5min_ave, na.rm = TRUE)), 
                     by = c("n_inds", "sex")) %>%
           group_by(n_inds, sex, seconds_total) %>%
           dplyr::summarize(speed = mean(speed, na.rm = TRUE),
                            speed_f5min_ave = mean(speed_f5min_ave, na.rm = TRUE)) %>%
           mutate(speed_normbyf5minave = speed / speed_f5min_ave),
         aes(x = seconds_total, 
             y = speed_normbyf5minave, 
             col = n_inds,
             group = n_inds)) +
  geom_line() +
  geom_vline(xintercept=300, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=315, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=330, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=345, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=360, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=375, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=390, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=405, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=420, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=435, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=450, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=465, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=480, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=495, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=510, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=525, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=540, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=555, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=570, linetype="dotted", alpha=.6, linewidth=.4) +
  geom_vline(xintercept=585, linetype="dotted", alpha=.6, linewidth=.4) +  scale_color_viridis_d(direction = -1) +
  xlab("Time (min)") +
  ylab("Moving speed relative to \n the average during the initial 5 min") +
  scale_x_continuous(breaks = seq(0, 600, 60), labels = seq(0, 10, 1)) +
  # scale_color_manual(values = rev(ggsci::pal_uchicago("default")(2))) +#values = viridis(3)[1:2]) +
  scale_color_manual(values = c("#aaaaa9", "#b3343a")) +
  facet_wrap(~ sex, nrow = 2) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.15, 0.65),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        strip.background = element_blank())
g_f10min_sg_speed_normbyf5minave_allmean

g_s5min_2995_5990_sg_speed_normbyf5minave_allmean <- 
  ggplot(df_s5min_2995_5990_speed_sg_normbyf5minave_strain %>%
           filter(strain != "norpA") %>%
           group_by(n_inds, sex, stim_time) %>%
           dplyr::summarize(speed_normbyf5minave = mean(speed_normbyf5minave, na.rm = TRUE)),
         aes(x = stim_time, 
             y = speed_normbyf5minave, 
             col = n_inds,
             group = n_inds)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  annotate("rect", xmin = 0, xmax = 0.5, ymin = -Inf, ymax = Inf, alpha = 0.4) +
  scale_x_continuous(breaks = seq(0, 15, 3)) +
  scale_color_manual(values = c("#aaaaa9", "#b3343a")) +
  coord_cartesian(xlim = c(-0.5, 15)) +
  xlab("Time after stimulus (s)") +
  ylab("Moving speed relative to \n the average during the initial 5 min") +
  facet_wrap(~ sex, nrow = 2) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.7, 0.65),
        legend.key.size = unit(0.4, "cm"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())
g_s5min_2995_5990_sg_speed_normbyf5minave_allmean

g_fig_sg_speed <- g_f10min_sg_speed_normbyf5minave_allmean +
  g_s5min_2995_5990_sg_speed_normbyf5minave_allmean +
  patchwork::plot_layout(ncol = 2,
                         width = c(0.8, 0.3))
ggsave("../figures/Figure1c.pdf", g_fig_sg_speed, w=6, h=3.5)


#### Figure 1d ####
##### load dataset #####
df_s5min_stim_vis <- read_parquet("../data/1_single_strain/df_s5min_stim_vis_0.5.parquet") %>%
  ungroup() 

df_s5min_stim_freez_vis <- df_s5min_stim_vis %>%
  dplyr::filter(n_inds == "Group") %>%
  transform(change_posture = factor(change_posture, levels=c("ww", "ws", "sw", "ss"))) %>%
  mutate(motion_cue_diff = lead(motion_cue) - motion_cue,
         motion_cue_diff2 = motion_cue - lead(motion_cue),
         motion_cue_diff3 = motion_cue - lag(motion_cue),
         motion_cue_next = lead(motion_cue))


##### make plot #####
gg_s5min_stim_motion_cue_exit <- ggplot(df_s5min_stim_freez_vis %>%
                                          filter(stimuli == "+10.0", change_posture %in% c("ss", "sw")) %>%
                                          mutate(change_posture = if_else(change_posture == "sw", 1, 0),
                                                 posture_num = if_else(posture == "walk", 1, 0)),
                                          aes(x = motion_cue_diff3, y = change_posture, col = sex)) +
  geom_smooth(aes(x = motion_cue_diff3, y = change_posture),
              method = "glm", method.args = list(family = "binomial")) +
  coord_cartesian(xlim = c(-40, 80), ylim = c(0, 1)) +
  # scale_color_manual(values = c("#c97586","#2a83a2")) +
  scale_color_manual(values = c("#c17181", "#68aac3")) +
  xlab("Increase in motion cue") +
  ylab("Probability of freezing exit") +
  theme_bw() +
  theme(legend.position = c(0.2, 0.85),
        legend.title = element_blank(),
        legend.background = element_blank())
gg_s5min_stim_motion_cue_exit
ggsave("../figures/Figure1d.pdf", gg_s5min_stim_motion_cue_exit, w=3, h=3)

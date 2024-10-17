#### load packages ####
library(tidyverse)
library(arrow)

#### Figure S11 ####
##### load dataset #####
###### diversity effect ######
dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial <- 
  read_parquet("../data/2_mixed_strain/dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial.parquet") %>%
  ungroup()

for_sec <- 3
a2 <- for_sec / 0.5
a <- (29 - a2) / a2

dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging <- 
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial %>%
  filter(n_inds == "Group", type == "Mixed", stim_time >= 0) %>%
  mutate(foraging_success = if_else(stim_time < for_sec, speed * (-a), speed)) %>%
  group_by(strain, sex, n_inds, trial, var, type, mixed, alpha) %>%
  dplyr::summarize(foraging_success = sum(foraging_success))

dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging2 <- 
  dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain_trial %>%
  filter(n_inds == "Group", stim_time >= 0) %>%
  mutate(foraging_success = if_else(stim_time < for_sec, speed * (-a), speed)) %>%
  group_by(strain, sex, n_inds, trial, var, type, mixed, alpha) %>%
  dplyr::summarize(foraging_success = sum(foraging_success))

dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging_tr <- 
  dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging2 %>%
  filter(n_inds == "Group") %>%
  group_by(strain, sex, n_inds, var, type, mixed, alpha) %>%
  dplyr::summarize(foraging_success = mean(foraging_success, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(strain1 = dplyr::lag(foraging_success, 3),
         strain2 = dplyr::lag(foraging_success, 2),
         mean_2strains =  dplyr::lag(foraging_success, 1)) %>%
  filter(var == "Group_Mixed") %>%
  mutate(overyield = foraging_success - mean_2strains,
         tr_overyield = if_else(strain1 > strain2, foraging_success - strain1, foraging_success - strain2),
         tr_underyield = if_else(strain1 > strain2, foraging_success - strain2, foraging_success - strain1),
         col = case_when(tr_overyield > 0 ~ "Transgressive overyielding",
                         tr_overyield <= 0 & overyield > 0 ~ "Overyielding",
                         tr_underyield >= 0 & overyield < 0 ~ "Underyielding",
                         tr_underyield < 0 ~ "Transgressive underyielding")) %>%
  transform(col = factor(col, levels = c("Transgressive overyielding", "Overyielding", "Underyielding", "Transgressive underyielding")))


###### motion cue exit diff ######
df_motion_cue_exit_coeff <- read_parquet("../data/1_single_strain/df_motion_cue_exit_coeff.parquet") %>%
  ungroup()

dfd_motion_cue_exit_coeff_diff <- data.frame()
for (i in unique(dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging_tr$strain)){
  strain1 <- str_split(i, "_")[[1]][1]
  strain2 <- str_split(i, "_")[[1]][2]
  dfd_motion_cue_exit_coeff_diff <- bind_rows(
    dfd_motion_cue_exit_coeff_diff,
    bind_rows(
      df_motion_cue_exit_coeff %>%
        filter(sex == "Female", strain == eval(strain1)),
      df_motion_cue_exit_coeff %>%
        filter(sex == "Female", strain == eval(strain2))
    ) %>%
      dplyr::select(!motion_cue_exit_coeff) %>%
      pivot_wider(id_cols = sex, names_from = strain, values_from = motion_cue_exit_intercept) %>%
      mutate(Mean = (get(strain1) + get(strain2)) / 2,
             Diff = abs(get(strain1) - get(strain2)),
             strain = i) %>%
      dplyr::rename(strain1_motion_cue_exit_intercept = strain1, 
                    strain2_motion_cue_exit_intercept = strain2) %>%
      mutate(strain1 = strain1,
             strain2 = strain2)
  )
}


##### make plot #####
g_gd_foraging_success_reaction_diff <-
  ggplot(inner_join(dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging_tr,
                    dfd_motion_cue_exit_coeff_diff, 
                    by = "strain"),
         aes(x = Diff, y = overyield)) +
  geom_point() +
  stat_smooth(linewidth = 2, color= "grey", method = "lm") +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(
                          after_stat(rr.label),
                          stat(p.value.label),
                          sep = "~~~")),
                        label.x = "left",
                        label.y = "top",
                        parse = TRUE, size = 4) + 
  xlab("Difference in visual reactivity") +
  ylab("Diversity effect on fitness") +
  theme_bw()
g_gd_foraging_success_reaction_diff
ggsave("../figures/FigureS11.pdf", g_gd_foraging_success_reaction_diff, w=3, h=3)

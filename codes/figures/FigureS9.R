#### load packages ####
library(tidyverse)
library(gtools)
library(arrow)

#### load functions ####
t_test_1sample_foraging_success <- function(dat){
  pval <-
    t.test(dat %>%
             filter(alpha == "Observed") %>%
             pull(foraging_success),
           mu = dat %>%
             filter(alpha == "Expected") %>%
             pull(foraging_success) %>%
             mean()
    )$p.value
  return(pval)
}

#### Figure S9 ####
##### load dataset #####
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

dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging_t_test <- 
  dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging %>% 
  group_by(strain) %>%
  nest() %>%
  dplyr::mutate(t.test.p.value = map(data, t_test_1sample_foraging_success)) %>%
  unnest(data) %>%
  mutate(t.test.p.value = as.numeric(t.test.p.value)) %>%
  dplyr::select(strain, t.test.p.value) %>%
  distinct() %>%
  mutate(col2 = if_else(t.test.p.value < 0.05, "sig", "insig"))

dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging3 <-
  dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging2 %>%
  group_by(strain, sex, n_inds, var, type, mixed, alpha) %>%
  dplyr::summarize(mean = mean(foraging_success, na.rm = T),
                   error = sd(foraging_success, na.rm = T)) %>%
  left_join(dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging_tr %>%
              dplyr::select(!c("var", "type", "mixed", "alpha"))) %>%
  left_join(dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging_t_test %>% 
              dplyr::select(strain, t.test.p.value, col2)) %>%
  transform(strain = factor(strain, levels = unique(.$strain) %>% gtools::mixedsort()))

##### make plot #####
g_d_s5min_stim_speed_raw_foraging_integrated2 <-  
  dfd_s5min_2995_5990_speed_sgd_raw_strain_trial_foraging3 %>%
  filter(var == "Group_Mixed") %>%
  ggplot(aes(x = reorder(strain, overyield, na.rm = TRUE), 
             y = overyield, col = col, alpha = col2)) +#col = t.test.p.value)) +
  geom_point(size = 4, shape = 16) +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_color_manual(values = c("#642426", "#CF7B2E", "#5A7C69", "#215A83")) +
  xlab("Strain") +
  ylab("Synergistic effects on fitness") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.35, 0.8),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.box = "horizontal")
g_d_s5min_stim_speed_raw_foraging_integrated2
ggsave("../figures/FigureS9.pdf", g_d_s5min_stim_speed_raw_foraging_integrated2, w=5, h=3)

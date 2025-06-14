#### load packages ####
targetPackages <- c('tidyverse','arrow','car','lme4','lmerTest','ggpmisc','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### load dataset ####
df_s5min_stim_vis <- read_parquet("../data/1_single_strain/parquet/df_s5min_stim_vis_0.5.parquet") %>%
  ungroup() 

df_s5min_stim_freez_vis <- df_s5min_stim_vis %>%
  dplyr::filter(n_inds == "Group") %>%
  transform(change_posture = factor(change_posture, levels=c("ww", "ws", "sw", "ss"))) %>%
  mutate(motion_cue_diff = lead(motion_cue) - motion_cue,
         motion_cue_diff2 = motion_cue - lead(motion_cue),
         motion_cue_diff3 = motion_cue - lag(motion_cue),
         motion_cue_next = lead(motion_cue))

df_motion_cue_exit_coeff <- read_parquet("../data/1_single_strain/parquet/df_motion_cue_exit_coeff.parquet") %>%
  ungroup()

df_motion_cue_exit_coeff_trial <- read_parquet("../data/1_single_strain/parquet/df_motion_cue_exit_coeff_trial.parquet") %>%
  ungroup()

#### analysis ####
var <- "motion_cue_exit_intercept"
r <- 20

##### LMM #####
# glmer_res <- lme4::glmer(change_posture ~ motion_cue_diff3*strain*sex +  
#                      (1|date) + (1|time) + (1|prefix) + (1|place), 
#                    df_s5min_stim_freez_vis %>% 
#                      filter(!str_detect(strain ,"norpA"),
#                             stimuli == "+10.0", change_posture %in% c("ss", "sw")) %>%
#                      mutate(change_posture = as.factor(change_posture)),
#                    family="binomial")
# anova_res <- car::Anova(glmer_res)
# anova_res

model_motion_cue_exit_intercept <- 
  lmerTest::lmer(get(var) ~ strain*sex +  
                   (1|date) + (1|time) + (1|prefix) + (1|place), 
                 df_motion_cue_exit_coeff_trial %>% 
                   filter(!str_detect(strain ,"norpA"))) #must be factor to use difflsmeans function
summary(model_motion_cue_exit_intercept)
car::Anova(model_motion_cue_exit_intercept)

##### male #####
###### LMM ######
model_motion_cue_exit_intercept_male <- 
  lmerTest::lmer(get(var) ~ strain +  
                (1|date) + (1|time) + (1|prefix) + (1|place), 
              df_motion_cue_exit_coeff_trial %>% 
                filter(!str_detect(strain ,"norpA"),
                       sex == "Male") %>%
                dplyr::mutate(strain = as.factor(strain))) #must be factor to use difflsmeans function
summary(model_motion_cue_exit_intercept_male)
car::Anova(model_motion_cue_exit_intercept_male)

###### calculate statistical difference ######
model_motion_cue_exit_intercept_male_all <- 
  lmerTest::lmer(get(var) ~ strain +  
                   (1|date) + (1|time) + (1|prefix) + (1|place), 
                 df_motion_cue_exit_coeff_trial %>% 
                   filter(sex == "Male") %>%
                   dplyr::mutate(strain = as.factor(strain))) #must be factor to use difflsmeans function

df_difflsmeans_motion_cue_exit_intercept_male <- lmerTest::difflsmeans(model_motion_cue_exit_intercept_male_all) %>% #multiple correction based on the Satterthwaite's approximation to degrees of freedom
  tibble::rownames_to_column() %>%
  as_tibble()

df_difflsmeans_motion_cue_exit_intercept_male2 <- df_difflsmeans_motion_cue_exit_intercept_male %>%
  dplyr::mutate(rowname = str_remove_all(rowname, "strain")) %>%
  filter(str_detect(rowname, "norpA")) %>%
  separate(rowname, sep = " - ", into = c("strain1", "strain2")) 

df_difflsmeans_motion_cue_exit_intercept_male2_norpA <- df_difflsmeans_motion_cue_exit_intercept_male2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "norpA", strain2, strain1),
                strain2 = if_else(strain1 == "norpA", "norpA", strain2),
                Estimate = if_else(strain1 == "norpA", Estimate, -Estimate),
                strain1 = if_else(strain1 == "norpA", strain_tmp, strain1)) %>%
  filter(strain2 == "norpA") %>%
  bind_rows(data.frame(strain1 = "norpA",
                       strain2 = "norpA"))

df_difflsmeans_motion_cue_exit_intercept_male2_208norpA <- df_difflsmeans_motion_cue_exit_intercept_male2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "DGRP208_norpA", strain2, strain1),
                strain2 = if_else(strain1 == "DGRP208_norpA", "DGRP208_norpA", strain2),
                Estimate = if_else(strain1 == "DGRP208_norpA", Estimate, -Estimate),
                strain1 = if_else(strain1 == "DGRP208_norpA", strain_tmp, strain1)) %>%
  filter(strain2 == "DGRP208_norpA") %>%
  bind_rows(data.frame(strain1 = "DGRP208_norpA",
                       strain2 = "DGRP208_norpA"))

df_s5min_motion_cue_exit_intercept_male_order <- df_motion_cue_exit_coeff_trial %>% 
  filter(sex == "Male",
         # !strain %in% c("norpA", "DGRP208_norpA")) %>%
         strain != "norpA") %>%
  group_by(strain) %>%
  dplyr::summarize(motion_cue_exit_intercept_mean = mean(motion_cue_exit_intercept, na.rm = T)) %>%
  ungroup() %>%
  arrange(motion_cue_exit_intercept_mean)

###### calculate heritability ######
res_male <- anova(lm(get(var) ~ strain, data = df_motion_cue_exit_coeff_trial %>% 
                       filter(!str_detect(strain ,"norpA"),
                              sex == "Male")))
sprintf("%.2e", res_male$`Pr(>F)`[1])
sg2_male <- (res_male$`Mean Sq`[1] - res_male$`Mean Sq`[2]) / r
se2_male <- res_male$`Mean Sq`[2]
sg2_male / (sg2_male + se2_male/r)


##### female #####
###### LMM ######
model_motion_cue_exit_intercept_female <- 
  lmerTest::lmer(get(var) ~ strain +  
                   (1|date) + (1|time) + (1|prefix) + (1|place), 
                 df_motion_cue_exit_coeff_trial %>% 
                   filter(!str_detect(strain ,"norpA"),
                          sex == "Female") %>%
                   dplyr::mutate(strain = as.factor(strain))) #must be factor to use difflsmeans function
summary(model_motion_cue_exit_intercept_female)
car::Anova(model_motion_cue_exit_intercept_female)

###### calculate statistical difference ######
model_motion_cue_exit_intercept_female_all <- 
  lmerTest::lmer(get(var) ~ strain +  
                   (1|date) + (1|time) + (1|prefix) + (1|place), 
                 df_motion_cue_exit_coeff_trial %>% 
                   filter(sex == "Female") %>%
                   dplyr::mutate(strain = as.factor(strain))) #must be factor to use difflsmeans function

df_difflsmeans_motion_cue_exit_intercept_female <- lmerTest::difflsmeans(model_motion_cue_exit_intercept_female_all) %>% #multiple correction based on the Satterthwaite's approximation to degrees of freedom
  tibble::rownames_to_column() %>%
  as_tibble()

df_difflsmeans_motion_cue_exit_intercept_female2 <- df_difflsmeans_motion_cue_exit_intercept_female %>%
  dplyr::mutate(rowname = str_remove_all(rowname, "strain")) %>%
  filter(str_detect(rowname, "norpA")) %>%
  separate(rowname, sep = " - ", into = c("strain1", "strain2")) 

df_difflsmeans_motion_cue_exit_intercept_female2_norpA <- df_difflsmeans_motion_cue_exit_intercept_female2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "norpA", strain2, strain1),
                strain2 = if_else(strain1 == "norpA", "norpA", strain2),
                Estimate = if_else(strain1 == "norpA", Estimate, -Estimate),
                strain1 = if_else(strain1 == "norpA", strain_tmp, strain1)) %>%
  filter(strain2 == "norpA") %>%
  bind_rows(data.frame(strain1 = "norpA",
                       strain2 = "norpA"))

df_difflsmeans_motion_cue_exit_intercept_female2_208norpA <- df_difflsmeans_motion_cue_exit_intercept_female2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "DGRP208_norpA", strain2, strain1),
                strain2 = if_else(strain1 == "DGRP208_norpA", "DGRP208_norpA", strain2),
                Estimate = if_else(strain1 == "DGRP208_norpA", Estimate, -Estimate),
                strain1 = if_else(strain1 == "DGRP208_norpA", strain_tmp, strain1)) %>%
  filter(strain2 == "DGRP208_norpA") %>%
  bind_rows(data.frame(strain1 = "DGRP208_norpA",
                       strain2 = "DGRP208_norpA"))

df_s5min_motion_cue_exit_intercept_female_order <- df_motion_cue_exit_coeff_trial %>% 
  filter(sex == "Female",
         # !strain %in% c("norpA", "DGRP208_norpA")) %>%
         strain != "norpA") %>%
  group_by(strain) %>%
  dplyr::summarize(motion_cue_exit_intercept_mean = mean(motion_cue_exit_intercept, na.rm = T)) %>%
  ungroup() %>%
  arrange(motion_cue_exit_intercept_mean)

###### calculate heritability ######
res_female <- anova(lm(get(var) ~ strain, data = df_motion_cue_exit_coeff_trial %>% 
                       filter(!strain %in% c("norpA", "DGRP208_norpA"),
                              sex == "Female")))
sprintf("%.2e", res_female$`Pr(>F)`[1])
sg2_female <- (res_female$`Mean Sq`[1] - res_female$`Mean Sq`[2]) / r
se2_female <- res_female$`Mean Sq`[2]
sg2_female / (sg2_female + se2_female/r)


#### make plot ####
##### Figure S5a #####
gg_s5min_stim_motion_cue_exit_strain <- ggplot(df_s5min_stim_freez_vis %>%
                                                 filter(!strain %in% c("norpA"),
                                                        stimuli == "+10.0", change_posture %in% c("ss", "sw")) %>%
                                                 mutate(change_posture = if_else(change_posture == "sw", 1, 0),
                                                        posture_num = if_else(posture == "walk", 1, 0)),
                                               aes(x = motion_cue_diff3, y = change_posture, col = sex)) +
  geom_point(alpha=0.1, size = 1) +
  geom_smooth(aes(x = motion_cue_diff3, y = change_posture),
              method = "glm", method.args = list(family = "binomial")) +
  coord_cartesian(xlim=c(-40,60), ylim=c(0,1)) +
  scale_x_continuous(breaks = c(0, 50)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_color_manual(values = c("#c17181", "#68aac3")) +
  xlab("Increase in motion cue") +
  ylab("Probability of freezing exit") +
  facet_wrap( ~ strain, nrow = 10) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        strip.background = element_blank())
gg_s5min_stim_motion_cue_exit_strain
ggsave("../figures/FigureS5a.pdf", gg_s5min_stim_motion_cue_exit_strain, w=10, h=10)

##### Figure S5b1 male #####
gg_s5min_motion_cue_exit_intercept_male <- 
  df_motion_cue_exit_coeff_trial %>% 
  mutate(plot_col = case_when(
    strain == "DGRP208" ~ "#833163",
    strain == "DGRP208_norpA" ~ "#e8c72d", #"#A8DADC", 
    TRUE ~ "black")) %>%
  filter(sex == "Male", strain != "norpA") %>%
  ggplot(aes(x = reorder(strain, motion_cue_exit_intercept, na.rm = TRUE),
             y = motion_cue_exit_intercept, col = plot_col)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0) +
  stat_summary(fun.data = "mean_se", geom = "point") +
  coord_cartesian(ylim = c(-80, 30)) +
  # scale_color_manual(values = c("black", "magenta")) + 
  scale_color_identity() +
  xlab("Strain") +
  ylab(expression("Reaction threshold" ~ (beta[0]))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        legend.position = 'none')
gg_s5min_motion_cue_exit_intercept_male

gg_s5min_motion_cue_exit_intercept_male_sig <- 
  ggplot(df_difflsmeans_motion_cue_exit_intercept_male2_208norpA %>%
           filter(strain1 != "norpA") %>%
           add_column(y = "value") %>%
           dplyr::rename(P = `Pr(>|t|)`) %>%
           dplyr::mutate(strain2 = str_replace(strain2, "208norpA", "vs. DGRP208_norpA")) %>%
           transform(strain1 = factor(strain1, 
                                      levels = df_s5min_motion_cue_exit_intercept_male_order$strain)), 
         aes(x = strain1, y = y, col = -Estimate, size = -log10(P))) +
  geom_point() +
  scale_colour_gradient2(low = "#1e50a2", mid = "#e5e4e6", high = "#c53d43",
                         limits = c(-50, 55), midpoint = 0) +
  scale_size_continuous(limits = c(0, 30)) +
  facet_wrap(~ strain2, nrow = 2) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.background = element_blank())#,

gg_s5min_motion_cue_exit_intercept_male_sig
# gg_s5min_motion_cue_exit_intercept_male_norpA_sig
# gg_s5min_motion_cue_exit_intercept_male_random_sig

gg_s5min_motion_cue_exit_intercept_male_sum <- 
  gg_s5min_motion_cue_exit_intercept_male_sig +
  gg_s5min_motion_cue_exit_intercept_male + 
  plot_layout(nrow = 2, ncol = 1, heights = c(0.3,1))
gg_s5min_motion_cue_exit_intercept_male_sum
ggsave("../figures/FigureS5b1_male.pdf", gg_s5min_motion_cue_exit_intercept_male_sum, w = 6, h = 2)


##### Figure S5b2 female #####
gg_s5min_motion_cue_exit_intercept_female <- 
  df_motion_cue_exit_coeff_trial %>% 
  mutate(plot_col = case_when(
    strain == "DGRP208" ~ "#833163",
    strain == "DGRP208_norpA" ~ "#e8c72d", #"#A8DADC", 
    TRUE ~ "black")) %>%
  filter(sex == "Female", strain != "norpA") %>%
  ggplot(aes(x = reorder(strain, motion_cue_exit_intercept, na.rm = TRUE),
             y = motion_cue_exit_intercept, col = plot_col)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0) +
  stat_summary(fun.data = "mean_se", geom = "point") +
  coord_cartesian(ylim = c(-80, 30)) +
  # scale_color_manual(values = c("black", "magenta")) + 
  scale_color_identity() +
  xlab("Strain") +
  ylab(expression("Reaction threshold" ~ (beta[0]))) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        legend.position = 'none')
gg_s5min_motion_cue_exit_intercept_female

gg_s5min_motion_cue_exit_intercept_female_sig <- 
  ggplot(df_difflsmeans_motion_cue_exit_intercept_female2_208norpA %>%
           filter(strain1 != "norpA") %>%
           add_column(y = "value") %>%
           dplyr::rename(P = `Pr(>|t|)`) %>%
           dplyr::mutate(strain2 = str_replace(strain2, "208norpA", "vs. DGRP208_norpA")) %>%
           transform(strain1 = factor(strain1, 
                                      levels = df_s5min_motion_cue_exit_intercept_female_order$strain)), 
         aes(x = strain1, y = y, col = -Estimate, size = -log10(P))) +
  geom_point() +
  scale_colour_gradient2(low = "#1e50a2", mid = "#e5e4e6", high = "#c53d43",
                         limits = c(-50, 55), midpoint = 0) +
  scale_size_continuous(limits = c(0, 30)) +
  facet_wrap(~ strain2, nrow = 2) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.background = element_blank())
gg_s5min_motion_cue_exit_intercept_female_sig


gg_s5min_motion_cue_exit_intercept_female_sum <- 
  gg_s5min_motion_cue_exit_intercept_female_sig +
  gg_s5min_motion_cue_exit_intercept_female + 
  plot_layout(nrow = 2, ncol = 1, heights = c(0.3,1))
gg_s5min_motion_cue_exit_intercept_female_sum
ggsave("../figures/FigureS5b2_female.pdf", gg_s5min_motion_cue_exit_intercept_female_sum, w = 6, h = 2)


##### Figure S5c #####
df_motion_cue_exit_coeff %>%
  filter(!strain %in% c("norpA", "DGRP208_norpA")) %>%
  group_by(strain, sex) %>%
  dplyr::summarize(motion_cue_exit_intercept = mean(motion_cue_exit_intercept, na.rm=T)) %>%
  pivot_wider(names_from = sex, values_from = motion_cue_exit_intercept, names_prefix = "motion_cue_exit_intercept_", names_sep = "_") %>%
  mutate(plot_col = if_else(motion_cue_exit_intercept_Female - motion_cue_exit_intercept_Male > 0, "over", "under")) %>%
  group_by(plot_col) %>%
  summarize(n = n())

# stat
df_figS5c <- df_motion_cue_exit_coeff %>%
  pivot_wider(id_cols = strain, names_from = sex, values_from = motion_cue_exit_intercept, names_prefix = "motion_cue_exit_intercept_", names_sep = "_") %>%
  mutate(plot_col = if_else(motion_cue_exit_intercept_Female - motion_cue_exit_intercept_Male > 0, "over", "under")) %>%
  filter(!strain %in% c("norpA", "DGRP208_norpA"))
df_figS5c_stats <- df_figS5c %>%
  ungroup() %>%
  dplyr::summarize(
    r = cor(motion_cue_exit_intercept_Male, motion_cue_exit_intercept_Female, use = "complete.obs"),
    cor_test = list(cor.test(motion_cue_exit_intercept_Male, motion_cue_exit_intercept_Female, use = "complete.obs"))
  ) %>%
  dplyr::mutate(
    p_value = map_dbl(cor_test, ~ .x$p.value),
    r_squared = r**2,
    p_formatted = sprintf("%.2e", p_value)
  ) %>%
  dplyr::select(r_squared, p_formatted) %>%
  print()


gg_s5min_stim_motion_cue_exit_vs_sex <- df_figS5c %>%
  ggplot(aes(x = motion_cue_exit_intercept_Male, y = motion_cue_exit_intercept_Female)) +
  geom_abline(slope = 1, linetype = "dashed") +
  stat_smooth(linewidth = 2, color= "grey", method = "lm") +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(
                          after_stat(rr.label),
                          stat(p.value.label),
                          sep = "~~~")),
                        label.x = "left",
                        label.y = "top",
                        parse = TRUE, size = 4) +
  geom_point(aes(col = plot_col), alpha = 0.5, size = 4, shape = 16) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_color_manual(values = c("#c17181", "#68aac3")) +
  labs(x = expression("Visual responsiveness to motion cue" ~ (beta[0]) ~ "in male"),
       y = expression("Visual responsiveness to motion cue" ~ (beta[0]) ~ "in female")) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        legend.position = 'none')

gg_s5min_stim_motion_cue_exit_vs_sex
ggsave("../figures/FigureS5c.pdf", gg_s5min_stim_motion_cue_exit_vs_sex, w=4, h=4)


#### load packages ####
targetPackages <- c('tidyverse','arrow','car','lmerTest','ggpmisc','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### load dataset ####
df_f5min_speed_wo_mutant <- read_parquet("../data/1_single_strain/parquet/df_f5min_speed.parquet") %>%
  filter(!str_detect(strain ,"norpA")) %>%
  group_by(strain, n_inds, sex, trial, prefix, date, time, place) %>%
  dplyr::summarize(speed = mean(speed, na.rm = TRUE))

#### analysis ####
var <- "speed"
r <- 20

##### LMM #####
model_speed <- lmerTest::lmer(get(var) ~ strain * sex * n_inds + 
                              (1|date) + (1|time) + (1|prefix) + (1|place), 
                              df_f5min_speed_wo_mutant)
summary(model_speed)
car::Anova(model_speed)

##### heritability single female #####
res_single_female <- anova(lm(get(var) ~ strain, data = df_f5min_speed_wo_mutant %>% 
                               filter(sex == "Female", n_inds == "Single"
                               )))
sg2_single_female <- (res_single_female$`Mean Sq`[1] - res_single_female$`Mean Sq`[2]) / r
se2_single_female <- res_single_female$`Mean Sq`[2]
sg2_single_female / (sg2_single_female + se2_single_female/r)

###### calculate heritability single female with Hoffmann and Parsons (1988) method ######

# n <- 104
# (n*sg2_single_female - se2_single_female) / (n*sg2_single_female + (n-1)*se2_single_female) #t_single_female


##### heritability group female #####
res_group_female <- anova(lm(get(var) ~ strain, data = df_f5min_speed_wo_mutant %>% 
                               filter(sex == "Female", n_inds == "Group"
                               )))
sg2_group_female <- (res_group_female$`Mean Sq`[1] - res_group_female$`Mean Sq`[2]) / r
se2_group_female <- res_group_female$`Mean Sq`[2]
sg2_group_female / (sg2_group_female + se2_group_female/r)

###### calculate heritability group female with Hoffmann and Parsons (1988) method ######

# n <- 104
# (n*sg2_group_female - se2_group_female) / (n*sg2_group_female + (n-1)*se2_group_female) #t_group_female


##### heritability single male #####
res_single_male <- anova(lm(get(var) ~ strain, data = df_f5min_speed_wo_mutant %>% 
                              filter(sex == "Male", n_inds == "Single"
                              )))
sg2_single_male <- (res_single_male$`Mean Sq`[1] - res_single_male$`Mean Sq`[2]) / r
se2_single_male <- res_single_male$`Mean Sq`[2]
sg2_single_male / (sg2_single_male + se2_single_male/r)

###### calculate heritability single male with Hoffmann and Parsons (1988) method ######

# n <- 104
# (n*sg2_single_male - se2_single_male) / (n*sg2_single_male + (n-1)*se2_single_male) #t_single_male


##### heritability group male #####
res_group_male <- anova(lm(get(var) ~ strain, data = df_f5min_speed_wo_mutant %>% 
                             filter(sex == "Male", n_inds == "Group"
                             )))
sg2_group_male <- (res_group_male$`Mean Sq`[1] - res_group_male$`Mean Sq`[2]) / r
se2_group_male <- res_group_male$`Mean Sq`[2]
sg2_group_male / (sg2_group_male + se2_group_male/r)

###### calculate heritability group male with Hoffmann and Parsons (1988) method ######

# n <- 104
# (n*sg2_group_male - se2_group_male) / (n*sg2_group_male + (n-1)*se2_group_male) #t_group_male


#### make plot ####
gig_f5min_speed_sex <- df_f5min_speed_wo_mutant %>% 
  ggplot(aes(x = reorder(strain, rep(.[.$n_inds=="Single",]$speed, each = 2), na.rm = TRUE), y = speed, col = n_inds)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0) +
  stat_summary(fun.data = "mean_se", geom = "point") +
  scale_y_continuous(breaks = seq(0, 10, 2), expand = c(0, 0)) +
  scale_color_manual(values = c("#aaaaa9", "#b3343a")) +
  xlab("Strain") +
  ylab("Moving speed (mm/s)") +
  coord_cartesian(ylim = c(0, 10)) +
  facet_wrap( ~ sex, nrow = 2, strip.position = "left") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.1,0.9),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        strip.text = element_text(size = 11),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.7, "lines"))
# gig_f5min_speed_sex

###### speed single vs group ######
df_f5min_speed_wo_mutant %>%
  group_by(strain, sex, n_inds) %>%
  dplyr::summarize(speed = mean(speed, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = n_inds, values_from = speed, names_prefix = "speed_", names_sep = "_") %>%
  mutate(plot_col = if_else(speed_Group / speed_Single > 1, "over", "under")) %>%
  group_by(sex, plot_col) %>%
  summarize(n = n())

gig_f5min_speed_vs_group <- df_f5min_speed_wo_mutant %>%
  group_by(strain, sex, n_inds) %>%
  dplyr::summarize(speed = mean(speed, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = n_inds, values_from = speed, names_prefix = "speed_", names_sep = "_") %>%
  mutate(plot_col = if_else(speed_Group / speed_Single > 1, "over", "under")) %>%
  ggplot(aes(x = speed_Single, y = speed_Group)) +
  geom_abline(slope = 1, linetype = "dashed") +
  stat_smooth(linewidth = 2, color= "grey", method = "lm") +
  geom_point(aes(col = plot_col, shape = sex), alpha=0.5, size=4) +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(
                          after_stat(rr.label),
                          after_stat(p.value.label),
                          sep = "~~~")),
                        label.x = "left",
                        label.y = "top",
                        parse = TRUE, size = 4) +
  scale_x_continuous(breaks = seq(0, 10, 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 10, 2), expand = c(0, 0)) +
  scale_color_manual(values = c("#b3343a", "#aaaaa9")) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
  xlab("Moving speed (mm/s)\nSingle flies") +
  ylab("Moving speed (mm/s)\nGroup flies") +
  facet_wrap(~ sex, nrow = 1) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.spacing = unit(0.7, "lines"))
# gig_f5min_speed_vs_group

###### speed male vs female ######
df_f5min_speed_wo_mutant %>%
  group_by(strain, sex, n_inds) %>%
  dplyr::summarize(speed = mean(speed, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = sex, values_from = speed, names_prefix = "speed_", names_sep = "_") %>%
  mutate(plot_col = if_else(speed_Female / speed_Male > 1, "over", "under")) %>%
  group_by(n_inds, plot_col) %>%
  summarize(n = n())

gig_f5min_speed_vs_sex <- df_f5min_speed_wo_mutant %>%
  group_by(strain, n_inds, sex) %>%
  dplyr::summarize(speed = mean(speed, na.rm=T)) %>%
  pivot_wider(names_from = sex, values_from = speed, names_prefix = "speed_", names_sep = "_") %>%
  mutate(plot_col = if_else(speed_Male / speed_Female < 1, "over", "under")) %>%
  ggplot(aes(x = speed_Male, y = speed_Female)) +
  geom_abline(slope = 1, linetype = "dashed") +
  stat_smooth(linewidth = 2, color= "grey", method = "lm") +
  geom_point(aes(col = plot_col), alpha = 0.5, size = 4, shape = 16) +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(
                          after_stat(rr.label),
                          after_stat(p.value.label),
                          sep = "~~~")),
                        label.x = "left",
                        label.y = "top",
                        parse = TRUE, size = 4) +
  scale_x_continuous(breaks = seq(0, 10, 2), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 10, 2), expand = c(0, 0)) +
  scale_color_manual(values = c("#c17181", "#68aac3")) +
  coord_cartesian(xlim = c(0, 10), ylim = c(0, 10)) +
  xlab("Moving speed (mm/s)\nMale") +
  ylab("Moving speed (mm/s)\nFemale") +
  facet_wrap(~ n_inds, ncol = 2) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.spacing = unit(0.7, "lines"))
# gig_f5min_speed_vs_sex


gig_f5min_speed <-
  gig_f5min_speed_sex /
  (gig_f5min_speed_vs_group | gig_f5min_speed_vs_sex)
# gig_f5min_speed

ggsave("../figures/FigureS1.pdf", gig_f5min_speed, w = 6, h = 6)

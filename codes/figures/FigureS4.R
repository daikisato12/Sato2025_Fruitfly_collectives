#### load packages ####
targetPackages <- c('tidyverse','arrow','car','lmerTest','ggpmisc','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### load dataset ####
df_s5min_freezing_duration <- read_parquet("../data/1_single_strain/parquet/df_s5min_2995_5990_freezing_duration.parquet") %>%
  ungroup()

#### analysis ####
var <- "freezing_duration"
r <- 20

##### LMM #####
model_freezing_duration <- lmerTest::lmer(get(var) ~ strain * sex * n_inds + 
                                (1|date) + (1|time) + (1|prefix) + (1|place), 
                              df_s5min_freezing_duration %>%
                                filter(!str_detect(strain ,"norpA")))
summary(model_freezing_duration)
car::Anova(model_freezing_duration)

##### heritability single female #####
res_single_female <- anova(lm(get(var) ~ strain, data = df_s5min_freezing_duration %>% 
                                filter(sex == "Female", n_inds == "Single", !str_detect(strain ,"norpA"))))
sprintf("%.2e", res_single_female$`Pr(>F)`[1])
sg2_single_female <- (res_single_female$`Mean Sq`[1] - res_single_female$`Mean Sq`[2]) / r
se2_single_female <- res_single_female$`Mean Sq`[2]
sg2_single_female / (sg2_single_female + se2_single_female/r)

##### heritability group female #####
res_group_female <- anova(lm(get(var) ~ strain, data = df_s5min_freezing_duration %>% 
                               filter(sex == "Female", n_inds == "Group", !str_detect(strain ,"norpA"))))
sprintf("%.2e", res_group_female$`Pr(>F)`[1])
sg2_group_female <- (res_group_female$`Mean Sq`[1] - res_group_female$`Mean Sq`[2]) / r
se2_group_female <- res_group_female$`Mean Sq`[2]
sg2_group_female / (sg2_group_female + se2_group_female/r)

##### heritability single male #####
res_single_male <- anova(lm(get(var) ~ strain, data = df_s5min_freezing_duration %>% 
                              filter(sex == "Male", n_inds == "Single", !str_detect(strain ,"norpA"))))
sprintf("%.2e", res_single_male$`Pr(>F)`[1])
sg2_single_male <- (res_single_male$`Mean Sq`[1] - res_single_male$`Mean Sq`[2]) / r
se2_single_male <- res_single_male$`Mean Sq`[2]
sg2_single_male / (sg2_single_male + se2_single_male/r)

##### heritability group male #####
res_group_male <- anova(lm(get(var) ~ strain, data = df_s5min_freezing_duration %>% 
                             filter(sex == "Male", n_inds == "Group", !str_detect(strain ,"norpA"))))
sprintf("%.2e", res_group_male$`Pr(>F)`[1])
sg2_group_male <- (res_group_male$`Mean Sq`[1] - res_group_male$`Mean Sq`[2]) / r
se2_group_male <- res_group_male$`Mean Sq`[2]
sg2_group_male / (sg2_group_male + se2_group_male/r)


#### make plot ####
##### Figure S4a #####
gig_s5min_freezing_duration_sex <- df_s5min_freezing_duration %>% 
  filter(!str_detect(strain ,"norpA")) %>%
  ggplot(aes(x = reorder(strain, rep(.[.$n_inds=="Single",]$freezing_duration, each=2), na.rm = TRUE), y = freezing_duration, col = n_inds)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0) +
  stat_summary(fun.data = "mean_se", geom = "point") +
  scale_y_continuous(breaks = seq(0, 15, 5), expand = c(0, 0)) +
  scale_color_manual(values = c("#aaaaa9", "#b3343a")) +
  xlab("Strain") +
  ylab("Freezing duration (s)") +
  coord_cartesian(ylim = c(0, 15)) +
  facet_wrap( ~ sex, nrow = 2, strip.position = "left") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.9,0.65),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        strip.text = element_text(size = 11),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.7, "lines"))
gig_s5min_freezing_duration_sex

##### Figure S4b #####
###### freezing_duration single vs group ######
df_s5min_freezing_duration %>%
  filter(!str_detect(strain ,"norpA")) %>%
  group_by(strain, sex, n_inds) %>%
  dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = n_inds, values_from = freezing_duration, names_prefix = "freezing_duration_", names_sep = "_") %>%
  mutate(plot_col = if_else(freezing_duration_Group / freezing_duration_Single > 1, "over", "under")) %>%
  group_by(sex, plot_col) %>%
  summarize(n = n())

# stat
df_figS4b1 <- df_s5min_freezing_duration %>%
  filter(!str_detect(strain ,"norpA")) %>%
  group_by(strain, sex, n_inds) %>%
  dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm = T)) %>%
  ungroup() %>%
  pivot_wider(names_from = n_inds, values_from = freezing_duration, names_prefix = "freezing_duration_", names_sep = "_") %>%
  mutate(plot_col = if_else(freezing_duration_Group / freezing_duration_Single > 1, "over", "under"))
df_figS4b1_stats <- df_figS4b1 %>%
  group_by(sex) %>%
  dplyr::summarize(
    r = cor(freezing_duration_Single, freezing_duration_Group, use = "complete.obs"),
    cor_test = list(cor.test(freezing_duration_Single, freezing_duration_Group, use = "complete.obs")),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_value = map_dbl(cor_test, ~ .x$p.value),
    r_squared = r**2,
    p_formatted = sprintf("%.2e", p_value)
  ) %>%
  dplyr::select(sex, r_squared, p_formatted) %>%
  print()

gig_s5min_freezing_duration_vs_group <- df_figS4b1 %>%
  ggplot(aes(x = freezing_duration_Single, y = freezing_duration_Group)) +
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
  scale_x_continuous(breaks = seq(0, 15, 5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 15, 5), expand = c(0, 0)) +
  scale_color_manual(values = c("#b3343a", "#aaaaa9")) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 15)) +
  xlab("Freezing duration (s)\nSingle flies") +
  ylab("Freezing duration (s)\nGroup flies") +
  facet_wrap(~ sex, nrow = 1) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.spacing = unit(0.7, "lines"))
gig_s5min_freezing_duration_vs_group

###### freezing_duration male vs female ######
df_s5min_freezing_duration %>%
  filter(!str_detect(strain ,"norpA")) %>%
  group_by(strain, n_inds, sex) %>%
  dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm=T)) %>%
  pivot_wider(names_from = sex, values_from = freezing_duration, names_prefix = "freezing_duration_", names_sep = "_") %>%
  mutate(plot_col = if_else(freezing_duration_Male / freezing_duration_Female < 1, "over", "under")) %>%
  group_by(n_inds, plot_col) %>%
  summarize(n = n())

# stat
df_figS4b2 <- df_s5min_freezing_duration %>%
  filter(!str_detect(strain ,"norpA")) %>%
  group_by(strain, n_inds, sex) %>%
  dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm=T)) %>%
  pivot_wider(names_from = sex, values_from = freezing_duration, names_prefix = "freezing_duration_", names_sep = "_") %>%
  mutate(plot_col = if_else(freezing_duration_Male / freezing_duration_Female < 1, "over", "under"))
df_figS4b2_stats <- df_figS4b2 %>%
  group_by(n_inds) %>%
  dplyr::summarize(
    r = cor(freezing_duration_Male, freezing_duration_Female, use = "complete.obs"),
    cor_test = list(cor.test(freezing_duration_Male, freezing_duration_Female, use = "complete.obs")),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_value = map_dbl(cor_test, ~ .x$p.value),
    r_squared = r**2,
    p_formatted = sprintf("%.2e", p_value)
  ) %>%
  dplyr::select(n_inds, r_squared, p_formatted) %>%
  print()

gig_s5min_freezing_duration_vs_sex <- df_figS4b2 %>%
  ggplot(aes(x = freezing_duration_Male, y = freezing_duration_Female)) +
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
  scale_x_continuous(breaks = seq(0, 15, 5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 15, 5), expand = c(0, 0)) +
  scale_color_manual(values = c("#c17181", "#68aac3")) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 15)) +
  xlab("Freezing duration (s)\nMale") +
  ylab("Freezing duration (s)\nFemale") +
  facet_wrap(~ n_inds, ncol = 2) +
  theme_bw() +
  theme(legend.position = 'none',
        strip.background = element_blank(),
        panel.spacing = unit(0.7, "lines"))
gig_s5min_freezing_duration_vs_sex


gig_s5min_freezing_duration <-
  gig_s5min_freezing_duration_sex /
  (gig_s5min_freezing_duration_vs_group | gig_s5min_freezing_duration_vs_sex) #+
gig_s5min_freezing_duration

ggsave("../figures/FigureS4.pdf", gig_s5min_freezing_duration, w = 6, h = 6)

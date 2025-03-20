#### load packages ####
targetPackages <- c('tidyverse','arrow','car','lmerTest','ggpmisc','patchwork','ggrepel')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### load dataset ####
df_f5min_group_male_seconds_nest_dist_chase_strain_mean <- read_parquet("../data/1_single_strain/parquet/df_f5min_group_male_seconds_nest_dist_chase_3mm_strain_mean.parquet") %>%
  ungroup()


#### analysis time_chase ####
var <- "time_chase"
r <- 20

##### LMM #####
model_time_chase <- lmerTest::lmer(get(var) ~ strain + 
                                     (1|date) + (1|time) + (1|prefix) + (1|place), 
                                   df_f5min_group_male_seconds_nest_dist_chase_strain_mean %>%
                                     mutate(date = str_sub(prefix, start=1, end=8),
                                            time = as.POSIXct(prefix, format = format('%Y%m%d%H%M%S')) %>% format("%H:%M")) %>%
                                     filter(!str_detect(strain ,"norpA")))
summary(model_time_chase)
car::Anova(model_time_chase)


##### heritability #####
res_time_chase <- anova(lm(get(var) ~ strain, data = df_f5min_group_male_seconds_nest_dist_chase_strain_mean %>% 
                             filter(!str_detect(strain, "norpA"))))
sg2_time_chase <- (res_time_chase$`Mean Sq`[1] - res_time_chase$`Mean Sq`[2]) / r
se2_time_chase <- res_time_chase$`Mean Sq`[2]
sg2_time_chase / (sg2_time_chase + se2_time_chase/r)


#### analysis time_chain ####
var <- "time_chain"

##### LMM #####
model_time_chain <- lmerTest::lmer(get(var) ~ strain + 
                                     (1|date) + (1|time) + (1|prefix) + (1|place), 
                                   df_f5min_group_male_seconds_nest_dist_chase_strain_mean %>%
                                     mutate(date = str_sub(prefix, start=1, end=8),
                                            time = as.POSIXct(prefix, format = format('%Y%m%d%H%M%S')) %>% format("%H:%M")) %>%
                                     filter(!str_detect(strain ,"norpA")))
summary(model_time_chain)
car::Anova(model_time_chain)


##### heritability #####
res_time_chain <- anova(lm(get(var) ~ strain, data = df_f5min_group_male_seconds_nest_dist_chase_strain_mean %>% 
                             filter(!str_detect(strain, "norpA"))))
sg2_time_chain <- (res_time_chain$`Mean Sq`[1] - res_time_chain$`Mean Sq`[2]) / r
se2_time_chain <- res_time_chain$`Mean Sq`[2]
sg2_time_chain / (sg2_time_chain + se2_time_chain/r)



#### make plot ####
##### gg_f5min_chasing_chaining_time #####
gg_f5min_chasing_chaining_time <- df_f5min_group_male_seconds_nest_dist_chase_strain_mean %>% 
  pivot_longer(cols = c(time_chase, time_chain), names_to = "var", values_to = "value") %>%
  filter(!str_detect(strain, "norpA")) %>%
  # filter(strain != "norpA") %>%
  # mutate(plot_col = case_when(strain == "DGRP208_norpA" ~ "#8c6c81", #"#A8DADC", 
  #                             TRUE ~ "black")) %>%
  mutate(var = case_when(var == "time_chase" ~ "Chasing duration (s)",
                         TRUE ~ "Chaining duration (s)")) %>%
  ggplot(aes(x = reorder(strain, value), y = value)) + #, col = plot_col
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0) +
  stat_summary(fun.data = "mean_se", geom = "point") +
  scale_y_continuous(breaks = seq(0, 30, 5), expand = c(0, 0)) +
  xlab("Strain") +
  # ylab("Duration (s)") +
  coord_cartesian(ylim = c(0, 30)) +
  # scale_color_identity() +
  facet_wrap( ~ var, nrow = 2, strip.position = "left") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        legend.background = element_blank(),
        legend.title = element_blank(),
        legend.position = c(0.9,0.65),
        legend.key = element_blank(),
        legend.key.size = unit(0.4, "cm"),
        strip.text = element_text(size = 11),
        strip.background = element_blank(),
        strip.placement = "outside",
        panel.spacing = unit(0.7, "lines"))
gg_f5min_chasing_chaining_time

##### gg_f5min_group_male_chase_chain_cor #####
# formula <- y ~ poly(x, 2, raw = TRUE)
formula <- y ~ x
gg_f5min_group_male_chase_chain_cor <- 
  df_f5min_group_male_seconds_nest_dist_chase_strain_mean %>% 
  filter(!str_detect(strain, "norpA")) %>%
  group_by(strain) %>%
  summarize(time_chase = mean(time_chase, na.rm = TRUE),
            time_chain = mean(time_chain, na.rm = TRUE)) %>%
  ggplot(aes(x = time_chase, y = time_chain, label = strain)) +
  ggpmisc::stat_poly_line(formula = formula,
                        color= "grey") +
  # ggpmisc::stat_poly_eq(aes(label = after_stat(eq.label)), 
  #                       formula = y ~ x, #poly(x, 2)
  #                       label.x = "left", label.y = 0.9, 
  #                       parse = TRUE, size = 2) +
  ggpmisc::stat_poly_eq(aes(label = paste(after_stat(rr.label),
                                          "*\", \"*", 
                                          after_stat(p.value.label))), 
                        formula = y ~ x, #poly(x, 2)
                        label.x = "left", label.y = 0.9, 
                        parse = TRUE, size = 2) +
  geom_point(shape = 16, alpha=0.5, size = 3) +
  geom_text_repel(size = 2) +
  xlab("Chasing duration (sec)") +
  ylab("Chaining duration (sec)") +
  theme_bw() +
  theme(
    legend.background = element_blank(),
    legend.title = element_blank(),
    legend.position = 'none',
    legend.key = element_blank())
gg_f5min_group_male_chase_chain_cor


gg_f5min_chasing_chaining <-
  gg_f5min_chasing_chaining_time +
  gg_f5min_group_male_chase_chain_cor +
  plot_layout(ncol = 2, widths = c(1, 0.5), guides = "collect")
gg_f5min_chasing_chaining

ggsave("../figures/FigureS3.pdf", gg_f5min_chasing_chaining, w = 6, h = 3)

#### load packages ####
targetPackages <- c('tidyverse','arrow','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### Figure S15a ####
##### load dataset #####
df_s5min_2995_5990_freezing_duration <- read_parquet("../data/1_single_strain/parquet/df_s5min_2995_5990_freezing_duration.parquet") %>%
  ungroup()

dfd_s5min_2995_5990_speed_d_normbyf5minave_strain_trial <- 
  read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_d_normbyf5minave_strain_trial.parquet") %>%
  ungroup()

list_strain_d <- read.csv("../data/2_mixed_strain/list_strain_d.txt") %>%
  pull(x)

dfd_s5min_2995_5990_gd_freezing_duration_dist <- 
  read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_gd_freezing_duration_dist.parquet") %>%
  ungroup()

##### stat #####
df_figS15a_stats <- dfd_s5min_2995_5990_gd_freezing_duration_dist %>%
  group_by(thr_sec) %>%
  dplyr::summarize(
    r = cor(freezing_duration_dist, performance, use = "complete.obs"),
    cor_test = list(cor.test(freezing_duration_dist, performance, use = "complete.obs")),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_value = map_dbl(cor_test, ~ .x$p.value),
    r_squared = r**2,
    p_formatted = sprintf("%.2e", p_value)
  ) %>%
  dplyr::select(thr_sec, r_squared, p_formatted) %>%
  print()

##### make plot #####
g_s5min_2995_5990_gd_freezing_duration_dist_facet <-
  ggplot(dfd_s5min_2995_5990_gd_freezing_duration_dist,
         aes(x = freezing_duration_dist, 
             y = performance)) +
  geom_point() +
  stat_smooth(linewidth = 1, color= "grey", method = "lm") +#, formula = y ~ poly(x, degree = 2, raw = TRUE) - 1) + #formula = y ~ log(x)) + #method = "lm") +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(#stat(eq.label),
                          after_stat(rr.label),
                          sep = "~~~")),
                        label.x = "right",
                        label.y = 0.95,
                        parse = TRUE, size = 2.5) +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(#stat(eq.label),
                          stat(p.value.label))),
                        label.x = "right",
                        label.y = 0.87,
                        parse = TRUE, size = 2.5) +
  xlab("Difference in freezing duration (s)") +
  ylab("Behavioral performance") +
  facet_wrap(~ thr_sec, nrow = 2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9),
        axis.text = element_text(color = "black"))
g_s5min_2995_5990_gd_freezing_duration_dist_facet


#### Figure S15b ####
##### load dataset #####
dfd_s5min_2995_5990_gd_motion_cue_exit_intercept_dist <- 
  read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_gd_motion_cue_exit_intercept_dist.parquet") %>%
  ungroup()

##### stat #####
df_figS15b_stats <- dfd_s5min_2995_5990_gd_motion_cue_exit_intercept_dist %>%
  group_by(thr_sec) %>%
  dplyr::summarize(
    r = cor(motion_cue_exit_intercept_dist, performance, use = "complete.obs"),
    cor_test = list(cor.test(motion_cue_exit_intercept_dist, performance, use = "complete.obs")),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_value = map_dbl(cor_test, ~ .x$p.value),
    r_squared = r**2,
    p_formatted = sprintf("%.2e", p_value)
  ) %>%
  dplyr::select(thr_sec, r_squared, p_formatted) %>%
  print()

##### make plot #####
g_s5min_2995_5990_gd_motion_cue_exit_intercept_dist_facet <-
  ggplot(dfd_s5min_2995_5990_gd_motion_cue_exit_intercept_dist,
         aes(x = motion_cue_exit_intercept_dist, 
             y = performance)) +
  geom_point() +
  stat_smooth(linewidth = 1, color= "grey", method = "lm") +#, formula = y ~ poly(x, degree = 2, raw = TRUE) - 1) + #formula = y ~ log(x)) + #method = "lm") +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(#stat(eq.label),
                          after_stat(rr.label),
                          sep = "~~~")),
                        label.x = "right",
                        label.y = 0.95,
                        parse = TRUE, size = 2.5) +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(#stat(eq.label),
                          stat(p.value.label))),
                        label.x = "right",
                        label.y = 0.87,
                        parse = TRUE, size = 2.5) +
  xlab("Difference in visual responsiveness") +
  ylab("Behavioral performance") +
  facet_wrap(~ thr_sec, nrow = 2) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 9),
        axis.text = element_text(color = "black"))
g_s5min_2995_5990_gd_motion_cue_exit_intercept_dist_facet


g_figS15 <- g_s5min_2995_5990_gd_freezing_duration_dist_facet +
  g_s5min_2995_5990_gd_motion_cue_exit_intercept_dist_facet

ggsave("../figures/FigureS15.pdf", g_figS15, w=10, h=5)

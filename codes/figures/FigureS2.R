#### load packages ####
library(tidyverse)
library(arrow)
library(lmerTest)
library(patchwork)


#### load dataset ####
df_f5min_nnd_rand <- read_parquet("../data/1_single_strain/df_f5min_nnd_rand.parquet") %>%
  ungroup()

#### analysis ####
model_nnd <- lmerTest::lmer(nnd ~ strain * sex +  
                              (1|date) + (1|time) + (1|prefix) + (1|place), 
                            df_f5min_nnd_rand %>%
                              filter(str_detect(strain ,"DGRP")) %>%
                              dplyr::mutate(strain = as.factor(strain), #must be factor to use difflsmeans function
                                     sex = as.factor(sex)))
summary(model_nnd)
car::Anova(model_nnd)

##### male #####
###### LMM ######
model_nnd_male <- lmerTest::lmer(nnd ~ strain +  
                                   (1|date) + (1|time) + (1|prefix) + (1|place), 
                                 df_f5min_nnd_rand %>%
                                   filter(sex == "Male") %>%
                                   dplyr::mutate(strain = as.factor(strain))) #must be factor to use difflsmeans function
summary(model_nnd_male)
car::Anova(model_nnd_male)

###### calculate statistical difference ######
df_difflsmeans_nnd_male <- lmerTest::difflsmeans(model_nnd_male) %>% #multiple correction based on the Satterthwaite's approximation to degrees of freedom
  tibble::rownames_to_column() %>%
  as_tibble()

df_difflsmeans_nnd_male2 <- df_difflsmeans_nnd_male %>%
  dplyr::mutate(rowname = str_remove_all(rowname, "strain")) %>%
  filter(str_detect(rowname, "norpA|random")) %>%
  separate(rowname, sep = " - ", into = c("strain1", "strain2")) 

df_difflsmeans_nnd_male2_norpA <- df_difflsmeans_nnd_male2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "norpA", strain2, strain1),
                strain2 = if_else(strain1 == "norpA", "norpA", strain2),
                Estimate = if_else(strain1 == "norpA", Estimate, -Estimate),
                strain1 = if_else(strain1 == "norpA", strain_tmp, strain1)) %>%
  filter(strain2 == "norpA") %>%
  bind_rows(data.frame(strain1 = "norpA",
                       strain2 = "norpA"))

df_difflsmeans_nnd_male2_random <- df_difflsmeans_nnd_male2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "random", strain2, strain1),
                strain2 = if_else(strain1 == "random", "random", strain2),
                Estimate = if_else(strain1 == "random", Estimate, -Estimate),
                strain1 = if_else(strain1 == "random", strain_tmp, strain1)) %>%
  filter(strain2 == "random") %>%
  bind_rows(data.frame(strain1 = "random",
                       strain2 = "random"))

df_f5min_nnd_male_order <- df_f5min_nnd_rand %>% 
  filter(sex == "Male") %>%
  group_by(strain) %>%
  dplyr::summarize(nnd_mean = mean(nnd, na.rm = T)) %>%
  ungroup() %>%
  arrange(nnd_mean)

###### calculate heritability ######
var <- "nnd"
r <- 20

res_single_female <- anova(lm(get(var) ~ strain, data = df_f5min_nnd_rand %>% 
                                filter(sex == "Female", str_detect(strain, "DGRP"))))
sg2_single_female <- (res_single_female$`Mean Sq`[1] - res_single_female$`Mean Sq`[2]) / r
se2_single_female <- res_single_female$`Mean Sq`[2]
sg2_single_female / (sg2_single_female + se2_single_female/r)


##### female #####
###### LMM ######
model_nnd_female <- lmerTest::lmer(nnd ~ strain +  
                                   (1|date) + (1|time) + (1|prefix) + (1|place), 
                                 df_f5min_nnd_rand %>%
                                   filter(sex == "Female") %>%
                                   dplyr::mutate(strain = as.factor(strain))) #must be factor to use difflsmeans function
summary(model_nnd_female)
car::Anova(model_nnd_female)

###### calculate statistical difference ######
df_difflsmeans_nnd_female <- lmerTest::difflsmeans(model_nnd_female) %>% #multiple correction based on the Satterthwaite's approximation to degrees of freedom
  tibble::rownames_to_column() %>%
  as_tibble()

df_difflsmeans_nnd_female2 <- df_difflsmeans_nnd_female %>%
  dplyr::mutate(rowname = str_remove_all(rowname, "strain")) %>%
  filter(str_detect(rowname, "norpA|random")) %>%
  separate(rowname, sep = " - ", into = c("strain1", "strain2")) 

df_difflsmeans_nnd_female2_norpA <- df_difflsmeans_nnd_female2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "norpA", strain2, strain1),
                strain2 = if_else(strain1 == "norpA", "norpA", strain2),
                Estimate = if_else(strain1 == "norpA", Estimate, -Estimate),
                strain1 = if_else(strain1 == "norpA", strain_tmp, strain1)) %>%
  filter(strain2 == "norpA") %>%
  bind_rows(data.frame(strain1 = "norpA",
                       strain2 = "norpA"))

df_difflsmeans_nnd_female2_random <- df_difflsmeans_nnd_female2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "random", strain2, strain1),
                strain2 = if_else(strain1 == "random", "random", strain2),
                Estimate = if_else(strain1 == "random", Estimate, -Estimate),
                strain1 = if_else(strain1 == "random", strain_tmp, strain1)) %>%
  filter(strain2 == "random") %>%
  bind_rows(data.frame(strain1 = "random",
                       strain2 = "random"))

df_f5min_nnd_female_order <- df_f5min_nnd_rand %>% 
  filter(sex == "Female") %>%
  group_by(strain) %>%
  dplyr::summarize(nnd_mean = mean(nnd, na.rm = T)) %>%
  ungroup() %>%
  arrange(nnd_mean)

###### calculate heritability ######
var <- "nnd"
r <- 20

res_single_male <- anova(lm(get(var) ~ strain, data = df_f5min_nnd_rand %>% 
                                filter(sex == "Male", str_detect(strain, "DGRP"))))
sg2_single_male <- (res_single_male$`Mean Sq`[1] - res_single_male$`Mean Sq`[2]) / r
se2_single_male <- res_single_male$`Mean Sq`[2]
sg2_single_male / (sg2_single_male + se2_single_male/r)


#### make plot ####
##### S2a1 male #####
gg_f5min_nnd_male <- df_f5min_nnd_rand %>% 
  mutate(plot_col = case_when(strain == "random" ~ "grey", 
                              strain == "norpA" ~ "magenta", 
                              TRUE ~ "black")) %>%
  filter(sex =="Male") %>%
  ggplot(aes(x = reorder(strain, nnd, na.rm = TRUE), y = nnd, col = plot_col)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0) +
  stat_summary(fun.data = "mean_se", geom = "point") +
  coord_cartesian(ylim = c(4, 8)) +
  scale_color_manual(values = c("black", "#25b7c0", "magenta")) +
  xlab("Strain") +
  ylab("NND (mm)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = 'none')
gg_f5min_nnd_male

gg_f5min_nnd_male_sig <- ggplot(bind_rows(df_difflsmeans_nnd_male2_norpA, 
                                          df_difflsmeans_nnd_male2_random) %>%
                                  add_column(y = "value") %>%
                                  rename(P = `Pr(>|t|)`) %>%
                                  mutate(strain2 = str_replace(strain2, "norpA", "vs. norpA"),
                                         strain2 = str_replace(strain2, "random", "vs. random")) %>%
                                  transform(strain1 = factor(strain1, 
                                                             levels = df_f5min_nnd_male_order$strain)),
                                aes(x = strain1, y = y, col = Estimate, size = -log10(P))) +
  geom_point() +
  scale_colour_gradient2(low = "#1e50a2", mid = "#e5e4e6", high = "#c53d43",
                         limits = c(-1.1, 2.4), midpoint = 0) + 
  scale_size_continuous(limits = c(0, 6)) +
  facet_wrap(~ strain2, nrow = 2) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.background = element_blank())
gg_f5min_nnd_male_sig

# gg_f5min_nnd_male_norpA_sig
# gg_f5min_nnd_male_random_sig

gg_f5min_nnd_male_sum <- 
  gg_f5min_nnd_male_sig +
  gg_f5min_nnd_male + 
  plot_layout(nrow = 2, ncol = 1, heights = c(0.5,1))
gg_f5min_nnd_male_sum
ggsave("../figures/FigureS2a_1.pdf", gg_f5min_nnd_male_sum, w = 6, h = 3)

##### S2a2 female #####
gg_f5min_nnd_female <- df_f5min_nnd_rand %>% 
  mutate(plot_col = case_when(strain == "random" ~ "grey", 
                              strain == "norpA" ~ "magenta", 
                              TRUE ~ "black"),
         type = if_else(plot_col == "black", "DGRP strains", "*")) %>%
  filter(sex =="Female") %>%
  ggplot(aes(x = reorder(strain, nnd, na.rm = TRUE), y = nnd, col = plot_col)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0) +
  stat_summary(fun.data = "mean_se", geom = "point") +
  coord_cartesian(ylim = c(4, 8)) +
  scale_color_manual(values = c("black", "#25b7c0", "magenta")) +
  xlab("Strain") +
  ylab("NND (mm)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = 'none')
gg_f5min_nnd_female

gg_f5min_nnd_female_sig <- ggplot(bind_rows(df_difflsmeans_nnd_female2_norpA, 
                                          df_difflsmeans_nnd_female2_random) %>%
                                  add_column(y = "value") %>%
                                  rename(P = `Pr(>|t|)`) %>%
                                  mutate(strain2 = str_replace(strain2, "norpA", "vs. norpA"),
                                         strain2 = str_replace(strain2, "random", "vs. random")) %>%
                                  transform(strain1 = factor(strain1, 
                                                             levels = df_f5min_nnd_female_order$strain)),
                                aes(x = strain1, y = y, col = Estimate, size = -log10(P))) +
  geom_point() +
  scale_colour_gradient2(low = "#1e50a2", mid = "#e5e4e6", high = "#c53d43",
                         limits = c(-1.1, 2.4), midpoint = 0) + 
  scale_size_continuous(limits = c(0, 6)) +
  facet_wrap(~ strain2, nrow = 2) +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(), 
        strip.background = element_blank())#,
gg_f5min_nnd_female_sig

# gg_f5min_nnd_female_norpA_sig
# gg_f5min_nnd_female_random_sig

gg_f5min_nnd_female_sum <- 
  gg_f5min_nnd_female_sig +
  gg_f5min_nnd_female + 
  patchwork::plot_layout(nrow = 2, ncol = 1, heights = c(0.5,1))
gg_f5min_nnd_female_sum
ggsave("../figures/FigureS2a_2.pdf", gg_f5min_nnd_female_sum, w = 6, h = 3)


##### S2b #####
gg_f5min_nnd_vs_sex <- df_f5min_nnd_rand %>%
  filter(str_detect(strain, "DGRP")) %>%
  group_by(strain, sex) %>%
  dplyr::summarize(nnd = mean(nnd, na.rm = T)) %>%
  pivot_wider(names_from = sex, values_from = nnd, names_prefix = "nnd_", names_sep = "_") %>%
  mutate(plot_col = if_else(nnd_Male / nnd_Female < 1, "over", "under")) %>%
  ggplot(aes(x = nnd_Male, y = nnd_Female)) +
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
  scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_color_manual(values = c("#c17181", "#68aac3")) +
  xlab("NND (mm)\nMale") +
  ylab("NND (mm)\nFemale") +
  theme_bw() +
  theme(legend.position = 'none')
gg_f5min_nnd_vs_sex

ggsave("../figures/FigureS2b.pdf", gg_f5min_nnd_vs_sex, w = 3, h = 3)

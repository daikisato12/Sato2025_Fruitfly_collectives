#### load packages ####
targetPackages <- c('tidyverse','arrow','lmerTest','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### load dataset ####
df_f5min_nnd_rand <- read_parquet("../data/1_single_strain/parquet/df_f5min_nnd_rand.parquet") %>%
  ungroup()

#### analysis ####
var <- "nnd"
r <- 20

model_nnd <- lmerTest::lmer(get(var) ~ strain * sex +  
                              (1|date) + (1|time) + (1|prefix) + (1|place), 
                            df_f5min_nnd_rand %>%
                              filter(!strain %in% c("random", "norpA", "DGRP208_norpA")) %>%
                              dplyr::mutate(strain = as.factor(strain), #must be factor to use difflsmeans function
                                     sex = as.factor(sex)))
summary(model_nnd)
car::Anova(model_nnd)

##### male #####
###### LMM ######
model_nnd_male <- lmerTest::lmer(get(var) ~ strain +  
                                   (1|date) + (1|time) + (1|prefix) + (1|place), 
                                 df_f5min_nnd_rand %>%
                                   filter(sex == "Male",
                                          !strain %in% c("random", "norpA", "DGRP208_norpA")) %>%
                                   dplyr::mutate(strain = as.factor(strain))) #must be factor to use difflsmeans function
summary(model_nnd_male)
car::Anova(model_nnd_male)

###### calculate statistical difference ######
model_nnd_male_all <- 
  lmerTest::lmer(get(var) ~ strain + (1|date) + (1|time) + (1|prefix) + (1|place), 
                 df_f5min_nnd_rand %>%
                   filter(sex == "Male") %>%
                   dplyr::mutate(strain = as.factor(strain))) #must be factor to use difflsmeans function

df_difflsmeans_nnd_male <- lmerTest::difflsmeans(model_nnd_male_all) %>% #multiple correction based on the Satterthwaite's approximation to degrees of freedom
  tibble::rownames_to_column() %>%
  as_tibble()

df_difflsmeans_nnd_male2 <- df_difflsmeans_nnd_male %>%
  dplyr::mutate(rowname = str_remove_all(rowname, "strain")) %>%
  filter(str_detect(rowname, "random|norpA")) %>%
  separate(rowname, sep = " - ", into = c("strain1", "strain2")) 

df_difflsmeans_nnd_male2_norpA <- df_difflsmeans_nnd_male2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "norpA", strain2, strain1),
                strain2 = if_else(strain1 == "norpA", "norpA", strain2),
                Estimate = if_else(strain1 == "norpA", Estimate, -Estimate),
                strain1 = if_else(strain1 == "norpA", strain_tmp, strain1)) %>%
  filter(strain2 == "norpA") %>%
  bind_rows(data.frame(strain1 = "norpA",
                       strain2 = "norpA"))

df_difflsmeans_nnd_male2_208norpA <- df_difflsmeans_nnd_male2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "DGRP208_norpA", strain2, strain1),
                strain2 = if_else(strain1 == "DGRP208_norpA", "DGRP208_norpA", strain2),
                Estimate = if_else(strain1 == "DGRP208_norpA", Estimate, -Estimate),
                strain1 = if_else(strain1 == "DGRP208_norpA", strain_tmp, strain1)) %>%
  filter(strain2 == "DGRP208_norpA") %>%
  bind_rows(data.frame(strain1 = "DGRP208_norpA",
                       strain2 = "DGRP208_norpA"))

df_difflsmeans_nnd_male2_random <- df_difflsmeans_nnd_male2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "random", strain2, strain1),
                strain2 = if_else(strain1 == "random", "random", strain2),
                Estimate = if_else(strain1 == "random", Estimate, -Estimate),
                strain1 = if_else(strain1 == "random", strain_tmp, strain1)) %>%
  filter(strain2 == "random") %>%
  bind_rows(data.frame(strain1 = "random",
                       strain2 = "random"))

df_f5min_nnd_male_order <- df_f5min_nnd_rand %>% 
  filter(sex == "Male", 
         # !strain %in% c("random", "norpA", "DGRP208_norpA")) %>%
         !strain %in% c("norpA")) %>%
  group_by(strain) %>%
  dplyr::summarize(nnd_mean = mean(nnd, na.rm = T)) %>%
  ungroup() %>%
  arrange(nnd_mean)

###### calculate heritability ######
res_male <- anova(lm(get(var) ~ strain, data = df_f5min_nnd_rand %>% 
                              filter(!strain %in% c("random", "norpA", "DGRP208_norpA"),
                                     sex == "Male")))
sg2_male <- (res_male$`Mean Sq`[1] - res_male$`Mean Sq`[2]) / r
se2_male <- res_male$`Mean Sq`[2]
sg2_male / (sg2_male + se2_male/r)


##### female #####
###### LMM ######
model_nnd_female <- lmerTest::lmer(get(var) ~ strain +  
                                   (1|date) + (1|time) + (1|prefix) + (1|place), 
                                 df_f5min_nnd_rand %>%
                                   filter(sex == "Female",
                                          !strain %in% c("random", "norpA", "DGRP208_norpA")) %>%
                                   dplyr::mutate(strain = as.factor(strain))) #must be factor to use difflsmeans function
summary(model_nnd_female)
car::Anova(model_nnd_female)

###### calculate statistical difference ######
model_nnd_female_all <- 
  lmerTest::lmer(get(var) ~ strain + (1|date) + (1|time) + (1|prefix) + (1|place), 
                 df_f5min_nnd_rand %>%
                   filter(sex == "Female") %>%
                   dplyr::mutate(strain = as.factor(strain))) #must be factor to use difflsmeans function

df_difflsmeans_nnd_female <- lmerTest::difflsmeans(model_nnd_female_all) %>% #multiple correction based on the Satterthwaite's approximation to degrees of freedom
  tibble::rownames_to_column() %>%
  as_tibble()

df_difflsmeans_nnd_female2 <- df_difflsmeans_nnd_female %>%
  dplyr::mutate(rowname = str_remove_all(rowname, "strain")) %>%
  filter(str_detect(rowname, "random|norpA")) %>%
  separate(rowname, sep = " - ", into = c("strain1", "strain2")) 

df_difflsmeans_nnd_female2_norpA <- df_difflsmeans_nnd_female2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "norpA", strain2, strain1),
                strain2 = if_else(strain1 == "norpA", "norpA", strain2),
                Estimate = if_else(strain1 == "norpA", Estimate, -Estimate),
                strain1 = if_else(strain1 == "norpA", strain_tmp, strain1)) %>%
  filter(strain2 == "norpA") %>%
  bind_rows(data.frame(strain1 = "norpA",
                       strain2 = "norpA"))

df_difflsmeans_nnd_female2_208norpA <- df_difflsmeans_nnd_female2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "DGRP208_norpA", strain2, strain1),
                strain2 = if_else(strain1 == "DGRP208_norpA", "DGRP208_norpA", strain2),
                Estimate = if_else(strain1 == "DGRP208_norpA", Estimate, -Estimate),
                strain1 = if_else(strain1 == "DGRP208_norpA", strain_tmp, strain1)) %>%
  filter(strain2 == "DGRP208_norpA") %>%
  bind_rows(data.frame(strain1 = "DGRP208_norpA",
                       strain2 = "DGRP208_norpA"))

df_difflsmeans_nnd_female2_random <- df_difflsmeans_nnd_female2 %>%
  dplyr::mutate(strain_tmp = if_else(strain1 == "random", strain2, strain1),
                strain2 = if_else(strain1 == "random", "random", strain2),
                Estimate = if_else(strain1 == "random", Estimate, -Estimate),
                strain1 = if_else(strain1 == "random", strain_tmp, strain1)) %>%
  filter(strain2 == "random") %>%
  bind_rows(data.frame(strain1 = "random",
                       strain2 = "random"))

df_f5min_nnd_female_order <- df_f5min_nnd_rand %>% 
  filter(sex == "Female",
         # !strain %in% c("random", "norpA", "DGRP208_norpA")) %>%
         !strain %in% c("norpA")) %>%
  group_by(strain) %>%
  dplyr::summarize(nnd_mean = mean(nnd, na.rm = T)) %>%
  ungroup() %>%
  arrange(nnd_mean)

###### calculate heritability ######
res_female <- anova(lm(get(var) ~ strain, data = df_f5min_nnd_rand %>% 
                                filter(!strain %in% c("random", "norpA", "DGRP208_norpA"),
                                       sex == "Female")))
sg2_female <- (res_female$`Mean Sq`[1] - res_female$`Mean Sq`[2]) / r
se2_female <- res_female$`Mean Sq`[2]
sg2_female / (sg2_female + se2_female/r)


#### make plot ####
##### S2a1 female #####
gg_f5min_nnd_female <- df_f5min_nnd_rand %>% 
  filter(strain != "norpA") %>%
  mutate(plot_col = case_when(strain == "random" ~ "darkgrey", 
                              # strain == "norpA" ~ "#FFABAB", 
                              strain == "DGRP208" ~ "#833163",
                              strain == "DGRP208_norpA" ~ "#e8c72d", #"#A8DADC", 
                              TRUE ~ "black")) %>%
  # transform(plot_col = factor(plot_col, levels = c("black", "darkgrey", "#FFABAB", "#A8DADC"))) %>%
  filter(sex =="Female") %>%
  ggplot(aes(x = reorder(strain, nnd, na.rm = TRUE), y = nnd, col = plot_col)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0) +
  stat_summary(fun.data = "mean_se", geom = "point") +
  coord_cartesian(ylim = c(4, 8)) +
  # scale_color_manual(values = c("black", "#25b7c0", "magenta")) +
  scale_color_identity() +
  xlab("Strain") +
  ylab("NND (mm)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = 'none')
gg_f5min_nnd_female

gg_f5min_nnd_female_sig <- ggplot(bind_rows(df_difflsmeans_nnd_female2_208norpA, 
                                            df_difflsmeans_nnd_female2_random) %>%
                                    # bind_rows(df_difflsmeans_nnd_female2_norpA) %>%
                                    filter(strain1 != "norpA") %>%
                                    add_column(y = "value") %>%
                                    rename(P = `Pr(>|t|)`) %>%
                                    mutate(strain2 = paste0("vs. ", strain2)) %>%
                                    transform(strain1 = factor(strain1, 
                                                               levels = df_f5min_nnd_female_order$strain)),
                                  aes(x = strain1, y = y, col = Estimate, size = -log10(P))) +
  geom_point() +
  scale_colour_gradient2(low = "#1e50a2", mid = "#e5e4e6", high = "#c53d43",
                         limits = c(-1.1, 3), midpoint = 0) + 
  scale_size_continuous(limits = c(0, 8)) +
  facet_wrap(~ strain2, nrow = 3) +
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
  patchwork::plot_layout(nrow = 2, ncol = 1, heights = c(0.4,1))
gg_f5min_nnd_female_sum
ggsave("../figures/FigureS2a_1.pdf", gg_f5min_nnd_female_sum, w = 6, h = 3)

##### S2a2 male #####
gg_f5min_nnd_male <- df_f5min_nnd_rand %>% 
  filter(strain != "norpA") %>%
  mutate(plot_col = case_when(strain == "random" ~ "darkgrey", 
                              # strain == "norpA" ~ "#FFABAB", 
                              strain == "DGRP208" ~ "#833163",
                              strain == "DGRP208_norpA" ~ "#e8c72d", #"#A8DADC", 
                              TRUE ~ "black")) %>%
  # transform(plot_col = factor(plot_col, levels = c("black", "darkgrey", "#FFABAB", "#A8DADC"))) %>%
  filter(sex == "Male") %>%
  ggplot(aes(x = reorder(strain, nnd, na.rm = TRUE), y = nnd, col = plot_col)) +
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0) +
  stat_summary(fun.data = "mean_se", geom = "point") +
  coord_cartesian(ylim = c(4, 8)) +
  # scale_color_manual(values = c("black", "darkgrey", "#804882", "#315587")) +
  scale_color_identity() +
  xlab("Strain") +
  ylab("NND (mm)") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        legend.position = 'none')
gg_f5min_nnd_male

gg_f5min_nnd_male_sig <- ggplot(bind_rows(df_difflsmeans_nnd_male2_208norpA, 
                                          df_difflsmeans_nnd_male2_random) %>%
                                  # bind_rows(df_difflsmeans_nnd_male2_norpA) %>%
                                  filter(strain1 != "norpA") %>%
                                  add_column(y = "value") %>%
                                  rename(P = `Pr(>|t|)`) %>%
                                  mutate(strain2 = paste0("vs. ", strain2)) %>%
                                  transform(strain1 = factor(strain1, 
                                                             levels = df_f5min_nnd_male_order$strain)),
                                aes(x = strain1, y = y, col = Estimate, size = -log10(P))) +
  geom_point() +
  scale_colour_gradient2(low = "#1e50a2", mid = "#e5e4e6", high = "#c53d43",
                         limits = c(-1.1, 3), midpoint = 0, na.value = NA) + 
  scale_size_continuous(limits = c(0, 8)) +
  facet_wrap(~ strain2, nrow = 3) +
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
  plot_layout(nrow = 2, ncol = 1, heights = c(0.4,1))
gg_f5min_nnd_male_sum
ggsave("../figures/FigureS2a_2.pdf", gg_f5min_nnd_male_sum, w = 6, h = 3)


##### S2b #####
df_f5min_nnd_rand %>%
  filter(!strain %in% c("random", "norpA", "DGRP208_norpA")) %>%
  group_by(strain, sex) %>%
  dplyr::summarize(nnd = mean(nnd, na.rm=T)) %>%
  pivot_wider(names_from = sex, values_from = nnd, names_prefix = "nnd_", names_sep = "_") %>%
  mutate(plot_col = if_else(nnd_Male / nnd_Female < 1, "over", "under")) %>%
  group_by(plot_col) %>%
  summarize(n = n())

gg_f5min_nnd_vs_sex <- df_f5min_nnd_rand %>%
  filter(!strain %in% c("random", "norpA", "DGRP208_norpA")) %>%
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

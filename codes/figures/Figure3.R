#### load packages ####
library(tidyverse)
library(ggpubr)
library(arrow)
library(lmerTest)
library(patchwork)

#### load functions ####
glm_exp <- function(dat, gene){

  tryCatch({
    dat_tmp <- dat %>% 
      group_by(genotype2) %>%
      dplyr::summarize(n = n())
    if(length(unique(dat$genotype2)) > 1 & length(unique(dat %>% pull(get(gene)))) > 1){
      if(nrow(dat_tmp[dat_tmp$n > 2,]) > 1){
        dat_tmp2 <- dat %>% 
          group_by(rep) %>%
          dplyr::summarize(n = n())
        if(nrow(dat_tmp2[dat_tmp2$n > 2,]) > 1){
          myformula <- as.formula(paste0(gene, " ~ genotype2 + (1|genotype) + (1|rep)"))
          p.val <- lme4::glmer(myformula,
                         data = dat,
                         family = "Gamma") %>%
            car::Anova() %>%
            pull(`Pr(>Chisq)`)
          method <- "glmer"
        }else{
          myformula <- as.formula(paste0(gene, " ~ genotype2 + (1|genotype)"))
          p.val <- lme4::glmer(myformula,
                         data = dat,
                         family = "Gamma") %>%
            car::Anova() %>%
            pull(`Pr(>Chisq)`)
          method <- "glmer"
        }
      }else{
        myformula <- as.formula(paste0(gene, " ~ genotype2"))
        p.val <- glm(myformula,
                     data = dat,
                     family = "Gamma") %>%
          car::Anova() %>%
          pull(`Pr(>Chisq)`)
        method <- "glm"
      }
    }else{
      p.val <- NA_real_
      method = NA_character_
    }
  }, 
  error = function(e) {
    message(e)           
    p.val <- NA_real_
    method = NA_character_
  })
  dat2 <- dat %>%
    mutate(pval = p.val,
           method = method)
  return(dat2)
}


makeStars <- function(x){
  stars <- c("****", "***", "**", "*", NA_character_)
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

#### Figure 3a ####
##### load dataset #####
df_early.ptp99a <- read.table("../data/4_scRNAseq/df_early_ptp99a.tsv", header = TRUE)
df_main.ptp99a <- read.table("../data/4_scRNAseq/df_main_ptp99a.tsv", header = TRUE)

df_integrated.ptp99A <-
  df_main.ptp99a %>%
  mutate(dataset = "main") %>%
  bind_rows(df_early.ptp99a %>%
              mutate(dataset = "early")) %>%
  ungroup() %>%
  group_nest(dataset, type, time) %>%
  mutate(data = map2(data, "ptp99a", glm_exp)) %>% 
  unnest(cols = c(data)) %>%
  mutate(genotype = case_when(str_detect(genotype, "_") ~ paste0("DGRP", parse_number(genotype)),
                              TRUE ~ genotype)) %>%
  transform(genotype = factor(genotype, levels = unique(.$genotype) %>% gtools::mixedsort()))

##### make plot #####
# g_early_ptp99a_normexp_tsne <- ggplot(df_integrated.ptp99A %>%
#                                         filter(dataset == "early") %>%
#                                         mutate(ptp99a = log(ptp99a)), 
#                                       aes(x = tSNE_1, y = tSNE_2, col = ptp99a)) +
#   geom_point(size = .5, alpha = .5, shape = 16) +
#   scale_color_gradientn(colours = pals::ocean.amp(100), guide = "colourbar") +
#   theme_bw() +
#   theme(legend.position = "none")
# ggsave("../figures/Figure3a1.png", g_early_ptp99a_normexp_tsne, w = 2.5, h = 3)

g_main_ptp99a_normexp_tsne <- ggplot(df_integrated.ptp99A %>%
                                       filter(dataset == "main") %>%
                                       mutate(ptp99a = log(ptp99a)), 
                                     aes(x = tSNE_1, y = tSNE_2, col = ptp99a)) +
  geom_point(size = .5, alpha = .5, shape = 16) +
  scale_color_gradientn(colours = pals::ocean.amp(100), guide = "colourbar") +
  theme_bw() +
  theme(legend.position = "none")
ggsave("../figures/Figure3a.png", g_main_ptp99a_normexp_tsne, w = 2.5, h = 3)


#### Figure 3b ####
##### load dataset #####
df_integrated.ptp99A_L <-
  df_integrated.ptp99A %>%
  filter(type %in% c("L1", "L2", "L3", "L4", "L5")) %>%
  mutate(ptp99a = log(ptp99a), #exprssion data analyzed by GLMM was already expr+1
         genotype2 = factor(genotype2, levels = c("T/T", "C/T"))) %>%
  group_by(dataset, type, time, genotype2) %>%
  summarize(mean = mean(ptp99a, na.rm = TRUE),
            sd = sd(ptp99a, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            pval = mean(pval, na.rm = TRUE)) %>%
  mutate(star = makeStars(pval))

##### make plot #####
g_main_normexp_ptp99a_L_meanse <- ggplot(df_integrated.ptp99A_L %>%
                                           filter(dataset == "main") %>%
                                           mutate(time = parse_number(time) %>% as.factor()),
                                         aes(x = time, y = mean, group = genotype2)) +
  geom_path(aes(color = genotype2), position = position_dodge(0.3)) +
  geom_point(aes(color = genotype2), shape = 16, size = 2, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd, color = genotype2), width = 0, position = position_dodge(0.3)) + 
  geom_text(aes(x = time,  y = 4.4, label = star), check_overlap = TRUE) +
  geom_text(x = 0.8, y = 4.6, label = "Error bar: ±SD", size = 2, hjust = 0,
            data = df_main.ptp99a %>%
              filter(type %in% c("L1")) %>%
              group_by(type, time, genotype2) %>%
              dplyr::summarize(mean = mean(ptp99a, na.rm = TRUE)),
            check_overlap = TRUE) +
  xlab("Time after pupal formation (h)") +
  ylab("Log(expr+1)") +
  scale_color_manual(values = c("#B8A3BA", "#401336")) +
  scale_alpha_manual(values = c(0, 0.01)) +
  coord_cartesian(ylim = c(0, 5)) +
  facet_grid(type ~ ., scales = "free_x", space = "free") +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = c(0.5, 0.98),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        strip.background = element_blank())#,
g_main_normexp_ptp99a_L_meanse

ggsave("../figures/Figure3b.pdf", g_main_normexp_ptp99a_L_meanse, w = 2.5, h = 4, limitsize = FALSE)


#### Figure 3c ####
##### load dataset #####
dfm_s5min_2995_5990_freezing_duration <- read_parquet("../data/3_mutants/dfm_s5min_2995_5990_freezing_duration.parquet") %>%
  ungroup()

dfm_motion_cue_exit_coeff_trial2 <- read_parquet("../data/3_mutants/dfm_motion_cue_exit_coeff_trial2.parquet") %>%
  ungroup()




##### make plot #####
###### freezing duration Ptp99ARNAi_R9B08 ######
my_comparisons_Ptp99ARNAi_R9B08 <- list( c("NA_R9B08", "Ptp99ARNAi_R9B08"), c("Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA") )
g_freezing_Ptp99ARNAi_R9B08 <- 
  ggplot(dfm_s5min_2995_5990_freezing_duration %>%
           filter(strain %in% c("NA_R9B08", "Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA"),
                  sex == "Female") %>%
           transform(strain = factor(strain, levels = c("NA_R9B08", "Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA"))),
         aes(x = strain, y = freezing_duration, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_Ptp99ARNAi_R9B08) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("R9B08-GAL4"), italic("R9B08>Ptp99A-RNAi"), italic("Ptp99A-RNAi"))) +
  scale_y_continuous(breaks = seq(0, 15, 5), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 18)) +
  ylab("Freezing duration (s)") +
  facet_grid(n_inds ~ .) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

g_freezing_Ptp99ARNAi_R9B08

###### freezing duration TNT_Ptp99A ######
my_comparisons_TNT_Ptp99A <- list( c("IMPTNT_Ptp99A", "TNT_Ptp99A") )
g_freezing_TNT_Ptp99A <- 
  ggplot(dfm_s5min_2995_5990_freezing_duration %>%
           filter(strain %in% c("IMPTNT_Ptp99A", "TNT_Ptp99A"),
                  sex == "Female") %>% #, n_inds == "Single"
           transform(strain = factor(strain, levels = c("IMPTNT_Ptp99A", "TNT_Ptp99A"))),
         aes(x = strain, y = freezing_duration, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_TNT_Ptp99A) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("Ptp99A>IMPTNT"), italic("Ptp99A>TNT"))) +
  scale_y_continuous(breaks = seq(0, 15, 5), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 16)) +
  ylab("Freezing duration (s)") +
  facet_grid(n_inds ~ .) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

g_freezing_TNT_Ptp99A

###### freezing duration TNT_R9B08 ######
my_comparisons_TNT_R9B08 <- list( c("IMPTNT_R9B08", "TNT_R9B08") )
g_freezing_TNT_R9B08 <- ggplot(dfm_s5min_2995_5990_freezing_duration %>%
                                  filter(strain %in% c("IMPTNT_R9B08", "TNT_R9B08"),
                                         sex == "Female") %>% #, n_inds == "Single"
                                  transform(strain = factor(strain, levels = c("IMPTNT_R9B08", "TNT_R9B08"))),
                               aes(x = strain, y = freezing_duration, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_TNT_R9B08) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("R9B08>IMPTNT"), italic("R9B08>TNT"))) +
  scale_y_continuous(breaks = seq(0, 15, 5), expand = c(0, 0)) +
  coord_cartesian(ylim = c(0, 16)) +
  ylab("Freezing duration (s)") +
  facet_grid(n_inds ~ .) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none",
        strip.background = element_blank())

g_freezing_TNT_R9B08


###### visual responsiveness Ptp99ARNAi_R9B08 ######
dfm_plot <- 
  dfm_motion_cue_exit_coeff_trial2 %>%
  pivot_longer(cols = contains("motion_cue"), names_to = "var", values_to = "value")
my_comparisons_Ptp99ARNAi_R9B08 <- list( c("NA_R9B08", "Ptp99ARNAi_R9B08"), c("Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA") )
g_visual_responsiveness_Ptp99ARNAi_R9B08 <- 
  ggplot(dfm_plot %>%
           filter(strain %in% c("NA_R9B08", "Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA"),
                  sex == "Female", var == "motion_cue_exit_intercept") %>%
           transform(strain = factor(strain, levels = c("NA_R9B08", "Ptp99ARNAi_R9B08", "Ptp99ARNAi_NA"))),
         aes(x = strain, y = value, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_Ptp99ARNAi_R9B08) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("R9B08-GAL4"), italic("R9B08>Ptp99A-RNAi"), italic("Ptp99A-RNAi"))) +
  coord_cartesian(ylim = c(-4.5, -2.2)) +
  ylab("Visual responsiveness") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")

g_visual_responsiveness_Ptp99ARNAi_R9B08


###### visual responsiveness TNT_Ptp99A ######
dfm_plot <- 
  dfm_motion_cue_exit_coeff_trial2 %>%
  pivot_longer(cols = contains("motion_cue"), names_to = "var", values_to = "value")
my_comparisons_TNT_Ptp99A <- list( c("IMPTNT_Ptp99A", "TNT_Ptp99A") )
g_visual_responsiveness_TNT_Ptp99A <- 
  ggplot(dfm_plot %>%
           filter(strain %in% c("IMPTNT_Ptp99A", "TNT_Ptp99A"),
                  sex == "Female", var == "motion_cue_exit_intercept") %>%
           transform(strain = factor(strain, levels = c("IMPTNT_Ptp99A", "TNT_Ptp99A"))),
         aes(x = strain, y = value, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_TNT_Ptp99A) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("Ptp99A>IMPTNT"), italic("Ptp99A>TNT"))) +
  coord_cartesian(ylim = c(-4.5, -2.2)) +
  ylab("Visual responsiveness") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")

g_visual_responsiveness_TNT_Ptp99A


###### visual responsiveness TNT_R9B08 ######
dfm_plot <- 
  dfm_motion_cue_exit_coeff_trial2 %>%
  pivot_longer(cols = contains("motion_cue"), names_to = "var", values_to = "value")
my_comparisons_TNT_R9B08 <- list( c("IMPTNT_R9B08", "TNT_R9B08") )
g_visual_responsiveness_TNT_R9B08 <- ggplot(dfm_plot %>%
                         filter(strain %in% c("IMPTNT_R9B08", "TNT_R9B08"),
                                sex == "Female", var == "motion_cue_exit_intercept") %>%
                         transform(strain = factor(strain, levels = c("IMPTNT_R9B08", "TNT_R9B08"))),
                       aes(x = strain, y = value, color = strain)) +
  gg.layers::geom_boxplot2(outlier.shape = NA, width.errorbar = 0) +
  geom_jitter(shape = 16, alpha=0.2, size = 2, width = 0.3, height = 0) +
  ggpubr::stat_compare_means(comparisons = my_comparisons_TNT_R9B08) +
  scale_color_manual(values = c("#5a5359", "#deb068", "#5a5359")) +
  scale_x_discrete(labels = expression(italic("R9B08>IMPTNT"), italic("R9B08>TNT"))) +
  coord_cartesian(ylim = c(-4.5, -2.2)) +
  ylab("Visual responsiveness") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")

g_visual_responsiveness_TNT_R9B08


g3c <- 
  g_freezing_TNT_R9B08 + g_freezing_TNT_Ptp99A + g_freezing_Ptp99ARNAi_R9B08 +
  g_visual_responsiveness_TNT_R9B08 + g_visual_responsiveness_TNT_Ptp99A + g_visual_responsiveness_Ptp99ARNAi_R9B08 +
  plot_layout(nrow = 2, ncol = 3, widths = c(0.66, 0.66, 1), heights = c(1.5, 1))
g3c
ggsave("../figures/Figure3c.pdf", g3c, w = 6, h = 6)



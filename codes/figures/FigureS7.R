#### load packages ####
library(tidyverse)


#### Figure S7a-c ####
##### load dataset #####
###### gwas_freezing_duration_single_scaleT_female ######
df_gwas_freezing_duration_single_scaleT_female <- read_parquet("../data/1_single_strain/gwas/df/df_gwas_freezing_duration_single_scaleT_female.parquet") %>%
  ungroup()

axisdf_gwas_freezing_duration_single_scaleT_female <- df_gwas_freezing_duration_single_scaleT_female %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_freezing_duration_single_scaleT_male ######
df_gwas_freezing_duration_single_scaleT_male <- read_parquet("../data/1_single_strain/gwas/df/df_gwas_freezing_duration_single_scaleT_male.parquet") %>%
  ungroup()

axisdf_gwas_freezing_duration_single_scaleT_male <- df_gwas_freezing_duration_single_scaleT_male %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_freezing_duration_group_scaleT_female ######
df_gwas_freezing_duration_group_scaleT_female <- read_parquet("../data/1_single_strain/gwas/df/df_gwas_freezing_duration_group_scaleT_female.parquet") %>%
  ungroup()

axisdf_gwas_freezing_duration_group_scaleT_female <- df_gwas_freezing_duration_group_scaleT_female %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_freezing_duration_group_scaleT_male ######
df_gwas_freezing_duration_group_scaleT_male <- read_parquet("../data/1_single_strain/gwas/df/df_gwas_freezing_duration_group_scaleT_male.parquet") %>%
  ungroup()

axisdf_gwas_freezing_duration_group_scaleT_male <- df_gwas_freezing_duration_group_scaleT_male %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_nnd_scaleT_female ######
df_gwas_nnd_scaleT_female <- read_parquet("../data/1_single_strain/gwas/df/df_gwas_nnd_scaleT_female.parquet") %>%
  ungroup()

axisdf_gwas_nnd_scaleT_female <- df_gwas_nnd_scaleT_female %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_nnd_scaleT_male ######
df_gwas_nnd_scaleT_male <- read_parquet("../data/1_single_strain/gwas/df/df_gwas_nnd_scaleT_male.parquet") %>%
  ungroup()

axisdf_gwas_nnd_scaleT_male <- df_gwas_nnd_scaleT_male %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_motion_cue_exit_intercept_scaleT_female ######
df_gwas_motion_cue_exit_intercept_scaleT_female <- read_parquet("../data/1_single_strain/gwas/df/df_gwas_motion_cue_exit_intercept_scaleT_female.parquet") %>%
  ungroup()

axisdf_gwas_motion_cue_exit_intercept_scaleT_female <- df_gwas_motion_cue_exit_intercept_scaleT_female %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_motion_cue_exit_intercept_scaleT_male ######
df_gwas_motion_cue_exit_intercept_scaleT_male <- read_parquet("../data/1_single_strain/gwas/df/df_gwas_motion_cue_exit_intercept_scaleT_male.parquet") %>%
  ungroup()

axisdf_gwas_motion_cue_exit_intercept_scaleT_male <- df_gwas_motion_cue_exit_intercept_scaleT_male %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

##### make plot #####
###### Figure S7a1 freezing_duration_single ######
g_gwas_freezing_duration_single_scaleT_female <- 
  ggplot(df_gwas_freezing_duration_single_scaleT_female, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c8c2be","#71686c"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_x_continuous(label = axisdf_gwas_freezing_duration_single_scaleT_female$CHR, breaks = axisdf_gwas_freezing_duration_single_scaleT_female$center) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 10)) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_freezing_duration_single_scaleT_female
ggsave("../figures/FigureS7a1.png", 
       g_gwas_freezing_duration_single_scaleT_female, 
       width = 6, height = 2)

g_gwas_freezing_duration_single_scaleT_male <- 
  ggplot(df_gwas_freezing_duration_single_scaleT_male, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c8c2be","#71686c"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_x_continuous(label = axisdf_gwas_freezing_duration_single_scaleT_male$CHR, breaks = axisdf_gwas_freezing_duration_single_scaleT_male$center) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 10)) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_freezing_duration_single_scaleT_male
ggsave("../figures/FigureS7a2.png", 
       g_gwas_freezing_duration_single_scaleT_male, 
       width = 6, height = 2)


###### Figure S7a2 freezing_duration_group ######
g_gwas_freezing_duration_group_scaleT_female <- 
  ggplot(df_gwas_freezing_duration_group_scaleT_female, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) + 
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_x_continuous(label = axisdf_gwas_freezing_duration_group_scaleT_female$CHR, breaks = axisdf_gwas_freezing_duration_group_scaleT_female$center) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 10)) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_freezing_duration_group_scaleT_female
ggsave("../figures/FigureS7a3.png", 
       g_gwas_freezing_duration_group_scaleT_female, 
       width = 6, height = 2)

g_gwas_freezing_duration_group_scaleT_male <- 
  ggplot(df_gwas_freezing_duration_group_scaleT_male, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_x_continuous(label = axisdf_gwas_freezing_duration_group_scaleT_male$CHR, breaks = axisdf_gwas_freezing_duration_group_scaleT_male$center) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 10)) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_freezing_duration_group_scaleT_male
ggsave("../figures/FigureS7a4.png", 
       g_gwas_freezing_duration_group_scaleT_male, 
       width = 6, height = 2)


###### Figure S7b nnd ######
g_gwas_nnd_scaleT_female <- 
  ggplot(df_gwas_nnd_scaleT_female, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_x_continuous(label = axisdf_gwas_nnd_scaleT_female$CHR, breaks = axisdf_gwas_nnd_scaleT_female$center) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 10)) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_nnd_scaleT_female
ggsave("../figures/FigureS7b1.png", 
       g_gwas_nnd_scaleT_female, 
       width = 6, height = 2)

g_gwas_nnd_scaleT_male <- 
  ggplot(df_gwas_nnd_scaleT_male, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) + 
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_x_continuous(label = axisdf_gwas_nnd_scaleT_male$CHR, breaks = axisdf_gwas_nnd_scaleT_male$center) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 10)) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_nnd_scaleT_male
ggsave("../figures/FigureS7b2.png", 
       g_gwas_nnd_scaleT_male, 
       width = 6, height = 2)

###### Figure S7c motion_cue_exit_intercept ######
g_gwas_motion_cue_exit_intercept_scaleT_female <- 
  ggplot(df_gwas_motion_cue_exit_intercept_scaleT_female, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) + 
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_x_continuous(label = axisdf_gwas_motion_cue_exit_intercept_scaleT_female$CHR, breaks = axisdf_gwas_motion_cue_exit_intercept_scaleT_female$center) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 10)) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_motion_cue_exit_intercept_scaleT_female
ggsave("../figures/FigureS7c1.png", 
       g_gwas_motion_cue_exit_intercept_scaleT_female, 
       width = 6, height = 2)

g_gwas_motion_cue_exit_intercept_scaleT_male <- 
  ggplot(df_gwas_motion_cue_exit_intercept_scaleT_male, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) + 
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_x_continuous(label = axisdf_gwas_motion_cue_exit_intercept_scaleT_male$CHR, breaks = axisdf_gwas_motion_cue_exit_intercept_scaleT_male$center) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 10)) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_motion_cue_exit_intercept_scaleT_male
ggsave("../figures/FigureS7c2.png", 
       g_gwas_motion_cue_exit_intercept_scaleT_male, 
       width = 6, height = 2)

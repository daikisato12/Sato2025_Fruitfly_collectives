#### load packages ####
targetPackages <- c('tidyverse','arrow','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### Figure S14a ####
##### load dataset #####
dfd_fig4 <- read_delim("../data/2_mixed_strain/dfd_fig4.tsv")

##### make plot #####
g_d_s5min_stim_speed_normbyf5minave_performance <-
  dfd_fig4 %>%
  # filter(!str_detect(strain, "norpA")) %>%
  dplyr::mutate(thr_sec = paste0("Threshold: ", thr_sec, " s")) %>%
  transform(thr_sec = factor(thr_sec, levels = unique(.$thr_sec) %>% gtools::mixedsort())) %>%
  ggplot(aes(x = alpha, y = performance,
             col = alpha)) +#col = t.test.p.value)) +
  geom_boxplot() +
  geom_path(aes(group = strain), color = "gray", linewidth = 0.4) +
  geom_point(aes(shape = alpha), show.legend = F, size = 2) +
  ggpubr::stat_compare_means(paired = TRUE, size = 1.8) +
  scale_color_manual(values = c("#9e8896", "#874c70")) +#values = viridis(3)[1:2]) +
  scale_shape_manual(values = c("Expected" = 17, "Observed" = 16)) +#values = viridis(3)[1:2]) +
  ylab("Behavioral performance") +
  facet_wrap( ~ thr_sec, nrow = 2) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9),
        legend.position = "bottom")
g_d_s5min_stim_speed_normbyf5minave_performance


#### Figure S14b ####
##### load dataset #####
dfd_fig4_2 <- read_delim("../data/2_mixed_strain/dfd_fig4_2.tsv")

##### make plot #####
dfd_figS14b <- dfd_fig4_2 %>%
  dplyr::mutate(thr_sec = paste0("Threshold: ", thr_sec, " s")) %>%
  transform(thr_sec = factor(thr_sec, levels = unique(.$thr_sec) %>% gtools::mixedsort()),
            time = factor(time, levels = c("Before", "After"))) %>%
  # group_by(thr_sec) %>%
  dplyr::mutate(xpos = as.numeric(factor(time)) + 
           ifelse(alpha == "Expected", -0.25, 0.25))  # Cごとにずらす

stat_d_s5min_stim_speed_normbyf5minave_performance_sep <-
  dfd_figS14b %>%
  group_by(thr_sec, time) %>%
  rstatix::wilcox_test(speed ~ alpha, paired = TRUE) %>%
  rstatix::adjust_pvalue(method = "bonferroni") %>%
  rstatix::add_significance("p.adj") %>%
  dplyr::mutate(xpos = if_else(time == "Before", 1, 2),
         ypos = if_else(p.adj.signif == "ns", 1.5, 1.47)) %>%
  print()


g_d_s5min_stim_speed_normbyf5minave_performance_sep <-
  dfd_figS14b %>%
  # filter(!str_detect(strain, "norpA")) %>%
  dplyr::mutate(group = paste(thr_sec, time, strain)) %>%
  ggplot(aes(x = xpos, y = speed,
             col = alpha)) +#col = t.test.p.value)) +
  # geom_boxplot() +
  geom_path(aes(group = group), color = "gray", linewidth = 0.4) +#aes(group = interaction(thr_sec, time, strain)), position = position_dodge(0.5), 
  geom_point(aes(shape = alpha), show.legend = F, size = 1) + #position = position_dodge(0.5), 
  # geom_jitter(aes(shape = mixed), height=0, width =0.2, size = 2) +
  # ggpubr::stat_compare_means(paired = TRUE, size = 1.8) +
  geom_text(data = stat_d_s5min_stim_speed_normbyf5minave_performance_sep, 
            aes(x = xpos, y = ypos, label = p.adj.signif), size = 3, inherit.aes = FALSE) +
  # annotate(geom = "text",
  #          x = c(1, 2),
  #          y = 1.2,
  #          label = stat_d_s5min_stim_speed_normbyf5minave_performance_sep %>%
  #            filter(thr_sec == thr_sec, ) %>%
  #            pull(p.adj.signif)) +
  coord_cartesian(xlim = c(0.5, 2.5), ylim = c(0.5, 1.6)) + #norm speed
  # coord_cartesian(xlim = c(0.5, 2.5), ylim = c(0, 10)) + #raw speed
  scale_x_continuous(breaks = c(1, 2), labels = c("Before", "After")) +
  scale_color_manual(values = c("#9e8896", "#874c70")) +#values = viridis(3)[1:2]) +
  scale_shape_manual(values = c("Expected" = 17, "Observed" = 16)) +#values = viridis(3)[1:2]) +
  ylab("Mean relative moving speed") +
  facet_wrap( ~ thr_sec, nrow = 2) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9),
        panel.grid.minor.x = element_blank())#,
g_d_s5min_stim_speed_normbyf5minave_performance_sep

g_figS14 <- g_d_s5min_stim_speed_normbyf5minave_performance +
     g_d_s5min_stim_speed_normbyf5minave_performance_sep

ggsave("../figures/FigureS14.pdf", g_figS14, w=10, h=5)


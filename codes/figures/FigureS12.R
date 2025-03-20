#### load packages ####
targetPackages <- c('tidyverse','arrow')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### Figure S12 ####
##### load dataset #####
dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain <- 
  read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain.parquet") %>%
  ungroup()

##### make plot #####
g_d_s5min_stim_speed_gd_onlymix_normbyf5minave_strain_boxplot_gd_alltime_facet <- 
  ggplot(dfd_s5min_2995_5990_speed_sgd_normbyf5minave_strain %>%
           filter(n_inds == "Group",
                  type == "Mixed",
                  !str_detect(strain, "norpA")) %>%
           mutate(stim_time = paste0(stim_time, " s")) %>%
           transform(stim_time = factor(stim_time, levels = unique(.$stim_time) %>% gtools::mixedsort())),
         aes(x = alpha, y = speed_norm_mean,
             color = alpha)) +
  geom_boxplot() + 
  geom_path(aes(group = strain), color = "gray", linewidth = 0.4) +
  geom_point(aes(shape = mixed), show.legend = F, size = 2) + 
  scale_color_manual(values = c("#9e8896", "#874c70")) +#values = viridis(3)[1:2]) +
  scale_shape_manual(values = c(17, 16)) +
  ggpubr::stat_compare_means(paired = TRUE, size = 2.5) +
  coord_cartesian(ylim = c(0.5, 1.3)) +
  ylab("Mean relative moving speed") +
  facet_wrap( ~ stim_time, nrow = 5) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        strip.background = element_blank(),
        legend.title = element_blank(),
        legend.position = "bottom")
g_d_s5min_stim_speed_gd_onlymix_normbyf5minave_strain_boxplot_gd_alltime_facet
ggsave("../figures/FigureS12.pdf", g_d_s5min_stim_speed_gd_onlymix_normbyf5minave_strain_boxplot_gd_alltime_facet, w=6, h=8)


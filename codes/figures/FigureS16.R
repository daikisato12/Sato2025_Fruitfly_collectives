#### load packages ####
targetPackages <- c('tidyverse','gtools','arrow')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### load functions ####
t_test_1sample_foraging_success <- function(dat){
  pval <-
    t.test(dat %>%
             filter(alpha == "Observed") %>%
             pull(foraging_success),
           mu = dat %>%
             filter(alpha == "Expected") %>%
             pull(foraging_success) %>%
             mean()
    )$p.value
  return(pval)
}

#### Figure S16 ####
##### load dataset #####
dfd_figS16 <- read_delim("../data/2_mixed_strain/dfd_figS16.tsv")

##### make plot #####
g_d_s5min_stim_speed_normbyf5minave_performance_integrated2 <-
  dfd_figS16 %>%
  # filter(!str_detect(strain, "norpA")) %>%
  mutate(thr_sec = paste0(thr_sec, " s"),
         col2 = if_else(col2 == "sig", "P ≦ 0.05", "P > 0.05")) %>%
  transform(thr_sec = factor(thr_sec, levels = unique(.$thr_sec) %>% gtools::mixedsort())) %>%
  filter(var == "Group_Mixed") %>%
  ggplot(aes(x = tidytext::reorder_within(strain, overyield, thr_sec, na.rm = TRUE),
             y = overyield, col = col, alpha = col2)) +#col = t.test.p.value)) +
  geom_hline(yintercept = 0, linetype = "dotted", linewidth = 0.4) +
  geom_point(size = 2, shape = 16) +
  scale_alpha_manual(values = c(0.4, 1)) +
  # scale_color_viridis_d(direction = -1) +
  scale_color_manual(values = c("Transgressive overyielding" = "#642426", 
                                "Overyielding"= "#CF7B2E", 
                                "Underyielding" = "#5A7C69", 
                                "Transgressive underyielding" = "#215A83")) +
  # scale_color_gradientn(colours = rev(tol.rainbow(100))) +
  xlab("Combination of strains") +
  ylab("Diversity effect on behavioral performance") +
  facet_wrap( ~ thr_sec, ncol = 1, scales = "free_x", strip.position="right") +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 9),
        legend.title = element_blank(),
        legend.text =  element_text(size = 6),
        legend.key = element_blank(),
        legend.key.size = unit(0.1, "cm"),
        # legend.position = c(0.88, 0.36),
        legend.position = "top",
        legend.background = element_blank(),
        legend.box = "vertical",
        legend.margin = margin(0, 0, unit(0.1, "cm"), 0), # レジェンドボックス間の距離
        legend.spacing.y =  unit(0.1, "cm") , # レジェンドボックス間の距離
  )
g_d_s5min_stim_speed_normbyf5minave_performance_integrated2
ggsave("../figures/FigureS16.pdf", g_d_s5min_stim_speed_normbyf5minave_performance_integrated2, w=6, h=4)


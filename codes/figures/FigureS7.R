
#### load packages ####
targetPackages <- c('tidyverse','arrow','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### load function ####
nls_get_ab_speed_normbystimminus05 <- function(dat){
  # dat <- test[[4]][[1]]
  res <- nls(formula = speed_normbystimminus05 ~ a * log(stim_time) + b, 
             start = c(a = 1, b = 1), trace = TRUE,
             data = dat) %>%
    summary()
  
  dat <- data.frame(dat,
                    a = res$coefficients[1,1],
                    b = res$coefficients[2,1])
  return(dat)
}

#### load dataset ####
dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial <- read_parquet("../data/2_mixed_strain/parquet/dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial.parquet") %>%
  ungroup()

#### nls for freezing curve ####
dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial_nls <-
  dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial %>%
  filter(stim_time >= 0.5, n_inds == "Group", type != "Mixed") %>%
  dplyr::mutate(strain = if_else(str_detect(var, "Strain1"), 
                                 str_split(strain, "_") %>% map_chr(1),
                                 str_split(strain, "_") %>% map_chr(2))) %>%
  distinct() %>%
  group_by(strain, stim_time) %>%
  dplyr::summarize(speed_normbystimminus05 = mean(speed_normbystimminus05, na.rm = TRUE)) %>%
  ungroup() %>%
  group_nest(strain) %>%
  mutate(data = map(data, nls_get_ab_speed_normbystimminus05)) %>%
  unnest() %>%
  group_by(strain) %>%
  dplyr::summarize(a = mean(a, na.rm = TRUE),
                   b = mean(b, na.rm = TRUE))

#### Figure S7a ####
df_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial <-
  dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial %>%
  filter(stim_time >= 0.5, n_inds == "Group", mixed == "Single strain", alpha == "Observed") %>%
  separate(strain, into = c("strain1", "strain2")) %>%
  dplyr::mutate(strain = if_else(type == "Strain1", strain1, strain2)) %>%
  group_by(strain, stim_time) %>%
  dplyr::summarize(speed_normbystimminus05 = mean(speed_normbystimminus05, na.rm = TRUE)) %>%
  ungroup() %>%
  transform(strain = factor(strain, levels = gtools::mixedsort(unique(.$strain))))
g_d_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial_nls <-
  ggplot(df_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial %>%
           filter(strain != "norpA"),
         aes(x = stim_time, y = speed_normbystimminus05)) +
  geom_point(shape = 16, alpha = 0.8) +
  stat_smooth(method = "nls", se = FALSE,
              formula = y ~ a * log(x) + b,#y ~ log(x),
              method.args = list(start = c(a = 1, b = 1)),
              linewidth = 1, color = "grey") +#, formula = y ~ poly(x, degree = 2, raw = TRUE) - 1) + #formula = y ~ log(x)) + #method = "lm") +
  scale_x_continuous(breaks = seq(0, 15, 3)) +
  coord_cartesian(xlim = c(0, 15)) +
  facet_wrap(~ strain, ncol = 5) +
  xlab("Time after stimulus (s)") +
  ylab("Moving speed relative to \n the average at 0.5 s before stimulus") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"),
        strip.background = element_blank())
g_d_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial_nls
ggsave(paste0("../figures/FigureS7a.pdf"), g_d_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial_nls, w=6, h=4)


#### Figure S7b ####
g_dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial_nls <-
  ggplot(dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial_nls %>%
           filter(!str_detect(strain, "norpA")) %>%
           dplyr::select(a, b) %>%
           distinct(),
         aes(x = a, y = b)) +
  geom_point(shape = 16, alpha = 0.8, size = 3) +
  stat_smooth(method = "lm", color = "black") +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(
                          stat(eq.label),
                          # after_stat(rr.label),
                          # stat(p.value.label),
                          sep = "~~~")),
                        label.x = "right",
                        label.y = "top",
                        parse = TRUE, size=4) +
  scale_x_continuous(breaks = seq(0, 1, 0.02)) +
  coord_cartesian(xlim = c(0, 0.098)) +
  xlab("Parameter a") +
  ylab("Parameter b") +
  theme_bw() +
  theme(axis.text = element_text(color = "black"))

g_dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial_nls
ggsave(paste0("../figures/FigureS7b.pdf"), 
       g_dfd_s5min_2995_5990_speed_sgd_normbystimminus05_strain_trial_nls, w=3, h=3)


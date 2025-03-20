#### load packages ####
targetPackages <- c('tidyverse','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)
# devtools::install_github("norment/normentR")

#### Figure S10abc ####
##### load dataset #####
###### gwas_freezing_duration_single_scaleT_female ######
df_gwas_freezing_duration_single_scaleT_female <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_single_scaleT_female.parquet") %>%
  ungroup()

axisdf_gwas_freezing_duration_single_scaleT_female <- df_gwas_freezing_duration_single_scaleT_female %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_freezing_duration_single_scaleT_male ######
df_gwas_freezing_duration_single_scaleT_male <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_single_scaleT_male.parquet") %>%
  ungroup()

axisdf_gwas_freezing_duration_single_scaleT_male <- df_gwas_freezing_duration_single_scaleT_male %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_freezing_duration_group_scaleT_female ######
df_gwas_freezing_duration_group_scaleT_female <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_group_scaleT_female.parquet") %>%
  ungroup()

axisdf_gwas_freezing_duration_group_scaleT_female <- df_gwas_freezing_duration_group_scaleT_female %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_freezing_duration_group_scaleT_male ######
df_gwas_freezing_duration_group_scaleT_male <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_group_scaleT_male.parquet") %>%
  ungroup()

axisdf_gwas_freezing_duration_group_scaleT_male <- df_gwas_freezing_duration_group_scaleT_male %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_nnd_scaleT_female ######
df_gwas_nnd_scaleT_female <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_nnd_scaleT_female.parquet") %>%
  ungroup()

axisdf_gwas_nnd_scaleT_female <- df_gwas_nnd_scaleT_female %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_nnd_scaleT_male ######
df_gwas_nnd_scaleT_male <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_nnd_scaleT_male.parquet") %>%
  ungroup()

axisdf_gwas_nnd_scaleT_male <- df_gwas_nnd_scaleT_male %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_motion_cue_exit_intercept_scaleT_female ######
df_gwas_motion_cue_exit_intercept_scaleT_female <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_motion_cue_exit_intercept_scaleT_female.parquet") %>%
  ungroup()

axisdf_gwas_motion_cue_exit_intercept_scaleT_female <- df_gwas_motion_cue_exit_intercept_scaleT_female %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

###### gwas_motion_cue_exit_intercept_scaleT_male ######
df_gwas_motion_cue_exit_intercept_scaleT_male <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_motion_cue_exit_intercept_scaleT_male.parquet") %>%
  ungroup()

axisdf_gwas_motion_cue_exit_intercept_scaleT_male <- df_gwas_motion_cue_exit_intercept_scaleT_male %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

##### make plot #####
###### Figure S10a1 freezing_duration_single ######
g_gwas_freezing_duration_single_scaleT_male <- 
  ggplot(df_gwas_freezing_duration_single_scaleT_male, 
         aes(y = BPcum, 
             x = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c8c2be","#71686c"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(label = axisdf_gwas_freezing_duration_single_scaleT_male$CHR, 
                     breaks = axisdf_gwas_freezing_duration_single_scaleT_male$center,
                     position = "right") +
  # scale_x_continuous(breaks = seq(0, 10, 2),
  #                    labels = seq(0, 10, 2)) +#scales::pretty_breaks(n = 6)) +
  scale_x_reverse(breaks = seq(0, 10, 2),
                  labels = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(10, 0)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  ggtitle("Male") +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
# g_gwas_freezing_duration_single_scaleT_male

nSNPs <- nrow(df_gwas_freezing_duration_single_scaleT_male)
ci <- 0.95
qqplot_df_gwas_freezing_duration_single_scaleT_male <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_freezing_duration_single_scaleT_male$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 9, color = "black"))

g_gwas_freezing_duration_single_scaleT_female <- 
  ggplot(df_gwas_freezing_duration_single_scaleT_female, 
         aes(y = BPcum, 
             x = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c8c2be","#71686c"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(label = axisdf_gwas_freezing_duration_single_scaleT_female$CHR, 
                     breaks = axisdf_gwas_freezing_duration_single_scaleT_female$center) +
  scale_x_continuous(breaks = seq(0, 10, 2),
                     labels = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(0, 10)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  ggtitle("Female") +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9, hjust = 0.5),
    axis.text = element_text(size = 9, color = "black"),
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
# g_gwas_freezing_duration_single_scaleT_female

nSNPs <- nrow(df_gwas_freezing_duration_single_scaleT_female)
ci <- 0.95
qqplot_df_gwas_freezing_duration_single_scaleT_female <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_freezing_duration_single_scaleT_female$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 9, color = "black"))

g_gwas_freezing_duration_single_scaleT <-
  g_gwas_freezing_duration_single_scaleT_male +
  g_gwas_freezing_duration_single_scaleT_female +
  plot_layout(ncol = 2, widths = c(0.95, 1)) & 
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 2))

g_gwas_freezing_duration_single_qq <-
  g_gwas_freezing_duration_single_scaleT /
  (qqplot_df_gwas_freezing_duration_single_scaleT_male + qqplot_df_gwas_freezing_duration_single_scaleT_female) +
  plot_layout(heights = c(1, 0.3))


ggsave("../figures/FigureS10a1.png",
       g_gwas_freezing_duration_single_qq,
       width = 4, height = 7)

###### Figure S10a2 freezing_duration_group ######
g_gwas_freezing_duration_group_scaleT_male <- 
  ggplot(df_gwas_freezing_duration_group_scaleT_male, 
         aes(y = BPcum, 
             x = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(label = axisdf_gwas_freezing_duration_group_scaleT_male$CHR, 
                     breaks = axisdf_gwas_freezing_duration_group_scaleT_male$center,
                     position = "right") +
  # scale_x_continuous(breaks = seq(0, 10, 2),
  #                    labels = seq(0, 10, 2)) +#scales::pretty_breaks(n = 6)) +
  scale_x_reverse(breaks = seq(0, 10, 2),
                  labels = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(10, 0)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  ggtitle("Male") +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
# g_gwas_freezing_duration_group_scaleT_male

nSNPs <- nrow(df_gwas_freezing_duration_group_scaleT_male)
ci <- 0.95
qqplot_df_gwas_freezing_duration_group_scaleT_male <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_freezing_duration_group_scaleT_male$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 9, color = "black"))

g_gwas_freezing_duration_group_scaleT_female <- 
  ggplot(df_gwas_freezing_duration_group_scaleT_female, 
         aes(y = BPcum, 
             x = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(label = axisdf_gwas_freezing_duration_group_scaleT_female$CHR, 
                     breaks = axisdf_gwas_freezing_duration_group_scaleT_female$center) +
  scale_x_continuous(breaks = seq(0, 10, 2),
                     labels = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(0, 10)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  ggtitle("Female") +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9, hjust = 0.5),
    axis.text = element_text(size = 9, color = "black"),
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
# g_gwas_freezing_duration_group_scaleT_female

nSNPs <- nrow(df_gwas_freezing_duration_group_scaleT_female)
ci <- 0.95
qqplot_df_gwas_freezing_duration_group_scaleT_female <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_freezing_duration_group_scaleT_female$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 9, color = "black"))

g_gwas_freezing_duration_group_scaleT <-
  g_gwas_freezing_duration_group_scaleT_male +
  g_gwas_freezing_duration_group_scaleT_female +
  plot_layout(ncol = 2, widths = c(0.95, 1)) & 
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 2))

g_gwas_freezing_duration_group_qq <-
  g_gwas_freezing_duration_group_scaleT /
  (qqplot_df_gwas_freezing_duration_group_scaleT_male + qqplot_df_gwas_freezing_duration_group_scaleT_female) +
  plot_layout(heights = c(1, 0.3))


ggsave("../figures/FigureS10a2.png",
       g_gwas_freezing_duration_group_qq,
       width = 4, height = 7)


###### Figure S10b nnd ######
g_gwas_nnd_scaleT_male <- 
  ggplot(df_gwas_nnd_scaleT_male, 
         aes(y = BPcum, 
             x = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(label = axisdf_gwas_nnd_scaleT_male$CHR, 
                     breaks = axisdf_gwas_nnd_scaleT_male$center,
                     position = "right") +
  # scale_x_continuous(breaks = seq(0, 10, 2),
  #                    labels = seq(0, 10, 2)) +#scales::pretty_breaks(n = 6)) +
  scale_x_reverse(breaks = seq(0, 10, 2),
                  labels = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(10, 0)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  ggtitle("Male") +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
# g_gwas_nnd_scaleT_male

nSNPs <- nrow(df_gwas_nnd_scaleT_male)
ci <- 0.95
qqplot_df_gwas_nnd_scaleT_male <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_nnd_scaleT_male$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 9, color = "black"))

g_gwas_nnd_scaleT_female <- 
  ggplot(df_gwas_nnd_scaleT_female, 
         aes(y = BPcum, 
             x = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(label = axisdf_gwas_nnd_scaleT_female$CHR, 
                     breaks = axisdf_gwas_nnd_scaleT_female$center) +
  scale_x_continuous(breaks = seq(0, 10, 2),
                     labels = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(0, 10)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  ggtitle("Female") +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9, hjust = 0.5),
    axis.text = element_text(size = 9, color = "black"),
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
# g_gwas_nnd_scaleT_female

nSNPs <- nrow(df_gwas_nnd_scaleT_female)
ci <- 0.95
qqplot_df_gwas_nnd_scaleT_female <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_nnd_scaleT_female$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 9, color = "black"))

g_gwas_nnd_scaleT <-
  g_gwas_nnd_scaleT_male +
  g_gwas_nnd_scaleT_female +
  plot_layout(ncol = 2, widths = c(0.95, 1)) & 
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 2))

g_gwas_nnd_qq <-
  g_gwas_nnd_scaleT /
  (qqplot_df_gwas_nnd_scaleT_male + qqplot_df_gwas_nnd_scaleT_female) +
  plot_layout(heights = c(1, 0.3))
  

ggsave("../figures/FigureS10b.png",
       g_gwas_nnd_qq,
       width = 4, height = 7)

###### Figure S10c motion_cue_exit_intercept ######
g_gwas_motion_cue_exit_intercept_scaleT_male <- 
  ggplot(df_gwas_motion_cue_exit_intercept_scaleT_male, 
         aes(y = BPcum, 
             x = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(label = axisdf_gwas_motion_cue_exit_intercept_scaleT_male$CHR, 
                     breaks = axisdf_gwas_motion_cue_exit_intercept_scaleT_male$center,
                     position = "right") +
  # scale_x_continuous(breaks = seq(0, 10, 2),
  #                    labels = seq(0, 10, 2)) +#scales::pretty_breaks(n = 6)) +
  scale_x_reverse(breaks = seq(0, 10, 2),
                  labels = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(10, 0)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  ggtitle("Male") +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    # axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 9, color = "black"),
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
# g_gwas_motion_cue_exit_intercept_scaleT_male

nSNPs <- nrow(df_gwas_motion_cue_exit_intercept_scaleT_male)
ci <- 0.95
qqplot_df_gwas_motion_cue_exit_intercept_scaleT_male <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_motion_cue_exit_intercept_scaleT_male$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 9, color = "black"))

g_gwas_motion_cue_exit_intercept_scaleT_female <- 
  ggplot(df_gwas_motion_cue_exit_intercept_scaleT_female, 
         aes(y = BPcum, 
             x = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(label = axisdf_gwas_motion_cue_exit_intercept_scaleT_female$CHR, 
                     breaks = axisdf_gwas_motion_cue_exit_intercept_scaleT_female$center) +
  scale_x_continuous(breaks = seq(0, 10, 2),
                     labels = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(0, 10)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  ggtitle("Female") +
  theme_bw() +
  theme(
    legend.position="none",
    plot.title = element_text(hjust = 0.5),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 9, hjust = 0.5),
    axis.text = element_text(size = 9, color = "black"),
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )
# g_gwas_motion_cue_exit_intercept_scaleT_female

nSNPs <- nrow(df_gwas_motion_cue_exit_intercept_scaleT_female)
ci <- 0.95
qqplot_df_gwas_motion_cue_exit_intercept_scaleT_female <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_motion_cue_exit_intercept_scaleT_female$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme(axis.text = element_text(size = 9, color = "black"))

g_gwas_motion_cue_exit_intercept_scaleT <-
  g_gwas_motion_cue_exit_intercept_scaleT_male +
  g_gwas_motion_cue_exit_intercept_scaleT_female +
  plot_layout(ncol = 2, widths = c(0.95, 1)) & 
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 2))

g_gwas_motion_cue_exit_intercept_qq <-
  g_gwas_motion_cue_exit_intercept_scaleT /
  (qqplot_df_gwas_motion_cue_exit_intercept_scaleT_male + qqplot_df_gwas_motion_cue_exit_intercept_scaleT_female) +
  plot_layout(heights = c(1, 0.3))


ggsave("../figures/FigureS10c.png",
       g_gwas_motion_cue_exit_intercept_qq,
       width = 4, height = 7)


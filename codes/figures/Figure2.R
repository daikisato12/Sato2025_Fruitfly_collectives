#### load packages ####
library(tidyverse)
library(tidytext)
library(arrow)

#### Figure 2a ####
##### load dataset #####
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
g_gwas_motion_cue_exit_intercept_scaleT_female <- 
  ggplot(df_gwas_motion_cue_exit_intercept_scaleT_female, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(is_sig),
             size = as.factor(is_sig))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) + #c("#c8c2be","#71686c") for single
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
ggsave("../figures/Figure2a1.png", 
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
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) + #c("#c8c2be","#71686c") for single
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
ggsave("../figures/Figure2a2.png", 
       g_gwas_motion_cue_exit_intercept_scaleT_male, 
       width = 6, height = 2)


#### Figure 2b ####
##### load dataset #####
df_go <- list.files(paste0("../data/1_single_strain/gwas/enrichment/"), 
                     pattern = "*.tsv", full.names = TRUE, recursive=T) %>% 
  map_df(~ data.table::fread(.)) %>%
  mutate(Description = str_to_sentence(Description),
         GeneRatio = str_split(GeneRatio, "/") %>% 
           map_chr(1) %>% 
           as.numeric() / str_split(GeneRatio, "/") %>% 
           map_chr(2) %>% 
           as.numeric()) %>%
  group_by(Type) %>%
  arrange(qvalue) %>%
  slice_head(n = 9) %>%
  mutate(type = "Female") %>%
  ungroup()


##### make plot #####
g_motion_cue_female <- 
  ggplot(df_go, 
         aes(y = tidytext::reorder_within(Description, -qvalue, type), 
             x = GeneRatio, col = qvalue)) +
  geom_point(aes(shape = Type, size = Count)) +#, size= 4) +
  tidytext::scale_y_reordered() +
  scale_color_viridis_c(name = expression(paste(italic(q), "-value"))) +
  scale_shape_manual(values = c(17, 16, 15), labels = c("Biological process", "Cellular component", "Molecular Function")) +
  xlab("Gene ratio") +
  ylab("GO term") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank())
g_motion_cue_female

ggplot2::ggsave("../figures/Figure2b.pdf", 
                g_motion_cue_female, width = 5.6, height = 4.5)


#### Figure 2d ####
##### load dataset #####
df_pheno_motion_cue <- read.table("../data/1_single_strain/gwas/pheno/df_out_motion_cue_exit_intercept_female_GREML.txt") %>%
  dplyr::select(!V1) %>%
  magrittr::set_colnames(c("strain", "value")) %>%
  mutate(sex = "Female") %>%
  bind_rows(read.table("../data/1_single_strain/gwas/pheno/df_out_motion_cue_exit_intercept_male_GREML.txt") %>%
              dplyr::select(!V1) %>%
              magrittr::set_colnames(c("strain", "value")) %>%
              mutate(sex = "Male"))

df_geno_ptp99a <- read.table("../data/1_single_strain/gwas/tophits/candidates/genotype/3R_25292746_SNP.txt", header = TRUE)

df_motion_cue_ptp99a <- inner_join(df_pheno_motion_cue,
                                   df_geno_ptp99a)

##### make plot #####
g_motion_cue_ptp99a <- ggplot(
  df_motion_cue_ptp99a %>%
    transform(genotype = factor(genotype, levels = c("T/T", "C/C"))), 
  aes(x = genotype,
      y = value,
      color = genotype)) +
  geom_boxplot(alpha = 0.8, outlier.colour = NA) +
  geom_jitter(shape = 16, alpha=0.4, size = 2, width = 0.3, height = 0) +
  scale_x_discrete(na.translate = FALSE) +
  scale_color_manual(values = c("#c0abbf", "#491c43")) +
  xlab("Genotype") +
  labs(y = expression("Visual responsiveness to motion cue" ~ (beta[0]))) +  
  facet_wrap(~ sex, ncol = 2) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        strip.background = element_blank())
g_motion_cue_ptp99a

ggsave("../figures/Figure2d.pdf", g_motion_cue_ptp99a, w = 3, h = 2)


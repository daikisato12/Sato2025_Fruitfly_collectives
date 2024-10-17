#### load packages ####
library(tidyverse)
library(plyr)
library(ggcorrplot)
library(corrr)

#### Figure S6a ####
##### load dataset #####
df_heritability <- read.table("../data/1_single_strain/gwas/heritability/1pheno_hG.txt",h=T) %>%
  mutate(trait = case_when(str_detect(trait, "freezing_duration") ~ str_replace(trait, "freezing_duration", "Freezing duration"),
                           str_detect(trait, "visual_reactivity") ~ str_replace(trait, "visual_reactivity", "Visual reactivity_group"),
                           str_detect(trait, "nnd") ~ str_replace(trait, "nnd", "NND_group"),
                           TRUE ~ trait),
         trait_tmp = trait) %>%
  separate(trait, sep = "_", into = c("trait", "n_inds", "sex")) %>% 
  tidyr::replace_na(., replace = list(hg = NA, hg_se = NA)) %>%
  mutate(n_inds = str_to_title(n_inds),
         sex = str_to_title(sex))

##### make plot #####
g_her <- ggplot(df_heritability, 
                aes(x = trait, y = hg, col = n_inds, shape = sex)) +
  geom_point(size = 4, position = position_dodge(width=0.45) ) +
  geom_errorbar(aes(ymin = hg - hg_se, ymax = hg + hg_se, width = 0), position = position_dodge(width=0.45) ) +
  scale_color_manual(values = c("#b3343a", "#aaaaa9")) +
  xlab("Trait") +
  ylab("SNP-based heritability") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(), 
        legend.title = element_blank(),
        legend.position = c(0.75, 0.85),
        legend.margin = unit(0, 'inch'),
        legend.key = element_blank(),
        legend.key.size = unit(0.2, 'inch'),
        legend.background = element_blank(),
        panel.background = element_rect(fill = "transparent",color = NA),
        plot.background = element_rect(fill = "transparent",color = NA))

g_her
ggsave("../figures/FigureS6a.pdf", g_her, width=2.2, height=4.5)


#### Figure S6b ####
##### load dataset #####
df_geno <- read.table("../data/1_single_strain/gwas/heritability/2pheno_rG.txt", h=T) %>%
  mutate(trait1 = case_when(str_detect(trait1, "freezing_duration") ~ str_replace(trait1, "freezing_duration", "Freezing duration"),
                            str_detect(trait1, "visual_reactivity") ~ str_replace(trait1, "visual_reactivity", "Visual reactivity_group"),
                            str_detect(trait1, "nnd") ~ str_replace(trait1, "nnd", "NND_group"),
                            TRUE ~ trait1),
         trait1_tmp = trait1) %>%
  mutate(trait2 = case_when(str_detect(trait2, "freezing_duration") ~ str_replace(trait2, "freezing_duration", "Freezing duration"),
                            str_detect(trait2, "visual_reactivity") ~ str_replace(trait2, "visual_reactivity", "Visual reactivity_group"),
                            str_detect(trait2, "nnd") ~ str_replace(trait2, "nnd", "NND_group"),
                            TRUE ~ trait2),
         trait2_tmp = trait2) %>%
  separate(trait1, sep = "_", into = c("trait1", "n_inds1", "sex1")) %>%
  separate(trait2, sep = "_", into = c("trait2", "n_inds2", "sex2")) %>%
  mutate(n_inds1 = str_to_title(n_inds1),
         sex1 = str_to_title(sex1),
         trait1 = paste0(trait1, " (", n_inds1, "/", sex1, ")"),
         n_inds2 = str_to_title(n_inds2),
         sex2 = str_to_title(sex2),
         trait2 = paste0(trait2, " (", n_inds2, "/", sex2, ")")) %>%
  mutate(compare = case_when(trait1 == trait2 & n_inds1 != n_inds2 & sex1 == sex2 ~ "Single-Group",
                             trait1 == trait2 & n_inds1 == n_inds2 & sex1 != sex2 ~ "Male-Female",
                             TRUE ~ "others")) %>% 
  tidyr::replace_na(., replace = list(cor = 0)) %>%
  mutate(cor = case_when(cor > 1 ~ 1,
                         cor < -1 ~ -1, 
                         TRUE ~ cor))


# cor all
df_genetic_cor <- df_geno %>%
  dplyr::select(trait1, trait2, cor) %>%
  pivot_wider(names_from = trait2, values_from = cor)

mat_genetic_cor <- df_genetic_cor %>%
  dplyr::select(-trait1) %>%
  as.matrix()
rownames(mat_genetic_cor) <- df_genetic_cor$trait1

##### make plot #####
g_geno_cor <- ggcorrplot(mat_genetic_cor,
                         outline.col = "white",
                         ggtheme = ggplot2::theme_bw,
                         colors = c("#522f60", "white", "#a25768"))
g_geno_cor
ggsave("../figures/FigureS6b.pdf", g_geno_cor, width=6, height=6)


#### Figure S6c ####
##### load dataset #####
df_pheno <- read_parquet("../data/1_single_strain/df_f5min_nnd_rand.parquet") %>%
  ungroup() %>%
  filter(!strain %in% c("norpA","random")) %>%
  mutate(strain = paste0("line_",parse_number(as.character(strain))),
         row = strain) %>%
  group_by(strain, sex) %>%
  dplyr::summarize(nnd = mean(nnd, na.rm=T)) %>%
  dplyr::select(strain, sex, nnd) %>%
  mutate(n_inds = "Group") %>%
  pivot_longer(cols = nnd, names_to = "trait", values_to = "value") %>%
  bind_rows(read_parquet("../data/1_single_strain/df_motion_cue_exit_coeff.parquet") %>%
              ungroup() %>%
              filter(!strain %in% c("norpA","random")) %>%
              mutate(strain = paste0("line_",parse_number(as.character(strain))),
                     row = strain) %>%
              group_by(strain, sex) %>%
              dplyr::summarize(motion_cue_exit_intercept = mean(motion_cue_exit_intercept, na.rm=T)) %>%
              dplyr::select(strain, sex, motion_cue_exit_intercept) %>%
              mutate(n_inds = "Group") %>%
              pivot_longer(cols = motion_cue_exit_intercept, names_to = "trait", values_to = "value")) %>%
  
  bind_rows(read_parquet("../data/1_single_strain/df_s5min_2995_5990_freezing_duration.parquet") %>%
              ungroup() %>%
              filter(!strain %in% c("norpA","random")) %>%
              mutate(strain = paste0("line_",parse_number(as.character(strain))),
                     row = strain) %>%
              group_by(strain, sex, n_inds) %>%
              dplyr::summarize(freezing_duration = mean(freezing_duration, na.rm=T)) %>%
              pivot_longer(cols = freezing_duration, names_to = "trait", values_to = "value")) %>%
  mutate(n_inds = tolower(n_inds),
         sex = tolower(sex),
         trait = case_when(trait == "nnd" ~ "NND",
                           trait == "freezing_duration" ~ "Freezing duration",
                           trait == "motion_cue_exit_intercept" ~ "Visual reactivity")) %>%
  arrange(trait, n_inds, sex) %>%
  mutate(n_inds = str_to_title(n_inds),
         sex = str_to_title(sex),
         trait2 = paste0(trait, " (", n_inds, "/", sex, ")")) %>%
  dplyr::select(!c(trait, n_inds, sex)) %>%
  pivot_wider(id_cols = strain, names_from = trait2, values_from = value)


df_phenotypic_cor <- df_pheno %>%
  dplyr::select(!strain) %>% 
  corrr::correlate(method = "pearson", use = "pairwise.complete.obs")

mat_phenotypic_cor <-  df_phenotypic_cor %>%
  dplyr::select(!term) %>%
  as.matrix()
rownames(mat_phenotypic_cor) <- df_phenotypic_cor$term


##### make plot #####
g_pheno_cor <- ggcorrplot(mat_phenotypic_cor, type = "upper",
                          outline.col = "white",
                          ggtheme = ggplot2::theme_bw,
                          colors = c("#522f60", "white", "#a25768"))
g_pheno_cor
ggsave("../figures/FigureS6c.pdf", g_pheno_cor, width=6, height=6)

#### Figure S6# corrleation between pheno vs. geno ####
##### load dataset #####
df_phenotypic_genetic_cor <-
  df_phenotypic_cor %>%
  pivot_longer(cols = !term, values_to = "phenotypic_cor") %>%
  right_join(df_genetic_cor %>%
               pivot_longer(cols = !trait1, 
                            values_to = "genetic_cor") %>%
               dplyr::rename(term = trait1))


##### make plot #####
g_cor_pheno_genet <- 
  ggplot(df_phenotypic_genetic_cor %>%
           mutate(plot_col = if_else(phenotypic_cor - genetic_cor < 0, "over", "under")),
         aes(x = genetic_cor, y = phenotypic_cor)) +
  geom_abline(slope = 1, linetype = "dashed") +
  stat_smooth(linewidth = 2, color= "grey", method = "lm") +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(
                          after_stat(rr.label),
                          stat(p.value.label),
                          sep = "~~~")),
                        label.x = "left",
                        label.y = "top",
                        parse = TRUE, size = 4) +
  geom_point(aes(color = plot_col), size = 4, shape = 16, alpha = 0.5) +
  scale_color_manual(values = c("#a25768","#71686c")) +
  xlab("SNP-based genetic correlation") +
  ylab("Phenotypic correlation") +
  theme_bw() +
  theme(legend.position = "none")
g_cor_pheno_genet

ggsave("../figures/FigureS6#.pdf", g_cor_pheno_genet, width = 2.5, height = 4)

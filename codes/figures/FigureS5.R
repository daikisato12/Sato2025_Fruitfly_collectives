#### load packages ####
library(tidyverse)
library(plyr)

#### Figure S5 ####
##### load dataset #####
strain_1pop <- read.table("../data/1_single_strain/gwas/pheno/df_out_motion_cue_exit_intercept_female_scaled_GREML.txt") %>%
  pull(V1)

strain_2pop <- c("line_88", "line_101", "line_136", "line_161",
                 "line_189", "line_301", "line_309", "line_324", 
                 "line_357", "line_358", "line_360", "line_399", 
                 "line_786", "line_852", "line_855")

df_pca <- read.table("../data/0_genome/dgrp2.eigenvec", header = F) %>%
  dplyr::select(!V2) %>%
  magrittr::set_colnames(c("strain", paste0("PC", seq(1, 10)))) %>%
  dplyr::mutate(color = case_when(strain %in% strain_2pop ~ "Used for GHAS",
                                  strain %in% strain_1pop ~ "Used for GWAS",
                                  TRUE ~ "All strains")) 

pc1_contri <- read.table("../data/0_genome/dgrp2.eigenval", header = F) %>%
  dplyr::pull(V1) %>%
  head(1)

pc2_contri <- read.table("../data/0_genome/dgrp2.eigenval", header = F) %>%
  dplyr::pull(V1) %>%
  head(2) %>%
  tail(1)

df_pca <- df_pca %>%
  bind_rows(df_pca %>%
              filter(color == "Used for GWAS") %>%
              dplyr::mutate(color = "All strains")) %>%
  bind_rows(df_pca %>%
              filter(color == "Used for GHAS") %>%
              dplyr::mutate(color = "All strains")) %>%
  bind_rows(df_pca %>%
              filter(color == "Used for GHAS") %>%
              dplyr::mutate(color = "Used for GWAS")) %>%  
  transform(color = factor(color, levels = c("Used for GWAS", "Used for GHAS", "All strains"))) %>%
  arrange(factor(color, levels = c("All strains", "Used for GWAS", "Used for GHAS")))

chulls <- plyr::ddply(df_pca, .(color), function(df_pca) df_pca[chull(df_pca$PC1, df_pca$PC2), ])

##### make plot #####
g_pca <- ggplot(df_pca,
                aes(x = PC1, y = PC2, color = color)) +
  geom_polygon(data=chulls, aes(x=PC1, y=PC2, fill=color), alpha=0.1) +
  geom_point(aes(alpha = color), shape = 16, size = 4) +
  # ggrepel::geom_label_repel(aes(label = strain)) +
  scale_color_manual(values = c("#b7282e", "#9d5b8b", "#afafb0"),
                     labels = c("Used for GWAS", "Used for GHAS", "All strains")) +
  scale_fill_manual(values = c("#b7282e", "#9d5b8b", "#afafb0"),
                    labels = c("Used for GWAS", "Used for GHAS", "All strains")) +
  scale_alpha_manual(values = c(0.8, 0.8, 0.5)) +
  xlab(paste0("PC1 (", round(pc1_contri, 2), " %)")) +
  ylab(paste0("PC2 (", round(pc2_contri, 2), " %)")) +
  theme_bw() +
  theme(legend.position = c(0.67, 0.88),
        legend.key.size = unit(0.2, "cm"),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank())

g_pca
ggsave("../figures/FigureS5.pdf", w = 3, h = 3)

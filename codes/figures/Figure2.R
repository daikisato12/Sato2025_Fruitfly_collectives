#### load packages ####
targetPackages <- c('tidyverse','tidytext','arrow','patchwork')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### Figure 2a ####
##### load dataset #####
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

###### gwas_metal_motion_cue_exit_intercept_scaleT ######
df_gwas_metal_motion_cue_exit_intercept_scaleT <- read_parquet("../data/1_single_strain/gwas/METAL/df/df_gwas_metal_motion_cue_exit_intercept_scaleT.parquet") %>%
  ungroup()

axisdf_gwas_metal_motion_cue_exit_intercept_scaleT <- df_gwas_metal_motion_cue_exit_intercept_scaleT %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

##### make plot #####
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
                     breaks = axisdf_gwas_motion_cue_exit_intercept_scaleT_female$center,
                     position = "right") +
  # scale_x_continuous(breaks = seq(0, 10, 2),
  #                    labels = seq(0, 10, 2)) +#scales::pretty_breaks(n = 6)) +
  scale_x_reverse(breaks = seq(0, 6, 2),
                  labels = seq(0, 6, 2)) +
  coord_cartesian(xlim = c(6.7, 0)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  ggtitle("Female") +
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
# g_gwas_motion_cue_exit_intercept_scaleT_female

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
                     breaks = axisdf_gwas_motion_cue_exit_intercept_scaleT_male$center) +
  scale_x_continuous(breaks = seq(0, 6, 2),
                     labels = seq(0, 6, 2)) +
  coord_cartesian(xlim = c(0, 6.7)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  ggtitle("Male") +
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
# g_gwas_motion_cue_exit_intercept_scaleT_male


g_gwas_motion_cue_exit_intercept_scaleT <-
  g_gwas_motion_cue_exit_intercept_scaleT_female +
  g_gwas_motion_cue_exit_intercept_scaleT_male +
  plot_layout(ncol = 2, widths = c(0.95, 1)) & 
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 2))


g_gwas_metal_motion_cue_exit_intercept_scaleT <- 
  ggplot(df_gwas_metal_motion_cue_exit_intercept_scaleT %>%
           dplyr::mutate(ptp99a = if_else(CHR == "3R" & BP == 25292746, "Y", "N")), 
         aes(y = BPcum, 
             x = -log10(P),
             color = as.factor(CHR),
             alpha = as.factor(ptp99a),
             size = as.factor(ptp99a))) +
  geom_point(shape = 16) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) +
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(label = axisdf_gwas_metal_motion_cue_exit_intercept_scaleT$CHR, 
                     breaks = axisdf_gwas_metal_motion_cue_exit_intercept_scaleT$center) +
  scale_x_continuous(breaks = seq(0, 10, 2),
                     labels = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(0, 11)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  ggtitle("Meta-analysis") +
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
# g_gwas_metal_motion_cue_exit_intercept_scaleT

g_gwas_motion_cue_exit_intercept_merge <-
  g_gwas_motion_cue_exit_intercept_scaleT +
  g_gwas_metal_motion_cue_exit_intercept_scaleT +
  plot_layout(ncol = 3, widths = c(0.8, 0.8, 1))


ggsave("../figures/Figure2a.png", 
       g_gwas_motion_cue_exit_intercept_merge, 
       width = 3.5, height = 6)

#### Figure 2b ####
##### load dataset #####
df_go <- list.files(paste0("../data/1_single_strain/gwas/METAL/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_metal/"), 
                    pattern = "*.tsv", full.names = TRUE, recursive=T) %>% 
  map_df(~ data.table::fread(., colClasses = c(ID = "character", 
                                               Description = "character",
                                               GeneRatio = "character",
                                               BgRatio = "character",
                                               geneID = "character",
                                               Type = "character"))) %>%
  dplyr::mutate(Description = str_to_sentence(Description),
                GeneRatio = str_split(GeneRatio, "/") %>% 
                  map_chr(1) %>% 
                  as.numeric() / str_split(GeneRatio, "/") %>% 
                  map_chr(2) %>% 
                  as.numeric()) %>%
  group_by(Type) %>%
  arrange(qvalue) %>%
  slice_head(n = 20) %>%
  dplyr::mutate(type = "Female") %>%
  filter(Type == "BP") %>%
  ungroup()


##### make plot #####
g_motion_cue_metal <- 
  ggplot(df_go, 
         aes(y = tidytext::reorder_within(Description, -qvalue, type), 
             x = GeneRatio, col = qvalue)) +
  geom_point(aes(size = Count), shape = 16) +#, size= 4) +
  tidytext::scale_y_reordered() +
  scale_color_viridis_c(name = expression(paste(italic(q), "-value"))) +
  # scale_shape_manual(values = c(17, 16, 15), labels = c("Biological process", "Cellular component", "Molecular Function")) +
  xlab("Gene ratio") +
  ylab("GO term") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(color = "black"),
        strip.background = element_blank())
g_motion_cue_metal

ggsave("../figures/Figure2b2.pdf", 
       g_motion_cue_metal, width = 5.5, height = 4.5)

#### Figure 2c ####
##### load dataset #####
gene_motion_cue_exit_female <- 
  read.table("../data/1_single_strain/gwas/result/tophits/tophits_motion_cue_exit_intercept_scaleT_female_geneid.txt") %>%
  pull(V1)

gene_motion_cue_exit_male <- 
  read.table("../data/1_single_strain/gwas/result/tophits/tophits_motion_cue_exit_intercept_scaleT_male_geneid.txt") %>%
  pull(V1)

mart <- biomaRt::useMart(biomart = "ensembl",  #useEnsembl or useMart
                         host = "https://www.ensembl.org", 
                         dataset = "dmelanogaster_gene_ensembl")

gene_shared <- intersect(gene_motion_cue_exit_female,
                         gene_motion_cue_exit_male)

gene_shared_name <- biomaRt::getBM(attributes = c("external_gene_name", 
                                                  "ensembl_gene_id"),
                                   filters = "ensembl_gene_id",
                                   values = gene_shared,
                                   mart = mart) %>%
  dplyr::rename(gene_name = external_gene_name)


gene_shared_name$gene_name %>% gtools::mixedsort()


#### Figure 2d ####
##### load dataset #####
df_pheno_motion_cue <- read.table("../data/1_single_strain/gwas/pheno/df_out_motion_cue_exit_intercept_female_scaled_GREML.txt") %>%
  dplyr::select(!V1) %>%
  magrittr::set_colnames(c("strain", "value")) %>%
  dplyr::mutate(sex = "Female") %>%
  bind_rows(read.table("../data/1_single_strain/gwas/pheno/df_out_motion_cue_exit_intercept_male_scaled_GREML.txt") %>%
              dplyr::select(!V1) %>%
              magrittr::set_colnames(c("strain", "value")) %>%
              dplyr::mutate(sex = "Male"))

df_geno_ptp99a <- read.table("../data/1_single_strain/gwas/result/tophits/candidates/genotype/3R_25292746_SNP.txt", header = TRUE)

df_motion_cue_ptp99a <- inner_join(df_pheno_motion_cue,
                                   df_geno_ptp99a)


##### stat #####
stats::wilcox.test(value ~ genotype, 
                   data = df_motion_cue_ptp99a %>% filter(sex == "Male"))
stats::wilcox.test(value ~ genotype, 
                   data = df_motion_cue_ptp99a %>% filter(sex == "Female"))

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
  ggpubr::stat_compare_means() +
  xlab("Genotype") +
  labs(y = expression("Visual responsiveness to motion cue" ~ (beta[0]))) +  
  facet_wrap(~ sex, ncol = 2) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text = element_text(color = "black"),
        strip.background = element_blank())
g_motion_cue_ptp99a

ggsave("../figures/Figure2d.pdf", g_motion_cue_ptp99a, w = 3, h = 2)


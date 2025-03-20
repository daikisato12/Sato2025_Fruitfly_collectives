#### load packages ####
targetPackages <- c('tidyverse','arrow','patchwork','Gviz','rtracklayer')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)
# BiocManager::install("Gviz")
# BiocManager::install("rtracklayer")
# BiocManager::install("ggbio")
# BiocManager::install("GenomicFeatures")

#### Figure 5a ####
##### load dataset #####
df_ghas <- read_parquet("../data/2_mixed_strain/gwas/result/rawdata/df_gwas_1kbp_performanceDE_3sec_regress_pi_grm.parquet") %>%
  ungroup()

axisdf_ghas <- df_ghas %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


##### make plot #####
g_ghas_manhattan <- 
  ggplot(df_ghas, 
         aes(y = BPcum, 
             x = -log10(P),
             color = as.factor(CHR))) +
  geom_point(shape = 16, size = 1.3, alpha = 0.2) +
  scale_color_manual(values = rep(c("#9e8896","#874c70"), 22)) +
  # scale_size_manual(values = c(1.3, 2.4)) +
  # scale_alpha_manual(values = c(0.2, 0.9)) +
  scale_fill_viridis_d() +
  scale_y_continuous(label = axisdf_ghas$CHR, 
                     breaks = axisdf_ghas$center) +
  scale_x_continuous(breaks = seq(0, 10, 2),
                     labels = seq(0, 10, 2)) +
  coord_cartesian(xlim = c(0, 11)) +
  ylab("Chromosome") +
  xlab(expression(paste("-Log"[10],"(",italic(P),")"))) +
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

g_ghas_manhattan
ggsave("../figures/Figure5a.png", g_ghas_manhattan, w = 2, h = 6)


#### Figure 5b ####
region <- GRanges("2R", IRanges(14665000, 14720000))
txdb <- GenomicFeatures::makeTxDbFromGFF("../data/0_genome/Drosophila_melanogaster.BDGP6.46.113.gtf", 
                                         format = "gtf")
p <- ggbio::autoplot(txdb, which = region) +
  coord_cartesian(xlim = c(14665000, 14720000)) +
  theme_bw()
pdf("../figures/Figure5b1.pdf", width=5, height=5)
p
dev.off()

g2 <- 
  ggplot(df_ghas %>%
           filter(CHR == "2R") %>%
           dplyr::mutate(plus = if_else(est > 0, "Y", "N")),
       aes(x = BP, y = -log10(P))) +
  geom_line(aes(color = plus)) +
  coord_cartesian(xlim = c(14665000, 14720000)) +
  ylab("-log(P)") +
  theme_bw() +
  theme(legend.position = "none")
g2
ggsave("../figures/Figure5b2.pdf", g2, w = 5, h = 3)


dfd_pi_1kbp2 <- read_parquet("../data/2_mixed_strain/gwas/pi/dfd_pi_1kbp2.parquet")
for_sec <- 3
dfd_pheno <- read_tsv("../data/2_mixed_strain/dfd_fig4.tsv") %>%
  filter(thr_sec == for_sec) %>%
  dplyr::mutate(DE = performance - lag(performance)) %>%
  filter(alpha == "Observed")

df_ghas %>% 
  arrange(P) %>% 
  head(1) %>%
  pull(BP)

dfd_pi_1kbp2_top <- dfd_pi_1kbp2 %>%
  filter(CHROM == "2R", BIN_START == 14689001) %>%
  left_join(dfd_pheno)

g_pi_de <-
  ggplot(dfd_pi_1kbp2_top,
       aes(x = PI, y = DE)) +
  geom_point(shape = 16, size = 3, alpha = 0.8) +
  stat_smooth(linewidth = 2, color= "grey", method = "lm") +
  ggpmisc::stat_poly_eq(formula = y ~ x,
                        aes(label = paste(
                          after_stat(rr.label),
                          stat(p.value.label),
                          sep = "~~~")),
                        label.x = "left",
                        label.y = "top",
                        parse = TRUE, size = 4) +
  xlab("pi") +
  ylab("Diversity effect on behavioral performance") +
  theme_bw()

g_pi_de
ggsave("../figures/Figure5b3.pdf", g_pi_de, w = 3, h = 3)


#### Figure 5c ####
##### load dataset #####
df_go <- list.files(paste0("../data/2_mixed_strain/gwas/result/enrichment/"),  #clusterProfiler_1kbp_foraging_success5_2000_regress_pi_mean_plus
                    pattern = "*.tsv", full.names = TRUE, recursive=T) %>% 
  map_df(~ data.table::fread(.) %>% mutate(Description = as.character(Description),
                                           ID = as.character(ID),
                                           GeneRatio = as.character(GeneRatio),
                                           BgRatio = as.character(BgRatio),
                                           geneID = as.character(geneID),
                                           Type = as.character(Type))) %>%
  mutate(Description = str_to_sentence(Description),
         GeneRatio = str_split(GeneRatio, "/") %>% 
           map_chr(1) %>% 
           as.numeric() / str_split(GeneRatio, "/") %>% 
           map_chr(2) %>% 
           as.numeric()) %>%
  group_by(Type) %>%
  arrange(qvalue) %>%
  slice_head(n = 20) %>%
  filter(Type == "BP") %>%
  ungroup()

##### make plot #####
g_1kbp_performanceDE_2000_regress_pi_grm_plus <- 
  ggplot(df_go, 
         aes(y = reorder(Description, GeneRatio), 
             x = GeneRatio, col = qvalue)) +
  geom_point(aes(size = Count), shape = 16) +
  tidytext::scale_y_reordered() +
  scale_color_viridis_c(name = expression(paste(italic(q), "-value"))) + #, limits = c(0.005, 0.02)
  # scale_shape_manual(values = c(17, 16, 15), labels = c("Biological process", "Cellular component", "Molecular Function")) +
  xlab("Gene ratio") +
  ylab("GO term") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text = element_text(color = "black"),
        strip.background = element_blank())


g_1kbp_performanceDE_2000_regress_pi_grm_plus
ggsave("../figures/Figure5c.pdf", g_1kbp_performanceDE_2000_regress_pi_grm_plus, w = 5.5, h = 3.5)

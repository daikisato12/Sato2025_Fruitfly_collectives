rm(list = ls(all.names = TRUE))

#### load libraries ####
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("pathview")
install.packages("wordcloud")

BiocManager::install("enrichplot")
install.packages("ggupset")
install.packages("ggnewscale")


library(clusterProfiler)
library(wordcloud)
library(ggupset)
library(ggnewscale)
library(enrichplot)
library(tidyverse)
library(tidytext)

organism = "org.Dm.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#### load dataset ####
genelist_background <- read.delim("../data/0_genome/dgrp2_1pop_SNPlist_geneid.txt", h=F)
genelist_background <- as.vector(genelist_background[,1])


#### motion_cue_exit_intercept_scaleT_female ####
genelist_motion_cue_exit_intercept_scaleT_female = read.delim("../data/1_single_strain/gwas/tophits/tophits_motion_cue_exit_intercept_scaleT_female_geneid.txt", header=F)
genelist_motion_cue_exit_intercept_scaleT_female <- as.vector(genelist_motion_cue_exit_intercept_scaleT_female[,1])
dir.create("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female", recursive = T)

##### BP #####

go_enrich_BP_motion_cue_exit_intercept_scaleT_female <- enrichGO(gene = genelist_motion_cue_exit_intercept_scaleT_female,
                                           universe = genelist_background,
                                           OrgDb = organism, 
                                           keyType = 'FLYBASE',
                                           readable = T,
                                           ont = "BP",
                                           pvalueCutoff = 0.05, 
                                           qvalueCutoff = 0.10)
head(go_enrich_BP_motion_cue_exit_intercept_scaleT_female)

df_enrich_BP_motion_cue_exit_intercept_scaleT_female <- go_enrich_BP_motion_cue_exit_intercept_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_enrich_BP_motion_cue_exit_intercept_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/motion_cue_exit_intercept_scaleT_female_BP.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_motion_cue_exit_intercept_scaleT_female <- upsetplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/UpsetPlot_BP.pdf", g1_BP_motion_cue_exit_intercept_scaleT_female, h=6, w=10)

g2_BP_motion_cue_exit_intercept_scaleT_female <- barplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_female, 
                                   drop = TRUE, 
                                   showCategory = 10, 
                                   title = "GO Biological Process",
                                   font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/BarPlot_BP.pdf", g2_BP_motion_cue_exit_intercept_scaleT_female, h=3, w=5)

g3_BP_motion_cue_exit_intercept_scaleT_female <- dotplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/DotPlot_BP.pdf", g3_BP_motion_cue_exit_intercept_scaleT_female, h=4, w=6)

g4_BP_motion_cue_exit_intercept_scaleT_female <- emapplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/EncrichmentMap_BP.pdf", g4_BP_motion_cue_exit_intercept_scaleT_female, h=6, w=8)

g5_BP_motion_cue_exit_intercept_scaleT_female <- goplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/EncrichedGOGraph_BP.pdf", g5_BP_motion_cue_exit_intercept_scaleT_female, h=6, w=8)

g6_BP_motion_cue_exit_intercept_scaleT_female <- cnetplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/CategoryNet_BP.pdf", g6_BP_motion_cue_exit_intercept_scaleT_female, h=12, w=10)

##### CC #####

go_enrich_CC_motion_cue_exit_intercept_scaleT_female <- enrichGO(gene = genelist_motion_cue_exit_intercept_scaleT_female,
                                           universe = genelist_background,
                                           OrgDb = organism, 
                                           keyType = 'FLYBASE',
                                           readable = T,
                                           ont = "CC",
                                           pvalueCutoff = 0.05, 
                                           qvalueCutoff = 0.10)
head(go_enrich_CC_motion_cue_exit_intercept_scaleT_female)

df_enrich_CC_motion_cue_exit_intercept_scaleT_female <- go_enrich_CC_motion_cue_exit_intercept_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_enrich_CC_motion_cue_exit_intercept_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/motion_cue_exit_intercept_scaleT_female_CC.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_motion_cue_exit_intercept_scaleT_female <- upsetplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/UpsetPlot_CC.pdf", g1_CC_motion_cue_exit_intercept_scaleT_female, h=6, w=10)

g2_CC_motion_cue_exit_intercept_scaleT_female <- barplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_female, 
                                   drop = TRUE, 
                                   showCategory = 10, 
                                   title = "GO Cellular Component",
                                   font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/BarPlot_CC.pdf", g2_CC_motion_cue_exit_intercept_scaleT_female, h=3, w=5)

g3_CC_motion_cue_exit_intercept_scaleT_female <- dotplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/DotPlot_CC.pdf", g3_CC_motion_cue_exit_intercept_scaleT_female, h=4, w=6)

g4_CC_motion_cue_exit_intercept_scaleT_female <- emapplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/EncrichmentMap_CC.pdf", g4_CC_motion_cue_exit_intercept_scaleT_female, h=6, w=8)

g5_CC_motion_cue_exit_intercept_scaleT_female <- goplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/EncrichedGOGraph_CC.pdf", g5_CC_motion_cue_exit_intercept_scaleT_female, h=6, w=8)

g6_CC_motion_cue_exit_intercept_scaleT_female <- cnetplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/CategoryNet_CC.pdf", g6_CC_motion_cue_exit_intercept_scaleT_female, h=12, w=10)

##### MF #####

go_enrich_MF_motion_cue_exit_intercept_scaleT_female <- enrichGO(gene = genelist_motion_cue_exit_intercept_scaleT_female,
                                           universe = genelist_background,
                                           OrgDb = organism, 
                                           keyType = 'FLYBASE',
                                           readable = T,
                                           ont = "MF",
                                           pvalueCutoff = 0.05, 
                                           qvalueCutoff = 0.10)
head(go_enrich_MF_motion_cue_exit_intercept_scaleT_female)

df_enrich_MF_motion_cue_exit_intercept_scaleT_female <- go_enrich_MF_motion_cue_exit_intercept_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_enrich_MF_motion_cue_exit_intercept_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/motion_cue_exit_intercept_scaleT_female_MF.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_motion_cue_exit_intercept_scaleT_female <- upsetplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/UpsetPlot_MF.pdf", g1_MF_motion_cue_exit_intercept_scaleT_female, h=6, w=10)

g2_MF_motion_cue_exit_intercept_scaleT_female <- barplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_female, 
                                   drop = TRUE, 
                                   showCategory = 10, 
                                   title = "GO Molecular Function",
                                   font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/BarPlot_MF.pdf", g2_MF_motion_cue_exit_intercept_scaleT_female, h=3, w=5)

g3_MF_motion_cue_exit_intercept_scaleT_female <- dotplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/DotPlot_MF.pdf", g3_MF_motion_cue_exit_intercept_scaleT_female, h=4, w=6)

g4_MF_motion_cue_exit_intercept_scaleT_female <- emapplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/EncrichmentMap_MF.pdf", g4_MF_motion_cue_exit_intercept_scaleT_female, h=6, w=8)

g5_MF_motion_cue_exit_intercept_scaleT_female <- goplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/EncrichedGOGraph_MF.pdf", g5_MF_motion_cue_exit_intercept_scaleT_female, h=6, w=8)

g6_MF_motion_cue_exit_intercept_scaleT_female <- cnetplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_female/CategoryNet_MF.pdf", g6_MF_motion_cue_exit_intercept_scaleT_female, h=12, w=10)


#### motion_cue_exit_intercept_scaleT_male ####
genelist_motion_cue_exit_intercept_scaleT_male = read.delim("../data/1_single_strain/gwas/tophits/tophits_motion_cue_exit_intercept_scaleT_male_geneid.txt", header=F)
genelist_motion_cue_exit_intercept_scaleT_male <- as.vector(genelist_motion_cue_exit_intercept_scaleT_male[,1])
dir.create("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male", recursive = T)

##### BP #####

go_enrich_BP_motion_cue_exit_intercept_scaleT_male <- enrichGO(gene = genelist_motion_cue_exit_intercept_scaleT_male,
                                                                 universe = genelist_background,
                                                                 OrgDb = organism, 
                                                                 keyType = 'FLYBASE',
                                                                 readable = T,
                                                                 ont = "BP",
                                                                 pvalueCutoff = 0.05, 
                                                                 qvalueCutoff = 0.10)
head(go_enrich_BP_motion_cue_exit_intercept_scaleT_male)

df_enrich_BP_motion_cue_exit_intercept_scaleT_male <- go_enrich_BP_motion_cue_exit_intercept_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_enrich_BP_motion_cue_exit_intercept_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/motion_cue_exit_intercept_scaleT_male_BP.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_motion_cue_exit_intercept_scaleT_male <- upsetplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/UpsetPlot_BP.pdf", g1_BP_motion_cue_exit_intercept_scaleT_male, h=6, w=10)

g2_BP_motion_cue_exit_intercept_scaleT_male <- barplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_male, 
                                                         drop = TRUE, 
                                                         showCategory = 10, 
                                                         title = "GO Biological Process",
                                                         font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/BarPlot_BP.pdf", g2_BP_motion_cue_exit_intercept_scaleT_male, h=3, w=5)

g3_BP_motion_cue_exit_intercept_scaleT_male <- dotplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/DotPlot_BP.pdf", g3_BP_motion_cue_exit_intercept_scaleT_male, h=4, w=6)

g4_BP_motion_cue_exit_intercept_scaleT_male <- emapplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/EncrichmentMap_BP.pdf", g4_BP_motion_cue_exit_intercept_scaleT_male, h=6, w=8)

g5_BP_motion_cue_exit_intercept_scaleT_male <- goplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/EncrichedGOGraph_BP.pdf", g5_BP_motion_cue_exit_intercept_scaleT_male, h=6, w=8)

g6_BP_motion_cue_exit_intercept_scaleT_male <- cnetplot(go_enrich_BP_motion_cue_exit_intercept_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/CategoryNet_BP.pdf", g6_BP_motion_cue_exit_intercept_scaleT_male, h=12, w=10)

##### CC #####

go_enrich_CC_motion_cue_exit_intercept_scaleT_male <- enrichGO(gene = genelist_motion_cue_exit_intercept_scaleT_male,
                                                                 universe = genelist_background,
                                                                 OrgDb = organism, 
                                                                 keyType = 'FLYBASE',
                                                                 readable = T,
                                                                 ont = "CC",
                                                                 pvalueCutoff = 0.05, 
                                                                 qvalueCutoff = 0.10)
head(go_enrich_CC_motion_cue_exit_intercept_scaleT_male)

df_enrich_CC_motion_cue_exit_intercept_scaleT_male <- go_enrich_CC_motion_cue_exit_intercept_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_enrich_CC_motion_cue_exit_intercept_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/motion_cue_exit_intercept_scaleT_male_CC.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_motion_cue_exit_intercept_scaleT_male <- upsetplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/UpsetPlot_CC.pdf", g1_CC_motion_cue_exit_intercept_scaleT_male, h=6, w=10)

g2_CC_motion_cue_exit_intercept_scaleT_male <- barplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_male, 
                                                         drop = TRUE, 
                                                         showCategory = 10, 
                                                         title = "GO Cellular Component",
                                                         font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/BarPlot_CC.pdf", g2_CC_motion_cue_exit_intercept_scaleT_male, h=3, w=5)

g3_CC_motion_cue_exit_intercept_scaleT_male <- dotplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/DotPlot_CC.pdf", g3_CC_motion_cue_exit_intercept_scaleT_male, h=4, w=6)

g4_CC_motion_cue_exit_intercept_scaleT_male <- emapplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/EncrichmentMap_CC.pdf", g4_CC_motion_cue_exit_intercept_scaleT_male, h=6, w=8)

g5_CC_motion_cue_exit_intercept_scaleT_male <- goplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/EncrichedGOGraph_CC.pdf", g5_CC_motion_cue_exit_intercept_scaleT_male, h=6, w=8)

g6_CC_motion_cue_exit_intercept_scaleT_male <- cnetplot(go_enrich_CC_motion_cue_exit_intercept_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/CategoryNet_CC.pdf", g6_CC_motion_cue_exit_intercept_scaleT_male, h=12, w=10)

##### MF #####

go_enrich_MF_motion_cue_exit_intercept_scaleT_male <- enrichGO(gene = genelist_motion_cue_exit_intercept_scaleT_male,
                                                                 universe = genelist_background,
                                                                 OrgDb = organism, 
                                                                 keyType = 'FLYBASE',
                                                                 readable = T,
                                                                 ont = "MF",
                                                                 pvalueCutoff = 0.05, 
                                                                 qvalueCutoff = 0.10)
head(go_enrich_MF_motion_cue_exit_intercept_scaleT_male)

df_enrich_MF_motion_cue_exit_intercept_scaleT_male <- go_enrich_MF_motion_cue_exit_intercept_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_enrich_MF_motion_cue_exit_intercept_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/motion_cue_exit_intercept_scaleT_male_MF.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_motion_cue_exit_intercept_scaleT_male <- upsetplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/UpsetPlot_MF.pdf", g1_MF_motion_cue_exit_intercept_scaleT_male, h=6, w=10)

g2_MF_motion_cue_exit_intercept_scaleT_male <- barplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_male, 
                                                         drop = TRUE, 
                                                         showCategory = 10, 
                                                         title = "GO Molecular Function",
                                                         font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/BarPlot_MF.pdf", g2_MF_motion_cue_exit_intercept_scaleT_male, h=3, w=5)

g3_MF_motion_cue_exit_intercept_scaleT_male <- dotplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/DotPlot_MF.pdf", g3_MF_motion_cue_exit_intercept_scaleT_male, h=4, w=6)

g4_MF_motion_cue_exit_intercept_scaleT_male <- emapplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/EncrichmentMap_MF.pdf", g4_MF_motion_cue_exit_intercept_scaleT_male, h=6, w=8)

g5_MF_motion_cue_exit_intercept_scaleT_male <- goplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/EncrichedGOGraph_MF.pdf", g5_MF_motion_cue_exit_intercept_scaleT_male, h=6, w=8)

g6_MF_motion_cue_exit_intercept_scaleT_male <- cnetplot(go_enrich_MF_motion_cue_exit_intercept_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_motion_cue_exit_intercept_scaleT_male/CategoryNet_MF.pdf", g6_MF_motion_cue_exit_intercept_scaleT_male, h=12, w=10)



#### nnd_scaleT_female ####
genelist_nnd_scaleT_female = read.delim("../data/1_single_strain/gwas/tophits/tophits_nnd_scaleT_female_geneid.txt", header=F)
genelist_nnd_scaleT_female <- as.vector(genelist_nnd_scaleT_female[,1])
dir.create("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female", recursive = T)

##### BP #####

go_enrich_BP_nnd_scaleT_female <- enrichGO(gene = genelist_nnd_scaleT_female,
                                                                 universe = genelist_background,
                                                                 OrgDb = organism, 
                                                                 keyType = 'FLYBASE',
                                                                 readable = T,
                                                                 ont = "BP",
                                                                 pvalueCutoff = 0.05, 
                                                                 qvalueCutoff = 0.10)
head(go_enrich_BP_nnd_scaleT_female)

df_enrich_BP_nnd_scaleT_female <- go_enrich_BP_nnd_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_enrich_BP_nnd_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/nnd_scaleT_female_BP.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_nnd_scaleT_female <- upsetplot(go_enrich_BP_nnd_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/UpsetPlot_BP.pdf", g1_BP_nnd_scaleT_female, h=6, w=10)

g2_BP_nnd_scaleT_female <- barplot(go_enrich_BP_nnd_scaleT_female, 
                                                         drop = TRUE, 
                                                         showCategory = 10, 
                                                         title = "GO Biological Process",
                                                         font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/BarPlot_BP.pdf", g2_BP_nnd_scaleT_female, h=3, w=5)

g3_BP_nnd_scaleT_female <- dotplot(go_enrich_BP_nnd_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/DotPlot_BP.pdf", g3_BP_nnd_scaleT_female, h=4, w=6)

g4_BP_nnd_scaleT_female <- emapplot(go_enrich_BP_nnd_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/EncrichmentMap_BP.pdf", g4_BP_nnd_scaleT_female, h=6, w=8)

g5_BP_nnd_scaleT_female <- goplot(go_enrich_BP_nnd_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/EncrichedGOGraph_BP.pdf", g5_BP_nnd_scaleT_female, h=6, w=8)

g6_BP_nnd_scaleT_female <- cnetplot(go_enrich_BP_nnd_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/CategoryNet_BP.pdf", g6_BP_nnd_scaleT_female, h=12, w=10)

##### CC #####

go_enrich_CC_nnd_scaleT_female <- enrichGO(gene = genelist_nnd_scaleT_female,
                                                                 universe = genelist_background,
                                                                 OrgDb = organism, 
                                                                 keyType = 'FLYBASE',
                                                                 readable = T,
                                                                 ont = "CC",
                                                                 pvalueCutoff = 0.05, 
                                                                 qvalueCutoff = 0.10)
head(go_enrich_CC_nnd_scaleT_female)

df_enrich_CC_nnd_scaleT_female <- go_enrich_CC_nnd_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_enrich_CC_nnd_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/nnd_scaleT_female_CC.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_nnd_scaleT_female <- upsetplot(go_enrich_CC_nnd_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/UpsetPlot_CC.pdf", g1_CC_nnd_scaleT_female, h=6, w=10)

g2_CC_nnd_scaleT_female <- barplot(go_enrich_CC_nnd_scaleT_female, 
                                                         drop = TRUE, 
                                                         showCategory = 10, 
                                                         title = "GO Cellular Component",
                                                         font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/BarPlot_CC.pdf", g2_CC_nnd_scaleT_female, h=3, w=5)

g3_CC_nnd_scaleT_female <- dotplot(go_enrich_CC_nnd_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/DotPlot_CC.pdf", g3_CC_nnd_scaleT_female, h=4, w=6)

g4_CC_nnd_scaleT_female <- emapplot(go_enrich_CC_nnd_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/EncrichmentMap_CC.pdf", g4_CC_nnd_scaleT_female, h=6, w=8)

g5_CC_nnd_scaleT_female <- goplot(go_enrich_CC_nnd_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/EncrichedGOGraph_CC.pdf", g5_CC_nnd_scaleT_female, h=6, w=8)

g6_CC_nnd_scaleT_female <- cnetplot(go_enrich_CC_nnd_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/CategoryNet_CC.pdf", g6_CC_nnd_scaleT_female, h=12, w=10)

##### MF #####

go_enrich_MF_nnd_scaleT_female <- enrichGO(gene = genelist_nnd_scaleT_female,
                                                                 universe = genelist_background,
                                                                 OrgDb = organism, 
                                                                 keyType = 'FLYBASE',
                                                                 readable = T,
                                                                 ont = "MF",
                                                                 pvalueCutoff = 0.05, 
                                                                 qvalueCutoff = 0.10)
head(go_enrich_MF_nnd_scaleT_female)

df_enrich_MF_nnd_scaleT_female <- go_enrich_MF_nnd_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_enrich_MF_nnd_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/nnd_scaleT_female_MF.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_nnd_scaleT_female <- upsetplot(go_enrich_MF_nnd_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/UpsetPlot_MF.pdf", g1_MF_nnd_scaleT_female, h=6, w=10)

g2_MF_nnd_scaleT_female <- barplot(go_enrich_MF_nnd_scaleT_female, 
                                                         drop = TRUE, 
                                                         showCategory = 10, 
                                                         title = "GO Molecular Function",
                                                         font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/BarPlot_MF.pdf", g2_MF_nnd_scaleT_female, h=3, w=5)

g3_MF_nnd_scaleT_female <- dotplot(go_enrich_MF_nnd_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/DotPlot_MF.pdf", g3_MF_nnd_scaleT_female, h=4, w=6)

g4_MF_nnd_scaleT_female <- emapplot(go_enrich_MF_nnd_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/EncrichmentMap_MF.pdf", g4_MF_nnd_scaleT_female, h=6, w=8)

g5_MF_nnd_scaleT_female <- goplot(go_enrich_MF_nnd_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/EncrichedGOGraph_MF.pdf", g5_MF_nnd_scaleT_female, h=6, w=8)

g6_MF_nnd_scaleT_female <- cnetplot(go_enrich_MF_nnd_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_female/CategoryNet_MF.pdf", g6_MF_nnd_scaleT_female, h=12, w=10)


#### nnd_scaleT_male ####
genelist_nnd_scaleT_male = read.delim("../data/1_single_strain/gwas/tophits/tophits_nnd_scaleT_male_geneid.txt", header=F)
genelist_nnd_scaleT_male <- as.vector(genelist_nnd_scaleT_male[,1])
dir.create("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male", recursive = T)

##### BP #####

go_enrich_BP_nnd_scaleT_male <- enrichGO(gene = genelist_nnd_scaleT_male,
                                                               universe = genelist_background,
                                                               OrgDb = organism, 
                                                               keyType = 'FLYBASE',
                                                               readable = T,
                                                               ont = "BP",
                                                               pvalueCutoff = 0.05, 
                                                               qvalueCutoff = 0.10)
head(go_enrich_BP_nnd_scaleT_male)

df_enrich_BP_nnd_scaleT_male <- go_enrich_BP_nnd_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_enrich_BP_nnd_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/nnd_scaleT_male_BP.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_nnd_scaleT_male <- upsetplot(go_enrich_BP_nnd_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/UpsetPlot_BP.pdf", g1_BP_nnd_scaleT_male, h=6, w=10)

g2_BP_nnd_scaleT_male <- barplot(go_enrich_BP_nnd_scaleT_male, 
                                                       drop = TRUE, 
                                                       showCategory = 10, 
                                                       title = "GO Biological Process",
                                                       font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/BarPlot_BP.pdf", g2_BP_nnd_scaleT_male, h=3, w=5)

g3_BP_nnd_scaleT_male <- dotplot(go_enrich_BP_nnd_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/DotPlot_BP.pdf", g3_BP_nnd_scaleT_male, h=4, w=6)

g4_BP_nnd_scaleT_male <- emapplot(go_enrich_BP_nnd_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/EncrichmentMap_BP.pdf", g4_BP_nnd_scaleT_male, h=6, w=8)

g5_BP_nnd_scaleT_male <- goplot(go_enrich_BP_nnd_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/EncrichedGOGraph_BP.pdf", g5_BP_nnd_scaleT_male, h=6, w=8)

g6_BP_nnd_scaleT_male <- cnetplot(go_enrich_BP_nnd_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/CategoryNet_BP.pdf", g6_BP_nnd_scaleT_male, h=12, w=10)

##### CC #####

go_enrich_CC_nnd_scaleT_male <- enrichGO(gene = genelist_nnd_scaleT_male,
                                                               universe = genelist_background,
                                                               OrgDb = organism, 
                                                               keyType = 'FLYBASE',
                                                               readable = T,
                                                               ont = "CC",
                                                               pvalueCutoff = 0.05, 
                                                               qvalueCutoff = 0.10)
head(go_enrich_CC_nnd_scaleT_male)

df_enrich_CC_nnd_scaleT_male <- go_enrich_CC_nnd_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_enrich_CC_nnd_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/nnd_scaleT_male_CC.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_nnd_scaleT_male <- upsetplot(go_enrich_CC_nnd_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/UpsetPlot_CC.pdf", g1_CC_nnd_scaleT_male, h=6, w=10)

g2_CC_nnd_scaleT_male <- barplot(go_enrich_CC_nnd_scaleT_male, 
                                                       drop = TRUE, 
                                                       showCategory = 10, 
                                                       title = "GO Cellular Component",
                                                       font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/BarPlot_CC.pdf", g2_CC_nnd_scaleT_male, h=3, w=5)

g3_CC_nnd_scaleT_male <- dotplot(go_enrich_CC_nnd_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/DotPlot_CC.pdf", g3_CC_nnd_scaleT_male, h=4, w=6)

g4_CC_nnd_scaleT_male <- emapplot(go_enrich_CC_nnd_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/EncrichmentMap_CC.pdf", g4_CC_nnd_scaleT_male, h=6, w=8)

g5_CC_nnd_scaleT_male <- goplot(go_enrich_CC_nnd_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/EncrichedGOGraph_CC.pdf", g5_CC_nnd_scaleT_male, h=6, w=8)

g6_CC_nnd_scaleT_male <- cnetplot(go_enrich_CC_nnd_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/CategoryNet_CC.pdf", g6_CC_nnd_scaleT_male, h=12, w=10)

##### MF #####

go_enrich_MF_nnd_scaleT_male <- enrichGO(gene = genelist_nnd_scaleT_male,
                                                               universe = genelist_background,
                                                               OrgDb = organism, 
                                                               keyType = 'FLYBASE',
                                                               readable = T,
                                                               ont = "MF",
                                                               pvalueCutoff = 0.05, 
                                                               qvalueCutoff = 0.10)
head(go_enrich_MF_nnd_scaleT_male)

df_enrich_MF_nnd_scaleT_male <- go_enrich_MF_nnd_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_enrich_MF_nnd_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/nnd_scaleT_male_MF.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_nnd_scaleT_male <- upsetplot(go_enrich_MF_nnd_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/UpsetPlot_MF.pdf", g1_MF_nnd_scaleT_male, h=6, w=10)

g2_MF_nnd_scaleT_male <- barplot(go_enrich_MF_nnd_scaleT_male, 
                                                       drop = TRUE, 
                                                       showCategory = 10, 
                                                       title = "GO Molecular Function",
                                                       font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/BarPlot_MF.pdf", g2_MF_nnd_scaleT_male, h=3, w=5)

g3_MF_nnd_scaleT_male <- dotplot(go_enrich_MF_nnd_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/DotPlot_MF.pdf", g3_MF_nnd_scaleT_male, h=4, w=6)

g4_MF_nnd_scaleT_male <- emapplot(go_enrich_MF_nnd_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/EncrichmentMap_MF.pdf", g4_MF_nnd_scaleT_male, h=6, w=8)

g5_MF_nnd_scaleT_male <- goplot(go_enrich_MF_nnd_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/EncrichedGOGraph_MF.pdf", g5_MF_nnd_scaleT_male, h=6, w=8)

g6_MF_nnd_scaleT_male <- cnetplot(go_enrich_MF_nnd_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_nnd_scaleT_male/CategoryNet_MF.pdf", g6_MF_nnd_scaleT_male, h=12, w=10)


#### freezing_duration_single_scaleT_female ####
genelist_freezing_duration_single_scaleT_female = read.delim("../data/1_single_strain/gwas/tophits/tophits_freezing_duration_single_scaleT_female_geneid.txt", header=F)
genelist_freezing_duration_single_scaleT_female <- as.vector(genelist_freezing_duration_single_scaleT_female[,1])
dir.create("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female", recursive = T)

##### BP #####

go_enrich_BP_freezing_duration_single_scaleT_female <- enrichGO(gene = genelist_freezing_duration_single_scaleT_female,
                                           universe = genelist_background,
                                           OrgDb = organism, 
                                           keyType = 'FLYBASE',
                                           readable = T,
                                           ont = "BP",
                                           pvalueCutoff = 0.05, 
                                           qvalueCutoff = 0.10)
head(go_enrich_BP_freezing_duration_single_scaleT_female)

df_enrich_BP_freezing_duration_single_scaleT_female <- go_enrich_BP_freezing_duration_single_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_enrich_BP_freezing_duration_single_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/freezing_duration_single_scaleT_female_BP.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_freezing_duration_single_scaleT_female <- upsetplot(go_enrich_BP_freezing_duration_single_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/UpsetPlot_BP.pdf", g1_BP_freezing_duration_single_scaleT_female, h=6, w=10)

g2_BP_freezing_duration_single_scaleT_female <- barplot(go_enrich_BP_freezing_duration_single_scaleT_female, 
                                   drop = TRUE, 
                                   showCategory = 10, 
                                   title = "GO Biological Process",
                                   font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/BarPlot_BP.pdf", g2_BP_freezing_duration_single_scaleT_female, h=3, w=5)

g3_BP_freezing_duration_single_scaleT_female <- dotplot(go_enrich_BP_freezing_duration_single_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/DotPlot_BP.pdf", g3_BP_freezing_duration_single_scaleT_female, h=4, w=6)

g4_BP_freezing_duration_single_scaleT_female <- emapplot(go_enrich_BP_freezing_duration_single_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/EncrichmentMap_BP.pdf", g4_BP_freezing_duration_single_scaleT_female, h=6, w=8)

g5_BP_freezing_duration_single_scaleT_female <- goplot(go_enrich_BP_freezing_duration_single_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/EncrichedGOGraph_BP.pdf", g5_BP_freezing_duration_single_scaleT_female, h=6, w=8)

g6_BP_freezing_duration_single_scaleT_female <- cnetplot(go_enrich_BP_freezing_duration_single_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/CategoryNet_BP.pdf", g6_BP_freezing_duration_single_scaleT_female, h=12, w=10)

##### CC #####

go_enrich_CC_freezing_duration_single_scaleT_female <- enrichGO(gene = genelist_freezing_duration_single_scaleT_female,
                                           universe = genelist_background,
                                           OrgDb = organism, 
                                           keyType = 'FLYBASE',
                                           readable = T,
                                           ont = "CC",
                                           pvalueCutoff = 0.05, 
                                           qvalueCutoff = 0.10)
head(go_enrich_CC_freezing_duration_single_scaleT_female)

df_enrich_CC_freezing_duration_single_scaleT_female <- go_enrich_CC_freezing_duration_single_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_enrich_CC_freezing_duration_single_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/freezing_duration_single_scaleT_female_CC.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_freezing_duration_single_scaleT_female <- upsetplot(go_enrich_CC_freezing_duration_single_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/UpsetPlot_CC.pdf", g1_CC_freezing_duration_single_scaleT_female, h=6, w=10)

g2_CC_freezing_duration_single_scaleT_female <- barplot(go_enrich_CC_freezing_duration_single_scaleT_female, 
                                   drop = TRUE, 
                                   showCategory = 10, 
                                   title = "GO Cellular Component",
                                   font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/BarPlot_CC.pdf", g2_CC_freezing_duration_single_scaleT_female, h=3, w=5)

g3_CC_freezing_duration_single_scaleT_female <- dotplot(go_enrich_CC_freezing_duration_single_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/DotPlot_CC.pdf", g3_CC_freezing_duration_single_scaleT_female, h=4, w=6)

g4_CC_freezing_duration_single_scaleT_female <- emapplot(go_enrich_CC_freezing_duration_single_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/EncrichmentMap_CC.pdf", g4_CC_freezing_duration_single_scaleT_female, h=6, w=8)

g5_CC_freezing_duration_single_scaleT_female <- goplot(go_enrich_CC_freezing_duration_single_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/EncrichedGOGraph_CC.pdf", g5_CC_freezing_duration_single_scaleT_female, h=6, w=8)

g6_CC_freezing_duration_single_scaleT_female <- cnetplot(go_enrich_CC_freezing_duration_single_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/CategoryNet_CC.pdf", g6_CC_freezing_duration_single_scaleT_female, h=12, w=10)

##### MF #####

go_enrich_MF_freezing_duration_single_scaleT_female <- enrichGO(gene = genelist_freezing_duration_single_scaleT_female,
                                           universe = genelist_background,
                                           OrgDb = organism, 
                                           keyType = 'FLYBASE',
                                           readable = T,
                                           ont = "MF",
                                           pvalueCutoff = 0.05, 
                                           qvalueCutoff = 0.10)
head(go_enrich_MF_freezing_duration_single_scaleT_female)

df_enrich_MF_freezing_duration_single_scaleT_female <- go_enrich_MF_freezing_duration_single_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_enrich_MF_freezing_duration_single_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/freezing_duration_single_scaleT_female_MF.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_freezing_duration_single_scaleT_female <- upsetplot(go_enrich_MF_freezing_duration_single_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/UpsetPlot_MF.pdf", g1_MF_freezing_duration_single_scaleT_female, h=6, w=10)

g2_MF_freezing_duration_single_scaleT_female <- barplot(go_enrich_MF_freezing_duration_single_scaleT_female, 
                                   drop = TRUE, 
                                   showCategory = 10, 
                                   title = "GO Molecular Function",
                                   font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/BarPlot_MF.pdf", g2_MF_freezing_duration_single_scaleT_female, h=3, w=5)

g3_MF_freezing_duration_single_scaleT_female <- dotplot(go_enrich_MF_freezing_duration_single_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/DotPlot_MF.pdf", g3_MF_freezing_duration_single_scaleT_female, h=4, w=6)

g4_MF_freezing_duration_single_scaleT_female <- emapplot(go_enrich_MF_freezing_duration_single_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/EncrichmentMap_MF.pdf", g4_MF_freezing_duration_single_scaleT_female, h=6, w=8)

g5_MF_freezing_duration_single_scaleT_female <- goplot(go_enrich_MF_freezing_duration_single_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/EncrichedGOGraph_MF.pdf", g5_MF_freezing_duration_single_scaleT_female, h=6, w=8)

g6_MF_freezing_duration_single_scaleT_female <- cnetplot(go_enrich_MF_freezing_duration_single_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_female/CategoryNet_MF.pdf", g6_MF_freezing_duration_single_scaleT_female, h=12, w=10)


#### freezing_duration_single_scaleT_male ####
genelist_freezing_duration_single_scaleT_male = read.delim("../data/1_single_strain/gwas/tophits/tophits_freezing_duration_single_scaleT_male_geneid.txt", header=F)
genelist_freezing_duration_single_scaleT_male <- as.vector(genelist_freezing_duration_single_scaleT_male[,1])
dir.create("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male", recursive = T)

##### BP #####

go_enrich_BP_freezing_duration_single_scaleT_male <- enrichGO(gene = genelist_freezing_duration_single_scaleT_male,
                                         universe = genelist_background,
                                         OrgDb = organism, 
                                         keyType = 'FLYBASE',
                                         readable = T,
                                         ont = "BP",
                                         pvalueCutoff = 0.05, 
                                         qvalueCutoff = 0.10)
head(go_enrich_BP_freezing_duration_single_scaleT_male)

df_enrich_BP_freezing_duration_single_scaleT_male <- go_enrich_BP_freezing_duration_single_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_enrich_BP_freezing_duration_single_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/freezing_duration_single_scaleT_male_BP.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_freezing_duration_single_scaleT_male <- upsetplot(go_enrich_BP_freezing_duration_single_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/UpsetPlot_BP.pdf", g1_BP_freezing_duration_single_scaleT_male, h=6, w=10)

g2_BP_freezing_duration_single_scaleT_male <- barplot(go_enrich_BP_freezing_duration_single_scaleT_male, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Biological Process",
                                 font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/BarPlot_BP.pdf", g2_BP_freezing_duration_single_scaleT_male, h=3, w=5)

g3_BP_freezing_duration_single_scaleT_male <- dotplot(go_enrich_BP_freezing_duration_single_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/DotPlot_BP.pdf", g3_BP_freezing_duration_single_scaleT_male, h=4, w=6)

g4_BP_freezing_duration_single_scaleT_male <- emapplot(go_enrich_BP_freezing_duration_single_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/EncrichmentMap_BP.pdf", g4_BP_freezing_duration_single_scaleT_male, h=6, w=8)

g5_BP_freezing_duration_single_scaleT_male <- goplot(go_enrich_BP_freezing_duration_single_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/EncrichedGOGraph_BP.pdf", g5_BP_freezing_duration_single_scaleT_male, h=6, w=8)

g6_BP_freezing_duration_single_scaleT_male <- cnetplot(go_enrich_BP_freezing_duration_single_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/CategoryNet_BP.pdf", g6_BP_freezing_duration_single_scaleT_male, h=12, w=10)

##### CC #####

go_enrich_CC_freezing_duration_single_scaleT_male <- enrichGO(gene = genelist_freezing_duration_single_scaleT_male,
                                         universe = genelist_background,
                                         OrgDb = organism, 
                                         keyType = 'FLYBASE',
                                         readable = T,
                                         ont = "CC",
                                         pvalueCutoff = 0.05, 
                                         qvalueCutoff = 0.10)
head(go_enrich_CC_freezing_duration_single_scaleT_male)

df_enrich_CC_freezing_duration_single_scaleT_male <- go_enrich_CC_freezing_duration_single_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_enrich_CC_freezing_duration_single_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/freezing_duration_single_scaleT_male_CC.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_freezing_duration_single_scaleT_male <- upsetplot(go_enrich_CC_freezing_duration_single_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/UpsetPlot_CC.pdf", g1_CC_freezing_duration_single_scaleT_male, h=6, w=10)

g2_CC_freezing_duration_single_scaleT_male <- barplot(go_enrich_CC_freezing_duration_single_scaleT_male, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Cellular Component",
                                 font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/BarPlot_CC.pdf", g2_CC_freezing_duration_single_scaleT_male, h=3, w=5)

g3_CC_freezing_duration_single_scaleT_male <- dotplot(go_enrich_CC_freezing_duration_single_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/DotPlot_CC.pdf", g3_CC_freezing_duration_single_scaleT_male, h=4, w=6)

g4_CC_freezing_duration_single_scaleT_male <- emapplot(go_enrich_CC_freezing_duration_single_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/EncrichmentMap_CC.pdf", g4_CC_freezing_duration_single_scaleT_male, h=6, w=8)

g5_CC_freezing_duration_single_scaleT_male <- goplot(go_enrich_CC_freezing_duration_single_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/EncrichedGOGraph_CC.pdf", g5_CC_freezing_duration_single_scaleT_male, h=6, w=8)

g6_CC_freezing_duration_single_scaleT_male <- cnetplot(go_enrich_CC_freezing_duration_single_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/CategoryNet_CC.pdf", g6_CC_freezing_duration_single_scaleT_male, h=12, w=10)

##### MF #####

go_enrich_MF_freezing_duration_single_scaleT_male <- enrichGO(gene = genelist_freezing_duration_single_scaleT_male,
                                         universe = genelist_background,
                                         OrgDb = organism, 
                                         keyType = 'FLYBASE',
                                         readable = T,
                                         ont = "MF",
                                         pvalueCutoff = 0.05, 
                                         qvalueCutoff = 0.10)
head(go_enrich_MF_freezing_duration_single_scaleT_male)

df_enrich_MF_freezing_duration_single_scaleT_male <- go_enrich_MF_freezing_duration_single_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_enrich_MF_freezing_duration_single_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/freezing_duration_single_scaleT_male_MF.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_freezing_duration_single_scaleT_male <- upsetplot(go_enrich_MF_freezing_duration_single_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/UpsetPlot_MF.pdf", g1_MF_freezing_duration_single_scaleT_male, h=6, w=10)

g2_MF_freezing_duration_single_scaleT_male <- barplot(go_enrich_MF_freezing_duration_single_scaleT_male, 
                                 drop = TRUE, 
                                 showCategory = 10, 
                                 title = "GO Molecular Function",
                                 font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/BarPlot_MF.pdf", g2_MF_freezing_duration_single_scaleT_male, h=3, w=5)

g3_MF_freezing_duration_single_scaleT_male <- dotplot(go_enrich_MF_freezing_duration_single_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/DotPlot_MF.pdf", g3_MF_freezing_duration_single_scaleT_male, h=4, w=6)

g4_MF_freezing_duration_single_scaleT_male <- emapplot(go_enrich_MF_freezing_duration_single_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/EncrichmentMap_MF.pdf", g4_MF_freezing_duration_single_scaleT_male, h=6, w=8)

g5_MF_freezing_duration_single_scaleT_male <- goplot(go_enrich_MF_freezing_duration_single_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/EncrichedGOGraph_MF.pdf", g5_MF_freezing_duration_single_scaleT_male, h=6, w=8)

g6_MF_freezing_duration_single_scaleT_male <- cnetplot(go_enrich_MF_freezing_duration_single_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_single_scaleT_male/CategoryNet_MF.pdf", g6_MF_freezing_duration_single_scaleT_male, h=12, w=10)




#### freezing_duration_group_scaleT_female ####
genelist_freezing_duration_group_scaleT_female = read.delim("../data/1_single_strain/gwas/tophits/tophits_freezing_duration_group_scaleT_female_geneid.txt", header=F)
genelist_freezing_duration_group_scaleT_female <- as.vector(genelist_freezing_duration_group_scaleT_female[,1])
dir.create("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female", recursive = T)

##### BP #####

go_enrich_BP_freezing_duration_group_scaleT_female <- enrichGO(gene = genelist_freezing_duration_group_scaleT_female,
                                                                 universe = genelist_background,
                                                                 OrgDb = organism, 
                                                                 keyType = 'FLYBASE',
                                                                 readable = T,
                                                                 ont = "BP",
                                                                 pvalueCutoff = 0.05, 
                                                                 qvalueCutoff = 0.10)
head(go_enrich_BP_freezing_duration_group_scaleT_female)

df_enrich_BP_freezing_duration_group_scaleT_female <- go_enrich_BP_freezing_duration_group_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_enrich_BP_freezing_duration_group_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/freezing_duration_group_scaleT_female_BP.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_freezing_duration_group_scaleT_female <- upsetplot(go_enrich_BP_freezing_duration_group_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/UpsetPlot_BP.pdf", g1_BP_freezing_duration_group_scaleT_female, h=6, w=10)

g2_BP_freezing_duration_group_scaleT_female <- barplot(go_enrich_BP_freezing_duration_group_scaleT_female, 
                                                         drop = TRUE, 
                                                         showCategory = 10, 
                                                         title = "GO Biological Process",
                                                         font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/BarPlot_BP.pdf", g2_BP_freezing_duration_group_scaleT_female, h=3, w=5)

g3_BP_freezing_duration_group_scaleT_female <- dotplot(go_enrich_BP_freezing_duration_group_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/DotPlot_BP.pdf", g3_BP_freezing_duration_group_scaleT_female, h=4, w=6)

g4_BP_freezing_duration_group_scaleT_female <- emapplot(go_enrich_BP_freezing_duration_group_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/EncrichmentMap_BP.pdf", g4_BP_freezing_duration_group_scaleT_female, h=6, w=8)

g5_BP_freezing_duration_group_scaleT_female <- goplot(go_enrich_BP_freezing_duration_group_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/EncrichedGOGraph_BP.pdf", g5_BP_freezing_duration_group_scaleT_female, h=6, w=8)

g6_BP_freezing_duration_group_scaleT_female <- cnetplot(go_enrich_BP_freezing_duration_group_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/CategoryNet_BP.pdf", g6_BP_freezing_duration_group_scaleT_female, h=12, w=10)

##### CC #####

go_enrich_CC_freezing_duration_group_scaleT_female <- enrichGO(gene = genelist_freezing_duration_group_scaleT_female,
                                                                 universe = genelist_background,
                                                                 OrgDb = organism, 
                                                                 keyType = 'FLYBASE',
                                                                 readable = T,
                                                                 ont = "CC",
                                                                 pvalueCutoff = 0.05, 
                                                                 qvalueCutoff = 0.10)
head(go_enrich_CC_freezing_duration_group_scaleT_female)

df_enrich_CC_freezing_duration_group_scaleT_female <- go_enrich_CC_freezing_duration_group_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_enrich_CC_freezing_duration_group_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/freezing_duration_group_scaleT_female_CC.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_freezing_duration_group_scaleT_female <- upsetplot(go_enrich_CC_freezing_duration_group_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/UpsetPlot_CC.pdf", g1_CC_freezing_duration_group_scaleT_female, h=6, w=10)

g2_CC_freezing_duration_group_scaleT_female <- barplot(go_enrich_CC_freezing_duration_group_scaleT_female, 
                                                         drop = TRUE, 
                                                         showCategory = 10, 
                                                         title = "GO Cellular Component",
                                                         font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/BarPlot_CC.pdf", g2_CC_freezing_duration_group_scaleT_female, h=3, w=5)

g3_CC_freezing_duration_group_scaleT_female <- dotplot(go_enrich_CC_freezing_duration_group_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/DotPlot_CC.pdf", g3_CC_freezing_duration_group_scaleT_female, h=4, w=6)

g4_CC_freezing_duration_group_scaleT_female <- emapplot(go_enrich_CC_freezing_duration_group_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/EncrichmentMap_CC.pdf", g4_CC_freezing_duration_group_scaleT_female, h=6, w=8)

g5_CC_freezing_duration_group_scaleT_female <- goplot(go_enrich_CC_freezing_duration_group_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/EncrichedGOGraph_CC.pdf", g5_CC_freezing_duration_group_scaleT_female, h=6, w=8)

g6_CC_freezing_duration_group_scaleT_female <- cnetplot(go_enrich_CC_freezing_duration_group_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/CategoryNet_CC.pdf", g6_CC_freezing_duration_group_scaleT_female, h=12, w=10)

##### MF #####

go_enrich_MF_freezing_duration_group_scaleT_female <- enrichGO(gene = genelist_freezing_duration_group_scaleT_female,
                                                                 universe = genelist_background,
                                                                 OrgDb = organism, 
                                                                 keyType = 'FLYBASE',
                                                                 readable = T,
                                                                 ont = "MF",
                                                                 pvalueCutoff = 0.05, 
                                                                 qvalueCutoff = 0.10)
head(go_enrich_MF_freezing_duration_group_scaleT_female)

df_enrich_MF_freezing_duration_group_scaleT_female <- go_enrich_MF_freezing_duration_group_scaleT_female@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_enrich_MF_freezing_duration_group_scaleT_female, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/freezing_duration_group_scaleT_female_MF.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_freezing_duration_group_scaleT_female <- upsetplot(go_enrich_MF_freezing_duration_group_scaleT_female)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/UpsetPlot_MF.pdf", g1_MF_freezing_duration_group_scaleT_female, h=6, w=10)

g2_MF_freezing_duration_group_scaleT_female <- barplot(go_enrich_MF_freezing_duration_group_scaleT_female, 
                                                         drop = TRUE, 
                                                         showCategory = 10, 
                                                         title = "GO Molecular Function",
                                                         font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/BarPlot_MF.pdf", g2_MF_freezing_duration_group_scaleT_female, h=3, w=5)

g3_MF_freezing_duration_group_scaleT_female <- dotplot(go_enrich_MF_freezing_duration_group_scaleT_female, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/DotPlot_MF.pdf", g3_MF_freezing_duration_group_scaleT_female, h=4, w=6)

g4_MF_freezing_duration_group_scaleT_female <- emapplot(go_enrich_MF_freezing_duration_group_scaleT_female %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/EncrichmentMap_MF.pdf", g4_MF_freezing_duration_group_scaleT_female, h=6, w=8)

g5_MF_freezing_duration_group_scaleT_female <- goplot(go_enrich_MF_freezing_duration_group_scaleT_female, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/EncrichedGOGraph_MF.pdf", g5_MF_freezing_duration_group_scaleT_female, h=6, w=8)

g6_MF_freezing_duration_group_scaleT_female <- cnetplot(go_enrich_MF_freezing_duration_group_scaleT_female, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_female/CategoryNet_MF.pdf", g6_MF_freezing_duration_group_scaleT_female, h=12, w=10)


#### freezing_duration_group_scaleT_male ####
genelist_freezing_duration_group_scaleT_male = read.delim("../data/1_single_strain/gwas/tophits/tophits_freezing_duration_group_scaleT_male_geneid.txt", header=F)
genelist_freezing_duration_group_scaleT_male <- as.vector(genelist_freezing_duration_group_scaleT_male[,1])
dir.create("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male", recursive = T)

##### BP #####

go_enrich_BP_freezing_duration_group_scaleT_male <- enrichGO(gene = genelist_freezing_duration_group_scaleT_male,
                                                               universe = genelist_background,
                                                               OrgDb = organism, 
                                                               keyType = 'FLYBASE',
                                                               readable = T,
                                                               ont = "BP",
                                                               pvalueCutoff = 0.05, 
                                                               qvalueCutoff = 0.10)
head(go_enrich_BP_freezing_duration_group_scaleT_male)

df_enrich_BP_freezing_duration_group_scaleT_male <- go_enrich_BP_freezing_duration_group_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_enrich_BP_freezing_duration_group_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/freezing_duration_group_scaleT_male_BP.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_BP_freezing_duration_group_scaleT_male <- upsetplot(go_enrich_BP_freezing_duration_group_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/UpsetPlot_BP.pdf", g1_BP_freezing_duration_group_scaleT_male, h=6, w=10)

g2_BP_freezing_duration_group_scaleT_male <- barplot(go_enrich_BP_freezing_duration_group_scaleT_male, 
                                                       drop = TRUE, 
                                                       showCategory = 10, 
                                                       title = "GO Biological Process",
                                                       font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/BarPlot_BP.pdf", g2_BP_freezing_duration_group_scaleT_male, h=3, w=5)

g3_BP_freezing_duration_group_scaleT_male <- dotplot(go_enrich_BP_freezing_duration_group_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/DotPlot_BP.pdf", g3_BP_freezing_duration_group_scaleT_male, h=4, w=6)

g4_BP_freezing_duration_group_scaleT_male <- emapplot(go_enrich_BP_freezing_duration_group_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/EncrichmentMap_BP.pdf", g4_BP_freezing_duration_group_scaleT_male, h=6, w=8)

g5_BP_freezing_duration_group_scaleT_male <- goplot(go_enrich_BP_freezing_duration_group_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/EncrichedGOGraph_BP.pdf", g5_BP_freezing_duration_group_scaleT_male, h=6, w=8)

g6_BP_freezing_duration_group_scaleT_male <- cnetplot(go_enrich_BP_freezing_duration_group_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/CategoryNet_BP.pdf", g6_BP_freezing_duration_group_scaleT_male, h=12, w=10)

##### CC #####

go_enrich_CC_freezing_duration_group_scaleT_male <- enrichGO(gene = genelist_freezing_duration_group_scaleT_male,
                                                               universe = genelist_background,
                                                               OrgDb = organism, 
                                                               keyType = 'FLYBASE',
                                                               readable = T,
                                                               ont = "CC",
                                                               pvalueCutoff = 0.05, 
                                                               qvalueCutoff = 0.10)
head(go_enrich_CC_freezing_duration_group_scaleT_male)

df_enrich_CC_freezing_duration_group_scaleT_male <- go_enrich_CC_freezing_duration_group_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_enrich_CC_freezing_duration_group_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/freezing_duration_group_scaleT_male_CC.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_freezing_duration_group_scaleT_male <- upsetplot(go_enrich_CC_freezing_duration_group_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/UpsetPlot_CC.pdf", g1_CC_freezing_duration_group_scaleT_male, h=6, w=10)

g2_CC_freezing_duration_group_scaleT_male <- barplot(go_enrich_CC_freezing_duration_group_scaleT_male, 
                                                       drop = TRUE, 
                                                       showCategory = 10, 
                                                       title = "GO Cellular Component",
                                                       font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/BarPlot_CC.pdf", g2_CC_freezing_duration_group_scaleT_male, h=3, w=5)

g3_CC_freezing_duration_group_scaleT_male <- dotplot(go_enrich_CC_freezing_duration_group_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/DotPlot_CC.pdf", g3_CC_freezing_duration_group_scaleT_male, h=4, w=6)

g4_CC_freezing_duration_group_scaleT_male <- emapplot(go_enrich_CC_freezing_duration_group_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/EncrichmentMap_CC.pdf", g4_CC_freezing_duration_group_scaleT_male, h=6, w=8)

g5_CC_freezing_duration_group_scaleT_male <- goplot(go_enrich_CC_freezing_duration_group_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/EncrichedGOGraph_CC.pdf", g5_CC_freezing_duration_group_scaleT_male, h=6, w=8)

g6_CC_freezing_duration_group_scaleT_male <- cnetplot(go_enrich_CC_freezing_duration_group_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/CategoryNet_CC.pdf", g6_CC_freezing_duration_group_scaleT_male, h=12, w=10)

##### MF #####

go_enrich_MF_freezing_duration_group_scaleT_male <- enrichGO(gene = genelist_freezing_duration_group_scaleT_male,
                                                               universe = genelist_background,
                                                               OrgDb = organism, 
                                                               keyType = 'FLYBASE',
                                                               readable = T,
                                                               ont = "MF",
                                                               pvalueCutoff = 0.05, 
                                                               qvalueCutoff = 0.10)
head(go_enrich_MF_freezing_duration_group_scaleT_male)

df_enrich_MF_freezing_duration_group_scaleT_male <- go_enrich_MF_freezing_duration_group_scaleT_male@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_enrich_MF_freezing_duration_group_scaleT_male, "../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/freezing_duration_group_scaleT_male_MF.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_freezing_duration_group_scaleT_male <- upsetplot(go_enrich_MF_freezing_duration_group_scaleT_male)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/UpsetPlot_MF.pdf", g1_MF_freezing_duration_group_scaleT_male, h=6, w=10)

g2_MF_freezing_duration_group_scaleT_male <- barplot(go_enrich_MF_freezing_duration_group_scaleT_male, 
                                                       drop = TRUE, 
                                                       showCategory = 10, 
                                                       title = "GO Molecular Function",
                                                       font.size = 8)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/BarPlot_MF.pdf", g2_MF_freezing_duration_group_scaleT_male, h=3, w=5)

g3_MF_freezing_duration_group_scaleT_male <- dotplot(go_enrich_MF_freezing_duration_group_scaleT_male, orderBy = "x")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/DotPlot_MF.pdf", g3_MF_freezing_duration_group_scaleT_male, h=4, w=6)

g4_MF_freezing_duration_group_scaleT_male <- emapplot(go_enrich_MF_freezing_duration_group_scaleT_male %>% pairwise_termsim())
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/EncrichmentMap_MF.pdf", g4_MF_freezing_duration_group_scaleT_male, h=6, w=8)

g5_MF_freezing_duration_group_scaleT_male <- goplot(go_enrich_MF_freezing_duration_group_scaleT_male, showCategory = 10)
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/EncrichedGOGraph_MF.pdf", g5_MF_freezing_duration_group_scaleT_male, h=6, w=8)

g6_MF_freezing_duration_group_scaleT_male <- cnetplot(go_enrich_MF_freezing_duration_group_scaleT_male, categorySize="pvalue")
ggplot2::ggsave("../data/1_single_strain/gwas/enrichment/clusterProfiler_freezing_duration_group_scaleT_male/CategoryNet_MF.pdf", g6_MF_freezing_duration_group_scaleT_male, h=12, w=10)




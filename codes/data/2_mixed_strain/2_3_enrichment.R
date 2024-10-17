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
library(ggplot2)
library(tidytext)

organism = "org.Dm.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#### load dataset ####
genelist_background <- read.delim("../data/genome/dgrp2_2pop_SNPlist_geneid.txt", h=F)
genelist_background <- as.vector(genelist_background[,1])

#### 1kbp_foraging_success5_2000_regress_pi_pca_dist_plus ####
genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus = read.delim("../data/2_mixed_strain/gwas/result/tophits/tophits_1kbp_foraging_success5_regress_pi_pca_dist_plus_geneid_2000.txt", header=F)
genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- as.vector(genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus[,1])
dir.create("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus", recursive = T)

##### BP #####

go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- enrichGO(gene = genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus,
                                                                         universe = genelist_background,
                                                                         OrgDb = organism, 
                                                                         keyType = 'FLYBASE',
                                                                         readable = T,
                                                                         ont = "BP",
                                                                         pvalueCutoff = 0.05, 
                                                                         qvalueCutoff = 0.10)
head(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus)
df_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, "../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/df_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus.tsv", quote = FALSE, sep = "\t", row.names = FALSE)


g1_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- upsetplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/UpsetPlot_BP.pdf", g1_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=6, w=10)

g2_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- barplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, 
                                                                 drop = TRUE, 
                                                                 showCategory = 10, 
                                                                 title = "GO Biological Process",
                                                                 font.size = 8)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/BarPlot_BP.pdf", g2_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=3, w=5)

g3_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- dotplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, orderBy = "x")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/DotPlot_BP.pdf", g3_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=4, w=6)

g4_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- emapplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus %>% pairwise_termsim())
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/EncrichmentMap_BP.pdf", g4_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=6, w=8)

g5_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- goplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, showCategory = 10)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/EncrichedGOGraph_BP.pdf", g5_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=6, w=8)

g6_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- cnetplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, categorySize="pvalue")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/CategoryNet_BP.pdf", g6_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=12, w=10)

##### CC #####

go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- enrichGO(gene = genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus,
                                                                         universe = genelist_background,
                                                                         OrgDb = organism, 
                                                                         keyType = 'FLYBASE',
                                                                         readable = T,
                                                                         ont = "CC",
                                                                         pvalueCutoff = 0.05, 
                                                                         qvalueCutoff = 0.10)
head(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus)
df_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, "../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/df_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- upsetplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/UpsetPlot_CC.pdf", g1_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=6, w=10)

g2_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- barplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, 
                                                                 drop = TRUE, 
                                                                 showCategory = 10, 
                                                                 title = "GO Cellular Component",
                                                                 font.size = 8)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/BarPlot_CC.pdf", g2_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=3, w=5)

g3_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- dotplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, orderBy = "x")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/DotPlot_CC.pdf", g3_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=4, w=6)

g4_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- emapplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus %>% pairwise_termsim())
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/EncrichmentMap_CC.pdf", g4_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=6, w=8)

g5_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- goplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, showCategory = 10)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/EncrichedGOGraph_CC.pdf", g5_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=6, w=8)

g6_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- cnetplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, categorySize="pvalue")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/CategoryNet_CC.pdf", g6_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=12, w=10)

##### MF #####

go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- enrichGO(gene = genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus,
                                                                         universe = genelist_background,
                                                                         OrgDb = organism, 
                                                                         keyType = 'FLYBASE',
                                                                         readable = T,
                                                                         ont = "MF",
                                                                         pvalueCutoff = 0.05, 
                                                                         qvalueCutoff = 0.10)
head(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus)
df_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, "../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/df_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- upsetplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/UpsetPlot_MF.pdf", g1_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=6, w=10)

g2_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- barplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, 
                                                                 drop = TRUE, 
                                                                 showCategory = 10, 
                                                                 title = "GO Molecular Function",
                                                                 font.size = 8)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/BarPlot_MF.pdf", g2_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=3, w=5)

g3_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- dotplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, orderBy = "x")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/DotPlot_MF.pdf", g3_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=4, w=6)

g4_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- emapplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus %>% pairwise_termsim())
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/EncrichmentMap_MF.pdf", g4_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=6, w=8)

g5_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- goplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, showCategory = 10)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/EncrichedGOGraph_MF.pdf", g5_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=6, w=8)

g6_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus <- cnetplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, categorySize="pvalue")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus/CategoryNet_MF.pdf", g6_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_plus, h=12, w=10)


#### 1kbp_foraging_success5_2000_regress_pi_pca_dist_minus ####
genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus = read.delim("../data/2_mixed_strain/gwas/result/tophits/tophits_1kbp_foraging_success5_regress_pi_pca_dist_minus_geneid_2000.txt", header=F)
genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- as.vector(genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus[,1])
dir.create("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus", recursive = T)

##### BP #####

go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- enrichGO(gene = genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus,
                                                                          universe = genelist_background,
                                                                          OrgDb = organism, 
                                                                          keyType = 'FLYBASE',
                                                                          readable = T,
                                                                          ont = "BP",
                                                                          pvalueCutoff = 0.05, 
                                                                          qvalueCutoff = 0.10)
head(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus)
df_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "BP")
write.table(df_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, "../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/df_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus.tsv", quote = FALSE, sep = "\t", row.names = FALSE)


g1_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- upsetplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/UpsetPlot_BP.pdf", g1_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=6, w=10)

g2_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- barplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, 
                                                                  drop = TRUE, 
                                                                  showCategory = 10, 
                                                                  title = "GO Biological Process",
                                                                  font.size = 8)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/BarPlot_BP.pdf", g2_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=3, w=5)

g3_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- dotplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, orderBy = "x")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/DotPlot_BP.pdf", g3_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=4, w=6)

g4_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- emapplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus %>% pairwise_termsim())
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/EncrichmentMap_BP.pdf", g4_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=6, w=8)

g5_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- goplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, showCategory = 10)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/EncrichedGOGraph_BP.pdf", g5_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=6, w=8)

g6_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- cnetplot(go_enrich_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, categorySize="pvalue")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/CategoryNet_BP.pdf", g6_BP_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=12, w=10)

##### CC #####

go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- enrichGO(gene = genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus,
                                                                          universe = genelist_background,
                                                                          OrgDb = organism, 
                                                                          keyType = 'FLYBASE',
                                                                          readable = T,
                                                                          ont = "CC",
                                                                          pvalueCutoff = 0.05, 
                                                                          qvalueCutoff = 0.10)
head(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus)
df_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "CC")
write.table(df_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, "../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/df_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- upsetplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/UpsetPlot_CC.pdf", g1_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=6, w=10)

g2_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- barplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, 
                                                                  drop = TRUE, 
                                                                  showCategory = 10, 
                                                                  title = "GO Cellular Component",
                                                                  font.size = 8)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/BarPlot_CC.pdf", g2_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=3, w=5)

g3_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- dotplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, orderBy = "x")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/DotPlot_CC.pdf", g3_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=4, w=6)

g4_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- emapplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus %>% pairwise_termsim())
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/EncrichmentMap_CC.pdf", g4_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=6, w=8)

g5_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- goplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, showCategory = 10)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/EncrichedGOGraph_CC.pdf", g5_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=6, w=8)

g6_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- cnetplot(go_enrich_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, categorySize="pvalue")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/CategoryNet_CC.pdf", g6_CC_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=12, w=10)

##### MF #####

go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- enrichGO(gene = genelist_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus,
                                                                          universe = genelist_background,
                                                                          OrgDb = organism, 
                                                                          keyType = 'FLYBASE',
                                                                          readable = T,
                                                                          ont = "MF",
                                                                          pvalueCutoff = 0.05, 
                                                                          qvalueCutoff = 0.10)
head(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus)
df_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus@result %>%
  filter(p.adjust < 0.05) %>%
  mutate(geneID = str_replace_all(geneID, "/", ", "),
         Type = "MF")
write.table(df_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, "../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/df_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

g1_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- upsetplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/UpsetPlot_MF.pdf", g1_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=6, w=10)

g2_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- barplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, 
                                                                  drop = TRUE, 
                                                                  showCategory = 10, 
                                                                  title = "GO Molecular Function",
                                                                  font.size = 8)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/BarPlot_MF.pdf", g2_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=3, w=5)

g3_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- dotplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, orderBy = "x")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/DotPlot_MF.pdf", g3_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=4, w=6)

g4_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- emapplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus %>% pairwise_termsim())
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/EncrichmentMap_MF.pdf", g4_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=6, w=8)

g5_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- goplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, showCategory = 10)
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/EncrichedGOGraph_MF.pdf", g5_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=6, w=8)

g6_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus <- cnetplot(go_enrich_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, categorySize="pvalue")
ggplot2::ggsave("../data/2_mixed_strain/gwas/result/enrichment/clusterProfiler_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus/CategoryNet_MF.pdf", g6_MF_1kbp_foraging_success5_2000_regress_pi_pca_dist_minus, h=12, w=10)





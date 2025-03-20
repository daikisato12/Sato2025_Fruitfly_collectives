#### load packages ####
targetPackages <- c('tidyverse','arrow','normentR','biomaRt')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)
# devtools::install_github("norment/normentR")
# BiocManager::install("biomaRt")

#### Figure S18a ####
##### load dataset #####
df_ghas <- read_parquet("../data/2_mixed_strain/gwas/result/rawdata/df_gwas_1kbp_performanceDE_3sec_regress_pi_grm.parquet") %>%
  ungroup()

axisdf_ghas <- df_ghas %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


##### make plot #####
nwindows <- nrow(df_ghas)
ci <- 0.95
qqplot_df_ghas <- ggplot(data.frame(
  observed = -log10(sort(df_ghas$P)),
  expected = -log10(ppoints(nwindows)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nwindows), shape2 = rev(seq(nwindows)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nwindows), shape2 = rev(seq(nwindows))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme(axis.text = element_text(color = "black"))
qqplot_df_ghas

ggsave("../figures/FigureS18a.pdf", qqplot_df_ghas, w=3, h=3)

#### Figure S18b ####
list_gwas <- read.delim("../data/1_single_strain/gwas/result/tophits/tophits_motion_cue_exit_intercept_scaleT_female_geneid.txt",
                                             header = FALSE) %>%
  pull(V1)
list_ghas <- read.delim("../data/2_mixed_strain/gwas/result/tophits/tophits_1kbp_performanceDE_3sec_regress_pi_grm_plus_geneid_2000.txt",
                                             header = FALSE) %>%
  pull(V1)


gene_shared <- intersect(list_gwas, list_ghas)

mart <- biomaRt::useMart(biomart = "ensembl",  #useEnsembl or useMart
                         host = "https://www.ensembl.org", 
                         dataset = "dmelanogaster_gene_ensembl")
gene_shared_name <- biomaRt::getBM(attributes = c("external_gene_name", 
                                                  "ensembl_gene_id"),
                           filters = "ensembl_gene_id",
                           values = gene_shared,
                           mart = mart) %>%
  dplyr::rename(gene_name = external_gene_name)

gene_shared_name

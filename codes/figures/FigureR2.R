#### load packages ####
targetPackages <- c('tidyverse','biomaRt')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

#### Figure R2 (Appendix Figure 2 in Peer Review File) ####
##### load dataset #####
df_genes_female <- read.table("../data/1_single_strain/gwas/result/tophits/tophits_motion_cue_exit_intercept_scaleT_female_geneid.txt") %>%
  pull(V1)
df_genes_male <- read.table("../data/1_single_strain/gwas/result/tophits/tophits_motion_cue_exit_intercept_scaleT_male_geneid.txt") %>%
  pull(V1)
df_genes_all <- c(df_genes_female, df_genes_male) %>%
  unique()
# biomaRt::listAttributes(mart)
mart <- biomaRt::useMart(biomart = "ensembl",
                         host = "https://www.ensembl.org", 
                         dataset = "dmelanogaster_gene_ensembl")
gene_name <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                  "external_gene_name"),
                   filters = "ensembl_gene_id",
                   values = df_genes_all,
                   mart = mart) #%>%
  dplyr::rename(gene_name = external_gene_name)
write.table(gene_name$external_gene_name,
            "../data/1_single_strain/gwas/result/tophits/tophits_motion_cue_exit_intercept_scaleT_bothsex_genename.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

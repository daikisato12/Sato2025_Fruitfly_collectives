#### load libraries ####
library(tidyverse)
library(plyr)
library(arrow)
library(patchwork)
# library(devtools)
# devtools::install_github("norment/normentR")
# library(normentR)

#### load dataset ####
##### pi #####
df_pi <- read.table("../data/0_genome/dgrp2_1kbp.windowed.pi", header = T) %>%
  dplyr::mutate(BP = BIN_END - 500)

df_TajimaD <- read.table("../data/0_genome/dgrp2_1kbp.Tajima.D", header = T) %>%
  dplyr::mutate(BP = BIN_START + 500)

##### GWAS loci frezing_duration #####
df_gwas_res_frezing_duration <- 
  read.table("../data/1_single_strain/gwas/result/scaled/rawdata/freezing_duration_group_female.fastGWA", header = TRUE) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  dplyr::mutate(BP = as.numeric(POS)) %>%
  filter(P < 1/(10**4)) %>%
  dplyr::rename(CHROM = CHR) %>%
  dplyr::mutate(BP = as.integer(BP/1000)*1000+500) %>%
  full_join(df_pi) %>%
  full_join(df_TajimaD %>%
              dplyr::select(CHROM, BP, TajimaD)) %>%
  dplyr::mutate(data = case_when(is.na(P) ~ "Not-associated",
                                 TRUE ~ "Associated"))

##### GWAS loci motioncue_intercept #####
df_gwas_res_motioncue_intercept <- 
  # read.table("../data/1_single_strain/gwas/result/scaled/tophits/tophits_motion_cue_exit_intercept_scaleT_female.txt", header = TRUE) %>%
  read.table("../data/1_single_strain/gwas/result/scaled/rawdata/motion_cue_exit_intercept_female.fastGWA", header = TRUE) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  dplyr::mutate(BP = as.numeric(POS)) %>%
  filter(P < 1/(10**4)) %>%
  dplyr::rename(CHROM = CHR) %>%
  dplyr::mutate(BP = as.integer(BP/1000)*1000+500) %>%
  full_join(df_pi) %>%
  full_join(df_TajimaD %>%
              dplyr::select(CHROM, BP, TajimaD)) %>%
  dplyr::mutate(data = case_when(is.na(P) ~ "Not-associated",
                                 TRUE ~ "Associated"))

##### DE effect loci #####
df_pi_res_diversity_plus <- read.table("../data/2_mixed_strain/gwas/tophits/tophits_1kbp_foraging_success5_regress_pi_pca_dist_plus.txt", header = TRUE) %>%
  dplyr::rename(CHROM = CHR) %>%
  right_join(df_pi) %>%
  right_join(df_TajimaD %>%
               dplyr::select(CHROM, BP, TajimaD)) %>%
  dplyr::mutate(data = case_when(is.na(P) ~ "others",
                                 TRUE ~ "DE involved"))

df_pi_res_diversity_minus <- read.table("../data/2_mixed_strain/gwas/tophits/tophits_1kbp_foraging_success5_regress_pi_pca_dist_minus.txt", header = TRUE) %>%
  dplyr::rename(CHROM = CHR) %>%
  right_join(df_pi) %>%
  right_join(df_TajimaD %>%
               dplyr::select(CHROM, BP, TajimaD)) %>%
  dplyr::mutate(data = case_when(is.na(P) ~ "others",
                                 TRUE ~ "DE involved"))


#### analysis ####
##### GWAS permutation test motion cue_intercept #####
numloci_gwas_motioncue_intercept <- df_gwas_res_motioncue_intercept %>%
  filter(data == "Associated") %>%
  nrow()
df_gwas_test_motioncue_intercept <- data.frame()
for (i in 1:100){
  set.seed(i)
  df_tmp <- 
    df_gwas_res_motioncue_intercept %>%
    filter(data == "Not-associated") %>%
    sample_n(size = as.integer(numloci_gwas_motioncue_intercept/100)*100) %>%
    bind_rows(df_gwas_res_motioncue_intercept %>%
                filter(data == "Associated") %>%
                sample_n(size = as.integer(numloci_gwas_motioncue_intercept/100)*100)) %>%
    group_by(data) %>%
    dplyr::summarize(PI = mean(PI, na.rm = TRUE),
                     TajimaD = mean(TajimaD, na.rm = TRUE)) %>%
    mutate(rep = i)
  
  df_gwas_test_motioncue_intercept <- df_gwas_test_motioncue_intercept %>%
    bind_rows(df_tmp)
}
write.table(df_gwas_test_motioncue_intercept, 
            "../data/1_single_strain/gwas/selection/df_gwas_test.tsv",
            row.names = F, quote = F, sep = "\t")


##### DE permutation test #####
numloci_plus <- df_pi_res_diversity_plus %>%
  filter(data == "DE involved") %>%
  nrow()
numloci_minus <- df_pi_res_diversity_minus %>%
  filter(data == "DE involved") %>%
  nrow()
df_DE_test <- data.frame()
for (i in 1:100){
  set.seed(i)
  df_tmp_plus <- 
    df_pi_res_diversity_plus %>%
    filter(data == "others") %>%
    sample_n(size = as.integer(numloci_plus/100)*100) %>%
    bind_rows(df_pi_res_diversity_plus %>%
                filter(data == "DE involved") %>%
                sample_n(size = as.integer(numloci_plus/100)*100)) %>%
    group_by(data) %>%
    dplyr::summarize(PI = mean(PI, na.rm = TRUE),
                     TajimaD = mean(TajimaD, na.rm = TRUE)) %>%
    mutate(rep = i)
  
  df_tmp_minus <- 
    df_pi_res_diversity_minus %>%
    filter(data == "others") %>%
    sample_n(size = as.integer(numloci_minus/100)*100) %>%
    bind_rows(df_pi_res_diversity_minus %>%
                filter(data == "DE involved") %>%
                sample_n(size = as.integer(numloci_minus/100)*100)) %>%
    group_by(data) %>%
    dplyr::summarize(PI = mean(PI, na.rm = TRUE),
                     TajimaD = mean(TajimaD, na.rm = TRUE)) %>%
    mutate(rep = i)
  
  df_DE_test <- df_DE_test %>%
    bind_rows(df_tmp_plus %>%
                dplyr::mutate(cor = "plus")) %>%
    bind_rows(df_tmp_minus %>%
                dplyr::mutate(cor = "minus"))
}
write.table(df_DE_test, "../data/2_mixed_strain/gwas/selection/df_DE_test.tsv",
            row.names = F, quote = F, sep = "\t")


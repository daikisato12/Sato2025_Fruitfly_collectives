#### load libraries ####
library(tidyverse)
library(arrow)
# devtools::install_github("norment/normentR")
library(normentR)
ci <- 0.95

alpha=1/(10**5)#0.05/1877723

#### GWAS results ####
##### GCTA GWA motion cue #####
###### gwas_motion_cue_exit_intercept_scaleT_female ######
gwas_motion_cue_exit_intercept_scaleT_female <- read.table("${dataset}/1_single_strain/gwas/motion_cue_exit_intercept_female.fastGWA",h=T) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(POS))

df_gwas_motion_cue_exit_intercept_scaleT_female <- gwas_motion_cue_exit_intercept_scaleT_female %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_motion_cue_exit_intercept_scaleT_female, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_motion_cue_exit_intercept_scaleT_female, "../data/_1_single_strain/gwas/df/df_gwas_motion_cue_exit_intercept_scaleT_female.parquet")
# df_gwas_motion_cue_exit_intercept_scaleT_female <- read_parquet("../data/_1_single_strain/gwas/df/df_gwas_motion_cue_exit_intercept_scaleT_female.parquet") %>%
#   ungroup()


df_gwas_motion_cue_exit_intercept_scaleT_female_top <- df_gwas_motion_cue_exit_intercept_scaleT_female %>%
  mutate(BETA_Z = BETA,
         SE_Z = SE) %>%
  as.data.frame() %>%
  mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_motion_cue_exit_intercept_scaleT_female_top %>%
              dplyr::select(c(POS, CHR, BP), everything()),
            "../data/_1_single_strain/gwas/tophits/tophits_motion_cue_exit_intercept_scaleT_female.txt", delim = "\t")


###### gwas_motion_cue_exit_intercept_scaleT_male ######
gwas_motion_cue_exit_intercept_scaleT_male <- read.table("${dataset}/1_single_strain/gwas/motion_cue_exit_intercept_male.fastGWA",h=T) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(POS))

df_gwas_motion_cue_exit_intercept_scaleT_male <- gwas_motion_cue_exit_intercept_scaleT_male %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_motion_cue_exit_intercept_scaleT_male, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_motion_cue_exit_intercept_scaleT_male, "../data/_1_single_strain/gwas/df/df_gwas_motion_cue_exit_intercept_scaleT_male.parquet")
# df_gwas_motion_cue_exit_intercept_scaleT_male <- read_parquet("../data/_1_single_strain/gwas/df/df_gwas_motion_cue_exit_intercept_scaleT_male.parquet") %>%
#   ungroup()


df_gwas_motion_cue_exit_intercept_scaleT_male_top <- df_gwas_motion_cue_exit_intercept_scaleT_male %>%
  mutate(BETA_Z = BETA,
         SE_Z = SE) %>%
  as.data.frame() %>%
  mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_motion_cue_exit_intercept_scaleT_male_top %>%
              dplyr::select(c(POS, CHR, BP), everything()),
            "../data/_1_single_strain/gwas/tophits/tophits_motion_cue_exit_intercept_scaleT_male.txt", delim = "\t")


##### GCTA GWA nnd #####
###### gwas_nnd_scaleT_female ######
gwas_nnd_scaleT_female <- read.table("${dataset}/1_single_strain/gwas/nnd_female.fastGWA",h=T) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(POS))

df_gwas_nnd_scaleT_female <- gwas_nnd_scaleT_female %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_nnd_scaleT_female, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_nnd_scaleT_female, "../data/_1_single_strain/gwas/df/df_gwas_nnd_scaleT_female.parquet")
# df_gwas_nnd_scaleT_female <- read_parquet("../data/_1_single_strain/gwas/df/df_gwas_nnd_scaleT_female.parquet") %>%
#   ungroup()


df_gwas_nnd_scaleT_female_top <- df_gwas_nnd_scaleT_female %>%
  mutate(BETA_Z = BETA,
         SE_Z = SE) %>%
  as.data.frame() %>%
  mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_nnd_scaleT_female_top %>%
              dplyr::select(c(POS, CHR, BP), everything()),
            "../data/_1_single_strain/gwas/tophits/tophits_nnd_scaleT_female.txt", delim = "\t")


###### gwas_nnd_scaleT_male ######
gwas_nnd_scaleT_male <- read.table("${dataset}/1_single_strain/gwas/nnd_male.fastGWA",h=T) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(POS))

df_gwas_nnd_scaleT_male <- gwas_nnd_scaleT_male %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_nnd_scaleT_male, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_nnd_scaleT_male, "../data/_1_single_strain/gwas/df/df_gwas_nnd_scaleT_male.parquet")
# df_gwas_nnd_scaleT_male <- read_parquet("../data/_1_single_strain/gwas/df/df_gwas_nnd_scaleT_male.parquet") %>%
#   ungroup()


df_gwas_nnd_scaleT_male_top <- df_gwas_nnd_scaleT_male %>%
  mutate(BETA_Z = BETA,
         SE_Z = SE) %>%
  as.data.frame() %>%
  mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_nnd_scaleT_male_top %>%
              dplyr::select(c(POS, CHR, BP), everything()),
            "../data/_1_single_strain/gwas/tophits/tophits_nnd_scaleT_male.txt", delim = "\t")


##### GCTA GWA freezing_duration_group #####
###### gwas_freezing_duration_group_scaleT_female ######
gwas_freezing_duration_group_scaleT_female <- read.table("${dataset}/1_single_strain/gwas/freezing_duration_group_female.fastGWA",h=T) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(POS))

df_gwas_freezing_duration_group_scaleT_female <- gwas_freezing_duration_group_scaleT_female %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_freezing_duration_group_scaleT_female, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_freezing_duration_group_scaleT_female, "../data/_1_single_strain/gwas/df/df_gwas_freezing_duration_group_scaleT_female.parquet")
# df_gwas_freezing_duration_group_scaleT_female <- read_parquet("../data/_1_single_strain/gwas/df/df_gwas_freezing_duration_group_scaleT_female.parquet") %>%
#   ungroup()


df_gwas_freezing_duration_group_scaleT_female_top <- df_gwas_freezing_duration_group_scaleT_female %>%
  mutate(BETA_Z = BETA,
         SE_Z = SE) %>%
  as.data.frame() %>%
  mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_freezing_duration_group_scaleT_female_top %>%
              dplyr::select(c(POS, CHR, BP), everything()),
            "../data/_1_single_strain/gwas/tophits/tophits_freezing_duration_group_scaleT_female.txt", delim = "\t")


###### gwas_freezing_duration_group_scaleT_male ######
gwas_freezing_duration_group_scaleT_male <- read.table("${dataset}/1_single_strain/gwas/freezing_duration_group_male.fastGWA",h=T) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(POS))

df_gwas_freezing_duration_group_scaleT_male <- gwas_freezing_duration_group_scaleT_male %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_freezing_duration_group_scaleT_male, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_freezing_duration_group_scaleT_male, "../data/_1_single_strain/gwas/df/df_gwas_freezing_duration_group_scaleT_male.parquet")
# df_gwas_freezing_duration_group_scaleT_male <- read_parquet("../data/_1_single_strain/gwas/df/df_gwas_freezing_duration_group_scaleT_male.parquet") %>%
#   ungroup()


df_gwas_freezing_duration_group_scaleT_male_top <- df_gwas_freezing_duration_group_scaleT_male %>%
  mutate(BETA_Z = BETA,
         SE_Z = SE) %>%
  as.data.frame() %>%
  mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_freezing_duration_group_scaleT_male_top %>%
              dplyr::select(c(POS, CHR, BP), everything()),
            "../data/_1_single_strain/gwas/tophits/tophits_freezing_duration_group_scaleT_male.txt", delim = "\t")


##### GCTA GWA freezing_duration_single #####
###### gwas_freezing_duration_single_scaleT_female ######
gwas_freezing_duration_single_scaleT_female <- read.table("${dataset}/1_single_strain/gwas/freezing_duration_single_female.fastGWA",h=T) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(POS))

df_gwas_freezing_duration_single_scaleT_female <- gwas_freezing_duration_single_scaleT_female %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_freezing_duration_single_scaleT_female, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_freezing_duration_single_scaleT_female, "../data/_1_single_strain/gwas/df/df_gwas_freezing_duration_single_scaleT_female.parquet")
# df_gwas_freezing_duration_single_scaleT_female <- read_parquet("../data/_1_single_strain/gwas/df/df_gwas_freezing_duration_single_scaleT_female.parquet") %>%
#   ungroup()


df_gwas_freezing_duration_single_scaleT_female_top <- df_gwas_freezing_duration_single_scaleT_female %>%
  mutate(BETA_Z = BETA,
         SE_Z = SE) %>%
  as.data.frame() %>%
  mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_freezing_duration_single_scaleT_female_top %>%
              dplyr::select(c(POS, CHR, BP), everything()),
            "../data/_1_single_strain/gwas/tophits/tophits_freezing_duration_single_scaleT_female.txt", delim = "\t")



###### gwas_freezing_duration_single_scaleT_male ######
gwas_freezing_duration_single_scaleT_male <- read.table("${dataset}/1_single_strain/gwas/freezing_duration_single_male.fastGWA",h=T) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(POS))

df_gwas_freezing_duration_single_scaleT_male <- gwas_freezing_duration_single_scaleT_male %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_freezing_duration_single_scaleT_male, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_freezing_duration_single_scaleT_male, "../data/_1_single_strain/gwas/df/df_gwas_freezing_duration_single_scaleT_male.parquet")
# df_gwas_freezing_duration_single_scaleT_male <- read_parquet("../data/_1_single_strain/gwas/df/df_gwas_freezing_duration_single_scaleT_male.parquet") %>%
#   ungroup()


df_gwas_freezing_duration_single_scaleT_male_top <- df_gwas_freezing_duration_single_scaleT_male %>%
  mutate(BETA_Z = BETA,
         SE_Z = SE) %>%
  as.data.frame() %>%
  mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_freezing_duration_single_scaleT_male_top %>%
              dplyr::select(c(POS, CHR, BP), everything()),
            "../data/_1_single_strain/gwas/tophits/tophits_freezing_duration_single_scaleT_male.txt", delim = "\t")



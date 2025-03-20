#### load packages ####
targetPackages <- c('tidyverse','arrow','normentR')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

# library(devtools)
# devtools::install_github("norment/normentR", force = TRUE)

ci <- 0.95
alpha=1/(10**5)#0.05/1877723

#### GWAS results ####
##### GCTA GWA freezing_duration_single #####
###### gwas_freezing_duration_single_scaleT_female ######
gwas_freezing_duration_single_scaleT_female <- read.table("../data/1_single_strain/gwas/result/rawdata/freezing_duration_single_female.fastGWA",h=T) %>%
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

write_parquet(df_gwas_freezing_duration_single_scaleT_female, "../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_single_scaleT_female.parquet")
# df_gwas_freezing_duration_single_scaleT_female <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_single_scaleT_female.parquet") %>%
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
            "../data/1_single_strain/gwas/result/tophits/tophits_freezing_duration_single_scaleT_female.txt", delim = "\t")



###### gwas_freezing_duration_single_scaleT_male ######
gwas_freezing_duration_single_scaleT_male <- read.table("../data/1_single_strain/gwas/result/rawdata/freezing_duration_single_male.fastGWA",h=T) %>%
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

write_parquet(df_gwas_freezing_duration_single_scaleT_male, "../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_single_scaleT_male.parquet")
# df_gwas_freezing_duration_single_scaleT_male <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_single_scaleT_male.parquet") %>%
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
            "../data/1_single_strain/gwas/result/tophits/tophits_freezing_duration_single_scaleT_male.txt", delim = "\t")


##### GCTA GWA freezing_duration_group #####
###### gwas_freezing_duration_group_scaleT_female ######
gwas_freezing_duration_group_scaleT_female <- read.table("../data/1_single_strain/gwas/result/rawdata/freezing_duration_group_female.fastGWA",h=T) %>%
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

write_parquet(df_gwas_freezing_duration_group_scaleT_female, "../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_group_scaleT_female.parquet")
# df_gwas_freezing_duration_group_scaleT_female <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_group_scaleT_female.parquet") %>%
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
            "../data/1_single_strain/gwas/result/tophits/tophits_freezing_duration_group_scaleT_female.txt", delim = "\t")


###### gwas_freezing_duration_group_scaleT_male ######
gwas_freezing_duration_group_scaleT_male <- read.table("../data/1_single_strain/gwas/result/rawdata/freezing_duration_group_male.fastGWA",h=T) %>%
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

write_parquet(df_gwas_freezing_duration_group_scaleT_male, "../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_group_scaleT_male.parquet")
# df_gwas_freezing_duration_group_scaleT_male <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_freezing_duration_group_scaleT_male.parquet") %>%
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
            "../data/1_single_strain/gwas/result/tophits/tophits_freezing_duration_group_scaleT_male.txt", delim = "\t")


##### GCTA GWA motion cue #####
###### gwas_motion_cue_exit_intercept_scaleT_female ######
gwas_motion_cue_exit_intercept_scaleT_female <- read.table("../data/1_single_strain/gwas/result/rawdata/motion_cue_exit_intercept_female.fastGWA",h=T) %>%
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

write_parquet(df_gwas_motion_cue_exit_intercept_scaleT_female, "../data/1_single_strain/gwas/result/df/df_gwas_motion_cue_exit_intercept_scaleT_female.parquet")
# df_gwas_motion_cue_exit_intercept_scaleT_female <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_motion_cue_exit_intercept_scaleT_female.parquet") %>%
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
            "../data/1_single_strain/gwas/result/tophits/tophits_motion_cue_exit_intercept_scaleT_female.txt", delim = "\t")


###### gwas_motion_cue_exit_intercept_scaleT_male ######
gwas_motion_cue_exit_intercept_scaleT_male <- read.table("../data/1_single_strain/gwas/result/rawdata/motion_cue_exit_intercept_male.fastGWA",h=T) %>%
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

write_parquet(df_gwas_motion_cue_exit_intercept_scaleT_male, "../data/1_single_strain/gwas/result/df/df_gwas_motion_cue_exit_intercept_scaleT_male.parquet")
# df_gwas_motion_cue_exit_intercept_scaleT_male <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_motion_cue_exit_intercept_scaleT_male.parquet") %>%
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
            "../data/1_single_strain/gwas/result/tophits/tophits_motion_cue_exit_intercept_scaleT_male.txt", delim = "\t")


###### gwas_motion_cue_stop_freezing_abs_scaleT_female ######
gwas_motion_cue_stop_freezing_abs_scaleT_female <- read.table("../data/1_single_strain/gwas/result/rawdata/motion_cue_stop_freezing_abs_female.fastGWA",h=T) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(POS))

df_gwas_motion_cue_stop_freezing_abs_scaleT_female <- gwas_motion_cue_stop_freezing_abs_scaleT_female %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_motion_cue_stop_freezing_abs_scaleT_female, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_motion_cue_stop_freezing_abs_scaleT_female, "../data/1_single_strain/gwas/result/df/df_gwas_motion_cue_stop_freezing_abs_scaleT_female.parquet")
# df_gwas_motion_cue_stop_freezing_abs_scaleT_female <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_motion_cue_stop_freezing_abs_scaleT_female.parquet") %>%
#   ungroup()


df_gwas_motion_cue_stop_freezing_abs_scaleT_female_top <- df_gwas_motion_cue_stop_freezing_abs_scaleT_female %>%
  mutate(BETA_Z = BETA,
         SE_Z = SE) %>%
  as.data.frame() %>%
  mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_motion_cue_stop_freezing_abs_scaleT_female_top %>%
              dplyr::select(c(POS, CHR, BP), everything()),
            "../data/1_single_strain/gwas/result/tophits/tophits_motion_cue_stop_freezing_abs_scaleT_female.txt", delim = "\t")


###### gwas_motion_cue_stop_freezing_abs_scaleT_male ######
gwas_motion_cue_stop_freezing_abs_scaleT_male <- read.table("../data/1_single_strain/gwas/result/rawdata/motion_cue_stop_freezing_abs_male.fastGWA",h=T) %>%
  separate(SNP, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(POS))

df_gwas_motion_cue_stop_freezing_abs_scaleT_male <- gwas_motion_cue_stop_freezing_abs_scaleT_male %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_motion_cue_stop_freezing_abs_scaleT_male, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_motion_cue_stop_freezing_abs_scaleT_male, "../data/1_single_strain/gwas/result/df/df_gwas_motion_cue_stop_freezing_abs_scaleT_male.parquet")
# df_gwas_motion_cue_stop_freezing_abs_scaleT_male <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_motion_cue_stop_freezing_abs_scaleT_male.parquet") %>%
#   ungroup()


df_gwas_motion_cue_stop_freezing_abs_scaleT_male_top <- df_gwas_motion_cue_stop_freezing_abs_scaleT_male %>%
  mutate(BETA_Z = BETA,
         SE_Z = SE) %>%
  as.data.frame() %>%
  mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_motion_cue_stop_freezing_abs_scaleT_male_top %>%
              dplyr::select(c(POS, CHR, BP), everything()),
            "../data/1_single_strain/gwas/result/tophits/tophits_motion_cue_stop_freezing_abs_scaleT_male.txt", delim = "\t")


##### GCTA GWA nnd #####
###### gwas_nnd_scaleT_female ######
gwas_nnd_scaleT_female <- read.table("../data/1_single_strain/gwas/result/rawdata/nnd_female.fastGWA",h=T) %>%
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

write_parquet(df_gwas_nnd_scaleT_female, "../data/1_single_strain/gwas/result/df/df_gwas_nnd_scaleT_female.parquet")
# df_gwas_nnd_scaleT_female <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_nnd_scaleT_female.parquet") %>%
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
            "../data/1_single_strain/gwas/result/tophits/tophits_nnd_scaleT_female.txt", delim = "\t")


###### gwas_nnd_scaleT_male ######
gwas_nnd_scaleT_male <- read.table("../data/1_single_strain/gwas/result/rawdata/nnd_male.fastGWA",h=T) %>%
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

write_parquet(df_gwas_nnd_scaleT_male, "../data/1_single_strain/gwas/result/df/df_gwas_nnd_scaleT_male.parquet")
# df_gwas_nnd_scaleT_male <- read_parquet("../data/1_single_strain/gwas/result/df/df_gwas_nnd_scaleT_male.parquet") %>%
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
            "../data/1_single_strain/gwas/result/tophits/tophits_nnd_scaleT_male.txt", delim = "\t")



#### METAL results ####
##### freezing_duration_single #####
###### gwas_metal_freezing_duration_single_scaleT ######
gwas_metal_freezing_duration_single_scaleT <- read.table("../data/1_single_strain/gwas/METAL/rawdata/metal_freezing_duration_single.tbl",h=T) %>%
  separate(MarkerName, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(BP))

# gwas_metal_freezing_duration_single_scaleT$FDR <-
#   p.adjust(gwas_metal_freezing_duration_single_scaleT$P.value, 
#            method = "fdr")

alpha <- 1/(10**6) #0.05/2672211
# alpha <- 0.05/(10**6)
df_gwas_metal_freezing_duration_single_scaleT <- gwas_metal_freezing_duration_single_scaleT %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_metal_freezing_duration_single_scaleT, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  dplyr::rename(P = P.value) %>%
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_metal_freezing_duration_single_scaleT, "../data/1_single_strain/gwas/METAL/df/df_gwas_metal_freezing_duration_single_scaleT.parquet")
# df_gwas_metal_freezing_duration_single_scaleT <- read_parquet("../data/1_single_strain/gwas/METAL/df/df_gwas_metal_freezing_duration_single_scaleT.parquet") %>%
#   ungroup()

df_gwas_metal_freezing_duration_single_scaleT_top <- df_gwas_metal_freezing_duration_single_scaleT %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_metal_freezing_duration_single_scaleT_top %>%
              dplyr::select(CHR, BP), 
            "../data/1_single_strain/gwas/METAL/tophits/tophits_freezing_duration_single_scaleT.bed", delim = " ", col_names = F)

write_delim(df_gwas_metal_freezing_duration_single_scaleT_top %>%
              dplyr::select(c(CHR, BP), everything()),
            "../data/1_single_strain/gwas/METAL/tophits/tophits_freezing_duration_single_scaleT.txt", delim = "\t")

axisdf_gwas_metal_freezing_duration_single_scaleT <- df_gwas_metal_freezing_duration_single_scaleT %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

g_gwas_metal_freezing_duration_single_scaleT <- 
  ggplot(df_gwas_metal_freezing_duration_single_scaleT, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR))) +#,
  # alpha = as.factor(is_sig),
  # size = as.factor(is_sig))) +
  
  # Show all points
  geom_point(shape = 16, alpha = 0.2) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) + #c("#c8c2be","#71686c") for single
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  #  geom_point( aes(fill=as.factor(is_col)), alpha=0.8, size=1.3) +
  scale_fill_viridis_d() +
  # scale_shape_manual(values = c(19,8)) +
  #  scale_color_manual(values = rep(c("grey", "#cd5c5c"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf_gwas_metal_freezing_duration_single_scaleT$CHR, breaks = axisdf_gwas_metal_freezing_duration_single_scaleT$center) +
  #  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 11)) +
  
  #  geom_hline(yintercept=-log10(alpha), linetype="dotted", color="#cd5c5c", size=1) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  # facet_wrap(~ sex, nrow = 2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_metal_freezing_duration_single_scaleT
ggsave("../data/1_single_strain/gwas/METAL/manhattan_plot/gwas_metal_freezing_duration_single_scaleT2.png", 
       g_gwas_metal_freezing_duration_single_scaleT, 
       width = 6, height = 2)


nSNPs <- nrow(df_gwas_metal_freezing_duration_single_scaleT)
ci <- 0.95
qqplot_df_gwas_metal_freezing_duration_single_scaleT <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_metal_freezing_duration_single_scaleT$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme()
ggsave("../data/1_single_strain/gwas/METAL/qqplot/qqplot_df_gwas_metal_freezing_duration_single_scaleT.png", 
       qqplot_df_gwas_metal_freezing_duration_single_scaleT, 
       width = 3, height = 3)


##### freezing_duration_group #####
###### gwas_metal_freezing_duration_group_scaleT ######
gwas_metal_freezing_duration_group_scaleT <- read.table("../data/1_single_strain/gwas/METAL/rawdata/metal_freezing_duration_group.tbl",h=T) %>%
  separate(MarkerName, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(BP))

# gwas_metal_freezing_duration_group_scaleT$FDR <-
#   p.adjust(gwas_metal_freezing_duration_group_scaleT$P.value, 
#            method = "fdr")

alpha <- 1/(10**6) #0.05/2672211
# alpha <- 0.05/(10**6)
df_gwas_metal_freezing_duration_group_scaleT <- gwas_metal_freezing_duration_group_scaleT %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_metal_freezing_duration_group_scaleT, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  dplyr::rename(P = P.value) %>%
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_metal_freezing_duration_group_scaleT, "../data/1_single_strain/gwas/METAL/df/df_gwas_metal_freezing_duration_group_scaleT.parquet")
# df_gwas_metal_freezing_duration_group_scaleT <- read_parquet("../data/1_single_strain/gwas/METAL/df/df_gwas_metal_freezing_duration_group_scaleT.parquet") %>%
#   ungroup()

df_gwas_metal_freezing_duration_group_scaleT_top <- df_gwas_metal_freezing_duration_group_scaleT %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_metal_freezing_duration_group_scaleT_top %>%
              dplyr::select(CHR, BP), 
            "../data/1_single_strain/gwas/METAL/tophits/tophits_freezing_duration_group_scaleT.bed", delim = " ", col_names = F)

write_delim(df_gwas_metal_freezing_duration_group_scaleT_top %>%
              dplyr::select(c(CHR, BP), everything()),
            "../data/1_single_strain/gwas/METAL/tophits/tophits_freezing_duration_group_scaleT.txt", delim = "\t")

axisdf_gwas_metal_freezing_duration_group_scaleT <- df_gwas_metal_freezing_duration_group_scaleT %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

g_gwas_metal_freezing_duration_group_scaleT <- 
  ggplot(df_gwas_metal_freezing_duration_group_scaleT, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR))) +#,
  # alpha = as.factor(is_sig),
  # size = as.factor(is_sig))) +
  
  # Show all points
  geom_point(shape = 16, alpha = 0.2) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) + #c("#c8c2be","#71686c") for single
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  #  geom_point( aes(fill=as.factor(is_col)), alpha=0.8, size=1.3) +
  scale_fill_viridis_d() +
  # scale_shape_manual(values = c(19,8)) +
  #  scale_color_manual(values = rep(c("grey", "#cd5c5c"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf_gwas_metal_freezing_duration_group_scaleT$CHR, breaks = axisdf_gwas_metal_freezing_duration_group_scaleT$center) +
  #  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 11)) +
  
  #  geom_hline(yintercept=-log10(alpha), linetype="dotted", color="#cd5c5c", size=1) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  # facet_wrap(~ sex, nrow = 2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_metal_freezing_duration_group_scaleT
ggsave("../data/1_single_strain/gwas/METAL/manhattan_plot/gwas_metal_freezing_duration_group_scaleT2.png", 
       g_gwas_metal_freezing_duration_group_scaleT, 
       width = 6, height = 2)


nSNPs <- nrow(df_gwas_metal_freezing_duration_group_scaleT)
ci <- 0.95
qqplot_df_gwas_metal_freezing_duration_group_scaleT <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_metal_freezing_duration_group_scaleT$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme()
ggsave("../data/1_single_strain/gwas/METAL/qqplot/qqplot_df_gwas_metal_freezing_duration_group_scaleT.png", 
       qqplot_df_gwas_metal_freezing_duration_group_scaleT, 
       width = 3, height = 3)


##### motion cue #####
###### gwas_metal_motion_cue_exit_intercept_scaleT ######
gwas_metal_motion_cue_exit_intercept_scaleT <- read.table("../data/1_single_strain/gwas/METAL/rawdata/metal_motion_cue_exit_intercept.tbl",h=T) %>%
  separate(MarkerName, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(BP))

# gwas_metal_motion_cue_exit_intercept_scaleT$FDR <-
#   p.adjust(gwas_metal_motion_cue_exit_intercept_scaleT$P.value, 
#            method = "fdr")

alpha <- 1/(10**6) #0.05/2672211
# alpha <- 0.05/(10**6)
df_gwas_metal_motion_cue_exit_intercept_scaleT <- gwas_metal_motion_cue_exit_intercept_scaleT %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_metal_motion_cue_exit_intercept_scaleT, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  dplyr::rename(P = P.value) %>%
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_metal_motion_cue_exit_intercept_scaleT, "../data/1_single_strain/gwas/METAL/df/df_gwas_metal_motion_cue_exit_intercept_scaleT.parquet")
# df_gwas_metal_motion_cue_exit_intercept_scaleT <- read_parquet("../data/1_single_strain/gwas/METAL/df/df_gwas_metal_motion_cue_exit_intercept_scaleT.parquet") %>%
#   ungroup()

df_gwas_metal_motion_cue_exit_intercept_scaleT_top <- df_gwas_metal_motion_cue_exit_intercept_scaleT %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_metal_motion_cue_exit_intercept_scaleT_top %>%
              dplyr::select(CHR, BP), 
            "../data/1_single_strain/gwas/METAL/tophits/tophits_motion_cue_exit_intercept_scaleT.bed", delim = " ", col_names = F)

write_delim(df_gwas_metal_motion_cue_exit_intercept_scaleT_top %>%
              dplyr::select(c(CHR, BP), everything()),
            "../data/1_single_strain/gwas/METAL/tophits/tophits_motion_cue_exit_intercept_scaleT.txt", delim = "\t")

axisdf_gwas_metal_motion_cue_exit_intercept_scaleT <- df_gwas_metal_motion_cue_exit_intercept_scaleT %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

g_gwas_metal_motion_cue_exit_intercept_scaleT <- 
  ggplot(df_gwas_metal_motion_cue_exit_intercept_scaleT, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR))) +#,
  # alpha = as.factor(is_sig),
  # size = as.factor(is_sig))) +
  
  # Show all points
  geom_point(shape = 16, alpha = 0.2) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) + #c("#c8c2be","#71686c") for single
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  #  geom_point( aes(fill=as.factor(is_col)), alpha=0.8, size=1.3) +
  scale_fill_viridis_d() +
  # scale_shape_manual(values = c(19,8)) +
  #  scale_color_manual(values = rep(c("grey", "#cd5c5c"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf_gwas_metal_motion_cue_exit_intercept_scaleT$CHR, breaks = axisdf_gwas_metal_motion_cue_exit_intercept_scaleT$center) +
  #  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 11)) +
  
  #  geom_hline(yintercept=-log10(alpha), linetype="dotted", color="#cd5c5c", size=1) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  # facet_wrap(~ sex, nrow = 2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_metal_motion_cue_exit_intercept_scaleT
ggsave("../data/1_single_strain/gwas/METAL/manhattan_plot/gwas_metal_motion_cue_exit_intercept_scaleT2.png", 
       g_gwas_metal_motion_cue_exit_intercept_scaleT, 
       width = 6, height = 2)


nSNPs <- nrow(df_gwas_metal_motion_cue_exit_intercept_scaleT)
ci <- 0.95
qqplot_df_gwas_metal_motion_cue_exit_intercept_scaleT <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_metal_motion_cue_exit_intercept_scaleT$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme()
ggsave("../data/1_single_strain/gwas/METAL/qqplot/qqplot_df_gwas_metal_motion_cue_exit_intercept_scaleT.png", 
       qqplot_df_gwas_metal_motion_cue_exit_intercept_scaleT, 
       width = 3, height = 3)



###### gwas_metal_motion_cue_stop_freezing_abs_scaleT ######
gwas_metal_motion_cue_stop_freezing_abs_scaleT <- read.table("../data/1_single_strain/gwas/METAL/rawdata/metal_motion_cue_stop_freezing_abs.tbl",h=T) %>%
  separate(MarkerName, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(BP))

# gwas_metal_motion_cue_stop_freezing_abs_scaleT$FDR <-
#   p.adjust(gwas_metal_motion_cue_stop_freezing_abs_scaleT$P.value, 
#            method = "fdr")

alpha <- 1/(10**6) #0.05/2672211
# alpha <- 0.05/(10**6)
df_gwas_metal_motion_cue_stop_freezing_abs_scaleT <- gwas_metal_motion_cue_stop_freezing_abs_scaleT %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_metal_motion_cue_stop_freezing_abs_scaleT, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  dplyr::rename(P = P.value) %>%
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_metal_motion_cue_stop_freezing_abs_scaleT, "../data/1_single_strain/gwas/METAL/df/df_gwas_metal_motion_cue_stop_freezing_abs_scaleT.parquet")
# df_gwas_metal_motion_cue_stop_freezing_abs_scaleT <- read_parquet("../data/1_single_strain/gwas/METAL/df/df_gwas_metal_motion_cue_stop_freezing_abs_scaleT.parquet") %>%
#   ungroup()

df_gwas_metal_motion_cue_stop_freezing_abs_scaleT_top <- df_gwas_metal_motion_cue_stop_freezing_abs_scaleT %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_metal_motion_cue_stop_freezing_abs_scaleT_top %>%
              dplyr::select(CHR, BP), 
            "../data/1_single_strain/gwas/METAL/tophits/tophits_motion_cue_stop_freezing_abs_scaleT.bed", delim = " ", col_names = F)

write_delim(df_gwas_metal_motion_cue_stop_freezing_abs_scaleT_top %>%
              dplyr::select(c(CHR, BP), everything()),
            "../data/1_single_strain/gwas/METAL/tophits/tophits_motion_cue_stop_freezing_abs_scaleT.txt", delim = "\t")

axisdf_gwas_metal_motion_cue_stop_freezing_abs_scaleT <- df_gwas_metal_motion_cue_stop_freezing_abs_scaleT %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

g_gwas_metal_motion_cue_stop_freezing_abs_scaleT <- 
  ggplot(df_gwas_metal_motion_cue_stop_freezing_abs_scaleT, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR))) +#,
  # alpha = as.factor(is_sig),
  # size = as.factor(is_sig))) +
  
  # Show all points
  geom_point(shape = 16, alpha = 0.2) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) + #c("#c8c2be","#71686c") for single
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  #  geom_point( aes(fill=as.factor(is_col)), alpha=0.8, size=1.3) +
  scale_fill_viridis_d() +
  # scale_shape_manual(values = c(19,8)) +
  #  scale_color_manual(values = rep(c("grey", "#cd5c5c"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf_gwas_metal_motion_cue_stop_freezing_abs_scaleT$CHR, breaks = axisdf_gwas_metal_motion_cue_stop_freezing_abs_scaleT$center) +
  #  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 11)) +
  
  #  geom_hline(yintercept=-log10(alpha), linetype="dotted", color="#cd5c5c", size=1) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  # facet_wrap(~ sex, nrow = 2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_metal_motion_cue_stop_freezing_abs_scaleT
ggsave("../data/1_single_strain/gwas/METAL/manhattan_plot/gwas_metal_motion_cue_stop_freezing_abs_scaleT2.png", 
       g_gwas_metal_motion_cue_stop_freezing_abs_scaleT, 
       width = 6, height = 2)


nSNPs <- nrow(df_gwas_metal_motion_cue_stop_freezing_abs_scaleT)
ci <- 0.95
qqplot_df_gwas_metal_motion_cue_stop_freezing_abs_scaleT <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_metal_motion_cue_stop_freezing_abs_scaleT$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme()
ggsave("../data/1_single_strain/gwas/METAL/qqplot/qqplot_df_gwas_metal_motion_cue_stop_freezing_abs_scaleT.png", 
       qqplot_df_gwas_metal_motion_cue_stop_freezing_abs_scaleT, 
       width = 3, height = 3)



##### nnd #####
###### gwas_metal_nnd_scaleT ######
gwas_metal_nnd_scaleT <- read.table("../data/1_single_strain/gwas/METAL/rawdata/metal_nnd.tbl",h=T) %>%
  separate(MarkerName, c("CHR", "BP", "TYPE"), sep = "_") %>%
  mutate(BP = as.numeric(BP))

# gwas_metal_nnd_scaleT$FDR <-
#   p.adjust(gwas_metal_nnd_scaleT$P.value, 
#            method = "fdr")

alpha <- 1/(10**6) #0.05/2672211
# alpha <- 0.05/(10**6)
df_gwas_metal_nnd_scaleT <- gwas_metal_nnd_scaleT %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len = max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas_metal_nnd_scaleT, ., by = c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  dplyr::rename(P = P.value) %>%
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas_metal_nnd_scaleT, "../data/1_single_strain/gwas/METAL/df/df_gwas_metal_nnd_scaleT.parquet")
# df_gwas_metal_nnd_scaleT <- read_parquet("../data/1_single_strain/gwas/METAL/df/df_gwas_metal_nnd_scaleT.parquet") %>%
#   ungroup()

df_gwas_metal_nnd_scaleT_top <- df_gwas_metal_nnd_scaleT %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_metal_nnd_scaleT_top %>%
              dplyr::select(CHR, BP), 
            "../data/1_single_strain/gwas/METAL/tophits/tophits_nnd_scaleT.bed", delim = " ", col_names = F)

write_delim(df_gwas_metal_nnd_scaleT_top %>%
              dplyr::select(c(CHR, BP), everything()),
            "../data/1_single_strain/gwas/METAL/tophits/tophits_nnd_scaleT.txt", delim = "\t")

axisdf_gwas_metal_nnd_scaleT <- df_gwas_metal_nnd_scaleT %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

g_gwas_metal_nnd_scaleT <- 
  ggplot(df_gwas_metal_nnd_scaleT, 
         aes(x = BPcum, 
             y = -log10(P),
             color = as.factor(CHR))) +#,
  # alpha = as.factor(is_sig),
  # size = as.factor(is_sig))) +
  
  # Show all points
  geom_point(shape = 16, alpha = 0.2) +
  scale_color_manual(values = rep(c("#c099a0","#a25768"), 22)) + #c("#c8c2be","#71686c") for single
  scale_size_manual(values = c(1.3, 2.4)) +
  scale_alpha_manual(values = c(0.2, 0.9)) +
  #  geom_point( aes(fill=as.factor(is_col)), alpha=0.8, size=1.3) +
  scale_fill_viridis_d() +
  # scale_shape_manual(values = c(19,8)) +
  #  scale_color_manual(values = rep(c("grey", "#cd5c5c"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf_gwas_metal_nnd_scaleT$CHR, breaks = axisdf_gwas_metal_nnd_scaleT$center) +
  #  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 11)) +
  
  #  geom_hline(yintercept=-log10(alpha), linetype="dotted", color="#cd5c5c", size=1) +
  xlab("Chromosome") +
  ylab(expression(paste("-Log"[10],"(",italic(P),")"))) +
  # facet_wrap(~ sex, nrow = 2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
# g_gwas_metal_nnd_scaleT
ggsave("../data/1_single_strain/gwas/METAL/manhattan_plot/gwas_metal_nnd_scaleT2.png", 
       g_gwas_metal_nnd_scaleT, 
       width = 6, height = 2)


nSNPs <- nrow(df_gwas_metal_nnd_scaleT)
ci <- 0.95
qqplot_df_gwas_metal_nnd_scaleT <- ggplot(data.frame(
  observed = -log10(sort(df_gwas_metal_nnd_scaleT$P)),
  expected = -log10(ppoints(nSNPs)),
  clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs)))),
  cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = seq(nSNPs), shape2 = rev(seq(nSNPs))))),
  aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymax = cupper, ymin = clower), fill = "grey30", alpha = 0.5) +
  geom_step(color = normentR::norment_colors[["purple"]], linewidth = 1.1, direction = "vh") +
  geom_segment(data = . %>% filter(expected == max(expected)), 
               aes(x = 0, xend = expected, y = 0, yend = expected),
               linewidth = 1.25, alpha = 0.5, color = "grey30", lineend = "round") +
  labs(x = expression(paste("Expected -Log"[10],"(", italic(P),")")),
       y = expression(paste("Observed -Log"[10],"(", italic(P),")"))) +
  theme_bw() +
  theme()
ggsave("../data/1_single_strain/gwas/METAL/qqplot/qqplot_df_gwas_metal_nnd_scaleT.png", 
       qqplot_df_gwas_metal_nnd_scaleT, 
       width = 3, height = 3)

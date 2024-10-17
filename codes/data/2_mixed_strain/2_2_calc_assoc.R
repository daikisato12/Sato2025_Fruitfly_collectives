#### load libraries ####
library(tidyverse)
library(arrow)
# library(devtools)
# devtools::install_github("norment/normentR")
library(normentR)

#### load functions ####
calc_regress_pi_pca_dist <- function(dat){

  dat <- dat %>%
    left_join(dfd_pheno %>% dplyr::select(strain, overyield), by = "strain") %>%
    left_join(dfd_pca_dist)
  
  if(nrow(dat) > 1){
    res.test <- lm(overyield ~ PI + dist, data = dat) %>% summary()
    dat.res <- data.frame(est = res.test$coefficients[2,1], p.value = res.test$coefficients[2,4])
  }else{
    dat.res <- data.frame(est = NA, p.value = NA)
  }
  
  return(dat.res)
  
}


#### load dataset ####
##### 1kbp #####
# dfd_pi_1kbp <- read.table("../data/2_mixed_strain/gwas/pi/pi_2pop_1000_all.windowed.pi", header = T) %>%
#   mutate(MIX = str_replace_all(MIX, "line_", "DGRP")) %>%
#   mutate(MIX = str_replace(MIX, "-", "_")) %>%
#   dplyr::rename(strain = MIX)
# 
# num_comb <- unique(dfd_pi_1kbp$strain) %>% length()
# 
# dfd_pi_1kbp_count_na <- dfd_pi_1kbp %>%
#   mutate(loci = paste0(CHROM, "_", BIN_START)) %>%
#   pivot_wider(id_cols = loci, names_from = strain, values_from = PI) %>%
#   tibble::column_to_rownames(var = "loci") %>%
#   mutate(count_na = rowSums(is.na(.))) %>%
#   filter(count_na != num_comb) %>%
#   replace(is.na(.), 0) %>%
#   dplyr::select(!count_na)
# 
# dfd_pi_1kbp2 <- dfd_pi_1kbp_count_na %>%
#   tibble::rownames_to_column(var = "loci") %>%
#   pivot_longer(cols = !loci, names_to = "strain", values_to = "PI") %>%
#   separate(loci, into = c("CHROM", "BIN_START")) %>%
#   mutate(BIN_START = as.integer(BIN_START)) %>%
#   left_join(dfd_pi_1kbp %>% dplyr::select(!PI)) %>%
#   mutate(BIN_END = BIN_START + 999,
#          N_VARIANTS = if_else(is.na(N_VARIANTS), 0, N_VARIANTS)) %>%
#   replace(is.na(.), 0)
# write_parquet(dfd_pi_1kbp2, "../data/2_mixed_strain/gwas/pi/dfd_pi_1kbp2.parquet")
# dfd_pi_1kbp2 <- read_parquet("../data/2_mixed_strain/gwas/pi/dfd_pi_1kbp2.parquet")


##### distance on pc1/2 dimension #####
strain_2pop <- c("line_88", "line_101", "line_136", "line_161",
                 "line_189", "line_301", "line_309", "line_324", 
                 "line_357", "line_358", "line_360", "line_399", 
                 "line_786", "line_852", "line_855")

df_pca <- read.table("../data/genome/dgrp2.eigenvec", header = F) %>%
  dplyr::select(!V2) %>%
  magrittr::set_colnames(c("strain", paste0("PC", seq(1, 10)))) %>%
  filter(strain %in% strain_2pop) %>%
  dplyr::select(strain, PC1, PC2)

dfd_pca_dist <- df_pca %>%
  sf::st_as_sf(., coords = c("PC1","PC2")) %>% 
  sf::st_distance() %>%
  magrittr::set_rownames(df_pca$strain) %>%
  magrittr::set_colnames(df_pca$strain) %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  pivot_longer(cols = !rowname, names_to = "Strain2", values_to = "Distance1") %>%
  dplyr::rename(Strain1 = rowname) %>%
  dplyr::mutate(Strain1 = parse_number(Strain1),
                Strain2 = parse_number(Strain2)) %>%
  filter(Strain2 > Strain1) %>%
  dplyr::mutate(Strain1 = paste0("DGRP", Strain1),
                Strain2 = paste0("DGRP", Strain2)) %>%
  unite(col = strain, Strain1, Strain2, sep = "_") %>%
  dplyr::rename(dist = Distance1)

write_parquet(dfd_pca_dist, "../data/2_mixed_strain/gwas/pi/dfd_pca_dist.parquet")
dfd_pca_dist <- read_parquet("../data/2_mixed_strain/gwas/pi/dfd_pca_dist.parquet")


#### 1kbp foraging_success5 regression cov dist pca ####
##### load dataset #####
dfd_pi_1kbp2 <- read_parquet("../data/2_mixed_strain/gwas/pi/dfd_pi_1kbp2.parquet")
dfd_pca_dist <- read_parquet("../data/2_mixed_strain/gwas/pi/dfd_pca_dist.parquet")
dfd_foraging_success5 <- read_tsv("../data/2_mixed_strain/gwas/pheno/dfd_s5min_2995_5990_speed_raw_foraging_sec3.tsv")

##### analysis #####
dfd_pheno <- dfd_foraging_success5
start <- Sys.time()
dfd_pi_1kbp.res <- dfd_pi_1kbp2 %>%
  group_nest(CHROM, BIN_START, BIN_END) %>%
  mutate(data = map(data, calc_regress_pi_pca_dist)) %>%
  unnest(data)
end <- Sys.time()
diff <- end - start

##### plot #####
gwas.res <- dfd_pi_1kbp.res %>%
  filter(!is.na(p.value)) %>%
  mutate(BP = (BIN_START - 1 + BIN_END)/2) %>%
  dplyr::rename(CHR = CHROM,
                est = est,
                P = p.value) %>%
  dplyr::select(!c(BIN_START, BIN_END)) %>%
  dplyr::select(CHR, BP, est, P)

alpha = 1/nrow(gwas.res)

df_gwas <- gwas.res %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas.res, by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum = BP + tot) %>%
  
  mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas, "../data/2_mixed_strain/gwas/df/df_gwas_1kbp_foraging_success5_regress_pi_pca_dist.parquet")
# df_gwas <- read_parquet("../data/2_mixed_strain/gwas/rawdata/df_gwas_1kbp_foraging_success5_regress_pi_pca_dist.parquet") %>%
#   ungroup()


df_gwas_1kbp_top <- df_gwas %>%
  mutate(est_Z = est) %>%
  as.data.frame() %>%
  mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")


write_delim(df_gwas_1kbp_top,
            "../data/2_mixed_strain/gwas/tophits/tophits_1kbp_foraging_success5_regress_pi_pca_dist.txt", delim = "\t")


write_delim(df_gwas_1kbp_top %>%
              filter(est > 0),
            "../data/2_mixed_strain/gwas/tophits/tophits_1kbp_foraging_success5_regress_pi_pca_dist_plus.txt", delim = "\t")

write_delim(df_gwas_1kbp_top %>%
              filter(est < 0),
            "../data/2_mixed_strain/gwas/tophits/tophits_1kbp_foraging_success5_regress_pi_pca_dist_minus.txt", delim = "\t")



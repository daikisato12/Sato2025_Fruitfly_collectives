#### load packages ####
targetPackages <- c('tidyverse','arrow','normentR')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

# library(devtools)
# devtools::install_github("norment/normentR")

#### load functions ####
# strain1とstrain2間の遺伝的相関を取得する関数
get_grm_value <- function(strain1, strain2, grm_matrix) {
  return(grm_matrix[as.character(strain1), as.character(strain2)])
}

calc_regress_pi_grm <- function(dat){
  # dat <- dfd_pi_1kbp2 %>%
  #   filter(CHROM == "2L", BIN_START == 5001) %>%
  #   left_join(dfd_pheno %>% dplyr::select(strain, DE), by = "strain") %>%
  #   separate(strain, into = c("strain1", "strain2"))
  
  dat <- dat %>%
    dplyr::mutate(strain2 = strain) %>%
    separate(strain2, into = c("strain1", "strain2"))
  dat$genetic_relatedness <- 
    mapply(get_grm_value, dat$strain1, dat$strain2, MoreArgs = list(grm_matrix))
  dat <- dat %>%
    left_join(dfd_pheno %>% dplyr::select(strain, DE), by = "strain")
  
  if(nrow(dat) > 1){
    res.test <- lm(DE ~ PI + genetic_relatedness, data = dat) %>% summary()
    dat.res <- data.frame(est = res.test$coefficients[2,1], p.value = res.test$coefficients[2,4])
  }else{
    dat.res <- data.frame(est = NA, p.value = NA)
  }
  
  return(dat.res)
  
}



#### load dataset ####
##### load GRM dataset #####
# Strain information
list_strain_d <- read.csv("../data/2_mixed_strain/list_strain_d.txt") %>%
  pull(x) 
list_strain_d_s <- list_strain_d %>%
  str_split("_") %>%
  unlist() %>%
  unique() %>%
  head(15)

# GRM dataset
grm_bin_file <- "../data/2_mixed_strain/gwas/grm/dgrp2.grm.bin"
grm_id_file <- "../data/2_mixed_strain/gwas/grm/dgrp2.grm.id"

# 1. 個体IDを読み込み、個体数を取得
id_data <- read.table(grm_id_file)
n <- nrow(id_data)  # 個体数

# 2. バイナリ形式のGRMデータを読み込み
grm_bin <- file(grm_bin_file, "rb")
grm_values <- readBin(grm_bin, what = numeric(), n = n * (n + 1) / 2, size = 4)
close(grm_bin)

# 3. 下三角行列の要素を完全な行列に変換
grm_matrix <- matrix(0, n, n)
k <- 1
for (i in 1:n) {
  for (j in 1:i) {
    grm_matrix[i, j] <- grm_values[k]
    grm_matrix[j, i] <- grm_values[k]  # 対称性を反映
    k <- k + 1
  }
}
grm_matrix <- 
  grm_matrix %>%
  as.data.frame() %>%
  # tibble::rowid_to_column() %>%
  magrittr::set_colnames(str_replace(id_data$V1, "line_", "DGRP")) %>%
  tibble::rowid_to_column() %>%
  dplyr::mutate(rowid = str_replace(id_data$V1, "line_", "DGRP")) %>%
  pivot_longer(cols = !rowid, names_to = "ID2", values_to = "dist") %>%
  dplyr::rename(rowname = rowid) %>%
  filter(rowname %in% list_strain_d_s, 
         ID2 %in% list_strain_d_s) %>%
  pivot_wider(id_cols = rowname, names_from = ID2, values_from = dist) %>%
  tibble::column_to_rownames() %>%
  as.matrix()



#### 1kbp performanceDE regression cov grm ####
dfd_pi_1kbp2 <- read_parquet("../data/2_mixed_strain/gwas/pi/dfd_pi_1kbp2.parquet")
for_sec <- 3
##### load dataset #####
dfd_pheno <- read_tsv("../data/2_mixed_strain/dfd_fig4.tsv") %>%
  filter(thr_sec == for_sec) %>%
  dplyr::mutate(DE = performance - lag(performance)) %>%
  filter(alpha == "Observed")

##### analysis #####
start <- Sys.time()
dfd_pi_1kbp.res <- dfd_pi_1kbp2 %>%
  group_nest(CHROM, BIN_START, BIN_END) %>%
  dplyr::mutate(data = map(data, calc_regress_pi_grm)) %>%
  unnest(data)
end <- Sys.time()
diff <- end - start
print(diff)

##### output #####
gwas.res <- dfd_pi_1kbp.res %>%
  filter(!is.na(p.value)) %>%
  dplyr::mutate(BP = (BIN_START - 1 + BIN_END)/2) %>%
  dplyr::rename(CHR = CHROM,
                est = est,
                P = p.value) %>%
  dplyr::select(!c(BIN_START, BIN_END)) %>%
  dplyr::select(CHR, BP, est, P)

alpha = 1/(10**5)#0.05/nrow(gwas.res)

df_gwas <- gwas.res %>%
  
  # Compute chromosome size
  group_by(CHR) %>% 
  dplyr::summarize(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  dplyr::mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas.res, by=c("CHR" = "CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  dplyr::mutate(BPcum = BP + tot) %>%
  
  dplyr::mutate(is_sig = if_else(P < alpha, "Y", "N"))

write_parquet(df_gwas, paste0("../data/2_mixed_strain/gwas/result/rawdata/df_gwas_1kbp_performanceDE_", for_sec, "sec_regress_pi_grm.parquet"))
# df_gwas <- read_parquet("../data/2_mixed_strain/gwas/result/rawdata/df_gwas_1kbp_performanceDE_regress_pi_grm.parquet") %>%
#   ungroup()

df_gwas_1kbp_top <- df_gwas %>%
  dplyr::mutate(est_Z = est) %>%
  as.data.frame() %>%
  dplyr::mutate_at(vars(ends_with("Z")), ~(scale(.) %>% as.vector)) %>%
  arrange(P) %>%
  filter(is_sig == "Y")

write_delim(df_gwas_1kbp_top %>%
              filter(est > 0) %>%
              dplyr::select(CHR, BP), 
            paste0("../data/2_mixed_strain/gwas/result/tophits/tophits_1kbp_performanceDE_", for_sec, "sec_regress_pi_grm_plus.bed"), delim = " ", col_names = F)

write_delim(df_gwas_1kbp_top %>%
              filter(est > 0),
            paste0("../data/2_mixed_strain/gwas/result/tophits/tophits_1kbp_performanceDE_", for_sec, "sec_regress_pi_grm_plus.txt"), delim = "\t")

axisdf_gwas <- df_gwas %>% 
  group_by(CHR) %>% 
  dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

g_gwas <- ggplot(df_gwas, 
                 aes(x = BPcum, 
                     y = -log10(P),
                     color = as.factor(CHR))) +#,
  # alpha = as.factor(is_sig))) +#,
  # size = as.factor(is_sig))) +
  
  # Show all points
  geom_point(shape = 16, alpha = 0.2) +
  scale_color_manual(values = rep(c("#9e8896","#874c70"), 22)) +
  # scale_size_manual(values = c(1.3, 2.4)) +
  # scale_alpha_manual(values = c(0.2, 0.9)) +
  #  geom_point( aes(fill=as.factor(is_col)), alpha=0.8, size=1.3) +
  scale_fill_viridis_d() +
  # scale_shape_manual(values = c(19,8)) +
  #  scale_color_manual(values = rep(c("grey", "#cd5c5c"), 22 )) +
  
  # custom X axis:
  scale_x_continuous(label = axisdf_gwas$CHR, breaks = axisdf_gwas$center) +
  #  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6)) +
  coord_cartesian(ylim = c(0, 10)) +
  
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
g_gwas
ggsave(paste0("../data/2_mixed_strain/gwas/result/manhattanplot/manhattanplot_1kbp_performanceDE_", for_sec, "sec_regress_pi_grm2.png"), g_gwas, width=6, height=2)


nwindows <- nrow(df_gwas)
ci <- 0.95
qqplot_df_gwas <- ggplot(data.frame(
  observed = -log10(sort(df_gwas$P)),
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
  theme()
qqplot_df_gwas
ggsave(paste0("../data/2_mixed_strain/gwas/result/qqplot/qqplot_df_gwas_1kbp_performanceDE_", for_sec, "sec_regress_pi_grm.png"), qqplot_df_gwas, width=3, height=3)


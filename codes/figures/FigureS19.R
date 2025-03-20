#### load packages ####
targetPackages <- c('sangerseqR','annotate')
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)
# BiocManager::install(c("sangerseqR","annotate")) 

#### load dataset ####
dir_path <- "../data/7_norpAsequence/ab1/"
file_list <- list.files(path = dir_path, pattern = "norpAseqF3-PREMIX\\.ab1$", full.names = TRUE, recursive = TRUE)
data_list <- lapply(file_list, sangerseqR::read.abif)
# list trim5 and trim3 number
trim5 <- c(185, 175, 178)
trim3 <- c(290, 287, 292)


#### plot chromatogram ####
for (i in seq_along(data_list)) {
  # i <- 8
  file <- basename(file_list[i]) %>% tools::file_path_sans_ext()
  pdf(paste0("../figures/FigureS19_", file, ".pdf"), w = 5, h = 2)
  par(mfrow = c(1, 1))
  sanger_seq <- sangerseq(data_list[[i]])
  chromatogram(sanger_seq, trim5 = trim5[i], trim3 = trim3[i], width = 28)  # 必要に応じてトリム値を調整
  dev.off()
}

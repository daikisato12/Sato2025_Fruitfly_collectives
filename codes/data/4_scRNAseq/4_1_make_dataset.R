#### load packages ####
library(tidyverse)
library(data.table)
library(magrittr)
library(pals)
library(ggrepel)
library(ggpubr)
# install.packages("Seurat")
library(Seurat)
library(car)
library(lme4)
Sys.setenv('R_MAX_VSIZE'=32000000000)

#### load functions ####
##### glmer function #####
glm_exp <- function(dat, gene){

  tryCatch({
    dat_tmp <- dat %>% 
      group_by(genotype2) %>%
      summarize(n = n())
    if(length(unique(dat$genotype2)) > 1 & length(unique(dat %>% pull(get(gene)))) > 1){
      if(nrow(dat_tmp[dat_tmp$n > 2,]) > 1){
        dat_tmp2 <- dat %>% 
          group_by(rep) %>%
          summarize(n = n())
        if(nrow(dat_tmp2[dat_tmp2$n > 2,]) > 1){
          myformula <- as.formula(paste0(gene, " ~ genotype2 + (1|genotype) + (1|rep)"))
          p.val <- glmer(myformula,
                         data = dat,
                         family = "Gamma") %>%
            Anova() %>%
            pull(`Pr(>Chisq)`)
          method <- "glmer"
        }else{
          myformula <- as.formula(paste0(gene, " ~ genotype2 + (1|genotype)"))
          p.val <- glmer(myformula,
                         data = dat,
                         family = "Gamma") %>%
            Anova() %>%
            pull(`Pr(>Chisq)`)
          method <- "glmer"
        }
      }else{
        myformula <- as.formula(paste0(gene, " ~ genotype2"))
        p.val <- glm(myformula,
                     data = dat,
                     family = "Gamma") %>%
          Anova() %>%
          pull(`Pr(>Chisq)`)
        method <- "glm"
      }
    }else{
      p.val <- NA_real_
      method = NA_character_
    }
  }, 
  error = function(e) {
    message(e)           
    p.val <- NA_real_
    method = NA_character_
  })
  dat2 <- dat %>%
    mutate(pval = p.val,
           method = method)
  return(dat2)
}

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", NA_character_)
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}


#### main ####
##### load dataset #####
df_main_matrix <- data.table::fread("../data/main/GSE156455_matrix_main.mtx", skip = 2, header = F, sep = " ") %>%
  magrittr::set_colnames(c("gene_rowname", "barcode_rowname", "count"))

df_main_feature <- data.table::fread("../data/main/GSE156455_features_main.tsv", header = F) %>%
  magrittr::set_colnames(c("GeneID", "Genename", "tmp")) %>%
  tibble::rownames_to_column() %>%
  rename(gene_rowname = rowname) %>%
  dplyr::select(!contains("tmp")) %>%
  mutate(gene_rowname = as.integer(gene_rowname))

df_main_meta <- data.table::fread("../data/main/GSE156455_metadata_main.tsv", header = T)

df_main_barcode <- data.table::fread("../data/main/GSE156455_barcodes_main.tsv", header = F) %>%
  magrittr::set_colnames(c("barcode")) %>%
  tibble::rownames_to_column() %>%
  rename(barcode_rowname = rowname) %>%
  inner_join(df_main_meta) %>%
  mutate(barcode_rowname = as.integer(barcode_rowname))

df_main_tsne <- data.table::fread("../data/main/GSE156455_tsne_main.tsv", header = T)


#### with seurat main ####
dros.main.data <-  Read10X(data.dir = "../data/main/seurat_format/")
dros.main <- CreateSeuratObject(counts = dros.main.data, project = "dros.main", min.cells = 3, min.features = 200)

##### QC ##### this is not necessary this time
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
dros.main[["percent.mt"]] <- PercentageFeatureSet(dros.main, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(dros.main, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(dros.main, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(dros.main, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

##### normalization #####
dros.main <- NormalizeData(dros.main, normalization.method = "LogNormalize", scale.factor = 10000)
# dros.main <- NormalizeData(dros.main) #same as avobe command


##### Ptp99A #####
###### visualize expression of Ptp99A ######
# dros.main.dpr13 <- dros.main@assays$RNA$data["dpr13",]
# dros.main.dpr13_name <- names(dros.main.dpr13)
# df.dros.main.dpr13 <- data.frame(dpr13 = dros.main.dpr13,
#                                   barcode = dros.main.dpr13_name)

dros.main.ptp99a <- dros.main@assays$RNA$data["Ptp99A",]
dros.main.ptp99a_name <- names(dros.main.ptp99a)

df.dros.main.ptp99a <- data.frame(ptp99a = dros.main.ptp99a,
                                  barcode = dros.main.ptp99a_name)
# write.table(df.dros.main.ptp99a, 
#             "../data/main/df.dros.main.ptp99a.tsv",
#             quote = F, row.names = F, sep = "\t")
# df.dros.main.ptp99a <- read.table("../data/main/df.dros.main.ptp99a.tsv", header = TRUE)

df_main.ptp99a <- inner_join(df_main_meta,
                             df_main_tsne) %>%
  inner_join(df.dros.main.ptp99a) %>%
  mutate(genotype2 = if_else(genotype %in% c(paste0("line_", c("21", "40", "129", "235", "304", "320", "395", "508", "748", "805", "819"))), "C/T", "T/T")) %>%
  mutate(ptp99a = exp(ptp99a)) %>%
  distinct()

write.table(df_main.ptp99a, 
            "../data/4_scRNAseq/df_main_ptp99a.tsv",
            quote = F, row.names = F, sep = "\t")
# df_main.ptp99a <- read.table("../data/main/df_main_ptp99a.tsv", header = TRUE)


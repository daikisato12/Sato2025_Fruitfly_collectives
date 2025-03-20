#### load packages ####
targetPackages <- c('tidyverse','Seurat','hdWGCNA','harmony',
                    'biomaRt','clusterProfiler') 
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)

# devtools::install_github('smorabit/hdWGCNA', ref='dev')
# devtools::install_github("immunogenomics/harmony")

#### load functions ####
glm_exp <- function(dat, gene){
  
  tryCatch({
    dat_tmp <- dat %>% 
      group_by(genotype2) %>%
      summarize(n = n())
    if(length(unique(dat$genotype2)) > 1 & length(unique(dat %>% pull(get(gene)))) > 1){
      direction <- if_else((dat[dat$genotype2=="C/T",] %>% pull(gene) %>% mean()) < (dat[dat$genotype2=="T/T",] %>% pull(gene) %>% mean()), "minus", "plus")
      if(nrow(dat_tmp[dat_tmp$n > 2,]) > 1){
        dat_tmp2 <- dat %>% 
          group_by(rep) %>%
          summarize(n = n())
        if(nrow(dat_tmp2[dat_tmp2$n > 2,]) > 1){
          tryCatch({
            myformula <<- as.formula(paste0(gene, " ~ genotype2 + (1|genotype) + (1|rep)"))
            p.val <<- lme4::glmer(myformula,
                                  data = dat,
                                  family = "Gamma") %>%
              car::Anova() %>%
              pull(`Pr(>Chisq)`)
            method <<- "lme4::glmer"
          }, error = function(e) {
            print("Error encountered, switching to glm.")
            myformula <<- as.formula(paste0(gene, " ~ genotype2"))
            method <<- "glm"
            p.val <<- glm(myformula,
                          data = dat,
                          family = "Gamma") %>%
              car::Anova() %>%
              pull(`Pr(>Chisq)`)
          })
        }else{
          myformula <- as.formula(paste0(gene, " ~ genotype2 + (1|genotype)"))
          p.val <- lme4::glmer(myformula,
                               data = dat,
                               family = "Gamma") %>%
            car::Anova() %>%
            pull(`Pr(>Chisq)`)
          method <- "lme4::glmer"
        }
      }else{
        myformula <- as.formula(paste0(gene, " ~ genotype2"))
        p.val <- glm(myformula,
                     data = dat,
                     family = "Gamma") %>%
          car::Anova() %>%
          pull(`Pr(>Chisq)`)
        method <- "glm"
      }
    }else{
      p.val <- NA_real_
      method = NA_character_
      direction = NA_character_
    }
  }, 
  error = function(e) {
    message(e)           
    p.val <- NA_real_
    method = NA_character_
    direction = NA_character_
  })
  dat2 <- dat %>%
    mutate(pval = p.val,
           direction = direction,
           method = method)
  return(dat2)
}

wilcox_test <- function(dat){
  # gene <- "expression"
  # dat_test <- ME_long %>%
  #   group_nest(type, time, Module)
  # dat <- dat_test$data[[1]]
  
  tryCatch({
    dat_tmp <- dat %>% 
      group_by(condition) %>%
      summarize(n = n())
    if(length(unique(dat$condition)) > 1 & length(unique(dat %>% pull(Eigengene))) > 1){
      if(nrow(dat_tmp[dat_tmp$n > 2,]) > 1){
        wilcox_res <- dat %>%
          rstatix::wilcox_test(formula = Eigengene ~ condition)
        p_val <- wilcox_res$p
        direction <- if_else(mean(dat[dat$condition=="C",]$Eigengene) < mean(dat[dat$condition=="T",]$Eigengene), "minus", "plus")
      }else{
        p_val <- NA_real_
        direction <- NA_character_
      }
    }else{
      p_val <- NA_real_
      direction <- NA_character_
    }
  }, 
  error = function(e) {
    message(e)           
    p_val <- NA_real_
    direction <- NA_character_
  })
  dat2 <- dat %>%
    mutate(pval = p_val,
           direction = direction,
           method = "wilcoxon test")
  return(dat2)
}

makeStars <- function(x){
  stars <- c("****", "***", "**", "*", NA_character_)
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}


#### Step0: 解析データの準備 ####
df_main_meta <- data.table::fread("../data/4_scRNAseq/GSE156455_metadata_main.tsv", header = T)
df_main_tsne <- data.table::fread("../data/4_scRNAseq/GSE156455_tsne_main.tsv", header = T)

#### with seurat main ####
# dataset can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156455
dros.main.data <-  Read10X(data.dir = "../../../../../Kurmangaliyev2020/data/main/seurat_format/")
dros.main <- CreateSeuratObject(counts = dros.main.data, project = "dros.main", min.cells = 3, min.features = 200)

##### QC #####
#this is not necessary this time
dros.main[["percent.mt"]] <- PercentageFeatureSet(dros.main, pattern = "^MT-")

##### normalization #####
dros.main <- NormalizeData(dros.main, normalization.method = "LogNormalize", scale.factor = 10000)
saveRDS(dros.main, "../data/4_scRNAseq/dros.main.rds")

##### prepare expression data of ptp99a #####
dros.main.ptp99a <- dros.main@assays$RNA$data["Ptp99A",]
dros.main.ptp99a_name <- names(dros.main.ptp99a)

df.dros.main.ptp99a <- data.frame(ptp99a = dros.main.ptp99a,
                                  barcode = dros.main.ptp99a_name)

df_main.ptp99a <- inner_join(df_main_meta,
                             df_main_tsne) %>%
  inner_join(df.dros.main.ptp99a) %>%
  filter(genotype != "W1118") %>%
  mutate(genotype2 = if_else(genotype %in% c(paste0("line_", c("21", "40", "129", "235", "304", "320", "395", "508", "748", "805", "819"))), "C/T", "T/T")) %>%
  mutate(expression = exp(ptp99a)) %>%
  distinct()

df_main.ptp99a_2 <-
  df_main.ptp99a %>%
  group_nest(type, time) %>%
  mutate(data = map2(data, "expression", glm_exp)) %>% 
  unnest() %>%
  mutate(expression = log(expression),
         genotype = case_when(str_detect(genotype, "_") ~ paste0("DGRP", parse_number(genotype)),
                              TRUE ~ genotype)) %>%
  transform(genotype = factor(genotype, levels = unique(.$genotype) %>% gtools::mixedsort()))
write.table(df_main.ptp99a_2, "../data/4_scRNAseq/df_main.ptp99a_2.tsv", quote = FALSE, sep = "\t", row.names = FALSE)


#### Step 1: Seurat オブジェクトを読み込む（すでに読み込んでいる場合はスキップ） ####
dros.main <- readRDS("../data/4_scRNAseq/dros.main.rds")
df_main_meta <- data.table::fread("../data/4_scRNAseq/GSE156455_metadata_main.tsv", header = T)


#### Step 2: メタデータを追加する ####
list_C <- c(paste0("line_", c("21", "181", "189", "437", "461", "508", "897")))
dros.main2 <- AddMetaData(object = dros.main, metadata = df_main_meta %>%
                            mutate(sample = paste(genotype, rep, trep, time, sep = "_"),
                                   rowname = barcode,
                                   condition = if_else(genotype %in% list_C,
                                                       "C",
                                                       "T"),
                                   condition_type = paste0(condition, "_", type)) %>%
                            tibble::column_to_rownames())
# list_T <- setdiff(unique(dros.main2@meta.data$genotype), list_C)


#### Step 3: 必要なデータのみ抽出したサブセットを作成 ####
cell_list <- c("R1.6", "L1", "L2", "L3", "L4", "L5", "Tm1", "Tm2", "Tm3", "Tm4", "Tm9", "T4.T5") # "Mi1", "Mi4", "Mi9", 
subset_obj <- subset(dros.main2, subset = type %in% cell_list & genotype != "W1118")


#### Step 4: 変動遺伝子を検出し、データをスケーリング ####
subset_obj <- FindVariableFeatures(subset_obj, selection.method = "vst", nfeatures = 2000)
print(length(VariableFeatures(subset_obj)))
subset_obj <- ScaleData(subset_obj, features = rownames(subset_obj))


#### Step 5: PCAによる次元削減を実行 ####
subset_obj <- RunPCA(subset_obj, features = VariableFeatures(object = subset_obj))


#### Step 6: バッチ効果の補正（Harmonyを使用）####
# バッチ効果を補正
subset_obj <- RunHarmony(subset_obj, group.by.vars = c("rep", "trep", "time"))
# saveRDS(subset_obj, "../data/4_scRNAseq/subset_obj.rds")


#### Step 7: 共通データセットで WGCNA ネットワークを構築 ####
# Seurat オブジェクト全体で WGCNA のセットアップを実行
seurat_obj <- SetupForWGCNA(
  subset_obj,
  gene_select = "variable",  # 最も変動の大きい遺伝子を選択
  wgcna_name = "all_conditions"
)

# WGCNA 用にデータをセットアップ
seurat_obj <- SetDatExpr(
  seurat_obj,
  assay = 'RNA',
  slot = 'data'
)

# スケーリングパワーを選択（全体で TestSoftPowers を実行）
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = "signed"
)

# 共通の共発現ネットワークを構築
seurat_obj <- ConstructNetwork(
  seurat_obj,
  minModuleSize = 30,          # 最小モジュールサイズ
  mergeCutHeight = 0.25,       # モジュールマージのカットオフ値
  networkType = "signed",      # signed ネットワーク
  TOMType = "signed",          # TOM のタイプ
  overwrite_tom = TRUE
)


#### Step 8: 共通の Module Eigengenes を計算 ####
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by = "condition"  # condition 列を基にモジュールごとに固有遺伝子を計算
)
saveRDS(seurat_obj, "../data/4_scRNAseq/hdWGCNA/seurat_obj.rds")

# モジュールの固有遺伝子（hMEs）を取得
MEs <- GetMEs(seurat_obj, harmonized = TRUE)  # Harmonized Eigengenes
mods <- colnames(MEs) %>% sort()

# get the module assignment table:
modules <- GetModules(seurat_obj)
mart <- biomaRt::useMart(biomart = "ensembl",  #useEnsembl or useMart
                         host = "https://www.ensembl.org", 
                         dataset = "dmelanogaster_gene_ensembl")
gene_ids <- biomaRt::getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                           filters = "external_gene_name",
                           values = modules,
                           mart = mart) %>%
  dplyr::rename(gene_name = external_gene_name)
modules <- modules %>%
  left_join(gene_ids) %>%
  dplyr::select(gene_name, ensembl_gene_id, everything())
write.table(modules, paste0("../data/4_scRNAseq/hdWGCNA/genelist_modules.tsv"), 
            quote = FALSE, sep = "\t", row.names = FALSE)


##### enrichment analysis #####
list_modules <- unique(modules$module)
organism = "org.Dm.eg.db"
# BiocManager::install(organism, character.only = TRUE)
library("org.Dm.eg.db")

for (module in list_modules){
  # module <- "green"
  dir.create(paste0("../data/4_scRNAseq/hdWGCNA/enrichment/module_", module, "/"), recursive = T)
  
  ###### BP ######
  go_enrich_BP_module <- 
    clusterProfiler::enrichGO(gene = modules[modules$module==module,]$ensembl_gene_id,
             universe = modules$ensembl_gene_id,
             OrgDb = organism, 
             keyType = 'FLYBASE',
             readable = T,
             ont = "BP",
             pvalueCutoff = 0.05, 
             qvalueCutoff = 0.10)
  
  head(go_enrich_BP_module)
  
  df_enrich_BP_module <- go_enrich_BP_module@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(geneID = str_replace_all(geneID, "/", ", "),
           Type = "BP")
  write.table(df_enrich_BP_module, paste0("../data/4_scRNAseq/hdWGCNA/enrichment/module_", module, "/module_", module, "_BP.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  ###### CC ######
  go_enrich_CC_module <- 
    clusterProfiler::enrichGO(gene = modules[modules$module==module,]$ensembl_gene_id,
             universe = modules$ensembl_gene_id,
             OrgDb = organism, 
             keyType = 'FLYBASE',
             readable = T,
             ont = "CC",
             pvalueCutoff = 0.05, 
             qvalueCutoff = 0.10)
  
  head(go_enrich_CC_module)
  
  df_enrich_CC_module <- go_enrich_CC_module@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(geneID = str_replace_all(geneID, "/", ", "),
           Type = "CC")
  write.table(df_enrich_CC_module, paste0("../data/4_scRNAseq/hdWGCNA/enrichment/module_", module, "/module_", module, "_CC.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  ###### MF ######
  go_enrich_MF_module <- 
    clusterProfiler::enrichGO(gene = modules[modules$module==module,]$ensembl_gene_id,
             universe = modules$ensembl_gene_id,
             OrgDb = organism, 
             keyType = 'FLYBASE',
             readable = T,
             ont = "MF",
             pvalueCutoff = 0.05, 
             qvalueCutoff = 0.10)
  
  head(go_enrich_MF_module)
  
  df_enrich_MF_module <- go_enrich_MF_module@result %>%
    filter(p.adjust < 0.05) %>%
    mutate(geneID = str_replace_all(geneID, "/", ", "),
           Type = "MF")
  write.table(df_enrich_MF_module, paste0("../data/4_scRNAseq/hdWGCNA/enrichment/module_", module, "/module_", module, "_MF.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  
}

#### Step 9: Module Eigengeneの発現量をPtp99Aの遺伝子型で比較 ####
# モジュールの固有遺伝子（hMEs）を取得
MEs <- GetMEs(seurat_obj, harmonized = TRUE)  # Harmonized Eigengenes
mods <- colnames(MEs) %>% sort()
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

ME_long <- bind_rows(seurat_obj@meta.data %>%
                       dplyr::select(sample, type, condition, all_of(mods))) %>% 
  tibble::rownames_to_column(var = "Sample") %>%
  pivot_longer(cols = all_of(mods), 
               names_to = "Module", values_to = "Eigengene") %>%
  left_join(data.frame(seurat_obj@meta.data) %>%
              dplyr::rename(Sample = barcode) %>%
              dplyr::select(Sample, time, rep))
# write.table(ME_long, "../data/4_scRNAseq/hdWGCNA/ME_long.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

wilcox_res_time <-
  ME_long %>%
  group_nest(Module, type, time) %>%
  mutate(data = map(data, wilcox_test)) %>% 
  unnest() %>%
  dplyr::select(Module, type, time, pval, direction) %>%
  distinct()
# write.table(wilcox_res_time, "../data/4_scRNAseq/hdWGCNA/ME_wilcox_res_time.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

df_ME_wilcox_time <-
  ME_long %>%
  left_join(wilcox_res_time) %>%
  filter(Module != "grey") %>%
  group_by(condition, type, time, Module, direction) %>%
  summarize(mean = mean(Eigengene, na.rm = TRUE),
            sd = sd(Eigengene, na.rm = TRUE),
            n = n(),
            se = sd / sqrt(n),
            pval = mean(pval, na.rm = TRUE)) %>%
  mutate(star = makeStars(pval),
         Module = str_to_sentence(Module)) %>%
  transform(type = factor(type, levels = cell_list))
write.table(df_ME_wilcox_time, "../data/4_scRNAseq/hdWGCNA/df_ME_wilcox_time.tsv", quote = FALSE, sep = "\t", row.names = FALSE)

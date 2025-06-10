#### load packages ####
targetPackages <- c('tidyverse','patchwork','WGCNA',"hdWGCNA")
newPackages <- targetPackages[!(targetPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = "http://cran.us.r-project.org")
for(package in targetPackages) library(package, character.only = T)
# devtools::install_github('smorabit/hdWGCNA', ref='dev')

#### load functions ####
makeStars <- function(x){
  stars <- c("****", "***", "**", "*", NA_character_)
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

#### load dataset ####
df_main_meta <- data.table::fread("../data/4_scRNAseq/GSE156455_metadata_main.tsv", header = T)
df_main_tsne <- data.table::fread("../data/4_scRNAseq/GSE156455_tsne_main.tsv", header = T)
cell_list <- c("R1.6", "L1", "L2", "L3", "L4", "L5", "Tm1", "Tm2", "Tm3", "Tm4", "Tm9", "T4.T5") # "Mi1", "Mi4", "Mi9", 
df_main.ptp99a_2 <- read.table("../data/4_scRNAseq/df_main.ptp99a_2.tsv", header = TRUE)

#### Figure S11a-b ####
##### Figure S11a plot tsne #####
# 色のリスト（カテゴリに対応する色）
group_colors <- c(pals::tol.rainbow(length(cell_list)), "grey")
# 名前付きベクトルを作成（名前がカテゴリ、値が色）
color_mapping <- setNames(group_colors, c(cell_list, "Others"))

tsneplot <- inner_join(df_main_meta, df_main_tsne) %>%
  mutate(color = if_else(type %in% cell_list, type, "Others")) %>%
  transform(color = factor(color, levels = c(cell_list, "Others"))) %>%
  # mutate(type = fct_shuffle(type)) %>%
  # arrange(type) %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2, color = color)) +
  geom_point(shape = 16, alpha = .2, size = .5) +
  # geom_text(stat = 'identity', position = 'identity', check_overlap = TRUE) +
  # scale_colour_manual(values = c(rep(pals::isol(20) %>% sample(), 9), pals::isol(16))) + #pals::kovesi.rainbow
  scale_colour_manual(values = color_mapping, name = "Cell type") + #pals::kovesi.rainbow
  # scale_color_gradientn(pals::kovesi.rainbow()) +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  theme_bw() +
  theme() +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
# tsneplot

##### Figure S11b Ptp99A expression plot #####
g_main_ptp99a_normexp_tsne <- ggplot(df_main.ptp99a_2, 
                                     aes(x = tSNE_1, y = tSNE_2, col = expression)) +
  geom_point(size = .5, alpha = .5, shape = 16) +
  scale_color_gradientn(colours = pals::ocean.amp(100), guide = "colourbar") +
  xlab("tSNE-1") +
  ylab("tSNE-2") +
  theme_bw()
# g_main_ptp99a_normexp_tsne

g_tsne_ptp99a <- tsneplot +
  g_main_ptp99a_normexp_tsne
ggsave("../figures/FigureS11ab.png", g_tsne_ptp99a, w = 9, h = 4)

#### Figure S11c ####
label_data <- data.frame(type = cell_list) %>%
  group_by(type) %>%
  dplyr::summarize(x = 0.5, y = 5) %>%  # 左上の位置を取得
  transform(type = factor(type, levels = cell_list))
g_main_normexp_ptp99a_meanse <- ggplot(df_main.ptp99a_2 %>%
                                         filter(type %in% cell_list) %>%
                                         transform(type = factor(type, levels = cell_list)) %>%
                                         group_by(type, time, genotype2, direction) %>%
                                         dplyr::summarize(mean = mean(expression, na.rm = TRUE),
                                                   sd = sd(expression, na.rm = TRUE),
                                                   n = n(),
                                                   se = sd / sqrt(n),
                                                   pval = mean(pval, na.rm = TRUE)) %>%
                                         dplyr::mutate(star = makeStars(pval)), 
                                       aes(x = time, y = mean, group = genotype2)) +
  geom_path(aes(color = genotype2)) +
  geom_point(aes(color = genotype2), shape = 16, size = 2, position = position_dodge(0.3)) +
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd, color = genotype2), width = 0, position = position_dodge(0.3)) + 
  geom_text(aes(x = time,  y = 4, label = star, color = direction), check_overlap = TRUE) +
  # geom_text(x = 0.8, y = 4, label = "Error bar: ±SE", size = 2, hjust = 0,
  #           data = df_main.ptp99a_2 %>%
  #             filter(type %in% c("L1")) %>%
  #             group_by(type, time, genotype2) %>%
  #             summarize(mean = mean(expression, na.rm = TRUE)),
  #           check_overlap = TRUE) +
  xlab("Time after pupal formation") +
  ylab("Normalized expression") +
  scale_color_manual(values = c("#AA95AC", "#37122D", "#AA95AC", "#37122D")) +
  # scale_fill_manual(values = c("C/C" = viridis(2)[1], "T/T" = viridis(2)[2], 
  #                              "no" = "white", "yes" = "#d3cbc6")) +
  scale_alpha_manual(values = c(0, 0.01)) +
  coord_cartesian(ylim = c(0, 5)) +
  facet_wrap( ~ type) + #nrow = 12, scales = "free"
  geom_text(data = label_data, aes(x = x, y = y, label = type), 
            inherit.aes = FALSE, hjust = -0.1, vjust = 1.1, size = 3) +  # 左上にラベルを配置
  theme_bw() +
  theme(legend.position = "none",
        # legend.title = element_blank(),
        # legend.position = c(0.86, 0.97),
        # # legend.text = element_text(size = 6),
        # legend.key = element_blank(),
        # legend.background = element_blank(),
        # legend.direction = "horizontal",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.text = element_blank(),
        strip.background = element_blank())  # ストリップをパネル内に配置
# g_main_normexp_ptp99a_meanse
ggsave("../figures/FigureS11c.pdf", g_main_normexp_ptp99a_meanse, w = 5, h = 4)


##### write Supplementary Data 9 #####
df_SuppData9_FigS11c <-
  df_main.ptp99a_2 %>%
  filter(type %in% cell_list) %>%
  transform(type = factor(type, levels = cell_list)) %>%
  group_by(type, time, genotype2, direction) %>%
  dplyr::summarize(mean = mean(expression, na.rm = TRUE),
                   sd = sd(expression, na.rm = TRUE),
                   n = n(),
                   se = sd / sqrt(n),
                   pval = mean(pval, na.rm = TRUE)) %>%
  ungroup() %>%
  dplyr::mutate(star = makeStars(pval)) %>%
  dplyr::rename(genotype = genotype2) %>%
  dplyr::mutate(across(everything(), ~ ifelse(is.na(.), "NA", as.character(.))))

openxlsx::write.xlsx(df_SuppData9_FigS11c, "../figures/SuppData9_FigS11c.xlsx")


#### Figure S11d ####
seurat_obj <- readRDS("../data/4_scRNAseq/hdWGCNA/seurat_obj.rds")

pdf(paste0("../figures/FigureS11d.pdf"), w = 6, h = 3)
hdWGCNA::PlotDendrogram(seurat_obj, main = paste0("hdWGCNA Dendrogram"))
dev.off()


#### Figure S11e ####
list_modules <- c("blue", "brown", "green", "turquoise", "yellow")
df_module_enrichment <- data.frame()
for (module in list_modules){
  # module <- "blue"
  tbl_fread <- 
    list.files(paste0("../data/4_scRNAseq/hdWGCNA/enrichment/module_", module, "/"), 
               pattern = "*.tsv", full.names = TRUE, recursive=T) %>% 
    map_df(~data.table::fread(., colClasses = c(ID = "character", 
                                                Description = "character",
                                                GeneRatio = "character",
                                                BgRatio = "character",
                                                geneID = "character",
                                                Type = "character"))) %>%
    mutate(module = module)
  
  df_module_enrichment <-
    bind_rows(df_module_enrichment,
              tbl_fread)
  
}

# enrichr dotplot  
p_enrichr <- df_module_enrichment %>%
  filter(module != "grey") %>%
  group_by(module) %>%
  arrange(p.adjust) %>%
  slice_head(n = 6) %>%
  ungroup() %>%
  # filter(db == "GO_Biological_Process_2021") %>%
  dplyr::mutate(padj_log_d = case_when(-log10(p.adjust) > -log10(0.0001) ~ "****", #Adjusted.P.value
                                       -log10(p.adjust) > -log10(0.001) ~ "***", #Adjusted.P.value
                                       -log10(p.adjust) > -log10(0.01) ~ "**",
                                       -log10(p.adjust) > -log10(0.05) ~ "*",
                                       TRUE ~ NA_character_),
                module = str_to_sentence(module)) %>%
  # transform(module = factor(module, levels = gtools::mixedsort(unique(.$module)))) %>%
  arrange(module) %>%
  transform(Description = factor(Description, levels = unique(.$Description) %>% rev())) %>%
  ggplot(., aes(y = Description,
                x = module)) + #col = padj_log_d
  geom_point(aes(size = padj_log_d, shape = Type, color = module)) +#, size= 4) +
  # scale_x_discrete(labels = scales::label_wrap(40)) +
  # tidytext::scale_y_reordered() +
  scale_color_manual(values = c("#4f5a97", "#82474a", "#7d9f56", "#56b7d9", "#e9bb30"),guide="none") +
  scale_size_discrete(name = expression(paste("Adjusted ", italic(P), "-value"))) +
  scale_shape_manual(name = "", values = c(16, 17, 15), labels = c("Biological process", "Cellular component", "Molecular Function")) +
  # scale_shape_manual(values = c(17, 16, 15), labels = c("Biological process", "Cellular component", "Molecular Function")) +
  xlab("Module") +
  ylab("GO term") +
  # facet_grid(type ~ ., scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank())
p_enrichr
ggsave(paste0("../figures/FigureS11e.pdf"), p_enrichr, w = 6, h = 5)


#### Figure S11f ####
df_ME_wilcox_time <- read.table("../data/4_scRNAseq/hdWGCNA/df_ME_wilcox_time.tsv",
                              header = TRUE)

g_wilcoxon_time <- 
  ggplot(df_ME_wilcox_time %>%
           transform(type = factor(type, levels = cell_list)),
         aes(x = time, y = mean)) +
  geom_path(aes(color = condition, group = condition)) +
  geom_point(aes(color = condition), shape = 16, size = 2, 
             position = position_dodge(0.3)) +
  geom_errorbar(aes(ymax = mean + sd, ymin = mean - sd, color = condition), 
                width = 0, position = position_dodge(0.3)) + 
  geom_text(aes(x = time,  y = 30, label = star, color = direction), check_overlap = TRUE) +
  scale_color_manual(values = c("#AA95AC", "#37122D", "#AA95AC", "#37122D")) +
  xlab("Time") +
  ylab("Eigengene expression") +
  facet_grid(type ~ Module) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "NA")#,
g_wilcoxon_time
ggsave("../figures/FigureS11f.pdf", g_wilcoxon_time, w = 7, h = 7)

##### write Supplementary Data 9 #####
df_SuppData9_FigS11f <-
  df_ME_wilcox_time %>%
  dplyr::mutate(genotype = if_else(condition == "C", "C/T", "T/T")) %>%
  dplyr::mutate(across(everything(), ~ ifelse(is.na(.), "NA", as.character(.)))) %>%
  dplyr::rename(module = Module) %>%
  transform(type = factor(type, levels = cell_list)) %>%
  arrange(module, type, time, condition) %>%
  dplyr::select(!condition) %>%
  dplyr::select(module, type, time, genotype, everything())

openxlsx::write.xlsx(df_SuppData9_FigS11f, "../figures/SuppData9_FigS11f.xlsx")

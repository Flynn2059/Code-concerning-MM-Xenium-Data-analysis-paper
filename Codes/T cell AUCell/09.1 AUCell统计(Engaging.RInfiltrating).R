# ==== 加载配置和工作环境 ====
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/")
getwd()
list.files()
library(qs)
library(scop)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# ==== 读取之前处理好的表达矩阵 ====
seurat_obj=read.csv('ALL_SAMPLES_aucell_results_with_secondary_annotation_version3_new_pathway.csv',row.names = 1)
seurat_obj=seurat_obj[,-1]
seurat_obj=CreateSeuratObject(t(seurat_obj))

# ==== 读取细胞注释和邻域信息 ====
Ctrl=read.csv('./Tcell_Plasma_GMM_clustering_annotation/Ctrl_4_samples_cell_niche_matrix.csv',
              row.names = 1)
MGUS=read.csv('./Tcell_Plasma_GMM_clustering_annotation/MGUS_2_samples_cell_niche_matrix.csv',
              row.names = 1)
SM=read.csv('./Tcell_Plasma_GMM_clustering_annotation/SM_5_samples_cell_niche_matrix.csv',
            row.names = 1)
MM=read.csv('./Tcell_Plasma_GMM_clustering_annotation/MM_10_samples_cell_niche_matrix.csv',
            row.names = 1)
annotation=rbind(Ctrl,MGUS,SM,MM)
rm(Ctrl,MGUS,SM,MM);gc()

# ==== 筛选细胞和注释信息 ====
annotation=annotation[colnames(seurat_obj),c("annotation","sample","niche")]
# 补充回去
seurat_obj@meta.data$sample=annotation$sample
seurat_obj@meta.data$celltype_annotation=annotation$annotation
seurat_obj@meta.data$niche=annotation$niche

# ==== 只保留有Infiltrating/Engaging的样本 ====
table(seurat_obj@meta.data$sample,seurat_obj@meta.data$niche)
seurat_obj=subset(seurat_obj,
                  seurat_obj@meta.data$sample%in%c("GSM9035026_Ctrl-4",
                                                    "GSM9035027_MGUS-1",
                                                    "GSM9035036_MM-8",
                                                    "GSM9035038_MM-10"))
print(seurat_obj)
seurat_obj@meta.data[seurat_obj@meta.data$niche%in%c('Engaging Zone','Infiltrating Zone'),'niche']=
  'Infiltrating/Engaging Zone'

# ==== 开始进行统计 ====
outdir <- "./aucell_score_stat_box_plot_infiltrating+_engaging"
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
for (feature in rownames(seurat_obj)) {
  p <- FeatureStatPlot(
    seurat_obj,
    stat.by = feature,
    group.by = "sample",
    split.by = "niche",
    layer = "counts",
    plot_type = "box"
  )+coord_flip()
  
  safe_feature <- gsub("[/\\\\:*?\"<>|]", "_", feature)
  outfile <- file.path(outdir, paste0(safe_feature, ".svg"))
  
  svglite::svglite(file = outfile, width = 11.76, height = 8.37)
  print(p)
  dev.off()
}

# ==== 对两个score进行统计 ====
library(dplyr)
library(tidyr)
library(purrr)
print(seurat_obj)
colnames(seurat_obj@meta.data)

Score=data.frame(
  row.names = colnames(seurat_obj),
  sample = seurat_obj@meta.data$sample,
  niche = seurat_obj@meta.data$niche
)
expr=GetAssayData(seurat_obj,layer = 'count') %>% as.data.frame() %>% t()
Score=cbind(Score,expr)
head(Score)
table(Score$sample)
table(Score$niche)

# 每个 sample 内做 niche 两两 Wilcoxon
# 需要检验的两个 score 列
score_cols <- c("AUCell-Exhaustion-Score",
                "AUCell-Inflammation.Cytotoxic-Score")

# 转成长表
Score_long <- Score %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "score_type",
    values_to = "score_value"
  )

# 先看每个 sample 内是否同时有两个 niche
niche_count <- Score_long %>%
  distinct(sample, niche) %>%
  count(sample, name = "n_niche")

print(niche_count)

# 分 sample、分 score 做双侧 Wilcoxon
wilcox_res <- Score_long %>%
  group_by(sample, score_type) %>%
  group_modify(~{
    dat <- .x
    
    # 只保留有至少两个 niche 组的数据
    if (length(unique(dat$niche)) < 2) {
      return(tibble(
        group1 = NA_character_,
        group2 = NA_character_,
        n1 = NA_integer_,
        n2 = NA_integer_,
        median1 = NA_real_,
        median2 = NA_real_,
        p_value = NA_real_
      ))
    }
    
    grp <- unique(dat$niche)
    g1 <- grp[1]
    g2 <- grp[2]
    
    x1 <- dat %>% filter(niche == g1) %>% pull(score_value)
    x2 <- dat %>% filter(niche == g2) %>% pull(score_value)
    
    wt <- wilcox.test(
      x1, x2,
      alternative = "two.sided",
      exact = FALSE
    )
    
    tibble(
      group1 = g1,
      group2 = g2,
      n1 = length(x1),
      n2 = length(x2),
      median1 = median(x1, na.rm = TRUE),
      median2 = median(x2, na.rm = TRUE),
      p_value = wt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  )

print(wilcox_res)

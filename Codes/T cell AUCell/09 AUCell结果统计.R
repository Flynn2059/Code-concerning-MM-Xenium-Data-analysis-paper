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
print(seurat_obj)
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

# ==== 只保留有T_dominant和Plasma_dominant的样本 ====
table(seurat_obj@meta.data$sample,seurat_obj@meta.data$niche)
seurat_obj=subset(seurat_obj,
                  !seurat_obj@meta.data$sample%in%c("GSM9035026_Ctrl-4",
                                                    "GSM9035027_MGUS-1",
                                                    "GSM9035036_MM-8",
                                                    "GSM9035038_MM-10"))
print(seurat_obj)

# ==== 开始进行统计 ====
seurat_sm_mm=subset(seurat_obj,
                    seurat_obj@meta.data$sample%in%c("GSM9035029_MM-1",
                    "GSM9035030_MM-2","GSM9035031_MM-3","GSM9035032_MM-4",
                    "GSM9035033_MM-5",
                    "GSM9035034_MM-6",
                    "GSM9035035_MM-7",
                    "GSM9035037_MM-9",
                    "GSM9035039_SM-1",
                    "GSM9035040_SM-2",
                    "GSM9035041_SM-3",
                    "GSM9035042_SM-4",
                    "GSM9035043_SM-5"))
seurat_ctrl_mgus=subset(seurat_obj,seurat_obj@meta.data$sample%in%c(
  "GSM9035023_Ctrl-1",
  "GSM9035024_Ctrl-2",
  "GSM9035025_Ctrl-3",
  "GSM9035028_MGUS-2"
))

outdir <- "./aucell_score_stat_box_plot_SM_MM"

for (feature in rownames(seurat_sm_mm)) {
  p <- FeatureStatPlot(
    seurat_sm_mm,
    stat.by = feature,
    group.by = "niche",
    split.by = "sample",
    layer = "counts",
    plot_type = "box"
  )
  
  safe_feature <- gsub("[/\\\\:*?\"<>|]", "_", feature)
  outfile <- file.path(outdir, paste0(safe_feature, ".png"))
  
  png(filename = outfile, width = 1176*2, height = 837*2, units = "px", res = 300)
  print(p)
  dev.off()
}

outdir <- "./aucell_score_stat_box_plot"
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}
for (feature in rownames(seurat_obj)) {
  p <- FeatureStatPlot(
    seurat_obj,
    stat.by = feature,
    group.by = "niche",
    split.by = "sample",
    layer = "counts",
    plot_type = "box"
  )
  
  safe_feature <- gsub("[/\\\\:*?\"<>|]", "_", feature)
  outfile <- file.path(outdir, paste0(safe_feature, ".png"))
  
  png(filename = outfile, width = 1176*3, height = 837*2, units = "px", res = 300)
  print(p)
  dev.off()
}

# ==== 统计学三组检验 ====
print(seurat_obj)

expr=GetAssayData(seurat_obj,layer = 'counts') %>% as.data.frame() %>% t()
metadata=seurat_obj@meta.data[,c('sample','niche')]

data=cbind(expr,metadata)
head(data)
table(data$sample)
table(data$niche)

library(dplyr)
library(tidyr)
library(purrr)

score_cols <- c(
  "AUCell-Exhaustion-Score",
  "AUCell-Inflammation.Cytotoxic-Score"
)

# 1. 转成长表
data_long <- data %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "score_type",
    values_to = "score_value"
  )

# 2. 分 sample、分 niche 统计两个 score
summary_res <- data_long %>%
  group_by(sample, niche, score_type) %>%
  summarise(
    n = n(),
    mean = mean(score_value, na.rm = TRUE),
    sd = sd(score_value, na.rm = TRUE),
    median = median(score_value, na.rm = TRUE),
    IQR = IQR(score_value, na.rm = TRUE),
    min = min(score_value, na.rm = TRUE),
    max = max(score_value, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_res)


# 3. 每个 sample 内，比较不同 niche 的 score 是否有差异
kw_res <- data_long %>%
  group_by(sample, score_type) %>%
  group_modify(~{
    dat <- .x %>%
      filter(!is.na(score_value), !is.na(niche))
    
    n_groups <- dplyr::n_distinct(dat$niche)
    
    # 如果某个 sample 只有 1 个 niche，无法做 Kruskal-Wallis
    if (n_groups < 2) {
      return(tibble(
        n_groups = n_groups,
        statistic = NA_real_,
        df = NA_real_,
        p_value = NA_real_
      ))
    }
    
    kt <- kruskal.test(score_value ~ niche, data = dat)
    
    tibble(
      n_groups = n_groups,
      statistic = unname(kt$statistic),
      df = unname(kt$parameter),
      p_value = kt$p.value
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj_BH = p.adjust(p_value, method = "BH")
  )

print(kw_res)

write.csv(kw_res,file='AUCell_T_2_feature_Kruskal_Wallis_Test.csv')

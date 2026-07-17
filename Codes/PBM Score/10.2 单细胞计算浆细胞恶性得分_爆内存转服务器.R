# ==== 加载配置和工作环境 ====
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("/mnt/bioSSD/XinYe_Li/")
getwd()
list.files()
library(qs)
library(scop)
library(dplyr)
library(purrr)
library(tidyr)
library(Seurat)
library(ggpubr)
library(ggthemes)
library(cowplot)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# ==== 读取之前处理好的数据 ====
seurat_obj=qread("all_sample_plasma_cells.qs")
print(seurat_obj)

# ==== 读取细胞注释和邻域信息 ====
Ctrl=read.csv('./Ctrl_4_samples_cell_niche_matrix.csv',
              row.names = 1)
MGUS=read.csv('./MGUS_2_samples_cell_niche_matrix.csv',
              row.names = 1)
SM=read.csv('./SM_5_samples_cell_niche_matrix.csv',
            row.names = 1)
MM=read.csv('./MM_10_samples_cell_niche_matrix.csv',
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

# ==== 读取plasma cell malignancy score的信息 ====
pathway=readxl::read_excel('./plasma_malignancy_gene_one_hot.xlsx')
pathway=tibble::column_to_rownames(pathway,var="Gene")
pathway=rownames(pathway)[pathway$`84_Plasma cells`=="1"]
common_genes=intersect(pathway,rownames(seurat_obj))

# ==== 类富集分析前准备 ====
colnames(seurat_obj@meta.data)
table(seurat_obj@meta.data$sample)
# 按照sample✖️niche去进行统计
seurat_obj@meta.data$group=paste0(seurat_obj@meta.data$sample,
                                  "_",
                                  seurat_obj@meta.data$niche)
# ==== 单细胞水平的类富集分析 ====
library(Seurat)
library(Matrix)

## -------------------------------------------------
## 单细胞版本的 unweighted running enrichment score
## -------------------------------------------------
calc_cell_es_unweighted <- function(cell_vec, gene_set) {
  cell_vec <- cell_vec[!is.na(cell_vec)]
  cell_vec <- sort(cell_vec, decreasing = TRUE)
  
  gene_set <- intersect(unique(gene_set), names(cell_vec))
  N <- length(cell_vec)
  Nh <- length(gene_set)
  
  if (Nh < 2 || Nh >= N) {
    return(NA_real_)
  }
  
  hits <- names(cell_vec) %in% gene_set
  
  Phit <- numeric(N)
  Pmiss <- numeric(N)
  
  ## 非加权 gene set
  Phit[hits] <- 1 / Nh
  Pmiss[!hits] <- 1 / (N - Nh)
  
  runningES <- cumsum(Phit) - cumsum(Pmiss)
  
  maxES <- max(runningES)
  minES <- min(runningES)
  
  if (abs(maxES) >= abs(minES)) maxES else minES
}

## -------------------------------------------------
## 单细胞版本：随机基因集置换归一化
## -------------------------------------------------
calc_cell_es_normalized <- function(cell_vec, gene_set, nperm = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  cell_vec <- cell_vec[!is.na(cell_vec)]
  gene_universe <- names(cell_vec)
  gene_set <- intersect(unique(gene_set), gene_universe)
  
  if (length(gene_set) < 2) {
    return(c(raw_es = NA_real_, norm_es = NA_real_))
  }
  
  raw_es <- calc_cell_es_unweighted(cell_vec, gene_set)
  
  perm_es <- replicate(nperm, {
    perm_set <- sample(gene_universe, length(gene_set), replace = FALSE)
    calc_cell_es_unweighted(cell_vec, perm_set)
  })
  
  if (is.na(raw_es)) {
    norm_es <- NA_real_
  } else if (raw_es >= 0) {
    denom <- mean(perm_es[perm_es >= 0], na.rm = TRUE)
    norm_es <- raw_es / denom
  } else {
    denom <- abs(mean(perm_es[perm_es < 0], na.rm = TRUE))
    norm_es <- raw_es / denom
  }
  
  c(raw_es = raw_es, norm_es = norm_es)
}

run_pbm_singlecell <- function(
    seu,
    pb_genes,
    assay = "RNA",
    slot = "data",
    center_by_gene = TRUE,
    nperm = 0,
    seed = 123,
    prefix = "PBM"
) {
  stopifnot(inherits(seu, "Seurat"))
  stopifnot(length(pb_genes) > 0)
  
  expr <- GetAssayData(seu, assay = assay, slot = slot)
  
  ## 转成普通 matrix，数据太大时要注意内存
  expr <- as.matrix(expr)
  
  ## 基因层面减中位数：把原文“跨样本”改成“跨细胞”
  if (center_by_gene) {
    gene_medians <- apply(expr, 1, median, na.rm = TRUE)
    expr <- sweep(expr, 1, gene_medians, FUN = "-")
  }
  
  pb_genes_use <- intersect(unique(pb_genes), rownames(expr))
  message("Signature genes found in matrix: ", length(pb_genes_use))
  
  if (length(pb_genes_use) < 10) {
    warning("可用 signature 基因太少，结果可能不稳定。")
  }
  
  nc <- ncol(expr)
  
  if (nperm <= 0) {
    raw_scores <- vapply(seq_len(nc), function(i) {
      v <- expr[, i]
      names(v) <- rownames(expr)
      calc_cell_es_unweighted(v, pb_genes_use)
    }, numeric(1))
    
    seu[[paste0(prefix, "_rank_raw")]] <- raw_scores
    return(seu)
  }
  
  score_mat <- t(vapply(seq_len(nc), function(i) {
    v <- expr[, i]
    names(v) <- rownames(expr)
    calc_cell_es_normalized(
      cell_vec = v,
      gene_set = pb_genes_use,
      nperm = nperm,
      seed = seed + i
    )
  }, numeric(2)))
  
  seu[[paste0(prefix, "_rank_raw")]] <- score_mat[, "raw_es"]
  seu[[paste0(prefix, "_rank_norm")]] <- score_mat[, "norm_es"]
  
  seu
}
seurat_obj=NormalizeData(seurat_obj,scale.factor = median(seurat_obj@meta.data$nCount_Xenium))
seurat_tmp <- run_pbm_singlecell(
  seu = seurat_obj,
  pb_genes = common_genes,
  assay = "Xenium",
  slot = "data",
  center_by_gene = TRUE,
  nperm = 1000,
  prefix = "PBM"
)
print(seurat_tmp)
qsave(seurat_tmp,file = "Plasma_Cells_All_Samples_single-cell-level_PBM_score_1000_permutation_norm.qs")
seurat_tmp=qread('Plasma_Cells_All_Samples_single-cell-level_PBM_score_1000_permutation_norm.qs')
colnames(seurat_tmp@meta.data)
hist(seurat_tmp@meta.data$PBM_rank_norm)
scop::FeatureStatPlot(seurat_tmp,
                      stat.by = "PBM_rank_norm",
                      plot_type = "box",
                      group.by = "sample",
                      pt.size = 0.01,
                      split.by = "niche")+
  labs(y='Normalized PBM Score')

# ==== 对PBM Score进行统计 ====
PBM_Score=data.frame(
  row.names = colnames(seurat_tmp),
  sample = seurat_tmp@meta.data$sample,
  niche = seurat_tmp@meta.data$niche,
  PBM_Score = seurat_tmp@meta.data$PBM_rank_norm
)
head(PBM_Score)
table(PBM_Score$sample)
table(PBM_Score$niche)

# 每个 sample 内做 niche 两两 Wilcoxon
wilcox_res <- PBM_Score %>%
  group_by(sample) %>%
  group_modify(~{
    dat <- .x
    
    comps <- combn(sort(unique(dat$niche)), 2, simplify = FALSE)
    
    res <- map_dfr(comps, function(cc) {
      x <- dat$PBM_Score[dat$niche == cc[1]]
      y <- dat$PBM_Score[dat$niche == cc[2]]
      
      wt <- wilcox.test(x, y, exact = FALSE)
      
      tibble(
        niche1 = cc[1],
        niche2 = cc[2],
        n1 = length(x),
        n2 = length(y),
        median1 = median(x, na.rm = TRUE),
        median2 = median(y, na.rm = TRUE),
        diff_median = median(x, na.rm = TRUE) - median(y, na.rm = TRUE),
        W = unname(wt$statistic),
        p_value = wt$p.value
      )
    })
    
    res %>%
      mutate(p_adj_BH = p.adjust(p_value, method = "BH"))
  }) %>%
  ungroup()

wilcox_res
write.csv(wilcox_res,file = 'PBM_Score_T_Plamsa_Dominant_Mixed_Wilcoxon_Test_result.csv')
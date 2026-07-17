# ==== 调整自己喜欢的系统设置 ====
Sys.setenv(LANGUAGE="en")
options(StringAsFactors=FALSE)
rm(list=ls());gc()
setwd("/Volumes/FlynnDisk/XinyeLi/GSE299207_RAW/")
getwd()
list.files()
# ==== 加载我们要用到的包 ====
library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(arrow)
library(qs)
library(RColorBrewer)
library(spacexr)

# ==== Control样本1 ====
ref=readRDS("./RCTD_ref/Healthy_Donor(Control).rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035023_Ctrl-1/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",layer = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_control1.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()


# ==== Control样本2 ====
ref=readRDS("./RCTD_ref/Healthy_Donor(Control).rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035024_Ctrl-2",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_control2.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()


# ==== Control样本3 ====
ref=readRDS("./RCTD_ref/Healthy_Donor(Control).rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035025_Ctrl-3/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",layer = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_control3.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()


# ==== Control样本4 ====
ref=readRDS("./RCTD_ref/Healthy_Donor(Control).rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035026_Ctrl-4/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",layer = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_control4.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()


# ==== MGUS样本1 ====
ref=readRDS("./RCTD_ref/MGUS.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035027_MGUS-1/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MGUS1.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()


# ==== MGUS样本2 ====
ref=readRDS("./RCTD_ref/MGUS.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035028_MGUS-2/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MGUS2.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== MM样本1 ====
ref=readRDS("./RCTD_ref/MM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035029_MM-1/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MM1.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== MM样本2 ====
ref=readRDS("./RCTD_ref/MM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035030_MM-2/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MM2.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== MM样本3 ====
ref=readRDS("./RCTD_ref/MM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035031_MM-3/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MM3.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== MM样本4 ====
ref=readRDS("./RCTD_ref/MM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035032_MM-4/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MM4.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== MM样本5 ====
ref=readRDS("./RCTD_ref/MM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035033_MM-5/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MM5.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== MM样本6 ====
ref=readRDS("./RCTD_ref/MM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035034_MM-6/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MM6.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== MM样本7 ====
ref=readRDS("./RCTD_ref/MM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035035_MM-7/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MM7.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== MM样本8 ====
ref=readRDS("./RCTD_ref/MM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035036_MM-8/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MM8.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== MM样本9 ====
ref=readRDS("./RCTD_ref/MM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035037_MM-9/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MM9.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== MM样本10 ====
ref=readRDS("./RCTD_ref/MM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035038_MM-10/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_MM10.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== SM样本1 ====
ref=readRDS("./RCTD_ref/SM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035039_SM-1/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_SM1.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()


# ==== SM样本2 ====
ref=readRDS("./RCTD_ref/SM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035040_SM-2/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_SM2.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()


# ==== SM样本3 ====
ref=readRDS("./RCTD_ref/SM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035040_SM-2/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_SM3.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()


# ==== SM样本4 ====
ref=readRDS("./RCTD_ref/SM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035042_SM-4/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_SM4.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

# ==== SM样本5 ====
ref=readRDS("./RCTD_ref/SM.rds")
print(ref)
table(ref@meta.data$celltype)
# 统计每个celltype的细胞数
celltype_counts <- table(ref@meta.data$celltype)
# 找到 >=25 个细胞的celltype
keep_celltypes <- names(celltype_counts[celltype_counts >= 25])
# subset Seurat对象
ref <- subset(ref, subset = celltype %in% keep_celltypes)
# 再检查一下结果
table(ref@meta.data$celltype)
query=LoadXenium("GSM9035043_SM-5/",flip.xy = T)
print(query)
# 首先分别创建参考和质询数据
# 参考数据
counts=GetAssayData(ref,assay="RNA",layer="counts")
cluster=ref@meta.data$celltype %>% as.factor()
names(cluster)=colnames(ref)
nUMI=colSums(GetAssayData(ref,layer = "counts")) %>% as.integer()
names(nUMI)=colnames(ref)
reference=spacexr::Reference(counts = counts,
                             cell_types = cluster,
                             nUMI = nUMI)
rm(counts,cluster,nUMI,ref);gc()
# 我们空间转录组数据
counts=GetAssayData(query,assay = "Xenium",slot = "counts")
coordinate=GetTissueCoordinates(query)
coordinate=coordinate[,c(1,2)]
rownames(coordinate)=colnames(query)
nUMI=colSums(counts)
names(nUMI)=colnames(query)
Spatial=spacexr::SpatialRNA(counts = counts,
                            coords = coordinate,
                            nUMI = nUMI)
rm(counts,coordinate,nUMI,query);gc()
# 可以开始runRCTD了
RCTD=spacexr::create.RCTD(Spatial,reference,max_cores = 4)
rm(reference,Spatial);gc()
RCTD=spacexr::run.RCTD(RCTD,doublet_mode = "doublet")

# 提取RCTD结果的信息
qsave(RCTD,file = "./RCTD/doublet_result_SM5.qs")
# 清空环境，运行下一个样本
rm(list=ls());gc()

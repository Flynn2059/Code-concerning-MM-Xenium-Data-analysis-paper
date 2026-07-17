# ==== 调整自己喜欢的系统设置 ====
Sys.setenv(LANGUAGE="en")
options(StringAsFactors=FALSE)
rm(list=ls());gc()
setwd("/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/")
getwd()
list.files()
# ==== 加载我们要用到的包 ====
library(Seurat)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyverse)
library(qs)
library(RColorBrewer)
library(spacexr)
library(clusterProfiler)
library(AUCell)

# ==== 首先处理细胞 ====
# 正常组
# 加载细胞注释数据并填回
Plasma_cells=read.csv("./4_control_major_annotation.csv",row.names = 1)
table(Plasma_cells$major_annotation)
Plasma_cells=Plasma_cells[Plasma_cells$major_annotation%in%c("Plasma"),,drop=F]
head(rownames(Plasma_cells),5)

# 创建数据框并处理
df <- tibble(original_id = rownames(Plasma_cells)) %>%
  mutate(
    group = str_extract(original_id, "Ctrl-\\d+"),
    cell_id = str_remove(original_id, "^.*Ctrl-\\d+_")
  ) %>%
  select(cell_id, group)

ctrl1=LoadXenium("GSM9035023_Ctrl-1/",flip.xy = T)
ctrl1=subset(ctrl1,cells = df$cell_id[df$group=="Ctrl-1"])
colnames(ctrl1)=paste0("GSM9035023_Ctrl-1_",colnames(ctrl1))

ctrl2=LoadXenium("GSM9035024_Ctrl-2/",flip.xy = T)
ctrl2=subset(ctrl2,cells = df$cell_id[df$group=="Ctrl-2"])
colnames(ctrl2)=paste0("GSM9035024_Ctrl-2_",colnames(ctrl2))

ctrl3=LoadXenium("GSM9035025_Ctrl-3/",flip.xy = T)
ctrl3=subset(ctrl3,cells = df$cell_id[df$group=="Ctrl-3"])
colnames(ctrl3)=paste0("GSM9035025_Ctrl-3_",colnames(ctrl3))

ctrl4=LoadXenium("GSM9035026_Ctrl-4/",flip.xy = T)
ctrl4=subset(ctrl4,cells = df$cell_id[df$group=="Ctrl-4"])
colnames(ctrl4)=paste0("GSM9035026_Ctrl-4_",colnames(ctrl4))

seurat_list=list(ctrl1,ctrl2,ctrl3,ctrl4)
seurat_ctrl=Reduce(merge,seurat_list)
print(seurat_ctrl)

rm(seurat_list,ctrl1,ctrl2,ctrl3,ctrl4,df);gc()
seurat_ctrl[["Xenium"]]=JoinLayers(seurat_ctrl[["Xenium"]])
seurat_ctrl=AddMetaData(seurat_ctrl,metadata = Plasma_cells)

# MGUS组
# 加载细胞注释数据并填回
Plasma_cells=read.csv("2_mgus_major_annotation.csv",row.names = 1)
table(Plasma_cells$major_annotation)
Plasma_cells=Plasma_cells[Plasma_cells$major_annotation%in%c("Plasma"),,drop=F]
head(rownames(Plasma_cells),5)

# 创建数据框并处理
df <- tibble(original_id = rownames(Plasma_cells)) %>%
  mutate(
    group = str_extract(original_id, "MGUS-\\d+"),
    cell_id = str_remove(original_id, "^.*MGUS-\\d+_")
  ) %>%
  select(cell_id, group)

mgus1=LoadXenium("GSM9035027_MGUS-1/",flip.xy = T)
mgus1=subset(mgus1,cells = df$cell_id[df$group=="MGUS-1"])
colnames(mgus1)=paste0("GSM9035027_MGUS-1_",colnames(mgus1))

mgus2=LoadXenium("GSM9035028_MGUS-2/",flip.xy = T)
mgus2=subset(mgus2,cells = df$cell_id[df$group=="MGUS-2"])
colnames(mgus2)=paste0("GSM9035028_MGUS-2_",colnames(mgus2))

seurat_list=list(mgus1,mgus2)
seurat_mgus=Reduce(merge,seurat_list)
print(seurat_mgus)

rm(seurat_list,mgus1,mgus2,df);gc()
seurat_mgus[["Xenium"]]=JoinLayers(seurat_mgus[["Xenium"]])
seurat_mgus=AddMetaData(seurat_mgus,metadata = Plasma_cells)

# SM组
# 加载细胞注释数据并填回
Plasma_cells=read.csv("5_sm_major_annotation.csv",row.names = 1)
table(Plasma_cells$major_annotation)
Plasma_cells=Plasma_cells[Plasma_cells$major_annotation%in%c("Plasma"),,drop=F]
head(rownames(Plasma_cells),5)

# 创建数据框并处理
df <- tibble(original_id = rownames(Plasma_cells)) %>%
  mutate(
    group = str_extract(original_id, "SM-\\d+"),
    cell_id = str_remove(original_id, "^.*SM-\\d+_")
  ) %>%
  select(cell_id, group)

sm1=LoadXenium("GSM9035039_SM-1/",flip.xy = T)
sm1=subset(sm1,cells = df$cell_id[df$group=="SM-1"])
colnames(sm1)=paste0("GSM9035039_SM-1_",colnames(sm1))

sm2=LoadXenium("GSM9035040_SM-2/",flip.xy = T)
sm2=subset(sm2,cells = df$cell_id[df$group=="SM-2"])
colnames(sm2)=paste0("GSM9035040_SM-2_",colnames(sm2))

sm3=LoadXenium("GSM9035041_SM-3/",flip.xy = T)
sm3=subset(sm3,cells = df$cell_id[df$group=="SM-3"])
colnames(sm3)=paste0("GSM9035041_SM-3_",colnames(sm3))

sm4=LoadXenium("GSM9035042_SM-4",flip.xy = T)
sm4=subset(sm4,cells = df$cell_id[df$group=="SM-4"])
colnames(sm4)=paste0("GSM9035042_SM-4_",colnames(sm4))

sm5=LoadXenium("GSM9035043_SM-5/",flip.xy = T)
sm5=subset(sm5,cells = df$cell_id[df$group=="SM-5"])
colnames(sm5)=paste0("GSM9035043_SM-5_",colnames(sm5))

seurat_list=list(sm1,sm2,sm3,sm4,sm5)
seurat_sm=Reduce(merge,seurat_list)
print(seurat_sm)

rm(seurat_list,sm1,sm2,sm3,sm4,sm5,df);gc()
seurat_sm[["Xenium"]]=JoinLayers(seurat_sm[["Xenium"]])
seurat_sm=AddMetaData(seurat_sm,metadata = Plasma_cells)

# MM组
# 加载细胞注释数据并填回
Plasma_cells=read.csv("10_MM_major_annotation.csv",row.names = 1)
Plasma_cells[Plasma_cells$major_annotation=="T Proliferation","major_annotation"]="T Proliferating"
table(Plasma_cells$major_annotation)
Plasma_cells=Plasma_cells[Plasma_cells$major_annotation%in%c("Malignant_Plasma"),,drop=F]
head(rownames(Plasma_cells),5)

# 创建数据框并处理
df <- tibble(original_id = rownames(Plasma_cells)) %>%
  mutate(
    group = str_extract(original_id, "MM-\\d+"),
    cell_id = str_remove(original_id, "^.*MM-\\d+_")
  ) %>%
  select(cell_id, group)

MM1=LoadXenium("GSM9035029_MM-1/",flip.xy = T)
MM1=subset(MM1,cells = df$cell_id[df$group=="MM-1"])
colnames(MM1)=paste0("GSM9035029_MM-1_",colnames(MM1))

MM2=LoadXenium("GSM9035030_MM-2/",flip.xy = T)
MM2=subset(MM2,cells = df$cell_id[df$group=="MM-2"])
colnames(MM2)=paste0("GSM9035030_MM-2_",colnames(MM2))

MM3=LoadXenium("GSM9035031_MM-3/",flip.xy = T)
MM3=subset(MM3,cells = df$cell_id[df$group=="MM-3"])
colnames(MM3)=paste0("GSM9035031_MM-3_",colnames(MM3))

MM4=LoadXenium("GSM9035032_MM-4/",flip.xy = T)
MM4=subset(MM4,cells = df$cell_id[df$group=="MM-4"])
colnames(MM4)=paste0("GSM9035032_MM-4_",colnames(MM4))

MM5=LoadXenium("GSM9035033_MM-5/",flip.xy = T)
MM5=subset(MM5,cells = df$cell_id[df$group=="MM-5"])
colnames(MM5)=paste0("GSM9035033_MM-5_",colnames(MM5))

MM6=LoadXenium("GSM9035034_MM-6/",flip.xy = T)
MM6=subset(MM6,cells = df$cell_id[df$group=="MM-6"])
colnames(MM6)=paste0("GSM9035034_MM-6_",colnames(MM6))

MM7=LoadXenium("GSM9035035_MM-7/",flip.xy = T)
MM7=subset(MM7,cells = df$cell_id[df$group=="MM-7"])
colnames(MM7)=paste0("GSM9035035_MM-7_",colnames(MM7))

MM8=LoadXenium("GSM9035036_MM-8/",flip.xy = T)
MM8=subset(MM8,cells = df$cell_id[df$group=="MM-8"])
colnames(MM8)=paste0("GSM9035036_MM-8_",colnames(MM8))

MM9=LoadXenium("GSM9035037_MM-9/",flip.xy = T)
MM9=subset(MM9,cells = df$cell_id[df$group=="MM-9"])
colnames(MM9)=paste0("GSM9035037_MM-9_",colnames(MM9))

MM10=LoadXenium("GSM9035038_MM-10/",flip.xy = T)
MM10=subset(MM10,cells = df$cell_id[df$group=="MM-10"])
colnames(MM10)=paste0("GSM9035038_MM-10_",colnames(MM10))


seurat_list=list(MM1,MM2,MM3,MM4,MM5,MM6,MM7,MM8,MM9,MM10)
seurat_mm=Reduce(merge,seurat_list)
print(seurat_mm)

rm(seurat_list,MM1,MM2,MM3,MM4,MM5,MM6,MM7,MM8,MM9,MM10,df);gc()
seurat_mm[["Xenium"]]=JoinLayers(seurat_mm[["Xenium"]])
seurat_mm=AddMetaData(seurat_mm,metadata = Plasma_cells)
print(seurat_mm)
# ==== 最后一轮整合 ====
seurat_list=list(seurat_ctrl,seurat_mgus,seurat_sm,seurat_mm)
seurat_obj=Reduce(merge,seurat_list)
print(seurat_obj)
rm(seurat_list,seurat_ctrl,seurat_mgus,seurat_sm,seurat_mm,Plasma_cells);gc()
seurat_obj[["Xenium"]]=JoinLayers(seurat_obj[["Xenium"]])
DefaultAssay(seurat_obj)
qsave(seurat_obj,file="all_sample_plasma_cells.qs")

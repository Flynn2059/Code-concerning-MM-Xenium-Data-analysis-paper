#%% imports
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata as sd
import spatialdata_plot
import matplotlib.pyplot as plt

# SpatialData stack
from spatialdata import SpatialData
from spatialdata.models import TableModel, ShapesModel, PointsModel, Image2DModel

import tifffile as tiff
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon

# some customized settings for better visulization
plt.rcParams["figure.figsize"] = (12, 12)
plt.rcParams["figure.dpi"] = 600
#%% 工作路径
os.chdir('/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/')
print(os.getcwd())
print("\n")
print(os.listdir())

#%% Control
adata=sc.read_h5ad('./4_control_major_annotation.h5ad')
print(adata)
print(adata.obs['sample'].value_counts())
# 我们读取之前的细胞注释结果
annotation_level1=pd.read_csv("./4_control_major_annotation.csv",index_col=0)
annotation_level2=pd.read_csv("./4_control_T_secondary_annotation.csv",index_col=0)
# 把 secondary_annotation 对应到 major_annotation 的同名细胞上
annotation_level1.loc[annotation_level2.index, "major_annotation"] = annotation_level2["secondary_annotation"]
annotation_level1['major_annotation'] = annotation_level1['major_annotation'].replace({
    'Mono/Macro': 'Myeloid',
})
# 检查一下
print(annotation_level1['major_annotation'].value_counts())
adata.obs['celltype_secondary']=annotation_level1.loc[adata.obs_names,"major_annotation"]
#sc.pl.umap(adata,color=['celltype',"celltype_secondary"])    
# 提取我们感兴趣的信息
data = pd.DataFrame(
    data={
        'X': adata.obsm['spatial'][:, 0],
        'Y': adata.obsm['spatial'][:, 1],
        'sample': adata.obs['sample'].values,
        'annotation': adata.obs['celltype_secondary'].values  # 建议加上 .values 转换为 numpy 数组
    },
    index=adata.obs_names
)
# 过滤掉 SOX2-OT+ Unknown 的行
data = data[data['annotation'] != 'SOX2-OT+ Unknown']
data.to_csv('4_ctrl_coordinate_annotation_matrix.csv')

#%% MGUS
adata=sc.read_h5ad('2_mgus_major_annotation.h5ad')
print(adata)
print(adata.obs['sample'].value_counts())
# 我们读取之前的细胞注释结果
annotation_level1=pd.read_csv("./2_mgus_major_annotation.csv",index_col=0)
annotation_level2=pd.read_csv("./2_mgus_T_secondary_annotation.csv",index_col=0)
# 把 secondary_annotation 对应到 major_annotation 的同名细胞上
annotation_level1.loc[annotation_level2.index, "major_annotation"] = annotation_level2["secondary_annotation"]
annotation_level1['major_annotation'] = annotation_level1['major_annotation'].replace({
    'Monocyte': 'Myeloid',
    'Macrophage': 'Myeloid',
})

# 检查一下
print(annotation_level1['major_annotation'].value_counts())
adata.obs['celltype_secondary']=annotation_level1.loc[adata.obs_names,"major_annotation"]
#sc.pl.umap(adata,color=['celltype',"celltype_secondary"])    
# 提取我们感兴趣的信息
data = pd.DataFrame(
    data={
        'X': adata.obsm['spatial'][:, 0],
        'Y': adata.obsm['spatial'][:, 1],
        'sample': adata.obs['sample'].values,
        'annotation': adata.obs['celltype_secondary'].values  # 建议加上 .values 转换为 numpy 数组
    },
    index=adata.obs_names
)
# 过滤掉 SOX2-OT+ Unknown 的行
data = data[data['annotation'] != 'SOX2-OT Unknown']
data.to_csv('2_mgus_coordinate_annotation_matrix.csv')
#%% SM
adata=sc.read_h5ad('5_sm_major_annotation.h5ad')
print(adata)
print(adata.obs['sample'].value_counts())
annotation_level1=pd.read_csv("./5_sm_major_annotation.csv",index_col=0)
annotation_level2=pd.read_csv("./5_sm_T_secondary_annotation.csv",index_col=0)
# 把 secondary_annotation 对应到 major_annotation 的同名细胞上
annotation_level1.loc[annotation_level2.index, "major_annotation"] = annotation_level2["secondary_annotation"]
# 细胞名有重复，重新调整一下
annotation_level1['major_annotation'] = annotation_level1['major_annotation'].replace({
    'SOX2-OT+ Unknown': 'SOX2-OT Unknown',
    'Mono/Macro': 'Myeloid'
})
# 检查一下
print(annotation_level1['major_annotation'].value_counts())
adata.obs['celltype_secondary']=annotation_level1.loc[adata.obs_names,"major_annotation"]
#sc.pl.umap(adata,color=['celltype',"celltype_secondary"])    
# 提取我们感兴趣的信息
data = pd.DataFrame(
    data={
        'X': adata.obsm['spatial'][:, 0],
        'Y': adata.obsm['spatial'][:, 1],
        'sample': adata.obs['sample'].values,
        'annotation': adata.obs['celltype_secondary'].values  # 建议加上 .values 转换为 numpy 数组
    },
    index=adata.obs_names
)
# 过滤掉 SOX2-OT+ Unknown 的行
data = data[data['annotation'] != 'SOX2-OT Unknown']
data.to_csv('5_sm_coordinate_annotation_matrix.csv')
#%% MM
adata=sc.read_h5ad('10_mm_major_annotation.h5ad')
print(adata)
print(adata.obs['sample'].value_counts())

annotation_level1=pd.read_csv("./10_mm_major_annotation.csv",index_col=0)
annotation_level2=pd.read_csv("./10_mm_T_secondary_annotation.csv",index_col=0)
# 把 secondary_annotation 对应到 major_annotation 的同名细胞上
annotation_level1.loc[annotation_level2.index, "major_annotation"] = annotation_level2["secondary_annotation"]
# 细胞名有重复，重新调整一下
annotation_level1['major_annotation'] = annotation_level1['major_annotation'].replace({
    'SOX2-OT+ Unknown': 'SOX2-OT Unknown',
    'Mono/Macro': 'Myeloid',
    'T Proliferation': 'T Proliferating'
})
# 检查一下
print(annotation_level1['major_annotation'].value_counts())

adata.obs['celltype_secondary']=annotation_level1.loc[adata.obs_names,"major_annotation"]
#sc.pl.umap(adata,color=['celltype',"celltype_secondary"])    
# 提取我们感兴趣的信息
data = pd.DataFrame(
    data={
        'X': adata.obsm['spatial'][:, 0],
        'Y': adata.obsm['spatial'][:, 1],
        'sample': adata.obs['sample'].values,
        'annotation': adata.obs['celltype_secondary'].values  # 建议加上 .values 转换为 numpy 数组
    },
    index=adata.obs_names
)
# 过滤掉 SOX2-OT+ Unknown 的行
data = data[data['annotation'] != 'SOX2-OT Unknown']
data.to_csv('10_mm_coordinate_annotation_matrix.csv')
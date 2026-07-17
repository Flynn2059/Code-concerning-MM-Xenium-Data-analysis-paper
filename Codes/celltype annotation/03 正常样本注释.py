#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 19 18:56:45 2026

@author: flynn
"""

#%% imports
import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import spatialdata as sd
import matplotlib.pyplot as plt

# SpatialData stack
from spatialdata import SpatialData
from spatialdata.models import TableModel, ShapesModel, PointsModel, Image2DModel

import tifffile as tiff
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon
from sklearn.neighbors import KernelDensity

# some customized settings for better visulization
plt.rcParams["figure.figsize"] = (12, 12)
plt.rcParams["figure.dpi"] = 300
#%% 工作路径
os.chdir('/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/')
print(os.getcwd())
print("\n")
print(os.listdir())

#%% 读取数据
# collect all zarr folders
zarr_files = [
    "GSM9035023_Ctrl-1_post_qc.zarr",
    "GSM9035024_Ctrl-2_post_qc.zarr",
    "GSM9035025_Ctrl-3_post_qc.zarr",
    "GSM9035026_Ctrl-4_post_qc.zarr",
]

for f in zarr_files:
    print(" -", f)

print("\nLoading spatialdata objects...\n")

# dictionaries for storing sdata and adata
sdata_dict = {}
adata_dict = {}


for zarr in zarr_files:
    # sample 名
    sample_name = zarr.replace("_post_qc.zarr", "").replace(".zarr", "")
    prefix = f"{sample_name}_"

    print(f"Reading {zarr} ...")

    # 读取 spatialdata
    sdata = sd.read_zarr(zarr)

    # ===== 1️⃣ 处理 AnnData table =====
    if "table" in sdata:
        adata = sdata["table"].copy()
        adata.obs["sample"]=sample_name

        # obs_names 加前缀
        adata.obs_names = prefix + adata.obs_names.astype(str)

        # 常见 cell_id 列
        for col in ["cell_id", "instance_id", "spot_id"]:
            if col in adata.obs.columns:
                adata.obs[col] = prefix + adata.obs[col].astype(str)

        sdata["table"] = adata
        print(f"  -> Updated AnnData IDs for {sample_name}")
    else:
        print(f"  -> Warning: no 'table' found in {zarr}")

    # ===== 2️⃣ 处理 shapes =====
    if "shapes" in sdata:
        for key in sdata["shapes"]:
            gdf = sdata["shapes"][key].copy()

            # index 加前缀
            gdf.index = prefix + gdf.index.astype(str)

            # 如果有 cell_id 列
            if "cell_id" in gdf.columns:
                gdf["cell_id"] = prefix + gdf["cell_id"].astype(str)

            sdata["shapes"][key] = gdf

        print(f"  -> Updated shapes IDs for {sample_name}")

    # ===== 3️⃣ 处理 points（如 transcripts）=====
    if "points" in sdata:
        for key in sdata["points"]:
            df = sdata["points"][key].copy()

            for col in ["cell_id", "instance_id"]:
                if col in df.columns:
                    df[col] = prefix + df[col].astype(str)

            sdata["points"][key] = df

        print(f"  -> Updated points IDs for {sample_name}")

    # ===== 保存到 dict =====
    sdata_key = f"sdata_{sample_name}"
    adata_key = f"adata_{sample_name}"

    sdata_dict[sdata_key] = sdata

    if "table" in sdata:
        adata_dict[adata_key] = sdata["table"]

print("\nAll files loaded and IDs prefixed.")
print(f"Total sdata objects: {len(sdata_dict)}")
print(f"Total adata objects: {len(adata_dict)}")

# 多样本整合
adata=sc.concat(adata_dict.values(), join="outer")
print(adata)
print(adata.obs['sample'].value_counts())

#%% 降维聚类
# 找高变基因
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000)
# Log Normalization
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()
# Scaling of data, we don't center the data here as the data centering is done in PCA and this keeps the X sparse for memory efficiency.
sc.pp.scale(adata, zero_center=False, max_value=10)
# PCA for dimension reduction further by UMAP for downstream analysis
sc.pp.pca(adata, n_comps=30)
sc.external.pp.harmony_integrate(adata, key="sample",max_iter_harmony=100)
sc.pp.neighbors(adata, metric="cosine",use_rep="X_pca_harmony")
sc.tl.umap(adata)
sc.pl.umap(adata,color="sample")
#%% 聚类与可视化
# clustering with leiden. We use igraph here as it scales better. The resolution is a hyper-parameter. A higher resolution will output more clusters.
sc.tl.leiden(adata, flavor="igraph", n_iterations=-1, resolution=0.2)
sc.pl.umap(adata, color="leiden",legend_loc="on data")

# FindAllMarkers
sc.tl.rank_genes_groups(adata, groupby="leiden", layer="lognorm", pts=True)
sc.tl.dendrogram(adata, groupby="leiden")


# Now we perform cell-type annotation based on marker genes
# We draw a dotplot for clear visulization
marker_genes = {
    "HSC": [
        "CD34",   # 经典 HSC / HSPC 表面标记
        "FLT3",
        "SELL",
        "SPN",
        "GATA2",
        "CD82",
        "CDK6",
        "SOCS2",
        "PBXIP1"
    ],
    "Erythroid": [
        "GYPA",
        "TFRC",   # CD71
        "EPOR"
    ],
    
    "Megakaryocyte": [
        "ITGA2B",
        "ITGB3",
        "GP1BA",
        "THBS1"
    ],
    
    "Eosinophil": [
        "IL5RA",
        "PRG2",
        "CCR3"
    ],
    
    "Basophil": [
        "ENPP3",
        "KIT",
        "HDC",
        "GATA2"
    ],
    
    # Bone Related
    "Osteoblast": [
        "ALPL",
        "RUNX2"
    ],
        
    "Osteoclast": [
        "CTSK",
        "ACP5",
        "CALCR"
    ],
    
    "Adipocyte": [
        "ADIPOQ",
        "PLIN1",
        "LPL"
    ],
    
    
    # T cells
    "general_T":["CD2","CD3E","CD3G","CD247"],
    "CD4_lineage":["CD4","CCR7"],
    "Regulatory":["FOXP3","IL2RA","IKZF2"],
    "CD8_lineage":["CD8A","CD8B","GZMB","GZMH","PRF1"],
    
    # NK cells
    "NK":["KLRB1","KLRC1","KLRD1","KLRK1"],
    
    # B cells
    "B":["CD19","MS4A1","CD79A","CD79B","BANK1"],
    
    # Plasma cells
    "Plsama":["MZB1","XBP1","PRDM1"],
    
    # Monocytes
    "Mono/Macro":["CD14","FCGR3A","CSF1R","CX3CR1","ITGAM"],
    
    # DCs
    "DC":["ITGAX","CD1C","CLEC9A","LILRA4","IRF4","IRF8"],
    
    # Neutrophils
    "Neutrophil":["MPO","ELANE","CXCR2"],
    
    "Mast":["KIT"],
    
    # Endothelium
    "Endothelium":["PECAM1","CDH5","KDR","ENG"],
    "Epithelium":["EPCAM","KRT8"],
    
    # Stromal
    "Fibroblast":["PDGFRA","PDGFRB","THY1","COL5A1"],
    "Pericyte":["CSPG4","RGS5"]
}

panel_genes = set(adata.var_names)

filtered_markers = {
    k: [g for g in v if g in panel_genes]
    for k, v in marker_genes.items()
}


sc.pl.dotplot(adata, 
              filtered_markers,
              groupby="leiden", 
              standard_scale="var",
              cmap='Blues'
             )
sc.pl.rank_genes_groups_dotplot(adata,n_genes=15,cmap="Blues")
#%% 整合RCTD结果
rctd_control_1=pd.read_csv("./RCTD/doublet_result_control1.csv",index_col=0)
rctd_control_2=pd.read_csv("./RCTD/doublet_result_control2.csv",index_col=0)
rctd_control_3=pd.read_csv("./RCTD/doublet_result_control3.csv",index_col=0)
rctd_control_4=pd.read_csv("./RCTD/doublet_result_control4.csv",index_col=0)
rctd_control_1.index = "GSM9035023_Ctrl-1_" + rctd_control_1.index.astype(str)
rctd_control_2.index = "GSM9035024_Ctrl-2_" + rctd_control_2.index.astype(str)
rctd_control_3.index = "GSM9035025_Ctrl-3_" + rctd_control_3.index.astype(str)
rctd_control_4.index = "GSM9035026_Ctrl-4_" + rctd_control_4.index.astype(str)
rctd=pd.concat([rctd_control_1,rctd_control_2,rctd_control_3,rctd_control_4],axis=0)
# 把RCTD的结果写回去
adata.obs["rctd_inferred_celltype"] = (
    adata.obs.index
    .map(rctd["celltype"])
    .fillna("NA")
)
sc.pl.umap(adata,color=['rctd_inferred_celltype'])

#%% 探索一下我们搞不清楚的cluster 7
meta = adata.obs.loc[adata.obs["leiden"] == "7", 
                     ["leiden", "rctd_inferred_celltype"]]
meta["rctd_inferred_celltype"].value_counts()
#%% 整合banksy结果
banksy_ctrl_pca=pd.read_csv("./Ctrl_4_samples_PCA.csv",index_col=0)
banksy_ctrl_umap=pd.read_csv("./Ctrl_4_samples_UMAP.csv",index_col=0)
banksy_ctrl_leiden=pd.read_csv("./Ctrl_4_samples_leiden_clustering.csv",index_col=0)

adata.obsm['X_banksy_pca'] = banksy_ctrl_pca.values.astype(np.float32) 
adata.obsm['X_banksy_umap'] = banksy_ctrl_umap.values.astype(np.float32)
adata.obs['banksy_leiden'] = pd.Categorical(banksy_ctrl_leiden.iloc[:, 0].astype(str))

sc.pl.embedding(adata, basis='banksy_umap', color='banksy_leiden',legend_loc="on data")

def banksy_enhance_umap(adata, color='celltype', figsize=(6,6)):
    fig, ax = plt.subplots(figsize=figsize)
    sc.pl.embedding(
        adata,
        basis='banksy_umap',
        color=color,
        legend_loc=None,
        frameon=False,
        size=15,
        ax=ax,
        add_outline=True,
        show=False
    )

    for cell_type in adata.obs[color].unique():
        mask = adata.obs[color] == cell_type
        x = adata.obsm["X_banksy_umap"][mask, 0]
        y = adata.obsm["X_banksy_umap"][mask, 1]

        x_min, x_max = x.min() - 0.5, x.max() + 0.5
        y_min, y_max = y.min() - 0.5, y.max() + 0.5
        xx, yy = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]

        xy_train = np.vstack([x, y]).T
        xy_test = np.vstack([xx.ravel(), yy.ravel()]).T

        kde = KernelDensity(bandwidth=0.5, metric='euclidean')
        kde.fit(xy_train)

        Z = np.exp(kde.score_samples(xy_test))
        Z = Z.reshape(xx.shape)

    cell_type_counts = adata.obs[color].value_counts(normalize=True) * 100
    for cell_type, percentage in cell_type_counts.items():
        mask = adata.obs[color] == cell_type
        x = np.median(adata.obsm["X_banksy_umap"][mask, 0])
        y = np.median(adata.obsm["X_banksy_umap"][mask, 1])
        ax.text(
            x, y, f"{cell_type}",
            fontsize=10, 
            ha="center", 
            color="black",
            bbox=dict(
                facecolor='white',         
                edgecolor='none', 
                alpha=0.65,
                boxstyle="round,pad=0.4"  
            )  
        )
    ax = plt.gca()
    ax.set_xlabel('UMAP1', fontsize=6)
    ax.set_ylabel('UMAP2', fontsize=6)
 
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()

    x_range = x_max - x_min
    y_range = y_max - y_min
    aspect_ratio = (fig.get_size_inches()[0] / fig.get_size_inches()[1]) * (y_range / x_range)  # 校正宽高比

    dx = x_range * 0.1 * aspect_ratio
    dy = y_range * 0.1 * aspect_ratio  

    head_width = dy * 0.1
    head_length = dx * 0.1

    ax.arrow(
        x_min, y_min,
        dx, 0,
        head_width=head_width,
        head_length=head_length,
        fc='black',
        ec='black'
    )

    ax.arrow(
        x_min, y_min,
        0, dy,
        head_width=head_width,
        head_length=head_length,
        fc='black',
        ec='black'
    )

    label_offset = dx * 0.1
    ax.text(
        x_min + dx/2, 
        y_min - label_offset,  
        'UMAP1',
        ha='center', 
        va='top', 
        fontsize=6
    )
    ax.text(
        x_min - label_offset, 
        y_min + dy/2,
        'UMAP2',
        rotation=90,
        ha='right', 
        va='center',
        fontsize=6
    )

banksy_enhance_umap(adata,color="banksy_leiden",)

#%% 看看banksy的结果能用吗
# FindAllMarkers
sc.tl.rank_genes_groups(adata, groupby="banksy_leiden", layer="lognorm", pts=True)
sc.tl.dendrogram(adata, groupby="banksy_leiden")

# Now we perform cell-type annotation based on marker genes
# We draw a dotplot for clear visulization
marker_genes = {
    "HSC": [
        "CD34",   # 经典 HSC / HSPC 表面标记
        "FLT3",
        "SELL",
        "SPN",
        "GATA2",
        "CD82",
        "CDK6",
        "SOCS2",
        "PBXIP1"
    ],
    "Erythroid": [
        "GYPA",
        "TFRC",   # CD71
        "EPOR"
    ],
    
    "Megakaryocyte": [
        "ITGA2B",
        "ITGB3",
        "GP1BA",
        "THBS1"
    ],
    
    "Eosinophil": [
        "IL5RA",
        "PRG2",
        "CCR3"
    ],
    
    "Basophil": [
        "ENPP3",
        "KIT",
        "HDC",
        "GATA2"
    ],
    
    # Bone Related
    "Osteoblast": [
        "ALPL",
        "RUNX2"
    ],
        
    "Osteoclast": [
        "CTSK",
        "ACP5",
        "CALCR"
    ],
    
    "Adipocyte": [
        "ADIPOQ",
        "PLIN1",
        "LPL"
    ],
    
    
    # T cells
    "general_T":["CD2","CD3E","CD3G","CD247"],
    "CD4_lineage":["CD4","CCR7"],
    "Regulatory":["FOXP3","IL2RA","IKZF2"],
    "CD8_lineage":["CD8A","CD8B","GZMB","GZMH","PRF1"],
    
    # NK cells
    "NK":["KLRB1","KLRC1","KLRD1","KLRK1"],
    
    # B cells
    "B":["CD19","MS4A1","CD79A","CD79B","BANK1"],
    
    # Plasma cells
    "Plsama":["MZB1","XBP1","PRDM1"],
    
    # Monocytes
    "Mono/Macro":["CD14","FCGR3A","CSF1R","CX3CR1","ITGAM"],
    
    # DCs
    "DC":["ITGAX","CD1C","CLEC9A","LILRA4","IRF4","IRF8"],
    
    # Neutrophils
    "Neutrophil":["MPO","ELANE","CXCR2"],
    
    "Mast":["KIT"],
    
    # Endothelium
    "Endothelium":["PECAM1","CDH5","KDR","ENG"],
    "Epithelium":["EPCAM","KRT8"],
    
    # Stromal
    "Fibroblast":["PDGFRA","PDGFRB","THY1","COL5A1"],
    "Pericyte":["CSPG4","RGS5"]
}

panel_genes = set(adata.var_names)

filtered_markers = {
    k: [g for g in v if g in panel_genes]
    for k, v in marker_genes.items()
}


sc.pl.dotplot(adata, 
              filtered_markers,
              groupby="banksy_leiden", 
              standard_scale="var",
              cmap='Reds'
             )
sc.pl.rank_genes_groups_dotplot(adata,groupby="banksy_leiden",n_genes=15,cmap="Reds")
#%% 细胞注释
# 我们对0.2分辨率的结果进行细胞注释
# Now we perform celltype annotation 
cluster_ids = sorted(adata.obs['leiden'].unique().tolist())
cluster_ids = [str(c) for c in cluster_ids]  # 转成字符串

# We create a dataframe first
celltype = pd.DataFrame({
    'ClusterID': cluster_ids,
    'celltype': ['Unknown'] * len(cluster_ids)
})

# We put in celltype information here
celltype.loc[celltype['ClusterID'].isin(['0']), 'celltype'] = 'Neutrophil'
celltype.loc[celltype['ClusterID'].isin(['1']), 'celltype'] = 'B'
celltype.loc[celltype['ClusterID'].isin(['2']), 'celltype'] = 'Mono/Macro'
celltype.loc[celltype['ClusterID'].isin(['3']), 'celltype'] = 'Erythrocyte'
celltype.loc[celltype['ClusterID'].isin(['4']), 'celltype'] = 'Endothelium'
celltype.loc[celltype['ClusterID'].isin(['5']), 'celltype'] = 'T'
celltype.loc[celltype['ClusterID'].isin(['6']), 'celltype'] = 'Fibroblast'
celltype.loc[celltype['ClusterID'].isin(['7']), 'celltype'] = 'Fibroblast'
celltype.loc[celltype['ClusterID'].isin(['8']), 'celltype'] = 'Mono/Macro'
celltype.loc[celltype['ClusterID'].isin(['9']), 'celltype'] = 'Plasma'
celltype.loc[celltype['ClusterID'].isin(['10']), 'celltype'] = 'Pericyte'
celltype.loc[celltype['ClusterID'].isin(['11']), 'celltype'] = 'Megakaryocyte'
celltype.loc[celltype['ClusterID'].isin(['12']), 'celltype'] = 'Adipocyte'

# Now we map the information back
adata.obs['leiden_str'] = adata.obs['leiden'].astype(str) # safety check
mapping = dict(zip(celltype['ClusterID'], celltype['celltype']))
adata.obs['celltype'] = adata.obs['leiden_str'].map(mapping)

print(adata.obs['celltype'].value_counts())

#%% Finally, we introduce a customized function to draw beautiful UMAP plot
from sklearn.neighbors import KernelDensity

def enhance_umap(adata, color='celltype', figsize=(6,6)):
    fig, ax = plt.subplots(figsize=figsize)
    sc.pl.umap(
        adata,
        color=color,
        legend_loc=None,
        frameon=False,
        size=15,
        ax=ax,
        add_outline=True,
        show=False
    )

    for cell_type in adata.obs[color].unique():
        mask = adata.obs[color] == cell_type
        x = adata.obsm["X_umap"][mask, 0]
        y = adata.obsm["X_umap"][mask, 1]

        x_min, x_max = x.min() - 0.5, x.max() + 0.5
        y_min, y_max = y.min() - 0.5, y.max() + 0.5
        xx, yy = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]

        xy_train = np.vstack([x, y]).T
        xy_test = np.vstack([xx.ravel(), yy.ravel()]).T

        kde = KernelDensity(bandwidth=0.5, metric='euclidean')
        kde.fit(xy_train)

        Z = np.exp(kde.score_samples(xy_test))
        Z = Z.reshape(xx.shape)

    cell_type_counts = adata.obs[color].value_counts(normalize=True) * 100
    for cell_type, percentage in cell_type_counts.items():
        mask = adata.obs[color] == cell_type
        x = np.median(adata.obsm["X_umap"][mask, 0])
        y = np.median(adata.obsm["X_umap"][mask, 1])
        ax.text(
            x, y, f"{cell_type}",
            fontsize=10, 
            ha="center", 
            color="black",
            bbox=dict(
                facecolor='white',         
                edgecolor='none', 
                alpha=0.65,
                boxstyle="round,pad=0.4"  
            )  
        )
    ax = plt.gca()
    ax.set_xlabel('UMAP1', fontsize=6)
    ax.set_ylabel('UMAP2', fontsize=6)
 
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()

    x_range = x_max - x_min
    y_range = y_max - y_min
    aspect_ratio = (fig.get_size_inches()[0] / fig.get_size_inches()[1]) * (y_range / x_range)  # 校正宽高比

    dx = x_range * 0.1 * aspect_ratio
    dy = y_range * 0.1 * aspect_ratio  

    head_width = dy * 0.1
    head_length = dx * 0.1

    ax.arrow(
        x_min, y_min,
        dx, 0,
        head_width=head_width,
        head_length=head_length,
        fc='black',
        ec='black'
    )

    ax.arrow(
        x_min, y_min,
        0, dy,
        head_width=head_width,
        head_length=head_length,
        fc='black',
        ec='black'
    )

    label_offset = dx * 0.1
    ax.text(
        x_min + dx/2, 
        y_min - label_offset,  
        'UMAP1',
        ha='center', 
        va='top', 
        fontsize=6
    )
    ax.text(
        x_min - label_offset, 
        y_min + dy/2,
        'UMAP2',
        rotation=90,
        ha='right', 
        va='center',
        fontsize=6
    )
enhance_umap(adata, color='leiden', figsize=(6,6))
enhance_umap(adata, color='celltype', figsize=(6,6))

#%% 保存h5ad
adata.write_h5ad("4_control_major_annotation.h5ad")

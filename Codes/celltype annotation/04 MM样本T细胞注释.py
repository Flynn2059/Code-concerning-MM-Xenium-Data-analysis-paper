#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  6 13:30:40 2026

@author: flynn
"""

#%% load packages
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

os.chdir('/Volumes/FlynnDisk/XinyeLi/GSE299207_RAW/')
print(os.getcwd())
print(os.listdir())
# some customized settings for better visulization
plt.rcParams["figure.figsize"] = (12, 12)
plt.rcParams["figure.dpi"] = 300
#%% load adata
# load the data that we annotated before
adata=sc.read_h5ad('./10_mm_major_annotation.h5ad')
print(adata)
print(adata.obs['celltype'].value_counts())
annotation = pd.DataFrame(
    {"major_annotation": adata.obs["celltype"]},
    index=adata.obs_names
)
annotation.to_csv("10_mm_major_annotation.csv")
print(adata.obs['celltype'].value_counts())
adata.X=adata.layers["counts"].copy()

# 提取T/NK亚群
# We extract T/NK cells, the major cells in tumor micro environment for further analysis.
adata_T_NK=adata[(adata.obs["celltype"]=="T/NK")]
print(adata_T_NK)
print(adata_T_NK.X)
#%% 重新降维聚类
sc.pp.normalize_total(adata_T_NK,target_sum=np.median(adata_T_NK.obs["total_counts"]))
sc.pp.log1p(adata_T_NK)
sc.pp.highly_variable_genes(adata_T_NK,n_top_genes=3000,subset=False)
sc.pp.scale(adata_T_NK,max_value=10)
sc.tl.pca(adata_T_NK)
sc.pl.pca_variance_ratio(adata_T_NK,log=True,n_pcs=50)
sc.external.pp.harmony_integrate(adata_T_NK,key="sample")
# Now we perform further dimentionality reduction and 
# unsupervised clustering based on harmony-corrected PCs
sc.pp.neighbors(adata_T_NK,n_pcs=20,use_rep="X_pca_harmony")
sc.tl.umap(adata_T_NK)
sc.pl.umap(adata_T_NK,color="sample")
#%% # We add a special cell for unsupervised clustering alone
# because this step is repeated many times until an ideal result is obtained
sc.tl.leiden(adata_T_NK,resolution=0.4,flavor="igraph",n_iterations=-1)
print(adata_T_NK.obs["leiden"].value_counts())
sc.pl.umap(adata_T_NK, color=["leiden"],legend_loc="on data")
# Now we perform cell-type annotation based on marker genes
# We draw a dotplot for clear visulization
marker_genes = {
    "common": ["CD2",'CD3E','CD4','CD8A',"CD8B"],
    "Naive/Memory":["TCF7","CCR7","SELL"],
    "Effect/Effect Memory":["CD44"],
    "NK":["GZMA","GZMB","GZMH","PRF1"],
    "Exhaustion":["PDCD1","CTLA4","TOX","HAVCR2","LAG3"],
    "Activation":["TIGIT"],
    "Proliferation":["MKI67","CCND1","CCNA2"],
    "Follicle Helper":["CXCL13","CD200"],
    "Tissue Residue Memory":["ITGAE","ZNF683","HOPX"],
    "Regulatory":["FOXP3","IL2RA"],
    "ILC":["KLRK1","THY1","IFNG","RORC","GATA6"],#"IL17A",
    "gdT":["KLRC2", "TRDC",  "KIR2DL4"]} # "TRDC1", "TRG-AS1",

sc.pl.dotplot(adata_T_NK, 
              marker_genes,
              groupby="leiden", 
              standard_scale="var",
              cmap='Blues'
             )

# FindAllMarkers
sc.tl.rank_genes_groups(adata_T_NK, groupby="leiden", layer="lognorm", pts=True)
sc.tl.dendrogram(adata_T_NK, groupby="leiden")
sc.pl.rank_genes_groups_dotplot(adata_T_NK,n_genes=15,cmap="Blues",swap_axes=True)
#%% 细胞注释
# Now we perform celltype annotation 
cluster_ids = sorted(adata_T_NK.obs['leiden'].unique().tolist())
cluster_ids = [str(c) for c in cluster_ids]  # 转成字符串

# We create a dataframe first
celltype = pd.DataFrame({
    'ClusterID': cluster_ids,
    'celltype': ['Unknown'] * len(cluster_ids)
})

# We put in celltype information here
celltype.loc[celltype['ClusterID'].isin(['0']), 'celltype'] = 'Myeloid'
celltype.loc[celltype['ClusterID'].isin(['1','2']), 'celltype'] = 'CD8+ T Effector'
celltype.loc[celltype['ClusterID'].isin(['3']), 'celltype'] = 'CD4+ T Naive'
celltype.loc[celltype['ClusterID'].isin(['4']), 'celltype'] = 'SOX2-OT Unknown'
celltype.loc[celltype['ClusterID'].isin(['5']), 'celltype'] = 'T Proliferating'
celltype.loc[celltype['ClusterID'].isin(['6']), 'celltype'] = 'CD8+ T Exhausted'
celltype.loc[celltype['ClusterID'].isin(['7','8']), 'celltype'] = 'Myeloid'


# Now we map the information back
adata_T_NK.obs['leiden_str'] = adata_T_NK.obs['leiden'].astype(str) # safety check
mapping = dict(zip(celltype['ClusterID'], celltype['celltype']))
adata_T_NK.obs['celltype_secondary'] = adata_T_NK.obs['leiden_str'].map(mapping)

print(adata_T_NK.obs['celltype_secondary'].value_counts())

#%% 细胞注释结束
# Finally, we introduce a customized function to draw beautiful UMAP plot
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
enhance_umap(adata_T_NK, color='leiden', figsize=(6,6))
enhance_umap(adata_T_NK, color='celltype_secondary', figsize=(6,6))
#%% 保存细胞注释结果
annotation = pd.DataFrame(
    {"secondary_annotation": adata_T_NK.obs["celltype_secondary"]},
    index=adata_T_NK.obs_names
)
annotation.to_csv("10_mm_T_secondary_annotation.csv")
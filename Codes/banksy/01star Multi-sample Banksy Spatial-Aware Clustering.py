#%%!/usr/bin/env python3
#%% -*- coding: utf-8 -*-
"""
Created on Wed Mar  4 11:16:07 2026

@author: flynn
"""


#%% Imports
import os
import random
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
import scanpy as sc
import anndata as ad
# BANKSY / Harmony / UMAP
from banksy.initialize_banksy import initialize_banksy
from banksy.embed_banksy import generate_banksy_matrix
from banksy_utils.umap_pca import pca_umap
from banksy.cluster_methods import run_Leiden_partition
from banksy.plot_banksy import plot_results

from harmony import harmonize
import umap

warnings.filterwarnings("ignore")

sc.logging.print_header()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 1  #%% errors (0), warnings (1), info (2), hints (3)
sns.set_style("white")
# Note that BANKSY itself is deterministic, here the seeds affect the umap clusters and leiden partition
seed = 1234
np.random.seed(seed)
random.seed(seed)
#%% 工作路径调整
os.chdir('/mnt/bioSSD/GSE299207_RAW/')
os.listdir('./banksy/')

#%% 加载我们要用到的数据
adata_SM_1 = sc.read_h5ad('./banksy/GSM9035039_SM-1.adata')
adata_SM_2 = sc.read_h5ad('./banksy/GSM9035040_SM-2.adata')
adata_SM_3 = sc.read_h5ad('./banksy/GSM9035041_SM-3.adata')
adata_SM_4 = sc.read_h5ad('./banksy/GSM9035042_SM-4.adata')
adata_SM_5 = sc.read_h5ad('./banksy/GSM9035043_SM-5.adata')

adatas = [adata_SM_1, adata_SM_2, adata_SM_3, adata_SM_4, adata_SM_5]
sample_names = ["SM-1", "SM-2", "SM-3", "SM-4", "SM-5"]


def add_spatial_obs(adata):
    if "spatial" not in adata.obsm_keys():
        raise KeyError('adata missing adata.obsm["spatial"].')

    spatial = adata.obsm["spatial"]

    if spatial.shape[1] < 2:
        raise ValueError(f"adata.obsm['spatial'] shape {spatial.shape}, need >=2 columns")

    adata.obs["x_slide_mm"] = spatial[:, 0]
    adata.obs["y_slide_mm"] = spatial[:, 1]

    return adata


adatas = [add_spatial_obs(a) for a in adatas]

for a in adatas:
    print(a.obs[["sample", "x_slide_mm", "y_slide_mm"]].head())

adata_all = ad.concat(adatas, label="sample", keys=sample_names, join="outer", merge="same")

# 保险起见：用 obs 里的 sample 字段作为 batch 信息
adata_all.obs["sample"] = adata_all.obs["sample"].astype("category")

print(adata_all)
print("samples:", adata_all.obs["sample"].value_counts())

#%% 样本根据坐标平移
coords_df = pd.DataFrame({
    "x": adata_all.obs["x_slide_mm"].values,
    "y": adata_all.obs["y_slide_mm"].values,
    "sample": adata_all.obs["sample"]
})

# 确保 sample 是 category
coords_df["sample"] = coords_df["sample"].astype("category")

# 直接得到 int 编码（0,1,2,3,4）
coords_df["sample_idx"] = coords_df["sample"].cat.codes.astype(int)

# 计算每个sample的范围
sample_ranges = coords_df.groupby("sample").agg(
    x_min=("x", "min"),
    x_max=("x", "max"),
    y_min=("y", "min"),
    y_max=("y", "max")
)
sample_ranges["x_range"] = sample_ranges["x_max"] - sample_ranges["x_min"]
sample_ranges["y_range"] = sample_ranges["y_max"] - sample_ranges["y_min"]

print(sample_ranges)

# spacing
x_spacing = float(sample_ranges["x_range"].max()) * 1.2

# 平移（不会再报 categorical * float）
coords_df["x_new"] = coords_df["x"] + coords_df["sample_idx"].astype(float) * x_spacing
coords_df["y_new"] = coords_df["y"]

print(coords_df)

# 写回
adata_all.obs["x_slide_mm"] = coords_df["x_new"].values
adata_all.obs["y_slide_mm"] = coords_df["y_new"].values

# banksy坐标
adata_all.obsm["coord_xy"] = np.c_[adata_all.obs["x_slide_mm"].values,
                                  adata_all.obs["y_slide_mm"].values]

#%% 检查平移结果
fig, ax = plt.subplots(figsize=(12, 6))

samples = adata_all.obs["sample"].unique()

for sample in samples:
    idx = adata_all.obs["sample"] == sample
    ax.scatter(
        adata_all.obs.loc[idx, "x_slide_mm"],
        adata_all.obs.loc[idx, "y_slide_mm"],
        s=2,
        label=sample
    )

ax.legend(markerscale=4, bbox_to_anchor=(1.02, 1), loc="upper left")

#%% 去掉所有坐标轴元素
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_title("")

for spine in ax.spines.values():
    spine.set_visible(False)
plt.tight_layout()
plt.show()

#%% banksy前数据预处理
adata_all.layers["counts"] = adata_all.X.copy()
sc.pp.normalize_total(adata_all, target_sum=np.median(adata_all.obs["total_counts"]))
sc.pp.log1p(adata_all)
sc.pp.highly_variable_genes(adata_all, n_top_genes=3000, flavor="seurat", inplace=True)
adata_all = adata_all[:, adata_all.var.highly_variable].copy()
print(adata_all)

#%% banksy步骤1: generate banksy matrix
coord_keys = ('x_slide_mm', 'y_slide_mm', 'coord_xy')

k_geom = 8
nbr_weight_decay = "scaled_gaussian"
max_m = 1

lambda_list = [0.1]
resolutions = [1]

banksy_dict = initialize_banksy(
    adata_all,
    coord_keys,
    k_geom,
    nbr_weight_decay=nbr_weight_decay,
    max_m=max_m,
    plt_edge_hist=False,
    plt_nbr_weights=False,
    plt_agf_angles=False,
    plt_theta=False
)

banksy_dict, banksy_matrix = generate_banksy_matrix(
    adata_all,
    banksy_dict,
    lambda_list,
    max_m,
    variance_balance=True
)



#%% banksy步骤2: PCA + Harmony批次效应矫正
# 这里使用的Harmony是Harmony-Pytorch，没有显卡会算得比较慢
pca_dims = [20]  #%% Dimensionality in which PCA reduces to

pca_umap(
    banksy_dict,
    pca_dims=[20],
    add_umap=False,
    plt_remaining_var=False
)

for pca_dim in pca_dims:
    # 获取 BANKSY embedding
    adata_banksy = banksy_dict[nbr_weight_decay][0.1]["adata"]

    Z = harmonize(
        adata_banksy.obsm[f"reduced_pc_{pca_dim}"],
        adata_banksy.obs,
        batch_key="sample",
    )

    print(f"Replacing reduced_pc_{pca_dim} with Harmony corrected embeddings")

    adata_banksy.obsm[f"reduced_pc_{pca_dim}"] = Z

    # 重新计算 UMAP
    reducer = umap.UMAP(
        n_neighbors=15,
        min_dist=0.3,
        random_state=42
    )

    umap_embedding = reducer.fit_transform(Z)
    adata_banksy.obsm[f"reduced_pc_{pca_dim}_umap"] = umap_embedding

#%% 保存harmony矫正的PCA和UMAP的结果

print(adata_banksy)

banksy_pca = adata_banksy.obsm['reduced_pc_20']
banksy_umap = adata_banksy.obsm['reduced_pc_20_umap']

pca_df = pd.DataFrame(
    banksy_pca,
    index=adata_banksy.obs_names
)

umap_df = pd.DataFrame(
    banksy_umap,
    columns=["UMAP1", "UMAP2"],
    index=adata_banksy.obs_names
)

pca_df.to_csv("SM_5_samples_PCA.csv")
umap_df.to_csv("SM_5_samples_UMAP.csv")
#%% 非常耗时的leiden聚类 5个Xenium样本花了50分钟
results_df, max_num_labels = run_Leiden_partition(
    banksy_dict,
    resolutions,
    num_nn=50,
    num_iterations=-1,
    partition_seed=1234,
    match_labels=True
)
#%% 提取leiden聚类的结果画图并保存
label_obj = results_df.iloc[0]["labels"]
clusters = label_obj.dense

adata_banksy.obs["banksy_leiden"] = clusters.astype(str)
adata_banksy.obsm["X_umap"] = adata_banksy.obsm["reduced_pc_20_umap"]

sc.pl.umap(
    adata_banksy,
    color="banksy_leiden",
    legend_loc="on data",
)

print(adata_banksy.obs['banksy_leiden'])

leiden_clustering = adata_banksy.obs['banksy_leiden']
leiden_df = pd.DataFrame(
    leiden_clustering,
    index=adata_banksy.obs_names
)
leiden_df.to_csv("SM_5_samples_leiden_clustering.csv")
#%% 一些banksy自带的画图，可以考虑跳过
coord_keys = ('x_slide_mm', 'y_slide_mm', 'coord_xy')
c_map = 'tab20'  #%% specify color map
weights_graph = banksy_dict['scaled_gaussian']['weights'][0]

plot_results(
    results_df,
    weights_graph,
    c_map,
    match_labels=True,
    coord_keys=coord_keys,
    max_num_labels=max_num_labels,
    save_path="./banksy_plot/SM.png",
    save_fig=True
)
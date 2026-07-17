#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 14:57:46 2026

@author: flynn
"""

#%% imports
import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# some customized settings for better visualization
plt.rcParams["figure.figsize"] = (12, 8)
plt.rcParams["figure.dpi"] = 300
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42

#%% 工作路径
os.chdir('/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/')
print(os.getcwd())

#%% 输出目录
outdir = "./Tcell_Dominant_niche_Tcell_composition"
os.makedirs(outdir, exist_ok=True)

#%% 读取前一步保存的 pooled 细胞级 niche 矩阵
niche_df = pd.read_csv(
    "./Tcell_Plasma_GMM_clustering_annotation/MM_10_samples_cell_niche_matrix.csv"
)

print("Loaded niche matrix:")
print(niche_df.shape)

#%% 定义 T 细胞类型
t_celltypes = [
    "CD4+ T Naive",
    "CD8+ T Effector",
    "CD8+ T Exhausted",
    "T Proliferating"
]

#%% 只保留有 T_Dominant niche 的样本
samples_with_t_dominant = (
    niche_df.loc[niche_df["niche"] == "T_Dominant", "sample"]
    .dropna()
    .unique()
    .tolist()
)
samples_with_t_dominant = sorted(samples_with_t_dominant)

print("\nSamples with T_Dominant niche:")
print(samples_with_t_dominant)

#%% 提取这些样本中，属于 T_Dominant niche 的所有细胞
td_niche_df = niche_df[
    (niche_df["sample"].isin(samples_with_t_dominant)) &
    (niche_df["niche"] == "T_Dominant")
].copy()

print("\nCells in T_Dominant niche:")
print(td_niche_df.shape)

#%% 统计每个样本中，不同 T cell type 在 T_Dominant niche 内的细胞数
# 分子：每个 T cell subtype 的细胞数
count_df = (
    td_niche_df[td_niche_df["annotation"].isin(t_celltypes)]
    .groupby(["sample", "annotation"])
    .size()
    .unstack(fill_value=0)
)

# 确保列顺序固定
count_df = count_df.reindex(columns=t_celltypes, fill_value=0)

# 分母：每个样本中 T_Dominant niche 的总细胞数（注意：是 niche 内所有细胞，而不只是 T 细胞）
total_td_niche_cells = (
    td_niche_df.groupby("sample")
    .size()
    .rename("total_t_dominant_niche_cells")
)

# 合并并计算比例
prop_df = count_df.div(total_td_niche_cells, axis=0).fillna(0)

# 补充总数信息，便于导出查看
summary_df = count_df.copy()
summary_df["total_t_dominant_niche_cells"] = total_td_niche_cells
for ct in t_celltypes:
    summary_df[f"{ct}_ratio"] = prop_df[ct]

# 按样本顺序整理
summary_df = summary_df.loc[samples_with_t_dominant]
prop_df = prop_df.loc[samples_with_t_dominant]

print("\nSummary table:")
print(summary_df)

#%% 保存统计结果
summary_df.to_csv(
    os.path.join(outdir, "T_Dominant_niche_Tcell_composition_summary.csv")
)
prop_df.to_csv(
    os.path.join(outdir, "T_Dominant_niche_Tcell_composition_ratio_only.csv")
)

#%% 绘制堆叠柱状图
color_map = {
    "CD4+ T Naive": "#4C72B0",
    "CD8+ T Effector": "#DD8452",
    "CD8+ T Exhausted": "#55A868",
    "T Proliferating": "#C44E52"
}

fig, ax = plt.subplots(figsize=(14, 8), dpi=300)

bottom = np.zeros(len(prop_df))
x = np.arange(len(prop_df.index))

for ct in t_celltypes:
    values = prop_df[ct].values
    ax.bar(
        x,
        values,
        bottom=bottom,
        label=ct,
        color=color_map[ct],
        width=0.8
    )
    bottom += values

ax.set_xticks(x)
ax.set_xticklabels(prop_df.index, rotation=45, ha="right")
ax.set_ylabel("Proportion of total T_Dominant niche cells")
ax.set_xlabel("Sample")
ax.set_title("T-cell composition within T_Dominant niche across samples")
ax.legend(frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")

ax.set_ylim(0, bottom.max() * 1.05)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(
    os.path.join(outdir, "T_Dominant_niche_Tcell_composition_stacked_bar.pdf"),
    bbox_inches="tight"
)
plt.savefig(
    os.path.join(outdir, "T_Dominant_niche_Tcell_composition_stacked_bar.svg"),
    bbox_inches="tight"
)
plt.close()

print("\nDone.")
print("Results saved to:", outdir)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 12:35:44 2026

@author: flynn
"""

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
os.chdir('/Volumes/FlynnDisk/XinyeLi/GSE299207_RAW/')
print(os.getcwd())
print("\n")
print(os.listdir())

#%% 读取数据
# collect all zarr folders
zarr_files = [
    "GSM9035039_SM-1_post_qc.zarr",
    "GSM9035040_SM-2_post_qc.zarr",
    "GSM9035041_SM-3_post_qc.zarr",
    "GSM9035042_SM-4_post_qc.zarr",
    "GSM9035043_SM-5_post_qc.zarr",
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

#%% 我们读取之前的细胞注释结果
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
annotation_level1['major_annotation'].value_counts()

celltype_order=[
    'CD4+ T Naive',
    'CD8+ T Effector',
    'T Proliferating',
    'B',
    'Plasma',
    'Myeloid',
    'Neutrophil',
    'Adipocyte',
    'Fibroblast',
    'Endothelium',
    'HSC',
    'Erythoid',
    'Megakaryocyte',
    'SOX2-OT Unknown'
    ]


annotation_split = {}

for sdata_key in sdata_dict.keys():
    # 从 sdata_key 提取 sample_name
    sample_name = sdata_key.replace("sdata_", "")
    prefix = sample_name + "_"

    # 按前缀筛选这个样本的注释
    annotation_sub = annotation_level1[annotation_level1.index.str.startswith(prefix)].copy()

    annotation_split[sample_name] = annotation_sub

    print(f"{sample_name}: {annotation_sub.shape[0]} cells")

for sdata_key, sdata in sdata_dict.items():
    sample_name = sdata_key.replace("sdata_", "")
    annotation_sub = annotation_split[sample_name]

    adata = sdata["table"].copy()

    adata.obs["annotation"] = pd.NA

    common_idx = adata.obs_names.intersection(annotation_sub.index)

    adata.obs.loc[common_idx, "annotation"] = annotation_sub.loc[
        common_idx, "major_annotation"
    ].astype(str)

    # 关键：这里重新指定 annotation 这一列的类别顺序
    adata.obs["annotation"] = pd.Categorical(
        adata.obs["annotation"],
        categories=celltype_order,
        ordered=True
    )

    # 可选：如果不想 legend 显示该样本中不存在的类别
    #adata.obs["annotation"] = adata.obs["annotation"].cat.remove_unused_categories()

    sdata["table"] = adata

    print(f"{sample_name}: matched {len(common_idx)} / {adata.n_obs} cells")

for sdata_key, sdata in sdata_dict.items():
    sample_name = sdata_key.replace("sdata_", "")
    print(f"\n{sample_name}")
    print(sdata["table"].obs["annotation"].value_counts(dropna=False))
    
#%% 补充图1:全局的细胞类型，每个样本分别画图
outdir_global = "./Global_celltype_plots"
os.makedirs(outdir_global, exist_ok=True)

palette=[
    '#1F77B4FF', #1
    '#FF7F0EFF', #2
#    '#8C564BFF', #3 
    '#D62728FF', #4
    '#9467BDFF', #5
    '#2CA02C55', #6
    '#E377C2FF', #7
    '#7F7F7FFF', #8
#    '#BCBD22FF', #9
    '#17BECFFF', #10
    '#AEC7E8FF', #11
#    '#FFBB78FF', #12
    '#98DF8AFF', #13
    '#FF9896FF', #14
    '#C5B0D5FF', #15
    '#C49C94FF', #16
#    '#F7B6D2FF', #17
    '#C7C7C7FF', #18
    ]

for sdata_key, sdata in sdata_dict.items():
    sample_name = sdata_key.replace("sdata_", "")
    prefix = sample_name + "_"

    adata = sdata.tables["table"].copy()

    # 恢复和 shapes 对应的原始 cell id
    adata.obs["cell_id"] = adata.obs_names.str.replace(f"^{prefix}", "", regex=True)
    adata.obs["region"] = "cell_boundaries"

    sdata["table"] = adata

    sdata.set_table_annotates_spatialelement(
        "table",
        region="cell_boundaries",
        region_key="region",
        instance_key="cell_id"
    )

    print(f"Plotting {sample_name} ...")

    plt.figure(figsize=(12, 12), dpi=600)

    sdata.pl.render_shapes(
        "cell_boundaries",
        groups=celltype_order,
        color="annotation",
        palette=palette).pl.show()

    ax = plt.gca()
    ax.set_title(sample_name, fontsize=16)
    ax.set_axis_off()

    plt.tight_layout()
    # 保存矢量图
    plt.savefig(
        os.path.join(outdir_global, f"{sample_name}_global_celltype.pdf"),
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(outdir_global, f"{sample_name}_global_celltype.svg"),
        bbox_inches="tight"
    )
    plt.show()
    plt.close()

#%% 补充图2：只绘制T细胞的（矢量图导出）
import os
import pandas as pd
import matplotlib.pyplot as plt

groups = ["Other", "CD4+ T Naive", "CD8+ T Effector","T Proliferating"]
palette = ["#BDBDBD", "#00A087FF", "#DC0000FF","#3C5488FF"]

outdir = "./Tcell_vector_plots"
os.makedirs(outdir, exist_ok=True)

# 让文字在 svg/pdf 中尽量保持可编辑
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42

for sdata_key, sdata in sdata_dict.items():
    sample_name = sdata_key.replace("sdata_", "")
    prefix = sample_name + "_"

    adata = sdata.tables["table"].copy()

    adata.obs["cell_id"] = adata.obs_names.str.replace(f"^{prefix}", "", regex=True)
    adata.obs["region"] = "cell_boundaries"

    adata.obs["Tcell_highlight"] = "Other"
    adata.obs.loc[adata.obs["annotation"] == "CD4+ T Naive", "Tcell_highlight"] = "CD4+ T Naive"
    adata.obs.loc[adata.obs["annotation"] == "CD8+ T Effector", "Tcell_highlight"] = "CD8+ T Effector"
    adata.obs.loc[adata.obs["annotation"] == "T Proliferating", "Tcell_highlight"] = "T Proliferating"
    adata.obs["Tcell_highlight"] = pd.Categorical(
        adata.obs["Tcell_highlight"],
        categories=groups
    )

    sdata["table"] = adata

    sdata.set_table_annotates_spatialelement(
        "table",
        region="cell_boundaries",
        region_key="region",
        instance_key="cell_id"
    )

    print(f"Plotting {sample_name} ...")

    fig = plt.figure(figsize=(12, 12),dpi=600)

    sdata.pl.render_shapes(
        "cell_boundaries",
        color="Tcell_highlight",
        groups=groups,
        palette=palette
    ).pl.show()

    ax = plt.gca()
    ax.set_title(sample_name, fontsize=16)
    ax.set_axis_off()

    plt.tight_layout()

    # 导出矢量图
    plt.savefig(
        os.path.join(outdir, f"{sample_name}_Tcell_highlight.pdf"),
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(outdir, f"{sample_name}_Tcell_highlight.svg"),
        bbox_inches="tight"
    )

    plt.show()
    plt.close()


#%% 统计每个样本的 CD4 和 CD8 组成占比并画图
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

outdir = "./Tcell_barplots"
os.makedirs(outdir, exist_ok=True)

sample_order = [
    "GSM9035039_SM-1",
    "GSM9035040_SM-2",
    "GSM9035041_SM-3",
    "GSM9035042_SM-4",
    "GSM9035043_SM-5",
]

sample_labels = ["SM-1", "SM-2", "SM-3", "SM-4", "SM-5"]

result_list = []

for sample_name, sample_label in zip(sample_order, sample_labels):
    sdata_key = f"sdata_{sample_name}"
    sdata = sdata_dict[sdata_key]
    adata = sdata.tables["table"]

    ann = adata.obs["annotation"].copy()

    total_cells = ann.notna().sum()

    CD4_label_mask = ann.isin([
        "CD4+ T Naive"
    ])

    CD8_label_mask = ann.isin([
        "CD8+ T Effector"
    ])

    cd4_count = CD4_label_mask.sum()
    cd8_count = CD8_label_mask.sum()

    cd4_pct_total = cd4_count / total_cells * 100
    cd8_pct_total = cd8_count / total_cells * 100

    result_list.append({
        "sample": sample_name,
        "sample_label": sample_label,
        "total_cells": total_cells,
        "CD4_count": cd4_count,
        "CD8_count": cd8_count,
        "CD4_pct_total": cd4_pct_total,
        "CD8_pct_total": cd8_pct_total
    })

df_plot = pd.DataFrame(result_list)

x = np.arange(len(df_plot))

cd4_values = df_plot["CD4_pct_total"].values
cd8_values = df_plot["CD8_pct_total"].values

fig, ax = plt.subplots(figsize=(6, 6), dpi=600)

ax.bar(
    x,
    cd4_values,
    width=0.7,
    label="CD4+ T Naive",
    color="#00A087FF"
)

ax.bar(
    x,
    cd8_values,
    width=0.7,
    bottom=cd4_values,
    label="CD8+ T Effector",
    color="#DC0000FF"
)

ax.set_xticks(x)
ax.set_xticklabels(df_plot["sample_label"], fontsize=11)
ax.set_ylabel("Percentage of total cells (%)", fontsize=12)
ax.set_xlabel("Sample", fontsize=12)

ax.set_title("")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.legend(
    frameon=False,
    loc="upper left",
    bbox_to_anchor=(1.02, 1),
    borderaxespad=0
)

plt.tight_layout()

plt.savefig(
    os.path.join(outdir, "SM_T_cell_barplot.pdf"),
    bbox_inches="tight"
)

plt.savefig(
    os.path.join(outdir, "SM_T_cell_barplot.svg"),
    bbox_inches="tight"
)

plt.close()


#%% 输出每个样本的 CD4/CD8 比例表格
outdir_ratio = "./Tcell_barplot"
os.makedirs(outdir_ratio, exist_ok=True)

df_plot["CD4_CD8_ratio"] = np.where(
    df_plot["CD8_count"] > 0,
    df_plot["CD4_count"] / df_plot["CD8_count"],
    np.nan
)

df_ratio = df_plot[["sample", "CD4_CD8_ratio"]].copy()
df_ratio = df_ratio.set_index("sample")

df_ratio["CD4_CD8_ratio"] = df_ratio["CD4_CD8_ratio"].round(4)

df_ratio.to_csv(
    os.path.join(outdir_ratio, "SM_CD4_CD8_ratio_table.csv")
)

print("CD4/CD8 ratio table:")
print(df_ratio)
#%% T细胞的Neighborhood空间密度图
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree
outdir = "./Tcell_spatial_density_plots"
os.makedirs(outdir, exist_ok=True)
radius = 200

for sdata_key, sdata in sdata_dict.items():
    sample_name = sdata_key.replace("sdata_", "")
    prefix = sample_name + "_"

    print(f"Processing {sample_name} ...")

    adata = sdata.tables["table"].copy()
    gdf = sdata.shapes["cell_boundaries"].copy()
    gdf.index = gdf.index.astype(str)

    # 临时 patch：把 obs_names 的样本前缀去掉，得到和 shape index 一致的 barcode
    adata.obs["cell_id"] = adata.obs_names.to_series().str.replace(f"^{prefix}", "", regex=True)
    adata.obs["region"] = "cell_boundaries"

    # 先把 annotation 映射到 shapes，方便筛 T 细胞
    ann_map = adata.obs.set_index("cell_id")["annotation"]
    gdf["annotation"] = gdf.index.map(ann_map)

    # 所有细胞 centroid
    all_centroids = gdf.geometry.centroid
    all_coords = np.column_stack([all_centroids.x.values, all_centroids.y.values])

    # T 细胞 centroid
    t_mask = gdf["annotation"].isin(["CD4+ T Naive", "CD8+ T Effector","T Proliferating"])
    t_gdf = gdf.loc[t_mask].copy()

    if t_gdf.shape[0] == 0:
        print(f"Skipping {sample_name}: no T cells found")
        continue

    t_centroids = t_gdf.geometry.centroid
    t_coords = np.column_stack([t_centroids.x.values, t_centroids.y.values])

    # 计算每个细胞周围 radius 内 T 细胞数
    tree_t = KDTree(t_coords)
    t_count = tree_t.query_radius(all_coords, r=radius, count_only=True)

    # 转成局部密度
    local_t_density = t_count / (np.pi * radius * radius)

    # 写回 table.obs
    density_series = pd.Series(local_t_density, index=gdf.index, name="local_t_density")
    adata.obs["local_t_density"] = adata.obs["cell_id"].map(density_series)

    # 放回 sdata
    sdata["table"] = adata

    # 重新建立 table 和 shapes 的关联
    sdata.set_table_annotates_spatialelement(
        "table",
        region="cell_boundaries",
        region_key="region",
        instance_key="cell_id"
    )

    print(f"Plotting {sample_name} ...")

    plt.figure(figsize=(12, 12), dpi=600)

    sdata.pl.render_shapes(
        "cell_boundaries",
        color="local_t_density",
        cmap="coolwarm"
    ).pl.show()

    ax = plt.gca()
    ax.set_title(sample_name, fontsize=16)
    ax.set_axis_off()

    plt.savefig(
    os.path.join(outdir, f"{sample_name}_T_cell_spatial_density_plot.pdf"),
    bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(outdir, f"{sample_name}_T_cell_spatial_density_plot.svg"),
        bbox_inches="tight"
    )
    plt.close()
#%% T细胞和浆细胞共定位的图
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree
from matplotlib.colors import to_hex

outdir = "./Tcell_Plasma_spatial_density_plots"
os.makedirs(outdir, exist_ok=True)

radius = 200

# 根据你的 annotation 名称自行补充
t_labels = ["CD4+ T Naive", "CD8+ T Effector", "T Proliferating"]
plasma_labels = ["Plasma"]

for sdata_key, sdata in sdata_dict.items():
    sample_name = sdata_key.replace("sdata_", "")
    prefix = sample_name + "_"

    print(f"Processing {sample_name} ...")

    adata = sdata.tables["table"].copy()
    gdf = sdata.shapes["cell_boundaries"].copy()
    gdf.index = gdf.index.astype(str)

    # 恢复和 shapes 对应的原始 cell_id
    adata.obs["cell_id"] = adata.obs_names.to_series().str.replace(f"^{prefix}", "", regex=True)

    # 把 annotation 映射到 shapes
    ann_map = adata.obs.set_index("cell_id")["annotation"]
    gdf["annotation"] = gdf.index.map(ann_map)

    # 所有细胞 centroid
    centroids = gdf.geometry.centroid
    coords = np.column_stack([centroids.x.values, centroids.y.values])

    # 建树
    tree = KDTree(coords)

    # 查询每个细胞在固定半径内的邻居
    neighbor_indices = tree.query_radius(coords, r=radius)

    # 准备 annotation 数组
    ann_array = gdf["annotation"].astype("object").values
    is_t = np.isin(ann_array, t_labels)
    is_plasma = np.isin(ann_array, plasma_labels)

    t_density = np.zeros(len(gdf), dtype=float)
    plasma_density = np.zeros(len(gdf), dtype=float)

    area = np.pi * radius * radius

    for i, nbr_idx in enumerate(neighbor_indices):
        t_count = is_t[nbr_idx].sum()
        plasma_count = is_plasma[nbr_idx].sum()

        t_density[i] = t_count / area
        plasma_density[i] = plasma_count / area

    # 写回表
    density_t_series = pd.Series(t_density, index=gdf.index, name="local_t_density")
    density_plasma_series = pd.Series(plasma_density, index=gdf.index, name="local_plasma_density")

    adata.obs["local_t_density"] = adata.obs["cell_id"].map(density_t_series)
    adata.obs["local_plasma_density"] = adata.obs["cell_id"].map(density_plasma_series)

    sdata["table"] = adata

    # 也写回 gdf，方便直接作图
    gdf["local_t_density"] = t_density
    gdf["local_plasma_density"] = plasma_density

    # 为了颜色混合，分别归一化到 0~1
    t_max = gdf["local_t_density"].max()
    p_max = gdf["local_plasma_density"].max()

    if t_max == 0:
        gdf["t_norm"] = 0.0
    else:
        gdf["t_norm"] = gdf["local_t_density"] / t_max

    if p_max == 0:
        gdf["p_norm"] = 0.0
    else:
        gdf["p_norm"] = gdf["local_plasma_density"] / p_max

    # 颜色逻辑
    # T-only 部分显示黄色
    # Plasma-only 部分显示蓝色
    # 两者共同高的部分显示绿色
    common = np.minimum(gdf["t_norm"].values, gdf["p_norm"].values)
    t_only = np.clip(gdf["t_norm"].values - gdf["p_norm"].values, 0, None)
    p_only = np.clip(gdf["p_norm"].values - gdf["t_norm"].values, 0, None)

    # overall signal strength
    signal = np.maximum(gdf["t_norm"].values, gdf["p_norm"].values)

    # color settings
    gray = np.array([0.88, 0.88, 0.88])   # light gray background
    red = np.array([1.0, 0.1, 0.1])   # T high
    blue = np.array([0.2, 0.45, 1.0])     # Plasma high
    green = np.array([0.2, 0.75, 0.2])    # both high

    # start from light gray, then blend toward target colors
    rgb = (1 - signal)[:, None] * gray
    rgb += t_only[:, None] * red
    rgb += p_only[:, None] * blue
    rgb += common[:, None] * green

    rgb = np.clip(rgb, 0, 1)
    gdf["mix_color"] = [to_hex(c) for c in rgb]

    print(f"Plotting {sample_name} ...")

    fig, ax = plt.subplots(figsize=(12, 12), dpi=600)
    gdf.plot(ax=ax, color=gdf["mix_color"], linewidth=0)

    ax.set_title(sample_name, fontsize=16)
    ax.set_axis_off()

    plt.tight_layout()
    plt.savefig(
        os.path.join(outdir, f"{sample_name}_T_Plasma_spatial_density_mix.pdf"),
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(outdir, f"{sample_name}_T_Plasma_spatial_density_mix.svg"),
        bbox_inches="tight"
    )
    plt.close()

print("Done.")

#%% 5个样本整合展示的 2D-KDE panel，仅用 contour 线条展示
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

outdir = "./Tcell_Plasma_2D_KDE"
os.makedirs(outdir, exist_ok=True)

sample_order = [
    "GSM9035039_SM-1",
    "GSM9035040_SM-2",
    "GSM9035041_SM-3",
    "GSM9035042_SM-4",
    "GSM9035043_SM-5",
]

fig, axes = plt.subplots(2, 3, figsize=(15, 10), dpi=600)
axes = axes.flatten()

for i, sample_name in enumerate(sample_order):
    ax = axes[i]
    sdata_key = f"sdata_{sample_name}"
    sdata = sdata_dict[sdata_key]
    adata = sdata.tables["table"].copy()

    df = adata.obs[["local_t_density", "local_plasma_density"]].copy()
    df = df.dropna()

    x = df["local_t_density"].values
    y = df["local_plasma_density"].values

    if len(df) < 10 or np.allclose(x, x[0]) or np.allclose(y, y[0]):
        ax.text(0.5, 0.5, "Insufficient data", ha="center", va="center", transform=ax.transAxes)
        ax.set_title(sample_name.split("_")[-1], fontsize=12)
        ax.set_xlabel("T density")
        ax.set_ylabel("Plasma density")
        continue

    values = np.vstack([x, y])
    kde = gaussian_kde(values)

    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()

    x_pad = (x_max - x_min) * 0.05 if x_max > x_min else 1e-6
    y_pad = (y_max - y_min) * 0.05 if y_max > y_min else 1e-6

    xx, yy = np.mgrid[
        (x_min - x_pad):(x_max + x_pad):200j,
        (y_min - y_pad):(y_max + y_pad):200j
    ]

    grid_coords = np.vstack([xx.ravel(), yy.ravel()])
    zz = kde(grid_coords).reshape(xx.shape)

    ax.contour(
        xx,
        yy,
        zz,
        levels=10,
        cmap="coolwarm",
        linewidths=1.0
    )

    ax.set_title(sample_name.split("_")[-1], fontsize=12)
    ax.set_xlabel("T density", fontsize=10)
    ax.tick_params(axis="x", labelrotation=30)
    for label in ax.get_xticklabels():
        label.set_ha("right")
    ax.set_ylabel("Plasma density", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

axes[-1].axis("off")

plt.tight_layout()
plt.savefig(
    os.path.join(outdir, "5_sm_T_vs_Plasma_2D_KDE_contour_panel.pdf"),
    bbox_inches="tight"
)
plt.savefig(
    os.path.join(outdir, "5_sm_T_vs_Plasma_2D_KDE_contour_panel.svg"),
    bbox_inches="tight"
)
plt.close()

print("Panel plot done.")


#%% 空间共现性分析
import squidpy as sq
import matplotlib.pyplot as plt
outdir = "./Tcell_spatial_co_localization_plots"
os.makedirs(outdir, exist_ok=True)

sample_order = [
    "sdata_GSM9035039_SM-1",
    "sdata_GSM9035040_SM-2",
    "sdata_GSM9035041_SM-3",
    "sdata_GSM9035042_SM-4",
    "sdata_GSM9035043_SM-5",
]

fig, axes = plt.subplots(2, 3, figsize=(38, 20), dpi=600)
axes = axes.flatten()

plt.subplots_adjust(
    left=0.06,
    right=0.98,
    bottom=0.06,
    top=0.92,
    wspace=0.18,
    hspace=0.22
)

for i, sdata_key in enumerate(sample_order):
    sdata = sdata_dict[sdata_key]
    sample_name = sdata_key.replace("sdata_", "")

    print(f"Processing {sample_name} ...")

    adata = sdata.tables["table"].copy()
    adata = adata[adata.obs["annotation"].notna()].copy()
    adata.obs["annotation"] = adata.obs["annotation"].astype("category")

    sq.gr.spatial_neighbors(
        adata,
        coord_type="generic",
        delaunay=True
    )

    sq.gr.nhood_enrichment(
        adata,
        cluster_key="annotation",
        n_jobs=1,
        show_progress_bar=False
    )

    sq.pl.nhood_enrichment(
        adata,
        cluster_key="annotation",
        ax=axes[i],
        cmap="coolwarm",
        title=sample_name,
        show=False
    )

for ax_ in fig.axes:
    if ax_.get_ylabel() == "annotation":
        ax_.set_ylabel("")
        ax_.yaxis.label.set_visible(False)

for j in range(len(sample_order), len(axes)):
    axes[j].axis("off")

plt.savefig(
    os.path.join(outdir, "SM_T_cell_spatial_co_localization_plot.pdf"),
    bbox_inches="tight"
)
plt.savefig(
    os.path.join(outdir, "SM_T_cell_spatial_co_localization_plot.svg"),
    bbox_inches="tight"
)
plt.close()
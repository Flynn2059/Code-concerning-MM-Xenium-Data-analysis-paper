#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 12:39:07 2026

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
from sklearn.neighbors import KDTree

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

    # 先创建一列，默认缺失
    adata.obs["annotation"] = pd.NA

    # 取共同细胞
    common_idx = adata.obs_names.intersection(annotation_sub.index)

    # 把注释写回去
    adata.obs.loc[common_idx, "annotation"] = annotation_sub.loc[common_idx, "major_annotation"]

    # 放回 sdata
    sdata["table"] = adata

    print(f"{sample_name}: matched {len(common_idx)} / {adata.n_obs} cells")

for sdata_key, sdata in sdata_dict.items():
    sample_name = sdata_key.replace("sdata_", "")
    print(f"\n{sample_name}")
    print(sdata["table"].obs["annotation"].value_counts(dropna=False))
#%% 查看 t_norm 和 p_norm 的峰形分布
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree

outdir = "./Tcell_Plasma_peak_inspection"
os.makedirs(outdir, exist_ok=True)

radius = 200

t_labels = ["CD4+ T Naive", "CD8+ T Effector", "T Proliferating"]
plasma_labels = ["Plasma"]

all_samples_df = []

for sdata_key, sdata in sdata_dict.items():
    sample_name = sdata_key.replace("sdata_", "")
    prefix = sample_name + "_"

    print(f"Processing {sample_name} ...")

    adata = sdata.tables["table"].copy()
    gdf = sdata.shapes["cell_boundaries"].copy()
    gdf.index = gdf.index.astype(str)

    # 恢复和 shapes 对应的原始 cell_id
    adata.obs["cell_id"] = adata.obs_names.to_series().str.replace(f"^{prefix}", "", regex=True)

    # annotation 映射到 shapes
    ann_map = adata.obs.set_index("cell_id")["annotation"]
    gdf["annotation"] = gdf.index.map(ann_map)

    # centroid 坐标
    centroids = gdf.geometry.centroid
    coords = np.column_stack([centroids.x.values, centroids.y.values])

    # KDTree
    tree = KDTree(coords)
    neighbor_indices = tree.query_radius(coords, r=radius)

    # annotation 数组
    ann_array = gdf["annotation"].astype("object").values
    is_t = np.isin(ann_array, t_labels)
    is_plasma = np.isin(ann_array, plasma_labels)

    # 计算 local density（沿用你前面面积定义）
    t_density = np.zeros(len(gdf), dtype=float)
    plasma_density = np.zeros(len(gdf), dtype=float)

    area = np.pi * radius * radius

    for i, nbr_idx in enumerate(neighbor_indices):
        t_count = is_t[nbr_idx].sum()
        plasma_count = is_plasma[nbr_idx].sum()

        t_density[i] = t_count / area
        plasma_density[i] = plasma_count / area

    gdf["local_t_density"] = t_density
    gdf["local_plasma_density"] = plasma_density

    # 归一化
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

    # 保存单样本表
    plot_df = gdf[["t_norm", "p_norm"]].copy()
    plot_df["sample"] = sample_name
    all_samples_df.append(plot_df)

    # -----------------------------
    # 单样本作图
    # -----------------------------
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), dpi=300)

    # 1. t_norm histogram
    axes[0].hist(gdf["t_norm"], bins=100)
    axes[0].set_title(f"{sample_name} - t_norm")
    axes[0].set_xlabel("t_norm")
    axes[0].set_ylabel("Cell count")

    # 2. p_norm histogram
    axes[1].hist(gdf["p_norm"], bins=100)
    axes[1].set_title(f"{sample_name} - p_norm")
    axes[1].set_xlabel("p_norm")
    axes[1].set_ylabel("Cell count")

    # 3. 2D histogram
    h = axes[2].hist2d(
        gdf["t_norm"],
        gdf["p_norm"],
        bins=100,
        range=[[0, 1], [0, 1]]
    )
    axes[2].set_title(f"{sample_name} - t_norm vs p_norm")
    axes[2].set_xlabel("t_norm")
    axes[2].set_ylabel("p_norm")
    plt.colorbar(h[3], ax=axes[2], label="Cell count")

    plt.tight_layout()
    plt.savefig(
        os.path.join(outdir, f"{sample_name}_t_norm_p_norm_peak_inspection.pdf"),
        bbox_inches="tight"
    )
    plt.savefig(
        os.path.join(outdir, f"{sample_name}_t_norm_p_norm_peak_inspection.svg"),
        bbox_inches="tight"
    )
    plt.close()
print("Done.")

#%% 基于 t_norm 和 p_norm 的 GMM 聚类 + 空间可视化（per-sample 多核版）
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.neighbors import KDTree
from sklearn.mixture import GaussianMixture
from joblib import Parallel, delayed

outdir = "./Tcell_Plasma_GMM_clustering"
os.makedirs(outdir, exist_ok=True)

radius = 200
area = np.pi * radius * radius

t_labels = ["CD4+ T Naive", "CD8+ T Effector", "T Proliferating"]
plasma_labels = ["Plasma"]

n_components_list = [4]

# 固定只用 6 个核心
N_JOBS = 6


def process_one_sample(sdata_key, sdata):
    sample_name = sdata_key.replace("sdata_", "")
    prefix = sample_name + "_"

    print(f"Processing {sample_name} ...")

    adata = sdata.tables["table"].copy()
    gdf = sdata.shapes["cell_boundaries"].copy()
    gdf.index = gdf.index.astype(str)

    # 恢复和 shapes 对应的原始 cell_id
    adata.obs["cell_id"] = adata.obs_names.to_series().str.replace(f"^{prefix}", "", regex=True)

    # annotation 映射到 shapes
    ann_map = adata.obs.set_index("cell_id")["annotation"]
    gdf["annotation"] = gdf.index.map(ann_map)

    # centroid
    centroids = gdf.geometry.centroid
    coords = np.column_stack([centroids.x.values, centroids.y.values])

    # KDTree
    tree = KDTree(coords)
    neighbor_indices = tree.query_radius(coords, r=radius)

    # annotation array
    ann_array = gdf["annotation"].astype("object").values
    is_t = np.isin(ann_array, t_labels)
    is_plasma = np.isin(ann_array, plasma_labels)

    # local density
    t_density = np.zeros(len(gdf), dtype=float)
    plasma_density = np.zeros(len(gdf), dtype=float)

    for i, nbr_idx in enumerate(neighbor_indices):
        t_count = is_t[nbr_idx].sum()
        plasma_count = is_plasma[nbr_idx].sum()

        t_density[i] = t_count / area
        plasma_density[i] = plasma_count / area

    gdf["local_t_density"] = t_density
    gdf["local_plasma_density"] = plasma_density

    # normalize
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

    # 写回 adata
    density_t_series = pd.Series(t_density, index=gdf.index, name="local_t_density")
    density_plasma_series = pd.Series(plasma_density, index=gdf.index, name="local_plasma_density")

    adata.obs["local_t_density"] = adata.obs["cell_id"].map(density_t_series)
    adata.obs["local_plasma_density"] = adata.obs["cell_id"].map(density_plasma_series)

    # GMM 输入
    X = gdf[["t_norm", "p_norm"]].values

    bic_records = []

    for n_components in n_components_list:
        print(f"  {sample_name}: Fitting GMM with n_components = {n_components}")

        # 每个 K 单独一个子文件夹
        k_outdir = os.path.join(outdir, f"K{n_components}")
        os.makedirs(k_outdir, exist_ok=True)

        gmm = GaussianMixture(
            n_components=n_components,
            covariance_type="full",
            random_state=0
        )
        gmm.fit(X)

        cluster = gmm.predict(X)
        prob = gmm.predict_proba(X)
        max_prob = prob.max(axis=1)

        bic = gmm.bic(X)
        aic = gmm.aic(X)
        bic_records.append({
            "sample": sample_name,
            "n_components": n_components,
            "BIC": bic,
            "AIC": aic
        })

        cluster_col = f"gmm_{n_components}cluster"
        conf_col = f"gmm_{n_components}cluster_confidence"

        gdf[cluster_col] = cluster.astype(str)
        gdf[conf_col] = max_prob

        # 写回 adata
        cluster_series = pd.Series(cluster.astype(str), index=gdf.index, name=cluster_col)
        conf_series = pd.Series(max_prob, index=gdf.index, name=conf_col)

        adata.obs[cluster_col] = adata.obs["cell_id"].map(cluster_series)
        adata.obs[conf_col] = adata.obs["cell_id"].map(conf_series)

        # cluster centers
        centers = pd.DataFrame(
            gmm.means_,
            columns=["t_norm_center", "p_norm_center"]
        )
        centers["cluster"] = centers.index.astype(str)
        centers = centers[["cluster", "t_norm_center", "p_norm_center"]]
        centers["sample"] = sample_name
        centers["n_components"] = n_components

        centers.to_csv(
            os.path.join(k_outdir, f"{sample_name}_gmm_{n_components}cluster_centers.csv"),
            index=False
        )

        print(centers)

        # -----------------------------
        # 图1：二维散点图（按 cluster 着色）
        # -----------------------------
        fig, ax = plt.subplots(figsize=(6, 6), dpi=300)

        ax.scatter(
            gdf["t_norm"],
            gdf["p_norm"],
            c=cluster,
            s=2
        )

        ax.scatter(
            gmm.means_[:, 0],
            gmm.means_[:, 1],
            marker="x",
            s=100
        )

        for i, (cx, cy) in enumerate(gmm.means_):
            ax.text(cx, cy, str(i), fontsize=10)

        ax.set_xlabel("t_norm")
        ax.set_ylabel("p_norm")
        ax.set_title(f"{sample_name} | GMM {n_components} clusters")

        plt.tight_layout()
        plt.savefig(
            os.path.join(k_outdir, f"{sample_name}_gmm_{n_components}cluster_scatter.pdf"),
            bbox_inches="tight"
        )
        plt.savefig(
            os.path.join(k_outdir, f"{sample_name}_gmm_{n_components}cluster_scatter.svg"),
            bbox_inches="tight"
        )
        plt.close()

        # -----------------------------
        # 图2：空间图（按 cluster 着色）
        # -----------------------------
        fig, ax = plt.subplots(figsize=(12, 12), dpi=600)
        gdf.plot(ax=ax, column=cluster_col, categorical=True, legend=True, linewidth=0)

        ax.set_title(f"{sample_name} | GMM {n_components} clusters", fontsize=16)
        ax.set_axis_off()

        plt.tight_layout()
        plt.savefig(
            os.path.join(k_outdir, f"{sample_name}_gmm_{n_components}cluster_spatial.pdf"),
            bbox_inches="tight"
        )
        plt.savefig(
            os.path.join(k_outdir, f"{sample_name}_gmm_{n_components}cluster_spatial.svg"),
            bbox_inches="tight"
        )
        plt.close()

    # 保存 BIC/AIC（总表仍放在主目录）
    bic_df = pd.DataFrame(bic_records)
    bic_df.to_csv(
        os.path.join(outdir, f"{sample_name}_gmm_model_selection.csv"),
        index=False
    )

    print(bic_df.sort_values("BIC"))

    # 返回结果，主进程再写回
    return sdata_key, adata


# -----------------------------
# per-sample 并行
# -----------------------------
results = Parallel(n_jobs=N_JOBS, backend="loky")(
    delayed(process_one_sample)(sdata_key, sdata)
    for sdata_key, sdata in sdata_dict.items()
)

# 主进程回写
for sdata_key, adata_res in results:
    sdata_dict[sdata_key]["table"] = adata_res

print("All samples finished.")
#%% 基于 GMM center 自动注释 cluster 类型 + 条件式人工确认 + Engaging Zone + 重绘 scatter/spatial 图
import os
import ast
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from sklearn.neighbors import KDTree
#from sklearn.mixture import GaussianMixture

outdir = "./Tcell_Plasma_GMM_clustering_annotation"
os.makedirs(outdir, exist_ok=True)

radius = 200
area = np.pi * radius * radius

t_labels = ["CD4+ T Naive", "CD8+ T Effector", "T Proliferating"]
plasma_labels = ["Plasma"]

n_components_list = [4]

# 颜色定义
label_color_map = {
    "T_Dominant": "red",
    "Plasma_Dominant": "blue",
    "Mixed": "green",
    "Engaging Zone": "darkorange"
}


def get_manual_override_if_needed(sample_name, centers, auto_t_cluster, auto_p_cluster):
    """
    只有当 auto_t_cluster == auto_p_cluster 时，才触发人工确认。
    返回:
        final_t_cluster, final_p_cluster, used_manual_override
    """
    if str(auto_t_cluster) != str(auto_p_cluster):
        return str(auto_t_cluster), str(auto_p_cluster), False

    print("\n" + "=" * 80)
    print(f"Sample: {sample_name}")
    print("Conflict detected: auto-assigned T_Dominant and Plasma_Dominant are the same cluster.")
    print("\nGMM cluster centers:")
    print(centers[["cluster", "t_norm_center", "p_norm_center"]].to_string(index=False))
    print("\nAuto assignment:")
    print(f"  T_Dominant      -> cluster {auto_t_cluster}")
    print(f"  Plasma_Dominant -> cluster {auto_p_cluster}")
    print("\nPress Enter to keep this result.")
    print("Or input a list like [T_cluster, Plasma_cluster], e.g. [2, 1]")

    while True:
        user_input = input("Your choice: ").strip()

        if user_input == "":
            print(f"Received. Using parameters: T_Dominant={auto_t_cluster}, Plasma_Dominant={auto_p_cluster}")
            return str(auto_t_cluster), str(auto_p_cluster), False

        try:
            parsed = ast.literal_eval(user_input)

            if isinstance(parsed, (list, tuple)) and len(parsed) == 2:
                t_manual = str(parsed[0])
                p_manual = str(parsed[1])

                valid_clusters = set(centers["cluster"].astype(str))
                if t_manual not in valid_clusters or p_manual not in valid_clusters:
                    print(f"Invalid cluster id. Valid cluster ids are: {sorted(valid_clusters)}")
                    continue

                print(f"Received. Using parameters: T_Dominant={t_manual}, Plasma_Dominant={p_manual}")
                return t_manual, p_manual, True
            else:
                print("Input format error. Please input like [2, 1], or press Enter.")
        except Exception:
            print("Cannot parse your input. Please input like [2, 1], or press Enter.")


def build_cluster_annotation(centers, final_t_cluster, final_p_cluster):
    """
    根据最终确认的 cluster 编号构建 cluster_annotation。
    """
    centers = centers.copy()
    centers["cluster_annotation"] = "Mixed"

    if str(final_t_cluster) == str(final_p_cluster):
        centers.loc[centers["cluster"] == str(final_t_cluster), "cluster_annotation"] = "Engaging Zone"
    else:
        centers.loc[centers["cluster"] == str(final_t_cluster), "cluster_annotation"] = "T_Dominant"
        centers.loc[centers["cluster"] == str(final_p_cluster), "cluster_annotation"] = "Plasma_Dominant"

    return centers


def make_legend_handles(categories_present):
    ordered_categories = ["T_Dominant", "Plasma_Dominant", "Mixed", "Engaging Zone"]
    handles = []

    for cat in ordered_categories:
        if cat in categories_present:
            handles.append(
                Patch(facecolor=label_color_map[cat], edgecolor="none", label=cat)
            )
    return handles


def process_one_sample_annotation(sdata_key, sdata):
    sample_name = sdata_key.replace("sdata_", "")
    prefix = sample_name + "_"

    print(f"\nProcessing {sample_name} ...")

    adata = sdata.tables["table"].copy()
    gdf = sdata.shapes["cell_boundaries"].copy()
    gdf.index = gdf.index.astype(str)

    # 恢复和 shapes 对应的原始 cell_id
    adata.obs["cell_id"] = adata.obs_names.to_series().str.replace(f"^{prefix}", "", regex=True)

    # annotation 映射到 shapes
    ann_map = adata.obs.set_index("cell_id")["annotation"]
    gdf["annotation"] = gdf.index.map(ann_map)

    # centroid
    centroids = gdf.geometry.centroid
    coords = np.column_stack([centroids.x.values, centroids.y.values])

    # KDTree
    tree = KDTree(coords)
    neighbor_indices = tree.query_radius(coords, r=radius)

    # annotation array
    ann_array = gdf["annotation"].astype("object").values
    is_t = np.isin(ann_array, t_labels)
    is_plasma = np.isin(ann_array, plasma_labels)

    # local density
    t_density = np.zeros(len(gdf), dtype=float)
    plasma_density = np.zeros(len(gdf), dtype=float)

    for i, nbr_idx in enumerate(neighbor_indices):
        t_count = is_t[nbr_idx].sum()
        plasma_count = is_plasma[nbr_idx].sum()

        t_density[i] = t_count / area
        plasma_density[i] = plasma_count / area

    gdf["local_t_density"] = t_density
    gdf["local_plasma_density"] = plasma_density

    # normalize
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

    # 写回 adata
    density_t_series = pd.Series(t_density, index=gdf.index, name="local_t_density")
    density_plasma_series = pd.Series(plasma_density, index=gdf.index, name="local_plasma_density")

    adata.obs["local_t_density"] = adata.obs["cell_id"].map(density_t_series)
    adata.obs["local_plasma_density"] = adata.obs["cell_id"].map(density_plasma_series)

    # GMM 输入
    X = gdf[["t_norm", "p_norm"]].values

    # 这里准备一个变量，最后导出细胞级矩阵
    cell_niche_df = None

    for n_components in n_components_list:
        print(f"  {sample_name}: Fitting GMM with n_components = {n_components}")

        k_outdir = os.path.join(outdir, f"K{n_components}")
        os.makedirs(k_outdir, exist_ok=True)

        gmm = GaussianMixture(
            n_components=n_components,
            covariance_type="full",
            random_state=0
        )
        gmm.fit(X)

        cluster = gmm.predict(X)
        prob = gmm.predict_proba(X)
        max_prob = prob.max(axis=1)

        cluster_col = f"gmm_{n_components}cluster"
        conf_col = f"gmm_{n_components}cluster_confidence"
        anno_col = f"gmm_{n_components}cluster_annotation"

        gdf[cluster_col] = cluster.astype(str)
        gdf[conf_col] = max_prob

        # cluster centers
        centers = pd.DataFrame(
            gmm.means_,
            columns=["t_norm_center", "p_norm_center"]
        )
        centers["cluster"] = centers.index.astype(str)

        # 自动命名
        auto_t_cluster = str(centers.loc[centers["t_norm_center"].idxmax(), "cluster"])
        auto_p_cluster = str(centers.loc[centers["p_norm_center"].idxmax(), "cluster"])

        # 仅在冲突时人工确认
        final_t_cluster, final_p_cluster, used_manual_override = get_manual_override_if_needed(
            sample_name=sample_name,
            centers=centers,
            auto_t_cluster=auto_t_cluster,
            auto_p_cluster=auto_p_cluster
        )

        # 构建最终注释
        centers = build_cluster_annotation(
            centers=centers,
            final_t_cluster=final_t_cluster,
            final_p_cluster=final_p_cluster
        )

        cluster_to_annotation = dict(zip(centers["cluster"], centers["cluster_annotation"]))
        gdf[anno_col] = gdf[cluster_col].map(cluster_to_annotation)

        centers["sample"] = sample_name
        centers["n_components"] = n_components
        centers["auto_t_cluster"] = auto_t_cluster
        centers["auto_p_cluster"] = auto_p_cluster
        centers["final_t_cluster"] = final_t_cluster
        centers["final_p_cluster"] = final_p_cluster
        centers["used_manual_override"] = used_manual_override

        centers.to_csv(
            os.path.join(k_outdir, f"{sample_name}_gmm_{n_components}cluster_centers_annotated.csv"),
            index=False
        )

        # 写回 adata
        cluster_series = pd.Series(cluster.astype(str), index=gdf.index, name=cluster_col)
        conf_series = pd.Series(max_prob, index=gdf.index, name=conf_col)
        anno_series = pd.Series(gdf[anno_col].values, index=gdf.index, name=anno_col)

        adata.obs[cluster_col] = adata.obs["cell_id"].map(cluster_series)
        adata.obs[conf_col] = adata.obs["cell_id"].map(conf_series)
        adata.obs[anno_col] = adata.obs["cell_id"].map(anno_series)

        # ========= 新增：生成这个样本的细胞级 niche 矩阵 =========
        cell_niche_df = adata.obs.copy()
        cell_niche_df["cell_full_id"] = cell_niche_df.index.astype(str)
        cell_niche_df["sample"] = sample_name
        cell_niche_df["niche"] = cell_niche_df[anno_col]

        # 这里只保留你最关心的列；如果你想保留更多列也可以继续加
        cell_niche_df = cell_niche_df[
            ["cell_full_id", "cell_id", "sample", "annotation", cluster_col, conf_col, "niche"]
        ].copy()

        # 当前样本实际出现的类别
        categories_present = pd.unique(gdf[anno_col])
        legend_handles = make_legend_handles(categories_present)

        # -----------------------------
        # 图1：scatter 图
        # -----------------------------
        fig, ax = plt.subplots(figsize=(6, 6), dpi=300)

        for label in ["T_Dominant", "Plasma_Dominant", "Mixed", "Engaging Zone"]:
            sub = gdf[gdf[anno_col] == label]
            if len(sub) > 0:
                ax.scatter(
                    sub["t_norm"],
                    sub["p_norm"],
                    c=label_color_map[label],
                    s=2,
                    alpha=0.8
                )

        ax.set_xlabel("t_norm")
        ax.set_ylabel("p_norm")
        ax.set_title(f"{sample_name} | GMM {n_components} annotated clusters")
        ax.legend(handles=legend_handles, frameon=False)

        plt.tight_layout()
        plt.savefig(
            os.path.join(k_outdir, f"{sample_name}_gmm_{n_components}cluster_scatter_annotated.pdf"),
            bbox_inches="tight"
        )
        plt.savefig(
            os.path.join(k_outdir, f"{sample_name}_gmm_{n_components}cluster_scatter_annotated.svg"),
            bbox_inches="tight"
        )
        plt.close()

        # -----------------------------
        # 图2：spatial 图
        # -----------------------------
        fig, ax = plt.subplots(figsize=(12, 12), dpi=600)

        for label in ["T_Dominant", "Plasma_Dominant", "Mixed", "Engaging Zone"]:
            sub = gdf[gdf[anno_col] == label]
            if len(sub) > 0:
                sub.plot(
                    ax=ax,
                    color=label_color_map[label],
                    linewidth=0
                )

        ax.set_title(f"{sample_name} | GMM {n_components} annotated clusters", fontsize=16)
        ax.set_axis_off()
        ax.legend(handles=legend_handles, frameon=False, loc="upper right")

        plt.tight_layout()
        plt.savefig(
            os.path.join(k_outdir, f"{sample_name}_gmm_{n_components}cluster_spatial_annotated.pdf"),
            bbox_inches="tight"
        )
        plt.savefig(
            os.path.join(k_outdir, f"{sample_name}_gmm_{n_components}cluster_spatial_annotated.svg"),
            bbox_inches="tight"
        )
        plt.close()

    return adata, cell_niche_df


# -----------------------------
# 串行处理（保留，因为仍有条件式终端交互）
# -----------------------------
all_cell_niche_df = []

for sdata_key, sdata in sdata_dict.items():
    adata_res, cell_niche_df = process_one_sample_annotation(sdata_key, sdata)
    sdata_dict[sdata_key]["table"] = adata_res
    all_cell_niche_df.append(cell_niche_df)

# 合并所有样本，一个总矩阵
all_cell_niche_df = pd.concat(all_cell_niche_df, axis=0)
all_cell_niche_df.index = range(len(all_cell_niche_df))

all_cell_niche_df.to_csv(
    os.path.join(outdir, "SM_5_samples_cell_niche_matrix.csv"),
    index=False
)

print("All samples finished.")
print("Cell-level niche matrix saved:")
print(os.path.join(outdir, "SM_5_samples_cell_niche_matrix.csv"))

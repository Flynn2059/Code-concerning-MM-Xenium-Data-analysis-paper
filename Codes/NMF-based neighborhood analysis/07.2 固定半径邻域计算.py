#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 16:03:50 2026

@author: flynn
"""

#%% imports
import os
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

# some customized settings for better visualization
plt.rcParams["figure.figsize"] = (12, 12)
plt.rcParams["figure.dpi"] = 600

#%% 工作路径
os.chdir('/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/')
print(os.getcwd())
print("\n")
print(os.listdir())

#%% 读取之前的邻域信息
ctrl = pd.read_csv("./4_ctrl_coordinate_annotation_matrix.csv", index_col=0)
mgus = pd.read_csv("./2_mgus_coordinate_annotation_matrix.csv", index_col=0)
sm = pd.read_csv("./5_sm_coordinate_annotation_matrix.csv", index_col=0)
mm = pd.read_csv("./10_mm_coordinate_annotation_matrix.csv", index_col=0)

#%% 定义半径
radii = [50, 75, 100]

#%% 按组存数据
group_dfs = {
    "Ctrl": ctrl.copy(),
    "MGUS": mgus.copy(),
    "SM": sm.copy(),
    "MM": mm.copy()
}

for group_name in group_dfs:
    group_dfs[group_name]["cell_id"] = group_dfs[group_name].index.astype(str)

#%% 计算单个样本在某个半径下的邻域结果 + 相对样本组成的倍数矩阵
def calculate_neighborhood_for_sample(sample_df, radius):
    sample_df = sample_df.copy()
    sample_df = sample_df.reset_index(drop=True)

    coords = sample_df[["X", "Y"]].to_numpy()
    annots = sample_df["annotation"].astype(str).to_numpy()
    sample_name = sample_df["sample"].iloc[0]

    tree = cKDTree(coords)
    neighbor_indices = tree.query_ball_point(coords, r=radius)

    all_annotations = sorted(sample_df["annotation"].astype(str).unique().tolist())

    sample_annot_counts = sample_df["annotation"].astype(str).value_counts()
    sample_annot_pct = {
        annot: sample_annot_counts.get(annot, 0) / len(sample_df)
        for annot in all_annotations
    }

    results = []

    for i, idxs in enumerate(neighbor_indices):
        idxs = [j for j in idxs if j != i]
        n_neighbors = len(idxs)
        row = {
            "cell_id": sample_df.loc[i, "cell_id"],
            "sample": sample_name,
            "celltype": sample_df.loc[i, "annotation"],
            "n_neighbors": n_neighbors
        }

        if len(idxs) == 0:
            for annot in all_annotations:
                row[f"pct_{annot}"] = 0
                row[f"fold_{annot}"] = 0
        else:
            neighbor_annots = annots[idxs]
            annot_counts = pd.Series(neighbor_annots).value_counts()

            for annot in all_annotations:
                count = int(annot_counts.get(annot, 0))
                neighborhood_pct = count / len(idxs)
                sample_pct = sample_annot_pct[annot]

                row[f"pct_{annot}"] = neighborhood_pct
                row[f"fold_{annot}"] = neighborhood_pct / sample_pct if sample_pct > 0 else 0

        results.append(row)

    result_df = pd.DataFrame(results)
    return result_df

#%% 对每个组分别计算
all_results = {}

for group_name, group_df in group_dfs.items():
    print(f"Processing group: {group_name}")

    samples = group_df["sample"].unique().tolist()
    all_results[group_name] = {}

    for radius in radii:
        print(f"  Radius: {radius} um")

        sample_results = []

        for sample_name in samples:
            print(f"    Sample: {sample_name}")

            sample_df = group_df[group_df["sample"] == sample_name].copy()
            one_result = calculate_neighborhood_for_sample(sample_df, radius)
            sample_results.append(one_result)

        group_radius_result = pd.concat(sample_results, axis=0, ignore_index=True)
        all_results[group_name][radius] = group_radius_result

#%% 保存百分比矩阵和倍数矩阵
output_dir = "./neighborhood_results"
os.makedirs(output_dir, exist_ok=True)

for group_name in all_results:
    for radius in all_results[group_name]:
        result_df = all_results[group_name][radius]

        meta_cols = ["cell_id", "sample", "celltype", "n_neighbors"]
        pct_cols = meta_cols + [
            col for col in result_df.columns if col.startswith("pct_")
        ]
        fold_cols = meta_cols + [
            col for col in result_df.columns if col.startswith("fold_")
        ]

        pct_matrix = result_df[pct_cols]
        fold_matrix = result_df[fold_cols]

        pct_out_path = os.path.join(
            output_dir,
            f"{group_name}_{radius}um_pct_matrix.csv"
        )
        fold_out_path = os.path.join(
            output_dir,
            f"{group_name}_{radius}um_fold_matrix.csv"
        )

        pct_matrix.to_csv(pct_out_path, index=False)
        fold_matrix.to_csv(fold_out_path, index=False)

        print(f"Saved: {pct_out_path}")
        print(f"Saved: {fold_out_path}")
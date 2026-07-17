#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 14:28:53 2026

@author: flynn
"""

#%% imports
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import seaborn as sns
from sklearn.decomposition import NMF

# some customized settings for better visualization
plt.rcParams["figure.figsize"] = (12, 12)
plt.rcParams["figure.dpi"] = 600
sc.settings.verbosity = 0


plot_dir="./neighborhood_results/plots/"
#%% helper
def format_resolution(resolution):
    return str(resolution).replace(".", "p")

#%% 工作路径
os.chdir('/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/')
print(os.getcwd())
print("\n")
print(os.listdir('./neighborhood_results'))

#%% 读取一个分组的邻域矩阵
mgus_2_samples = pd.read_csv('./neighborhood_results/MGUS_50um_pct_matrix.csv', index_col=0)
mgus_2_samples = mgus_2_samples[
    mgus_2_samples["celltype"].isin(["CD4+ T Central Memory", "pct_CD8+ T Effector",
                                     'pct_CD8+ T Exhausted','pct_T Proliferating']) &
    (mgus_2_samples["n_neighbors"] > 5)
].copy()

print(mgus_2_samples["sample"].value_counts())

# 提取 pct 列
pct_cols = [col for col in mgus_2_samples.columns if col.startswith("pct_")]

# 获取样本列表
sample_names = sorted(mgus_2_samples["sample"].unique().tolist())
print(sample_names)
#%% 绘制特征分布的直方图
# 每个样本画一张图：每个 pct 特征一个子图，使用直方图
import math

output_dir = "./pct_feature_distribution_hist_plots"
os.makedirs(output_dir, exist_ok=True)

for sample_name in sample_names:
    print(f"Plotting pct feature histograms for sample: {sample_name}")
    
    sample_df = mgus_2_samples[mgus_2_samples["sample"] == sample_name].copy()
    
    n_features = len(pct_cols)
    n_cols = 4
    n_rows = math.ceil(n_features / n_cols)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    axes = axes.flatten()
    
    for i, pct_col in enumerate(pct_cols):
        ax = axes[i]
        
        values = sample_df[pct_col].dropna()
        ax.hist(values, bins=30)
        ax.set_title(pct_col, fontsize=10)
        ax.set_xlabel("")
        ax.set_ylabel("Cell count")
        ax.tick_params(axis="x", labelrotation=45)
    
    for j in range(n_features, len(axes)):
        fig.delaxes(axes[j])
    
    fig.suptitle(f"{sample_name} pct feature distributions", fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    
    save_path = os.path.join(output_dir, f"{sample_name}_pct_feature_histograms.svg")
    fig.savefig(save_path, bbox_inches="tight")
    plt.close(fig)
    
    print(f"Saved: {save_path}")
#%% 参数设置
k_values = [3, 4, 5]
# first-pass resolution grid
resolutions = [0.05,0.1,0.2]

# 对每个样本分别跑 NMF 和 Leiden
nmf_results_by_sample = {}
ctrl_labeled_by_k_res = {}
all_label_tables = []
all_cluster_summaries = []

for sample_name in sample_names:
    print(f"\nProcessing sample: {sample_name}")
    
    sample_df = mgus_2_samples[mgus_2_samples["sample"] == sample_name].copy()
    X_raw = sample_df[pct_cols].copy()
    X = np.log1p(X_raw * 100)

    X = pd.DataFrame(
        X,
        index=X_raw.index,
        columns=X_raw.columns
    )
    
    nmf_results_by_sample[sample_name] = {}
    
    for k in k_values:
        print(f"Running NMF for {sample_name}, n_components={k} ...")
        
        model = NMF(
            n_components=k,
            init="nndsvda",
            random_state=42,
            max_iter=1000
        )
        
        W = model.fit_transform(X)
        H = model.components_
        
        W_df = pd.DataFrame(
            W,
            index=X.index.astype(str),
            columns=[f"NMF{k}_component_{i+1}" for i in range(k)]
        )
        
        H_df = pd.DataFrame(
            H,
            index=[f"NMF{k}_component_{i+1}" for i in range(k)],
            columns=X.columns
        )
        
        sample_df_base = sample_df.copy()
        sample_df_base["cell_id"] = sample_df_base.index.astype(str)
        sample_df_base["k"] = k
        
        sample_df_wide = sample_df_base.copy()
        
        # Leiden on W space
        adata = sc.AnnData(W_df.copy())
        adata.obs_names = W_df.index.astype(str)
        adata.var_names = W_df.columns.astype(str)
        
        sc.pp.neighbors(adata)
        
        leiden_summary = {}
        
        for res in resolutions:
            res_tag = format_resolution(res)
            leiden_key = f"leiden_r{res_tag}"
            
            print(
                f"Running Leiden for {sample_name}, "
                f"k={k}, resolution={res} ..."
            )
            
            sc.tl.leiden(adata, resolution=res, key_added=leiden_key)
            
            # 原始 Leiden cluster 编号
            labels = adata.obs[leiden_key].astype(str).copy()
            label_map_raw = labels.to_dict()
            
            # 加上 sample/k/resolution 前缀，避免不同样本 cluster 编号冲突
            labels_unique = pd.Series(
                [
                    f"{sample_name}_k{k}_r{res_tag}_cluster{x}"
                    for x in labels.tolist()
                ],
                index=labels.index
            )
            label_map_unique = labels_unique.to_dict()
            
            # 宽表里同时保留原始 cluster id 和唯一标签
            sample_df_wide[f"{leiden_key}_cluster_id"] = sample_df_wide["cell_id"].map(label_map_raw)
            sample_df_wide[leiden_key] = sample_df_wide["cell_id"].map(label_map_unique)
            
            n_clusters = labels.nunique()
            cluster_sizes = labels.value_counts().sort_values(ascending=False)
            
            leiden_summary[res] = {
                "n_clusters": int(n_clusters),
                "cluster_sizes": cluster_sizes.to_dict()
            }
            
            # 长表
            label_df = sample_df_base[["cell_id", "sample", "k"]].copy()
            label_df["resolution"] = res
            label_df["cluster_id"] = label_df["cell_id"].map(label_map_raw)
            label_df["neighborhood_label"] = label_df["cell_id"].map(label_map_unique)
            all_label_tables.append(label_df)
            
            # 按 k + resolution 收集完整结果
            if k not in ctrl_labeled_by_k_res:
                ctrl_labeled_by_k_res[k] = {}
            if res not in ctrl_labeled_by_k_res[k]:
                ctrl_labeled_by_k_res[k][res] = []
            
            labeled_res_df = sample_df_base.copy()
            labeled_res_df["resolution"] = res
            labeled_res_df["cluster_id"] = labeled_res_df["cell_id"].map(label_map_raw)
            labeled_res_df["neighborhood_label"] = labeled_res_df["cell_id"].map(label_map_unique)
            ctrl_labeled_by_k_res[k][res].append(labeled_res_df)
            
            # 生成 cluster summary
            cluster_summary_df = (
                labels.value_counts()
                .rename_axis("cluster_id")
                .reset_index(name="n_cells")
                .sort_values("n_cells", ascending=False)
                .reset_index(drop=True)
            )
            cluster_summary_df["sample"] = sample_name
            cluster_summary_df["k"] = k
            cluster_summary_df["resolution"] = res
            cluster_summary_df["unique_cluster_label"] = cluster_summary_df["cluster_id"].apply(
                lambda x: f"{sample_name}_k{k}_r{res_tag}_cluster{x}"
            )
            all_cluster_summaries.append(cluster_summary_df)
            
            print(
                f"Leiden done: sample={sample_name}, k={k}, resolution={res}, "
                f"clusters={n_clusters}, largest={int(cluster_sizes.iloc[0])}, "
                f"smallest={int(cluster_sizes.iloc[-1])}"
            )
        
        nmf_results_by_sample[sample_name][k] = {
            "model": model,
            "W": W_df,
            "H": H_df,
            "labeled_df_wide": sample_df_wide,
            "leiden_summary": leiden_summary,
            "reconstruction_err": model.reconstruction_err_
        }
        
        print(
            f"Done: sample={sample_name}, k={k}, "
            f"reconstruction error={model.reconstruction_err_:.4f}"
        )

# 合并每个 k + resolution 下所有样本的结果
for k in ctrl_labeled_by_k_res:
    for res in ctrl_labeled_by_k_res[k]:
        if isinstance(ctrl_labeled_by_k_res[k][res], list):
            ctrl_labeled_by_k_res[k][res] = pd.concat(
                ctrl_labeled_by_k_res[k][res],
                axis=0,
                ignore_index=True
            )

# 合并所有样本所有 k + resolution 的标签长表
ctrl_labels_long = pd.concat(all_label_tables, axis=0, ignore_index=True)

# 合并所有 cluster summary
cluster_summary_long = pd.concat(all_cluster_summaries, axis=0, ignore_index=True)

#%% 画每个 sample 的 NMF H loading heatmap
selected_k = sorted(k_values)

for k in selected_k:
    fig, axes = plt.subplots(2, 1, figsize=(11, 18), dpi=600)
    axes = axes.flatten()

    for ax, sample_name in zip(axes, sample_names):
        if sample_name not in nmf_results_by_sample:
            ax.axis("off")
            ax.set_title(f"{sample_name}\nNo NMF result")
            continue

        if k not in nmf_results_by_sample[sample_name]:
            ax.axis("off")
            ax.set_title(f"{sample_name}\nNo k={k} result")
            continue

        H_df = nmf_results_by_sample[sample_name][k]["H"].copy()

        H_plot = H_df.copy()
        H_plot.columns = [col.replace("pct_", "") for col in H_plot.columns]

        sns.heatmap(
            H_plot,
            ax=ax,
            cmap="coolwarm",
            cbar=True,
            linewidths=0.2,
            linecolor="white"
        )

        ax.set_title(f"{sample_name} | NMF H loadings | k={k}")
        ax.set_xlabel("Celltype")
        ax.set_ylabel("NMF_component")
        ax.tick_params(axis="x", rotation=45)
        ax.tick_params(axis="y", rotation=0)

    for idx in range(len(sample_names), len(axes)):
        axes[idx].axis("off")

    fig.suptitle(
        f"NMF component loading heatmap by sample | k={k}",
        fontsize=18,
        y=0.98
    )
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    out_path = os.path.join(
        plot_dir,
        f"MGUS_50um_nmf_H_heatmap_k{k}_by_sample.svg"
    )
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved figure: {out_path}")
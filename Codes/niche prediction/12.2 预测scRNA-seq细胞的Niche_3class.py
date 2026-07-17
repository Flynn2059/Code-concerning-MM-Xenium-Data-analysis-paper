#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12.2 预测 scRNA-seq 细胞的空间 Niche（3 分类器版 + PyTorch MPS 加速）

使用 12.1 训练的 disease-specific 3 分类器：
  SM 模型 → 预测 SM scRNA-seq
  MM 模型 → 预测 MM scRNA-seq

与 11.2 的区别：
  - 加载 ./niche_classifier_3class/ 中的模型
  - 软聚类：主要输出是 3 列连续概率，不做硬切分
  - predicted_niche 仅保留 argmax 硬标签供参考，不用于下游统计
  - 概率列：prob_Mixed, prob_Plasma_Dominant, prob_T_Dominant（加和=1，软成分）

输出（./niche_classifier_3class/predictions/）：
  {disease}_T_cell_niche_pred.csv
  {disease}_Plasma_cell_niche_pred.csv
  列：barcode, sample_id, celltype_secondary, clinical_grade,
       predicted_niche,        argmax 硬标签（仅参考）
       prob_Mixed, prob_Plasma_Dominant, prob_T_Dominant   ← 主要输出
"""

#%% ── 依赖 ──────────────────────────────────────────────────────────────────
import os
import numpy as np
import pandas as pd
import anndata
import scanpy as sc
import joblib
import scipy.sparse as sp

import torch
import torch.nn as nn

#%% ── 配置 ──────────────────────────────────────────────────────────────────
os.chdir('/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/')

SCRNA_T_LABELS = [
    "CD4_T_memory", "CD4_T_naive", "CD4_Tcm",
    "CD8_T_memory", "CD8_T_naive",
    "CTL",
    "T_prolifer",
    "MAIT_gdT",
]
SCRNA_PLASMA_LABELS = ["Plasma_cell"]

SCRNA_DATASETS = [
    ("SM", "../external/data_for_scAnalysis/SM_scAnalysis.h5ad"),
    ("MM", "../external/data_for_scAnalysis/MM_scAnalysis.h5ad"),
]

CLASSIFIER_DIR = "./niche_classifier_3class"
PRED_OUTDIR    = "./niche_classifier_3class/predictions"
os.makedirs(PRED_OUTDIR, exist_ok=True)

# ── 设备 ──────────────────────────────────────────────────────────────────────
if torch.backends.mps.is_available():
    DEVICE = torch.device("mps")
elif torch.cuda.is_available():
    DEVICE = torch.device("cuda")
else:
    DEVICE = torch.device("cpu")
print(f"Using device: {DEVICE}")

#%% ── 模型定义（与 12.1 一致）───────────────────────────────────────────────

class MLPClassifier(nn.Module):
    def __init__(self, n_features, n_classes, dropout=0.3):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_features, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(512, 256),
            nn.BatchNorm1d(256),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(256, 128),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.Linear(128, n_classes)
        )

    def forward(self, x):
        return self.net(x)


def load_classifier(disease, cell_type_name, classifier_dir, device):
    meta = joblib.load(
        os.path.join(classifier_dir, f"{disease}_{cell_type_name}_niche_meta.pkl")
    )
    le         = meta["label_encoder"]
    genes      = meta["selected_genes"]
    n_features = meta["n_features"]
    n_classes  = meta["n_classes"]
    dropout    = meta.get("dropout", 0.3)

    model = MLPClassifier(n_features, n_classes, dropout)
    state = torch.load(
        os.path.join(classifier_dir, f"{disease}_{cell_type_name}_niche_classifier.pt"),
        map_location="cpu",
        weights_only=True,
    )
    model.load_state_dict(state)
    model = model.to(device)
    model.eval()
    return model, le, genes


def predict_niche(adata_sub, model, le, model_genes, device,
                  cell_type_tag, batch_size=4096):
    """
    对已归一化的单细胞子集做 niche 预测。
    缺失基因补 0。

    主输出：3 列软概率（prob_Mixed / prob_Plasma_Dominant / prob_T_Dominant），
    每行加和 = 1，供下游软聚类分析使用。
    predicted_niche 为 argmax 硬标签，仅供参考。
    """
    genes_present = [g for g in model_genes if g in adata_sub.var_names]
    genes_missing  = [g for g in model_genes if g not in adata_sub.var_names]
    if genes_missing:
        print(f"  [{cell_type_tag}] Missing genes: {len(genes_missing)}/{len(model_genes)}")

    X_present = adata_sub[:, genes_present].X
    X_present  = X_present.toarray() if sp.issparse(X_present) else np.array(X_present)

    gene_idx = {g: i for i, g in enumerate(model_genes)}
    X_full = np.zeros((adata_sub.n_obs, len(model_genes)), dtype=np.float32)
    for j, g in enumerate(genes_present):
        X_full[:, gene_idx[g]] = X_present[:, j]

    all_probs = []
    with torch.no_grad():
        for i in range(0, len(X_full), batch_size):
            xb    = torch.tensor(X_full[i:i+batch_size], dtype=torch.float32).to(device)
            probs = torch.softmax(model(xb), dim=1).cpu().numpy()
            all_probs.append(probs)
    probs_all = np.vstack(all_probs)

    classes   = le.classes_                   # ["Mixed", "Plasma_Dominant", "T_Dominant"]
    pred_hard = le.inverse_transform(probs_all.argmax(axis=1))   # argmax，仅参考

    prob_cols = {f"prob_{c}": probs_all[:, i] for i, c in enumerate(classes)}

    return pd.DataFrame({
        "barcode":         adata_sub.obs_names,
        "predicted_niche": pred_hard,          # argmax 硬标签，仅参考
        **prob_cols                             # 软概率：主要输出
    }, index=adata_sub.obs_names)


#%% ── 主流程 ─────────────────────────────────────────────────────────────────

for disease, h5ad_path in SCRNA_DATASETS:
    print("\n" + "=" * 60)
    print(f"Disease: {disease}")
    print("=" * 60)

    print(f"  Loading {h5ad_path} ...")
    adata = sc.read_h5ad(h5ad_path)

    x_sample  = adata.X[0]
    x_sample  = x_sample.toarray().flatten() if sp.issparse(x_sample) else np.array(x_sample).flatten()
    x_nonzero = x_sample[x_sample > 0]
    already_log = len(x_nonzero) > 0 and float(x_nonzero.max()) < 20

    if already_log:
        print("  X appears already log-normalized, using as-is.")
    else:
        print("  Normalizing (normalize_total + log1p) ...")
        sc.pp.normalize_total(adata, target_sum=10000)
        sc.pp.log1p(adata)

    id_col = "sample_id" if "sample_id" in adata.obs.columns else "donor_id"

    for cell_tag, cell_labels in [
        ("T_cell",      SCRNA_T_LABELS),
        ("Plasma_cell", SCRNA_PLASMA_LABELS),
    ]:
        mask = adata.obs["celltype_secondary"].isin(cell_labels)
        print(f"\n  {cell_tag}: {mask.sum()} cells")
        if mask.sum() == 0:
            continue

        try:
            model, le, model_genes = load_classifier(
                disease, cell_tag, CLASSIFIER_DIR, DEVICE
            )
        except FileNotFoundError as e:
            print(f"  → Classifier not found: {e}, skip.")
            continue

        print(f"  Classes: {le.classes_}")

        adata_sub = adata[mask].copy()
        pred_df   = predict_niche(adata_sub, model, le, model_genes, DEVICE, cell_tag)

        pred_df["sample_id"]         = adata_sub.obs[id_col].values
        pred_df["celltype_secondary"] = adata_sub.obs["celltype_secondary"].values
        pred_df["disease"]            = disease
        if "clinical_grade" in adata_sub.obs.columns:
            pred_df["clinical_grade"] = adata_sub.obs["clinical_grade"].values

        print(f"  predicted_niche dist:\n{pred_df['predicted_niche'].value_counts().to_dict()}")
        for col in ["prob_Mixed", "prob_Plasma_Dominant", "prob_T_Dominant"]:
            print(f"  {col}: mean={pred_df[col].mean():.3f}  median={pred_df[col].median():.3f}")

        out_path = os.path.join(PRED_OUTDIR, f"{disease}_{cell_tag}_niche_pred.csv")
        pred_df.to_csv(out_path, index=False)
        print(f"  → Saved: {out_path}")

    del adata

print("\n" + "=" * 60)
print("Done. Predictions in:", PRED_OUTDIR)
print("=" * 60)

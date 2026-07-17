#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
12.1 构建 Niche 3 分类器（T_Dominant / Plasma_Dominant / Mixed + PyTorch MPS 加速）

与 11.1 的区别：
  - 训练目标改为 3 类（新增 Mixed），丢弃 Engaging Zone
  - LabelEncoder 固定顺序：Mixed=0, Plasma_Dominant=1, T_Dominant=2
  - 类权重改为 3 类平衡
  - AUC 改为 multi-class OvR macro

流程：
  1. 检查 Xenium ∩ scRNA-seq 基因交集（复用 11.1 输出的 common_genes.csv）
  2. 对每个 disease（MM / SM）独立处理：
     a. 从 Xenium h5ad 提取 T 细胞 / Plasma 细胞 + niche 标签
        （保留 T_Dominant、Plasma_Dominant、Mixed，丢弃 Engaging Zone）
     b. ANOVA F-test 筛选 top K 基因
     c. MLP（3 隐藏层）5-fold 交叉验证 + 全量训练
     d. 保存模型、选中基因列表、特征重要性

输出（./niche_classifier_3class/）：
  {disease}_{cell_type}_selected_genes.csv
  {disease}_{cell_type}_niche_classifier.pt
  {disease}_{cell_type}_niche_meta.pkl
  {disease}_{cell_type}_cv_results.csv
  {disease}_{cell_type}_feature_importance.csv
"""

#%% ── 依赖 ──────────────────────────────────────────────────────────────────
import os
import numpy as np
import pandas as pd
import anndata
import joblib
import scipy.sparse as sp

import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

from sklearn.model_selection import StratifiedGroupKFold
from sklearn.preprocessing import LabelEncoder
from sklearn.metrics import classification_report, f1_score, roc_auc_score
from sklearn.feature_selection import SelectKBest, f_classif

#%% ── 配置 ──────────────────────────────────────────────────────────────────
os.chdir('/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/')

T_LABELS      = ["CD4+ T Naive", "CD8+ T Effector", "CD8+ T Exhausted", "T Proliferating"]
PLASMA_LABELS = ["Plasma", "Malignant_Plasma"]

NICHE_CLASSES = ["T_Dominant", "Plasma_Dominant", "Mixed"]   # 丢弃 Engaging Zone

XENIUM_DATASETS = [
    ("MM", "10_mm_major_annotation.h5ad",
           "Tcell_Plasma_GMM_clustering_annotation/MM_10_samples_cell_niche_matrix.csv"),
    ("SM", "5_sm_major_annotation.h5ad",
           "Tcell_Plasma_GMM_clustering_annotation/SM_5_samples_cell_niche_matrix.csv"),
]

SCRNA_REF = "../external/data_for_scAnalysis/MM_scAnalysis.h5ad"

OUTDIR = "./niche_classifier_3class"
os.makedirs(OUTDIR, exist_ok=True)

EXCLUDE_GENES = ["SOX2-OT"]

# ── 超参数（与 11.1 一致）────────────────────────────────────────────────────
TOP_K    = 500
BATCH    = 512
EPOCHS   = 300
LR       = 1e-3
WD       = 1e-4
DROPOUT  = 0.3
N_SPLITS = 5

# ── 设备选择 ──────────────────────────────────────────────────────────────────
if torch.backends.mps.is_available():
    DEVICE = torch.device("mps")
elif torch.cuda.is_available():
    DEVICE = torch.device("cuda")
else:
    DEVICE = torch.device("cpu")
print(f"Using device: {DEVICE}")

#%% ── MLP 模型定义（与 11.1 相同）──────────────────────────────────────────

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


def fit(model, X, y, device, n_classes, epochs=EPOCHS, batch=BATCH, lr=LR, wd=WD,
        verbose=False):
    w = torch.tensor(
        len(y) / (n_classes * np.bincount(y, minlength=n_classes).clip(1).astype(float)),
        dtype=torch.float32
    ).to(device)
    loss_fn   = nn.CrossEntropyLoss(weight=w)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)
    loader = DataLoader(
        TensorDataset(torch.tensor(X, dtype=torch.float32),
                      torch.tensor(y, dtype=torch.long)),
        batch_size=batch, shuffle=True
    )
    model.train()
    for epoch in range(epochs):
        total_loss = 0.0
        for xb, yb in loader:
            xb, yb = xb.to(device), yb.to(device)
            optimizer.zero_grad()
            loss = loss_fn(model(xb), yb)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        scheduler.step()
        if verbose and (epoch + 1) % 100 == 0:
            print(f"    epoch {epoch+1}/{epochs}  loss={total_loss/len(loader):.4f}")
    return model


def predict(model, X, device, batch_size=4096):
    model.eval()
    parts = []
    with torch.no_grad():
        for i in range(0, len(X), batch_size):
            xb    = torch.tensor(X[i:i+batch_size], dtype=torch.float32).to(device)
            probs = torch.softmax(model(xb), dim=1).cpu().numpy()
            parts.append(probs)
    return np.vstack(parts)


#%% ── Step 1：基因交集（复用 11.1 的 common_genes.csv）────────────────────
print("=" * 60)
print("Step 1: Load common genes from 11.1")
print("=" * 60)

common_genes_path = "./niche_classifier/common_genes.csv"
common_genes = pd.read_csv(common_genes_path, header=None)[0].tolist()
print(f"Common genes loaded: {len(common_genes)}")

#%% ── Step 2：按 disease 分别训练 ────────────────────────────────────────────

def run_disease(disease, h5ad_file, niche_csv):
    print("\n" + "=" * 60)
    print(f"Disease: {disease}  |  3-class: T_Dominant / Plasma_Dominant / Mixed")
    print("=" * 60)

    niche_df       = pd.read_csv(niche_csv, index_col="cell_full_id")
    three_cls_mask = niche_df["niche"].isin(NICHE_CLASSES)

    t_ids      = niche_df.index[niche_df["annotation"].isin(T_LABELS)      & three_cls_mask].tolist()
    plasma_ids = niche_df.index[niche_df["annotation"].isin(PLASMA_LABELS) & three_cls_mask].tolist()
    print(f"  T cells    : {len(t_ids)}")
    print(f"  Plasma cells: {len(plasma_ids)}")

    print(f"  Loading {h5ad_file} (backed) ...")
    adata_full    = anndata.read_h5ad(h5ad_file, backed='r')
    genes_in_h5ad = [g for g in common_genes
                     if g in adata_full.var_names and g not in EXCLUDE_GENES]

    for cell_type_name, cell_ids in [("T_cell", t_ids), ("Plasma_cell", plasma_ids)]:
        if not cell_ids:
            print(f"  → {cell_type_name}: no cells, skip")
            continue

        valid_ids = [cid for cid in cell_ids if cid in adata_full.obs_names]
        if not valid_ids:
            print(f"  → {cell_type_name}: no valid IDs found in h5ad, skip")
            continue

        print(f"\n  ── {cell_type_name} ({len(valid_ids)} cells) ──────────────")

        # ── 加载表达矩阵 ──────────────────────────────────────────────────
        sub   = adata_full[valid_ids, genes_in_h5ad].to_memory()
        X_raw = sub.layers["lognorm"] if "lognorm" in sub.layers else sub.X
        X_all = X_raw.toarray() if sp.issparse(X_raw) else np.array(X_raw, dtype=np.float32)

        y_raw  = niche_df.loc[valid_ids, "niche"].values
        groups = niche_df.loc[valid_ids, "sample"].values

        # 固定 LabelEncoder 顺序：Mixed=0, Plasma_Dominant=1, T_Dominant=2
        le = LabelEncoder()
        le.fit(["Mixed", "Plasma_Dominant", "T_Dominant"])
        y  = le.transform(y_raw)
        for cls, idx in zip(le.classes_, range(len(le.classes_))):
            print(f"  {cls}: {(y == idx).sum()}")

        n_classes = len(le.classes_)

        # ── ANOVA F-test 特征筛选 ──────────────────────────────────────────
        print(f"  Feature selection: ANOVA F-test, top {TOP_K} genes ...")
        sel = SelectKBest(f_classif, k=min(TOP_K, len(genes_in_h5ad)))
        sel.fit(X_all, y)
        idx       = sel.get_support(indices=True)
        sel_genes = [genes_in_h5ad[i] for i in idx]
        X         = X_all[:, idx].astype(np.float32)

        top_genes = pd.DataFrame({
            "gene":    sel_genes,
            "f_score": sel.scores_[idx]
        }).sort_values("f_score", ascending=False)
        top_genes.to_csv(
            os.path.join(OUTDIR, f"{disease}_{cell_type_name}_selected_genes.csv"),
            index=False
        )
        print(f"  Selected {len(sel_genes)} genes → saved.")
        print(f"  Top 10 discriminating genes:\n{top_genes.head(10).to_string(index=False)}")

        n_features = X.shape[1]

        # ── 5-fold CV ─────────────────────────────────────────────────────
        cv = StratifiedGroupKFold(n_splits=N_SPLITS, shuffle=True, random_state=42)
        fold_f1s, fold_aucs = [], []

        for fold, (tr, te) in enumerate(cv.split(X, y, groups)):
            m = MLPClassifier(n_features, n_classes, DROPOUT).to(DEVICE)
            fit(m, X[tr], y[tr], DEVICE, n_classes, verbose=False)
            probs  = predict(m, X[te], DEVICE)
            y_pred = probs.argmax(axis=1)
            f1  = f1_score(y[te], y_pred, average="macro", zero_division=0)
            auc = roc_auc_score(y[te], probs, multi_class='ovr', average='macro')
            fold_f1s.append(f1)
            fold_aucs.append(auc)
            print(f"  Fold {fold}: macro-F1 = {f1:.3f}  AUC(ovr) = {auc:.3f}")
            print(classification_report(
                y[te], y_pred, target_names=le.classes_, zero_division=0
            ))

        print(f"\n  ── CV Summary ──")
        print(f"  macro-F1 : {np.mean(fold_f1s):.3f} ± {np.std(fold_f1s):.3f}")
        print(f"  AUC-ROC  : {np.mean(fold_aucs):.3f} ± {np.std(fold_aucs):.3f}")

        pd.DataFrame({
            "fold":     range(N_SPLITS),
            "macro_f1": fold_f1s,
            "auc_roc":  fold_aucs,
        }).to_csv(
            os.path.join(OUTDIR, f"{disease}_{cell_type_name}_cv_results.csv"),
            index=False
        )

        # ── 全量训练最终模型 ───────────────────────────────────────────────
        print("  Training final model on all data ...")
        final_model = MLPClassifier(n_features, n_classes, DROPOUT).to(DEVICE)
        fit(final_model, X, y, DEVICE, n_classes, verbose=True)

        # ── 保存模型 + meta ────────────────────────────────────────────────
        torch.save(
            final_model.cpu().state_dict(),
            os.path.join(OUTDIR, f"{disease}_{cell_type_name}_niche_classifier.pt")
        )
        joblib.dump({
            "label_encoder":  le,
            "selected_genes": sel_genes,
            "n_features":     n_features,
            "n_classes":      n_classes,
            "dropout":        DROPOUT,
        }, os.path.join(OUTDIR, f"{disease}_{cell_type_name}_niche_meta.pkl"))
        print(f"  → Saved: {disease}_{cell_type_name}_niche_classifier.pt + meta.pkl")

        # ── 特征重要性：输入层权重 L2 范数 ────────────────────────────────
        w = final_model.net[0].weight.detach().numpy()   # shape: 512 × n_features
        gene_importance = np.linalg.norm(w, axis=0)
        imp_df = pd.DataFrame({
            "gene":       sel_genes,
            "importance": gene_importance
        }).sort_values("importance", ascending=False)
        imp_df.to_csv(
            os.path.join(OUTDIR, f"{disease}_{cell_type_name}_feature_importance.csv"),
            index=False
        )
        print(f"  → Feature importance saved.")

    adata_full.file.close()


for disease, h5ad_file, niche_csv in XENIUM_DATASETS:
    run_disease(disease, h5ad_file, niche_csv)

print("\n" + "=" * 60)
print("Done. Output in:", OUTDIR)
print("=" * 60)

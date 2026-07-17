#%% imports
import os
from pathlib import Path
import numpy as np
import scanpy as sc
import spatialdata as sd

#%% 工作路径
os.chdir('/Volumes/FlynnBio/XinyeLi/GSE299207_RAW/')
print(os.getcwd())
print("\n")

#%% 创建输出目录
outdir = Path("banksy")
outdir.mkdir(exist_ok=True)

#%% 收集所有 *_post_qc.zarr
zarr_files = sorted([
    f for f in os.listdir()
    if f.endswith("_post_qc.zarr")
])

print("Found zarr files:")
for f in zarr_files:
    print(" -", f)

print("\nProcessing samples...\n")

#%% 遍历读取
for zarr in zarr_files:

    # sample 名（复用之前逻辑）
    sample_name = zarr.replace("_post_qc.zarr", "").replace(".zarr", "")
    prefix = f"{sample_name}_"

    print(f"Reading {zarr} ...")

    # 读取 SpatialData
    sdata = sd.read_zarr(zarr)

    if "table" not in sdata:
        print(f"  -> Warning: no 'table' found in {zarr}")
        continue

    # 提取 AnnData
    adata = sdata["table"].copy()

    # ========================
    # 样本前缀处理
    # ========================

    adata.obs["sample"] = sample_name

    # obs_names 加 prefix
    adata.obs_names = prefix + adata.obs_names.astype(str)

    # 常见 ID 列
    for col in ["cell_id", "instance_id", "spot_id"]:
        if col in adata.obs.columns:
            adata.obs[col] = prefix + adata.obs[col].astype(str)

    # ========================
    # 预处理步骤
    # ========================

    adata.layers["counts"] = adata.X.copy()
    
    sc.pp.normalize_total(
    adata,
    target_sum=np.median(adata.obs["total_counts"])
    )

    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=3000,
        subset=False
    )

    print(adata)

    # ========================
    # 保存 adata
    # ========================

    outfile = outdir / f"{sample_name}.adata"
    adata.write(outfile)

    print(f"  -> saved to {outfile}\n")

print("All samples processed.")
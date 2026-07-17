#!/usr/bin/env python
# coding: utf-8

"""
Batch Xenium QC pipeline (with resume & error tracking)

New features
------------
1. 自动跳过已经完成 QC（存在 *_post_qc.zarr）的样本
2. 记录失败样本并在最后统一报告
3. 支持断点续跑：只处理未完成样本
4. 生成 metrics_summary.csv（仅包含成功样本）
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# SpatialData stack
from spatialdata import SpatialData
from spatialdata.models import TableModel, ShapesModel, PointsModel, Image2DModel

import tifffile as tiff
import geopandas as gpd
from shapely.geometry import Polygon, MultiPolygon


# ==============================
# Unified QC thresholds
# ==============================
MIN_COUNTS = 20
MIN_GENES = 10
plt.rcParams["figure.figsize"] = (6, 6)


# ==============================
# Xenium reader
# ==============================

def read_xenium_geo_sample(sample_dir: str | Path) -> SpatialData:
    sample_dir = Path(sample_dir)

    matrix_h5 = sample_dir / "cell_feature_matrix.h5"
    cells_csv = sample_dir / "cells.csv.gz"
    boundaries_csv = sample_dir / "cell_boundaries.csv.gz"
    transcripts_parquet = sample_dir / "transcripts.parquet"
    morphology_tif = sample_dir / "morphology.ome.tif"

    if not matrix_h5.exists():
        raise FileNotFoundError(f"Missing {matrix_h5}")

    adata = sc.read_10x_h5(matrix_h5)
    adata.var_names_make_unique()

    adata.layers["counts"] = adata.X.copy()
    adata.X = adata.layers["counts"]

    # ---------- cells ----------
    cells = pd.read_csv(cells_csv, compression="gzip").set_index("cell_id")
    cells = cells.reindex(adata.obs_names)

    adata.obs = cells.drop(columns=["x_centroid", "y_centroid"], errors="ignore")
    adata.obsm["spatial"] = cells[["x_centroid", "y_centroid"]].to_numpy()

    # ---------- images ----------
    images = {}
    if morphology_tif.exists():
        img = tiff.imread(morphology_tif)
        images["morphology"] = Image2DModel.parse(img, dims=("c", "y", "x"))

    # ---------- transcripts ----------
    points = {}
    if transcripts_parquet.exists():
        tx = pd.read_parquet(transcripts_parquet)

        tx = tx.rename(columns={"x_location": "x", "y_location": "y", "z_location": "z"})
        tx["cell_id"] = tx["cell_id"].astype("category")
        tx["feature_name"] = tx["feature_name"].astype("category")

        coord_map = {"x": "x", "y": "y"}
        if "z" in tx.columns:
            coord_map["z"] = "z"

        points["transcripts"] = PointsModel.parse(
            tx.reset_index(drop=True),
            coordinates=coord_map,
            feature_key="feature_name",
            instance_key="cell_id",
            npartitions=8,
        )

    # ---------- boundaries ----------
    shapes = {}
    if boundaries_csv.exists():
        b = pd.read_csv(boundaries_csv, compression="gzip")

        def make_polygon(df):
            coords = df[["vertex_x", "vertex_y"]].to_numpy()
            if len(coords) < 3:
                return None
            if not np.allclose(coords[0], coords[-1]):
                coords = np.vstack([coords, coords[0]])
            poly = Polygon(coords).buffer(0)
            return None if poly.is_empty else poly

        polys = (
            b.groupby(["cell_id", "label_id"], observed=True)
            .apply(make_polygon)
            .dropna()
        )

        cell_to_geom = {}
        for (cid, _), geom in polys.items():
            cell_to_geom.setdefault(cid, []).append(geom)

        geoms = [g[0] if len(g) == 1 else MultiPolygon(g) for g in cell_to_geom.values()]

        gdf = gpd.GeoDataFrame({"cell_id": list(cell_to_geom.keys())}, geometry=geoms).set_index("cell_id")
        gdf = gdf.loc[gdf.index.intersection(adata.obs_names)]

        shapes["cell_boundaries"] = ShapesModel.parse(gdf)

    table = TableModel.parse(adata)

    sdata = SpatialData(images=images, points=points, shapes=shapes, tables={"table": table})

    return sdata


# ==============================
# QC plotting
# ==============================

def save_qc_plots(adata, sample_name, qc_dir):
    qc_dir.mkdir(exist_ok=True)

    plt.figure()
    plt.hist(adata.obs["total_counts"], bins=100)
    plt.xlabel("Transcripts per cell")
    plt.ylabel("Number of cells")
    plt.title(sample_name)
    plt.savefig(qc_dir / f"{sample_name}_total_counts.png", dpi=150)
    plt.close()

    plt.figure()
    plt.hist(adata.obs["n_genes_by_counts"], bins=100)
    plt.xlabel("Genes detected per cell")
    plt.ylabel("Number of cells")
    plt.title(sample_name)
    plt.savefig(qc_dir / f"{sample_name}_n_genes.png", dpi=150)
    plt.close()


# ==============================
# Main batch loop (with resume)
# ==============================

def main(workdir: str | Path):
    workdir = Path(workdir)
    qc_dir = workdir / "QC"
    qc_dir.mkdir(exist_ok=True)

    summary_rows = []
    failed_samples = []

    for sample in sorted(workdir.iterdir()):
        import gc
        gc.collect()

        if not sample.is_dir():
            continue
        if not sample.name.startswith("GSM"):
            continue

        zarr_path = workdir / f"{sample.name}_post_qc.zarr"

        # ---------- skip finished samples ----------
        if zarr_path.exists():
            print(f"[SKIP] {sample.name} already processed.")
            continue

        print(f"Processing {sample.name} ...")

        try:
            sdata = read_xenium_geo_sample(sample)
            adata = sdata["table"]

            sc.pp.calculate_qc_metrics(adata, inplace=True)

            raw_cells = adata.n_obs

            sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
            sc.pp.filter_cells(adata, min_genes=MIN_GENES)

            filtered_cells = adata.n_obs

            min_gene_val = adata.obs["n_genes_by_counts"].min()
            min_gene_cell_count = (adata.obs["n_genes_by_counts"] == min_gene_val).sum()

            save_qc_plots(adata, sample.name, qc_dir)

            sdata.tables["table"] = adata
            sdata.write(zarr_path, overwrite=True)

            summary_rows.append(
                {
                    "sample": sample.name,
                    "raw_cells": raw_cells,
                    "post_qc_cells": filtered_cells,
                    "min_n_genes": min_gene_val,
                    "cells_at_min_n_genes": min_gene_cell_count,
                }
            )

        except Exception as e:
            print(f"[ERROR] {sample.name}: {e}")
            failed_samples.append({"sample": sample.name, "error": str(e)})
            continue

    import gc
    gc.collect()

    # ---------- save success summary ----------
    if summary_rows:
        summary_df = pd.DataFrame(summary_rows)
        summary_df.to_csv(qc_dir / "metrics_summary.csv", index=False)

    # ---------- save failure log ----------
    if failed_samples:
        fail_df = pd.DataFrame(failed_samples)
        fail_df.to_csv(qc_dir / "failed_samples.csv", index=False)

    print("\nAll samples processed.")
    print(f"Success: {len(summary_rows)} | Failed: {len(failed_samples)}")


if __name__ == "__main__":
    main(".")

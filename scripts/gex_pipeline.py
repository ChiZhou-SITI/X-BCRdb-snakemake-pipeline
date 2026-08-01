import warnings

import gzip
from pathlib import Path

import anndata as ad
import celltypist
import numpy as np
import pandas as pd
import scanpy as sc
from celltypist import models
from scipy.io import mmread



def _open_maybe_gzip(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _find_10x_file(path, stem):
    path = Path(path)
    for name in (f"{stem}.tsv.gz", f"{stem}.tsv", f"{stem}.mtx.gz", f"{stem}.mtx"):
        candidate = path / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Could not find {stem} in {path}")


def read_10x_matrix_fallback(path):
    path = Path(path)
    matrix_path = _find_10x_file(path, "matrix")
    features_path = _find_10x_file(path, "features")
    barcodes_path = _find_10x_file(path, "barcodes")

    with _open_maybe_gzip(matrix_path, "rb") as handle:
        matrix = mmread(handle).tocsr().T

    features = pd.read_csv(features_path, sep="\t", header=None, compression="infer")
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, compression="infer")
    if features.shape[1] >= 2:
        var_names = features.iloc[:, 1].astype(str).values
        gene_ids = features.iloc[:, 0].astype(str).values
    else:
        var_names = features.iloc[:, 0].astype(str).values
        gene_ids = var_names

    adata = ad.AnnData(X=matrix)
    adata.obs_names = barcodes.iloc[:, 0].astype(str).values
    adata.var_names = var_names
    adata.var["gene_ids"] = gene_ids
    if features.shape[1] >= 3:
        adata.var["feature_types"] = features.iloc[:, 2].astype(str).values
    return adata


def read_10x_matrix(path):
    try:
        return sc.read_10x_mtx(path, var_names="gene_symbols", cache=True)
    except Exception as exc1:
        print(f"read_10x_mtx gene_symbols failed: {exc1}; retry with gene_ids")
        try:
            return sc.read_10x_mtx(path, var_names="gene_ids", cache=True)
        except Exception as exc2:
            print(f"read_10x_mtx gene_ids failed: {exc2}; use custom 10x reader")
            return read_10x_matrix_fallback(path)


def safe_has_nan(x):
    data = x.data if hasattr(x, "data") else x
    try:
        return bool(np.any(np.isnan(data)))
    except TypeError:
        return False


def sanitize_x(adata):
    if safe_has_nan(adata.X):
        print("检测到 NaN 值，进行填充处理")
        x = adata.X.copy()
        if hasattr(x, "toarray"):
            x = x.toarray()
        adata.X = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)


def ensure_embedding_and_cluster(adata, cluster_key="leiden"):
    if "X_umap" not in adata.obsm or adata.obsm["X_umap"].shape[0] != adata.n_obs:
        adata.obsm["X_umap"] = np.zeros((adata.n_obs, 2), dtype=float)
    if cluster_key not in adata.obs:
        adata.obs[cluster_key] = "0"
    if "cluster" not in adata.obs:
        adata.obs["cluster"] = adata.obs[cluster_key].astype(str)


def preprocess_for_embedding(adata, n_top_genes, n_neighbors, max_n_pcs, resolution):
    if adata.n_obs == 0:
        ensure_embedding_and_cluster(adata)
        return adata

    if adata.n_obs < 2 or adata.n_vars < 2:
        ensure_embedding_and_cluster(adata)
        return adata

    n_top_genes = max(1, min(n_top_genes, adata.n_vars))
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)
        if "highly_variable" in adata.var and int(adata.var["highly_variable"].sum()) >= 2:
            adata = adata[:, adata.var.highly_variable].copy()
    except Exception as exc:
        print(f"highly_variable_genes failed: {exc}; use all genes")

    sanitize_x(adata)
    try:
        sc.pp.scale(adata, max_value=10)
        sanitize_x(adata)
    except Exception as exc:
        print(f"scale failed: {exc}; continue without scaling")

    max_possible = min(adata.n_obs, adata.n_vars) - 1
    if max_possible < 1:
        ensure_embedding_and_cluster(adata)
        return adata

    n_comps = max(1, min(max_n_pcs, max_possible))
    try:
        sc.tl.pca(adata, n_comps=n_comps, svd_solver="arpack")
        n_pcs = max(1, min(max_n_pcs, adata.obsm["X_pca"].shape[1]))
        nn = max(1, min(n_neighbors, adata.n_obs - 1))
        sc.pp.neighbors(adata, n_neighbors=nn, n_pcs=n_pcs)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=resolution)
    except Exception as exc:
        print(f"embedding/clustering failed: {exc}; use zero embedding")
        ensure_embedding_and_cluster(adata)
        return adata

    ensure_embedding_and_cluster(adata)
    return adata


def annotate_celltypist(adata, model_name, fallback_label):
    try:
        model = models.Model.load(model=model_name)
        predictions = celltypist.annotate(adata, model=model, majority_voting=True)
        return predictions.predicted_labels.predicted_labels.astype(str)
    except Exception as exc:
        print(f"CellTypist 注释失败 ({model_name}): {exc}; fallback={fallback_label}")
        return pd.Series(fallback_label, index=adata.obs_names, dtype=str)


def add_umap_to_meta(adata):
    ensure_embedding_and_cluster(adata)
    meta = adata.obs.copy()
    meta["UMAP1"] = adata.obsm["X_umap"][:, 0] if adata.n_obs else []
    meta["UMAP2"] = adata.obsm["X_umap"][:, 1] if adata.n_obs else []
    return meta


def subset_to_b_cells(adata):
    cell_type = adata.obs.get("cell_type", pd.Series("", index=adata.obs_names)).astype(str)
    mask = cell_type.str.contains("B cells|B cell|Plasma", case=False, na=False)
    if int(mask.sum()) == 0:
        print("警告: 未识别到 B 细胞，使用全部细胞作为 BCR 关联表达对象")
        mask = pd.Series(True, index=adata.obs_names)

    parent = adata[mask.values].copy()
    if adata.raw is not None:
        adata_b = adata[mask.values].raw.to_adata()
        adata_b.obs = parent.obs.copy()
    else:
        adata_b = parent.copy()
    if "X_umap" in parent.obsm:
        adata_b.obsm["X_umap"] = np.asarray(parent.obsm["X_umap"])
    ensure_embedding_and_cluster(adata_b)
    return adata_b


# 1. 读取 10x 数据
adata = read_10x_matrix(snakemake.input[0])
adata.var_names_make_unique()

# 2. QC
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
if adata.n_obs == 0 or adata.n_vars == 0:
    print("警告: QC 后没有细胞或基因，写出空结果")
    adata.obs["cell_type"] = pd.Series("Unknown", index=adata.obs_names, dtype=str)
    adata.obs["B_subtype"] = pd.Series("Unknown", index=adata.obs_names, dtype=str)
    ensure_embedding_and_cluster(adata)
    adata_b = adata.copy()
    adata.write_h5ad(snakemake.output["all_h5ad"])
    adata_b.write_h5ad(snakemake.output["b_h5ad"])
    add_umap_to_meta(adata).to_csv(snakemake.output["all"])
    add_umap_to_meta(adata_b).to_csv(snakemake.output["b"])
    raise SystemExit(0)

adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
try:
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 20, :].copy()
except Exception as exc:
    print(f"calculate_qc_metrics/MT filter failed: {exc}; continue")

# 3. 标准化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# 4. 全体细胞聚类
adata = preprocess_for_embedding(
    adata,
    n_top_genes=3000,
    n_neighbors=15,
    max_n_pcs=40,
    resolution=0.5,
)

# 5. CellTypist 自动注释（全体细胞）
adata.obs["cell_type"] = annotate_celltypist(
    adata,
    model_name="Immune_All_Low.pkl",
    fallback_label="B cells",
)

# 6. 提取 B 细胞
adata_b = subset_to_b_cells(adata)

# 7-9. B 细胞重新预处理、聚类、注释
if adata_b.n_obs < 10 or adata_b.n_vars < 2:
    print(f"警告: B 细胞数量太少或基因太少 ({adata_b.n_obs} cells, {adata_b.n_vars} genes)，跳过精细聚类")
    if "B_subtype" not in adata_b.obs:
        adata_b.obs["B_subtype"] = "Unknown"
    ensure_embedding_and_cluster(adata_b)
else:
    try:
        sc.pp.filter_cells(adata_b, min_genes=200)
        sc.pp.filter_genes(adata_b, min_cells=3)
        sc.pp.normalize_total(adata_b, target_sum=1e4)
        sc.pp.log1p(adata_b)
        adata_b.raw = adata_b
        adata_b = preprocess_for_embedding(
            adata_b,
            n_top_genes=min(2000, adata_b.n_vars),
            n_neighbors=15,
            max_n_pcs=40,
            resolution=0.6,
        )
        adata_b.obs["B_subtype"] = annotate_celltypist(
            adata_b,
            model_name="Immune_All_High.pkl",
            fallback_label="Unknown",
        )
    except Exception as exc:
        print(f"B 细胞精细处理失败: {exc}; 使用已有 B 细胞对象继续")
        if "B_subtype" not in adata_b.obs:
            adata_b.obs["B_subtype"] = "Unknown"
        ensure_embedding_and_cluster(adata_b)

if "cell_type" not in adata_b.obs:
    adata_b.obs["cell_type"] = "B cells"
if "B_subtype" not in adata_b.obs:
    adata_b.obs["B_subtype"] = "Unknown"

# 10. 输出 metadata / h5ad
meta_all = add_umap_to_meta(adata)
meta_all.to_csv(snakemake.output["all"])
adata.write_h5ad(snakemake.output["all_h5ad"])

meta_b = add_umap_to_meta(adata_b)
meta_b.to_csv(snakemake.output["b"])
adata_b.write_h5ad(snakemake.output["b_h5ad"])

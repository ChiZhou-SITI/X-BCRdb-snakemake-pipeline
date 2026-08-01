#!/usr/bin/env python3
"""
Precompute clone-level gene-expression summaries for fast BCR-transcriptome
queries.

The output is intentionally a compact summary table rather than all per-cell
expression values.  It supports fast clone-clone transcriptome similarity and
interactive differential-expression screening between expanded clonotypes.
"""

from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


DEFAULT_ALWAYS_INCLUDE_GENES = [
    "MS4A1",
    "CD79A",
    "CD79B",
    "CD74",
    "HLA-DRA",
    "HLA-DRB1",
    "CD27",
    "CD38",
    "IGHM",
    "IGHD",
    "IGHG1",
    "IGHA1",
    "IGKC",
    "IGLC1",
    "MZB1",
    "JCHAIN",
    "XBP1",
    "PRDM1",
    "AICDA",
    "BCL6",
    "MKI67",
    "TOP2A",
    "FCRL5",
    "ITGAX",
]


def normalize_barcode(value: str) -> str:
    """Normalize Cell Ranger VDJ contig IDs to GEX barcodes."""
    value = str(value or "")
    value = re.sub(r"_contig_.*$", "", value, flags=re.IGNORECASE)
    return value


def get_expression_source(adata, prefer_raw: bool = True):
    if prefer_raw and adata.raw is not None:
        return adata.raw, "raw_log1p"
    return adata, "X"


def matrix_column_stats(matrix):
    """Return mean, variance and detected counts for a cell x gene matrix."""
    n = matrix.shape[0]
    if n == 0:
        raise ValueError("Cannot summarize an empty expression matrix")
    if sparse.issparse(matrix):
        mean = np.asarray(matrix.mean(axis=0)).ravel()
        sq_mean = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel()
        detected = np.asarray((matrix > 0).sum(axis=0)).ravel()
    else:
        arr = np.asarray(matrix)
        mean = arr.mean(axis=0)
        sq_mean = np.square(arr).mean(axis=0)
        detected = (arr > 0).sum(axis=0)
    variance = np.maximum(sq_mean - np.square(mean), 0.0)
    if n > 1:
        variance = variance * n / (n - 1)
    return mean.astype(float), variance.astype(float), detected.astype(int)


def load_clone_cells(clone_path: Path, min_cells: int, max_clones: int) -> pd.DataFrame:
    clone_df = pd.read_csv(clone_path, sep="\t", dtype=str)
    required = {"clone_id", "cell_barcode"}
    missing = required - set(clone_df.columns)
    if missing:
        raise ValueError(f"{clone_path} missing required columns: {', '.join(sorted(missing))}")

    clone_df = clone_df[["clone_id", "cell_barcode"]].copy()
    clone_df["clone_id"] = clone_df["clone_id"].fillna("").astype(str)
    clone_df["cell_barcode"] = clone_df["cell_barcode"].fillna("").astype(str).map(normalize_barcode)
    clone_df = clone_df[(clone_df["clone_id"] != "") & (clone_df["cell_barcode"] != "")]
    clone_df = clone_df.drop_duplicates()

    counts = (
        clone_df.groupby("clone_id", as_index=False)["cell_barcode"]
        .nunique()
        .rename(columns={"cell_barcode": "n_cells"})
    )
    counts = counts[counts["n_cells"] >= min_cells].copy()
    if counts.empty:
        return clone_df.iloc[0:0].assign(n_cells=pd.Series(dtype=int))

    def clone_number(clone_id: str) -> int:
        match = re.search(r"(\d+)$", str(clone_id))
        return int(match.group(1)) if match else 2_147_483_647

    counts["clone_rank_number"] = counts["clone_id"].map(clone_number)
    counts = counts.sort_values(["n_cells", "clone_rank_number", "clone_id"], ascending=[False, True, True])
    counts = counts.head(max_clones)
    return clone_df.merge(counts[["clone_id", "n_cells"]], on="clone_id", how="inner")


def select_genes(expression_source, max_genes: int, min_detected_cells: int, always_include: list[str]):
    var_names = pd.Index(expression_source.var_names.astype(str))
    X = expression_source.X
    mean, variance, detected = matrix_column_stats(X)
    eligible = detected >= min_detected_cells
    if not np.any(eligible):
        eligible = detected > 0
    scores = np.where(eligible, variance, -1.0)
    ranked = np.argsort(scores)[::-1]
    selected = []
    selected_set = set()

    for gene in always_include:
        if gene in var_names and gene not in selected_set:
            idx = int(var_names.get_loc(gene))
            if detected[idx] > 0:
                selected.append(idx)
                selected_set.add(gene)

    for idx in ranked:
        if len(selected) >= max_genes:
            break
        if scores[idx] < 0:
            break
        gene = str(var_names[idx])
        if gene in selected_set:
            continue
        selected.append(int(idx))
        selected_set.add(gene)

    selected = selected[:max_genes]
    return selected, [str(var_names[idx]) for idx in selected]


def precompute_clone_gene_expression(
    h5ad_path: Path,
    clone_path: Path,
    output_path: Path,
    study: str,
    sample_id: str,
    min_cells: int = 6,
    max_clones: int = 40,
    max_genes: int = 1500,
    min_detected_cells: int = 10,
    prefer_raw: bool = True,
):
    logger.info("Reading B-cell h5ad: %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    expression_source, expression_layer = get_expression_source(adata, prefer_raw=prefer_raw)

    clone_cells = load_clone_cells(clone_path, min_cells=min_cells, max_clones=max_clones)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    header = [
        "study",
        "sample_id",
        "clone_id",
        "gene",
        "expression_layer",
        "n_cells",
        "mean_expression",
        "variance_expression",
        "detected_cells",
        "detection_fraction",
    ]

    if clone_cells.empty:
        logger.warning("No clonotypes with >= %s cells in %s; writing header only", min_cells, sample_id)
        pd.DataFrame(columns=header).to_csv(output_path, sep="\t", index=False)
        return

    obs_names = pd.Index(adata.obs_names.astype(str))
    barcode_to_idx = {normalize_barcode(barcode): idx for idx, barcode in enumerate(obs_names)}
    clone_cells["obs_idx"] = clone_cells["cell_barcode"].map(barcode_to_idx)
    clone_cells = clone_cells.dropna(subset=["obs_idx"]).copy()
    clone_cells["obs_idx"] = clone_cells["obs_idx"].astype(int)

    if clone_cells.empty:
        logger.warning("No clone cells from %s matched h5ad obs names; writing header only", sample_id)
        pd.DataFrame(columns=header).to_csv(output_path, sep="\t", index=False)
        return

    gene_indices, gene_names = select_genes(
        expression_source,
        max_genes=max_genes,
        min_detected_cells=min_detected_cells,
        always_include=DEFAULT_ALWAYS_INCLUDE_GENES,
    )
    if not gene_indices:
        logger.warning("No eligible genes found in %s; writing header only", sample_id)
        pd.DataFrame(columns=header).to_csv(output_path, sep="\t", index=False)
        return

    rows = []
    X = expression_source.X
    for clone_id, group in clone_cells.groupby("clone_id", sort=False):
        obs_idx = np.asarray(sorted(group["obs_idx"].unique()), dtype=int)
        if obs_idx.size < min_cells:
            continue
        sub = X[obs_idx, :][:, gene_indices]
        mean, variance, detected = matrix_column_stats(sub)
        n_cells = int(obs_idx.size)
        for idx, gene in enumerate(gene_names):
            rows.append(
                {
                    "study": study,
                    "sample_id": sample_id,
                    "clone_id": clone_id,
                    "gene": gene,
                    "expression_layer": expression_layer,
                    "n_cells": n_cells,
                    "mean_expression": mean[idx],
                    "variance_expression": variance[idx],
                    "detected_cells": int(detected[idx]),
                    "detection_fraction": float(detected[idx]) / n_cells if n_cells else 0.0,
                }
            )

    pd.DataFrame(rows, columns=header).to_csv(
        output_path,
        sep="\t",
        index=False,
        float_format="%.6g",
    )
    logger.info(
        "Wrote %s clone-gene summary rows for %s clones and %s genes to %s",
        len(rows),
        clone_cells["clone_id"].nunique(),
        len(gene_names),
        output_path,
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Precompute clone-level B-cell gene-expression summaries")
    parser.add_argument("--h5ad", required=True, help="Input B-cell h5ad file")
    parser.add_argument("--clone", required=True, help="Input clone TSV with clone_id and cell_barcode")
    parser.add_argument("--output", required=True, help="Output clone-gene summary TSV")
    parser.add_argument("--study", required=True, help="Study ID")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--min-cells", type=int, default=6, help="Minimum cells per clonotype")
    parser.add_argument("--max-clones", type=int, default=40, help="Maximum expanded clonotypes to summarize")
    parser.add_argument("--max-genes", type=int, default=1500, help="Maximum genes retained per clone")
    parser.add_argument("--min-detected-cells", type=int, default=10, help="Minimum detected cells for gene selection")
    parser.add_argument("--use-x", action="store_true", help="Use adata.X instead of adata.raw.X")
    return parser.parse_args()


def main():
    if "snakemake" in globals():
        precompute_clone_gene_expression(
            h5ad_path=Path(snakemake.input.h5ad),
            clone_path=Path(snakemake.input.clone),
            output_path=Path(snakemake.output[0]),
            study=snakemake.params.study,
            sample_id=snakemake.params.sample,
            min_cells=getattr(snakemake.params, "min_cells", 6),
            max_clones=getattr(snakemake.params, "max_clones", 40),
            max_genes=getattr(snakemake.params, "max_genes", 1500),
            min_detected_cells=getattr(snakemake.params, "min_detected_cells", 10),
            prefer_raw=getattr(snakemake.params, "prefer_raw", True),
        )
        return

    args = parse_args()
    precompute_clone_gene_expression(
        h5ad_path=Path(args.h5ad),
        clone_path=Path(args.clone),
        output_path=Path(args.output),
        study=args.study,
        sample_id=args.sample,
        min_cells=args.min_cells,
        max_clones=args.max_clones,
        max_genes=args.max_genes,
        min_detected_cells=args.min_detected_cells,
        prefer_raw=not args.use_x,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logger.error("Failed to precompute clone gene expression: %s", exc)
        sys.exit(1)

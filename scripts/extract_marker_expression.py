#!/usr/bin/env python3
"""
Extract marker-gene expression from each sample's B-cell h5ad for fast
featureplot queries.
"""

from pathlib import Path
import argparse
import logging
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


DEFAULT_B_CELL_MARKERS = [
    "CD19",
    "MS4A1",
    "CD79A",
    "CD79B",
    "CD74",
    "HLA-DRA",
    "HLA-DRB1",
    "CD37",
    "CD22",
    "CD24",
    "CD27",
    "CD38",
    "CD40",
    "CR2",
    "TCL1A",
    "IGHD",
    "IGHM",
    "IGHG1",
    "IGHG2",
    "IGHG3",
    "IGHG4",
    "IGHA1",
    "IGHA2",
    "IGKC",
    "IGLC1",
    "IGLC2",
    "JCHAIN",
    "MZB1",
    "XBP1",
    "PRDM1",
    "SDC1",
    "AICDA",
    "BCL6",
    "MKI67",
    "TOP2A",
    "CD69",
    "ITGAX",
    "FCRL5",
    "FCGR2B",
    "TNFRSF13B",
    "TNFRSF13C",
    "BANK1",
    "PAX5",
]


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def read_marker_genes(marker_file=None):
    """Return an ordered, de-duplicated marker gene list."""
    if marker_file:
        marker_path = Path(marker_file)
        with open(marker_path, "r", encoding="utf-8") as f:
            markers = [
                line.strip().split()[0]
                for line in f
                if line.strip() and not line.lstrip().startswith("#")
            ]
    else:
        markers = DEFAULT_B_CELL_MARKERS

    seen = set()
    ordered = []
    for gene in markers:
        if gene not in seen:
            seen.add(gene)
            ordered.append(gene)
    return ordered


def get_expression_source(adata, prefer_raw=True):
    """Choose the expression matrix used for marker expression export."""
    if prefer_raw and adata.raw is not None:
        return adata.raw, "raw_log1p"
    return adata, "X"


def dense_column(matrix):
    """Convert one gene column into a 1D numpy array."""
    if sparse.issparse(matrix):
        return matrix.toarray().ravel()
    return np.asarray(matrix).ravel()


def extract_marker_expression(
    h5ad_path,
    output_path,
    study,
    sample_id,
    marker_file=None,
    cell_scope="B_cells",
    prefer_raw=True,
):
    """Write long-format marker expression cache TSV for one sample."""
    h5ad_path = Path(h5ad_path)
    output_path = Path(output_path)

    if not h5ad_path.exists():
        raise FileNotFoundError(f"h5ad file not found: {h5ad_path}")

    logger.info("Reading %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    expression_source, expression_layer = get_expression_source(adata, prefer_raw=prefer_raw)

    marker_genes = read_marker_genes(marker_file)
    var_names = pd.Index(expression_source.var_names.astype(str))
    available_genes = [gene for gene in marker_genes if gene in var_names]
    missing_genes = [gene for gene in marker_genes if gene not in var_names]

    if missing_genes:
        logger.warning(
            "Sample %s missing %s marker genes: %s",
            sample_id,
            len(missing_genes),
            ", ".join(missing_genes),
        )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cell_barcodes = pd.Index(adata.obs_names.astype(str))
    with open(output_path, "w", encoding="utf-8") as out:
        out.write(
            "study\tsample_id\tcell_barcode\tcell_scope\tgene\t"
            "expression_layer\texpression\tdetected\n"
        )
        if not available_genes:
            logger.warning("No marker genes were found in %s; wrote header-only marker table", h5ad_path)
            return
        for gene in available_genes:
            gene_idx = int(var_names.get_loc(gene))
            values = dense_column(expression_source.X[:, gene_idx]).astype(float)

            if values.shape[0] != len(cell_barcodes):
                raise ValueError(f"{gene} expression length does not match obs length")

            gene_df = pd.DataFrame(
                {
                    "study": study,
                    "sample_id": sample_id,
                    "cell_barcode": cell_barcodes,
                    "cell_scope": cell_scope,
                    "gene": gene,
                    "expression_layer": expression_layer,
                    "expression": values,
                    "detected": values > 0,
                }
            )
            gene_df.to_csv(
                out,
                sep="\t",
                index=False,
                header=False,
                float_format="%.6g",
            )

    logger.info(
        "Wrote %s marker genes for %s cells to %s",
        len(available_genes),
        adata.n_obs,
        output_path,
    )


def parse_args():
    parser = argparse.ArgumentParser(description="Extract B-cell marker expression cache")
    parser.add_argument("--input", required=True, help="Input B-cell h5ad file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--study", required=True, help="Study ID")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--markers", help="Optional marker-gene text file")
    parser.add_argument("--cell-scope", default="B_cells", help="Cell group represented by this h5ad")
    parser.add_argument("--use-x", action="store_true", help="Use adata.X instead of adata.raw.X")
    return parser.parse_args()


def main():
    if "snakemake" in globals():
        extract_marker_expression(
            h5ad_path=snakemake.input.h5ad,
            output_path=snakemake.output[0],
            study=snakemake.params.study,
            sample_id=snakemake.params.sample,
            marker_file=getattr(snakemake.params, "markers", None),
            cell_scope=getattr(snakemake.params, "cell_scope", "B_cells"),
            prefer_raw=getattr(snakemake.params, "prefer_raw", True),
        )
        return

    args = parse_args()
    extract_marker_expression(
        h5ad_path=args.input,
        output_path=args.output,
        study=args.study,
        sample_id=args.sample,
        marker_file=args.markers,
        cell_scope=args.cell_scope,
        prefer_raw=not args.use_x,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logger.error("Failed to extract marker expression: %s", exc)
        sys.exit(1)

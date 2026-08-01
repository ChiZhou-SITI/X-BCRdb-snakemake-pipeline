#!/usr/bin/env python3
"""
Convert each sample's B-cell h5ad to zarr and write a database manifest.
"""

from pathlib import Path
import argparse
import logging
import shutil
import sys
from datetime import datetime

import anndata as ad
import pandas as pd


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def convert_h5ad_to_zarr(
    h5ad_path,
    zarr_path,
    manifest_path,
    study,
    sample_id,
    cell_scope="B_cells",
    chunks=(1024, 256),
):
    """Write zarr store plus one-row manifest for database loading."""
    h5ad_path = Path(h5ad_path)
    zarr_path = Path(zarr_path)
    manifest_path = Path(manifest_path)

    if not h5ad_path.exists():
        raise FileNotFoundError(f"h5ad file not found: {h5ad_path}")

    logger.info("Reading %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)

    zarr_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    if zarr_path.exists():
        shutil.rmtree(zarr_path)

    logger.info("Writing zarr store to %s", zarr_path)
    adata.write_zarr(zarr_path, chunks=chunks)

    if adata.raw is not None:
        expression_layer = "raw_log1p"
        n_genes = adata.raw.n_vars
    else:
        expression_layer = "X"
        n_genes = adata.n_vars

    manifest = pd.DataFrame(
        [
            {
                "study": study,
                "sample_id": sample_id,
                "cell_scope": cell_scope,
                "store_format": "zarr",
                "zarr_path": str(zarr_path),
                "expression_layer": expression_layer,
                "n_cells": int(adata.n_obs),
                "n_genes": int(n_genes),
                "source_h5ad_path": str(h5ad_path),
                "created_at": datetime.now().isoformat(timespec="seconds"),
            }
        ]
    )
    manifest.to_csv(manifest_path, sep="\t", index=False)
    logger.info("Wrote manifest to %s", manifest_path)


def parse_args():
    parser = argparse.ArgumentParser(description="Convert B-cell h5ad to zarr")
    parser.add_argument("--input", required=True, help="Input B-cell h5ad file")
    parser.add_argument("--zarr", required=True, help="Output zarr directory")
    parser.add_argument("--manifest", required=True, help="Output manifest TSV")
    parser.add_argument("--study", required=True, help="Study ID")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--cell-scope", default="B_cells", help="Cell group represented by this h5ad")
    return parser.parse_args()


def main():
    if "snakemake" in globals():
        convert_h5ad_to_zarr(
            h5ad_path=snakemake.input.h5ad,
            zarr_path=snakemake.output.zarr,
            manifest_path=snakemake.output.manifest,
            study=snakemake.params.study,
            sample_id=snakemake.params.sample,
            cell_scope=getattr(snakemake.params, "cell_scope", "B_cells"),
        )
        return

    args = parse_args()
    convert_h5ad_to_zarr(
        h5ad_path=args.input,
        zarr_path=args.zarr,
        manifest_path=args.manifest,
        study=args.study,
        sample_id=args.sample,
        cell_scope=args.cell_scope,
    )


if __name__ == "__main__":
    try:
        main()
    except ModuleNotFoundError as exc:
        if exc.name == "zarr":
            logger.error("zarr is required for h5ad->zarr conversion. Install it in the panbcr environment.")
        else:
            logger.error("Missing Python module: %s", exc)
        sys.exit(1)
    except Exception as exc:
        logger.error("Failed to convert h5ad to zarr: %s", exc)
        sys.exit(1)

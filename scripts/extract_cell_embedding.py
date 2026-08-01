#!/usr/bin/env python3
"""
Extract per-cell embedding coordinates and annotations from a sample h5ad.

In the current Snakemake workflow this script is wired to each sample's
B-cell h5ad file, so the generated TSV contains B-cell UMAP coordinates and
B-cell annotations for database loading.
"""

from pathlib import Path
import argparse
import logging
import sys

import anndata as ad
import numpy as np
import pandas as pd


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def extract_cell_embedding(
    h5ad_path,
    output_path,
    study,
    sample_id,
    embedding_key="X_umap",
    cell_scope="B_cells",
):
    """Write a clean TSV for database loading from a B-cell AnnData embedding."""
    h5ad_path = Path(h5ad_path)
    output_path = Path(output_path)

    if not h5ad_path.exists():
        raise FileNotFoundError(f"h5ad file not found: {h5ad_path}")

    logger.info("Reading %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)

    if embedding_key not in adata.obsm:
        available = ", ".join(adata.obsm.keys())
        raise KeyError(f"{embedding_key} not found in obsm. Available embeddings: {available}")

    embedding = np.asarray(adata.obsm[embedding_key])
    if embedding.ndim != 2 or embedding.shape[1] < 2:
        raise ValueError(f"{embedding_key} must be a 2D embedding with at least two columns")

    obs = adata.obs.copy()
    raw_cell_type = obs["cell_type"].astype(str) if "cell_type" in obs else pd.Series("", index=obs.index)
    raw_cell_subtype = obs["B_subtype"].astype(str) if "B_subtype" in obs else pd.Series("", index=obs.index)
    cell_type = raw_cell_type.copy()
    broad_b_mask = cell_type.str.strip().str.lower().isin(["b cells", "b cell"])
    cell_type.loc[broad_b_mask] = "Other B cells"

    if "cluster" in obs:
        cluster_clean = obs["cluster"].astype(str)
    elif "leiden" in obs:
        cluster_clean = obs["leiden"].astype(str)
    else:
        cluster_clean = pd.Series("", index=obs.index)
    cluster_clean = cluster_clean.replace({"nan": "", "None": ""})

    df = pd.DataFrame(
        {
            "study": study,
            "sample_id": sample_id,
            "cell_barcode": obs.index.astype(str),
            "cell_scope": cell_scope,
            "embedding_key": embedding_key,
            "umap_1": embedding[:, 0],
            "umap_2": embedding[:, 1],
            "cell_type": cell_type,
            "cell_subtype": raw_cell_subtype,
            "cluster": cluster_clean,
        }
    )

    df = df.replace({"nan": "", "None": ""})
    df = df[np.isfinite(df["umap_1"]) & np.isfinite(df["umap_2"])]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False, na_rep="")
    logger.info("Wrote %s rows to %s", len(df), output_path)


def parse_args():
    parser = argparse.ArgumentParser(description="Extract B-cell embedding from h5ad")
    parser.add_argument("--input", required=True, help="Input B-cell h5ad file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--study", required=True, help="Study ID")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--embedding-key", default="X_umap", help="AnnData obsm key")
    parser.add_argument("--cell-scope", default="B_cells", help="Cell group represented by this embedding")
    return parser.parse_args()


def main():
    if "snakemake" in globals():
        extract_cell_embedding(
            h5ad_path=snakemake.input.h5ad,
            output_path=snakemake.output[0],
            study=snakemake.params.study,
            sample_id=snakemake.params.sample,
            embedding_key=getattr(snakemake.params, "embedding_key", "X_umap"),
            cell_scope=getattr(snakemake.params, "cell_scope", "B_cells"),
        )
        return

    args = parse_args()
    extract_cell_embedding(
        h5ad_path=args.input,
        output_path=args.output,
        study=args.study,
        sample_id=args.sample,
        embedding_key=args.embedding_key,
        cell_scope=args.cell_scope,
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        logger.error("Failed to extract cell embedding: %s", exc)
        sys.exit(1)

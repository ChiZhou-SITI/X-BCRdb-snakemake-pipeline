#!/usr/bin/env python3
"""Fill Change-O c_call values from Cell Ranger filtered_contig_annotations.csv."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

MISSING = {"", "na", "nan", "none", "null", "unknown"}


def is_missing(value) -> bool:
    return str(value or "").strip().lower() in MISSING


def load_c_gene_map(annotation_path: Path) -> dict[str, str]:
    if not annotation_path.exists():
        return {}
    annot = pd.read_csv(annotation_path, dtype=str).fillna("")
    if "contig_id" not in annot.columns or "c_gene" not in annot.columns:
        return {}
    annot = annot[["contig_id", "c_gene"]].copy()
    annot["contig_id"] = annot["contig_id"].astype(str).str.strip()
    annot["c_gene"] = annot["c_gene"].astype(str).str.strip()
    annot = annot[(annot["contig_id"] != "") & (annot["c_gene"] != "")]
    return dict(zip(annot["contig_id"], annot["c_gene"]))


def fill_plain_table(df: pd.DataFrame, c_map: dict[str, str]) -> int:
    if "sequence_id" not in df.columns:
        return 0
    if "c_call" not in df.columns:
        df["c_call"] = ""
    before = df["c_call"].map(is_missing)
    mapped = df["sequence_id"].astype(str).map(c_map).fillna("")
    use = before & (mapped != "")
    df.loc[use, "c_call"] = mapped[use]
    return int(use.sum())


def fill_paired_table(df: pd.DataFrame, c_map: dict[str, str]) -> int:
    changed = 0
    for suffix in ["_H", "_L"]:
        seq_col = f"sequence_id{suffix}"
        c_col = f"c_call{suffix}"
        if seq_col not in df.columns:
            continue
        if c_col not in df.columns:
            df[c_col] = ""
        before = df[c_col].map(is_missing)
        mapped = df[seq_col].astype(str).map(c_map).fillna("")
        use = before & (mapped != "")
        df.loc[use, c_col] = mapped[use]
        changed += int(use.sum())
    return changed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--annotations", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--paired", action="store_true")
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    c_map = load_c_gene_map(Path(args.annotations))

    df = pd.read_csv(input_path, sep="\t", dtype=str).fillna("")
    changed = fill_paired_table(df, c_map) if args.paired else fill_plain_table(df, c_map)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    print(f"filled_c_call={changed} rows using {len(c_map)} contig annotations")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Export paired or IGH-only BCR records to a consistent wide JSON schema."""

import argparse
import json
from pathlib import Path

import pandas as pd


IDENTITY_COLUMNS = {"cell_barcode", "sample_id"}


def read_tsv(path):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep="\t", low_memory=False)


def clean_value(value):
    if value is None:
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        value = value.item()
    return value


def clean_records(df):
    return [
        {key: clean_value(value) for key, value in row.items()}
        for row in df.to_dict(orient="records")
    ]


def build_clone_lookup(clone_df):
    if clone_df.empty or "cell_barcode" not in clone_df.columns:
        return pd.DataFrame()
    columns = [
        column
        for column in ("cell_barcode", "clone_id", "cdr3_aa")
        if column in clone_df.columns
    ]
    return clone_df[columns].drop_duplicates("cell_barcode")


def enrich_paired_records(paired_df, clone_df):
    lookup = build_clone_lookup(clone_df)
    if lookup.empty:
        return paired_df
    rename = {
        column: f"{column}_H"
        for column in lookup.columns
        if column != "cell_barcode"
    }
    return paired_df.merge(lookup.rename(columns=rename), on="cell_barcode", how="left")


def build_igh_only_records(heavy_df, clone_df, paired_columns, sample_id):
    source = clone_df if not clone_df.empty else heavy_df
    if source.empty:
        return pd.DataFrame(columns=paired_columns)

    source = source.copy()
    if "sample_id" not in source.columns:
        source["sample_id"] = sample_id
    if "cell_barcode" not in source.columns and "sequence_id" in source.columns:
        source["cell_barcode"] = source["sequence_id"].astype(str).str.replace(
            r"_contig_.*$", "", regex=True
        )

    output_data = {"cell_barcode": source.get("cell_barcode")}
    for column in source.columns:
        if column in IDENTITY_COLUMNS:
            continue
        output_data[f"{column}_H"] = source[column]

    for column in paired_columns:
        if column not in output_data and column not in {"paired_status", "sample_id"}:
            output_data[column] = None
    output_data["paired_status"] = False
    output_data["sample_id"] = source["sample_id"].fillna(sample_id)
    output = pd.DataFrame(output_data, index=source.index)

    ordered = [column for column in paired_columns if column in output.columns]
    extra = [column for column in output.columns if column not in ordered]
    return output[ordered + extra]


def export_bcr_json(paired_path, heavy_path, clone_path, output_path, sample_id=""):
    paired_df = read_tsv(paired_path)
    heavy_df = read_tsv(heavy_path)
    clone_df = read_tsv(clone_path)

    if not paired_df.empty:
        output_df = enrich_paired_records(paired_df, clone_df)
    else:
        output_df = build_igh_only_records(
            heavy_df, clone_df, list(paired_df.columns), sample_id
        )

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(
            clean_records(output_df),
            handle,
            ensure_ascii=False,
            indent=2,
            allow_nan=False,
        )
    return len(output_df)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paired", required=True)
    parser.add_argument("--heavy", required=True)
    parser.add_argument("--clone", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample-id", default="")
    return parser.parse_args()


def main():
    args = parse_args()
    count = export_bcr_json(
        args.paired, args.heavy, args.clone, args.output, args.sample_id
    )
    print(f"Exported {count} BCR records to {args.output}")


if __name__ == "__main__":
    main()

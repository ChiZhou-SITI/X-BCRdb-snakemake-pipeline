import re
from pathlib import Path

import pandas as pd

input_file = snakemake.input[0]
study = snakemake.wildcards.study
sample = snakemake.wildcards.sample


def empty_paired_columns(base_columns):
    cols = ["cell_barcode"]
    cols += [f"{c}_H" for c in base_columns if c != "cell_barcode"]
    cols += [f"{c}_L" for c in base_columns if c != "cell_barcode"]
    cols += ["paired_status", "sample_id"]
    return cols


def extract_cell_barcode(sequence_id):
    """Extract the cell-level barcode from Cell Ranger or custom contig IDs."""
    sid = str(sequence_id)
    if "_contig_" in sid:
        return sid.rsplit("_contig_", 1)[0]

    match = re.match(r"^(.+)_(IGH|IGK|IGL|TRA|TRB|TRD|TRG)_\d+$", sid)
    if match:
        return match.group(1)

    return sid.split("_")[0]


def sample_requires_paired_bcr(study, sample):
    meta_path = Path("metadata/sample.tsv")
    if not meta_path.exists():
        return True

    meta = pd.read_csv(meta_path, sep="\t", dtype=str).fillna("")
    if "paired_BCR" not in meta.columns:
        return True

    rows = meta[(meta["study"] == study) & (meta["sample_id"] == sample)]
    if rows.empty:
        return True

    value = str(rows.iloc[0]["paired_BCR"]).strip().lower()
    return value not in {"false", "0", "f", "no", "n"}


df = pd.read_table(input_file)
if "sequence_id" not in df.columns:
    raise ValueError(f"{input_file} missing sequence_id column")
if "locus" not in df.columns:
    raise ValueError(f"{input_file} missing locus column")

# Extract cell barcode and keep productive contigs. Some Change-O files store
# booleans as strings, so normalize before filtering.
df["cell_barcode"] = df["sequence_id"].map(extract_cell_barcode)
productive = df.get("productive", True)
if hasattr(productive, "astype"):
    df = df[productive.astype(str).str.lower().isin(["true", "1", "t", "yes"])]

heavy = df[df["locus"] == "IGH"].copy()
paired = pd.DataFrame(columns=empty_paired_columns(df.columns))

if not df.empty and sample_requires_paired_bcr(study, sample):
    locus_counts = df.groupby(["cell_barcode", "locus"]).size().unstack(fill_value=0)
    valid_cells = locus_counts[
        (locus_counts.get("IGH", 0) == 1)
        & ((locus_counts.get("IGK", 0) + locus_counts.get("IGL", 0)) == 1)
    ]
    paired_df = df[df["cell_barcode"].isin(valid_cells.index)].copy()
    paired_heavy = paired_df[paired_df["locus"] == "IGH"].copy()
    paired_light = paired_df[paired_df["locus"].isin(["IGK", "IGL"])].copy()

    if not paired_heavy.empty and not paired_light.empty:
        paired = pd.merge(paired_heavy, paired_light, on="cell_barcode", suffixes=("_H", "_L"))
        paired["paired_status"] = True
        paired["sample_id"] = sample

heavy["sample_id"] = sample
paired.to_csv(snakemake.output[0], index=False, sep="\t")
heavy.to_csv(snakemake.output[1], index=False, sep="\t")

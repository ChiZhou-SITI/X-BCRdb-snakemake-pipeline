#!/usr/bin/env python3
"""Relabel clone_id by clone size (descending): clone_1, clone_2, ..."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def clone_sort_key(clone_id: str) -> tuple[int, int | str]:
    """Stable tie-breaker: numeric ids first (ascending), then lexical."""
    if clone_id.isdigit():
        return (0, int(clone_id))
    return (1, clone_id)


def relabel_clone_ids(df: pd.DataFrame) -> pd.DataFrame:
    if "clone_id" not in df.columns:
        raise ValueError("Input table does not contain 'clone_id' column.")

    clone_raw = df["clone_id"]
    clone_clean = clone_raw.astype(str).str.strip()
    valid_mask = clone_raw.notna() & clone_clean.ne("")

    if not valid_mask.any():
        return df

    counts = clone_clean[valid_mask].value_counts()
    ordered_ids = sorted(
        counts.index.tolist(),
        key=lambda cid: (-int(counts[cid]), clone_sort_key(cid)),
    )
    id_map = {old_id: f"clone_{idx + 1}" for idx, old_id in enumerate(ordered_ids)}

    out = df.copy()
    out.loc[valid_mask, "clone_id"] = clone_clean[valid_mask].map(id_map)
    return out


def translate_nt_to_aa(seq: str) -> str:
    clean = "".join(ch for ch in seq.upper() if ch in {"A", "C", "G", "T"})
    if len(clean) < 3:
        return ""
    aa = []
    for i in range(0, len(clean) - 2, 3):
        codon = clean[i:i + 3]
        aa.append(CODON_TABLE.get(codon, "X"))
    return "".join(aa).rstrip("*")


def add_cdr3_aa(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "cdr3_aa" not in out.columns:
        out["cdr3_aa"] = ""

    junction_aa = out["junction_aa"] if "junction_aa" in out.columns else pd.Series([""] * len(out))
    cdr3_nt = out["cdr3"] if "cdr3" in out.columns else pd.Series([""] * len(out))

    def derive(row_idx: int) -> str:
        aa_from_junction = str(junction_aa.iloc[row_idx]) if pd.notna(junction_aa.iloc[row_idx]) else ""
        aa_from_junction = aa_from_junction.strip()
        if aa_from_junction:
            return aa_from_junction
        nt = str(cdr3_nt.iloc[row_idx]) if pd.notna(cdr3_nt.iloc[row_idx]) else ""
        nt = nt.strip()
        if nt:
            return translate_nt_to_aa(nt)
        return ""

    out["cdr3_aa"] = [derive(i) for i in range(len(out))]
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rename clone_id by clone size descending with clone_ prefix."
    )
    parser.add_argument("--input", required=True, help="Input clone TSV file")
    parser.add_argument(
        "--output",
        required=False,
        help="Output clone TSV file (default: overwrite input)",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output) if args.output else input_path

    df = pd.read_csv(input_path, sep="\t", dtype={"clone_id": "object"})
    relabeled = relabel_clone_ids(df)
    with_cdr3_aa = add_cdr3_aa(relabeled)
    with_cdr3_aa.to_csv(output_path, sep="\t", index=False)


if __name__ == "__main__":
    main()

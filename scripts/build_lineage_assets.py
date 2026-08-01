#!/usr/bin/env python3
"""Build lineage edge table + lineage JSON for web visualization."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd


def parse_isotype(value: str | float | None) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "Unknown"
    raw = str(value).strip()
    if not raw:
        return "Unknown"
    first = raw.split(",")[0].strip()
    if "*" in first:
        first = first.split("*")[0].strip()
    return first or "Unknown"


def clean_seq(seq: str | float | None) -> str:
    if seq is None or (isinstance(seq, float) and math.isnan(seq)):
        return ""
    s = str(seq).upper().replace(".", "-")
    return "".join(ch if ch in {"A", "C", "G", "T", "N", "-"} else "N" for ch in s)


def hamming_distance(a: str, b: str) -> int:
    if not a and not b:
        return 0
    max_len = max(len(a), len(b))
    a = a.ljust(max_len, "-")
    b = b.ljust(max_len, "-")
    dist = 0
    for ca, cb in zip(a, b):
        if ca == cb:
            continue
        if ca == "N" or cb == "N":
            continue
        dist += 1
    return dist


def maybe_run_buildtrees(
    clone_tsv: Path,
    outdir: Path,
    outname: str,
    min_seq: int,
    nproc: int,
    igphyml_path: str | None = None,
) -> tuple[str, str]:
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "BuildTrees.py",
        "-d",
        str(clone_tsv),
        "--outdir",
        str(outdir),
        "--outname",
        outname,
        "--collapse",
        "--minseq",
        str(max(1, min_seq)),
    ]
    env = os.environ.copy()
    igphyml_bin = None
    if igphyml_path:
        p = Path(igphyml_path).expanduser().resolve()
        if p.is_file() and os.access(p, os.X_OK):
            igphyml_bin = str(p)
            env["PATH"] = f"{p.parent}:{env.get('PATH', '')}"
    if igphyml_bin is None:
        igphyml_bin = shutil.which("igphyml")

    method = "buildtrees_prep_only"
    if igphyml_bin:
        cmd.extend(["--igphyml", "--nproc", str(max(1, nproc)), "--clean", "none"])
        method = "buildtrees_igphyml"

    try:
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True, env=env)
        if proc.returncode == 0:
            return method, (proc.stdout or "").strip()
        return "buildtrees_failed", (proc.stderr or proc.stdout or "").strip()
    except FileNotFoundError:
        return "buildtrees_missing", "BuildTrees.py not found in PATH"


def build_clone_tree(clone_df: pd.DataFrame) -> tuple[list[dict], list[dict], list[dict], list[dict]]:
    clone_df = clone_df.copy()
    clone_df["mu_freq_num"] = pd.to_numeric(clone_df.get("mu_freq"), errors="coerce")
    clone_df["clean_sequence_alignment"] = clone_df["sequence_alignment"].map(clean_seq)
    clone_df["isotype_clean"] = clone_df.get("c_call").map(parse_isotype)
    if "locus" in clone_df.columns:
        locus = clone_df["locus"].fillna("").astype(str).str.upper().str.strip()
        unresolved_heavy = (clone_df["isotype_clean"] == "Unknown") & (locus == "IGH")
        clone_df.loc[unresolved_heavy, "isotype_clean"] = "IGH"
    clone_df = clone_df[clone_df["clean_sequence_alignment"] != ""].copy()
    clone_df = clone_df.sort_values(
        by=["mu_freq_num", "sequence_id"], ascending=[True, True], na_position="last"
    )

    germline = ""
    for col in ("germline_alignment_d_mask", "germline_alignment"):
        if col in clone_df.columns:
            vals = clone_df[col].dropna().astype(str).str.strip()
            vals = vals[vals != ""]
            if not vals.empty:
                germline = clean_seq(vals.iloc[0])
                break

    if not germline:
        first_seq = (
            clean_seq(clone_df["clean_sequence_alignment"].iloc[0]) if len(clone_df) else ""
        )
        germline = "N" * len(first_seq)

    clone_id = str(clone_df["clone_id"].iloc[0])
    germline_node = f"{clone_id}__germline"

    node_seq: dict[str, str] = {germline_node: germline}
    node_meta: dict[str, dict] = {
        germline_node: {"label": "Germline", "type": "germline", "mu_freq": 0.0}
    }
    order: list[str] = []
    edges: list[dict] = []
    edge_json: list[dict] = []
    grouped_nodes = []
    for idx, (_, group_df) in enumerate(
        clone_df.groupby("clean_sequence_alignment", dropna=False, sort=False), start=1
    ):
        seq = str(group_df["clean_sequence_alignment"].iloc[0])
        if not seq:
            continue

        node_id = f"{clone_id}__node_{idx}"
        count = int(len(group_df))
        mu_series = pd.to_numeric(group_df.get("mu_freq"), errors="coerce").dropna()
        mu_val = float(mu_series.mean()) if not mu_series.empty else None
        isotype = (
            group_df["isotype_clean"].fillna("Unknown").astype(str).value_counts().index[0]
            if "isotype_clean" in group_df.columns and len(group_df)
            else "Unknown"
        )
        representative = str(group_df["sequence_id"].iloc[0])
        grouped_nodes.append(
            {
                "node_id": node_id,
                "sequence": seq,
                "count": count,
                "mu_freq": mu_val,
                "isotype": isotype,
                "rep_id": representative,
            }
        )

    grouped_nodes.sort(
        key=lambda x: (
            x["mu_freq"] if x["mu_freq"] is not None else float("inf"),
            -x["count"],
            x["rep_id"],
        )
    )

    for node in grouped_nodes:
        node_id = node["node_id"]
        seq = node["sequence"]
        node_seq[node_id] = seq
        node_meta[node_id] = {
            "label": node["rep_id"],
            "type": "sequence",
            "mu_freq": node["mu_freq"],
            "cell_count": node["count"],
            "isotype": node["isotype"],
        }

        candidates = [germline_node] + order
        parent = min(candidates, key=lambda nid: hamming_distance(node_seq[nid], seq))
        dist = hamming_distance(node_seq[parent], seq)

        parent_type = "germline" if parent == germline_node else "sequence"
        child_type = "sequence"
        weight = 1.0 / (1.0 + float(dist))
        edges.append(
            {
                "parent_node": parent,
                "child_node": node_id,
                "parent_type": parent_type,
                "child_type": child_type,
                "distance": int(dist),
                "edge_weight": weight,
                "is_germline_edge": parent_type == "germline",
            }
        )
        edge_json.append(
            {
                "data": {
                    "id": f"{parent}->{node_id}",
                    "source": parent,
                    "target": node_id,
                    "distance": int(dist),
                    "weight": weight,
                }
            }
        )

        order.append(node_id)

    nodes_json = []
    for node_id, meta in node_meta.items():
        nodes_json.append(
            {
                "data": {
                    "id": node_id,
                    "label": meta["label"],
                    "node_type": meta["type"],
                    "mu_freq": meta["mu_freq"],
                    "cell_count": meta.get("cell_count", 0),
                    "isotype": meta.get("isotype", "Germline" if meta["type"] == "germline" else "Unknown"),
                }
            }
        )

    return nodes_json, edge_json, edges, clone_df.to_dict(orient="records")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build lineage edge table and lineage JSON from clone TSV."
    )
    parser.add_argument("--input", required=True, help="Input clone TSV")
    parser.add_argument("--study", required=True, help="Study ID")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--edge-out", required=True, help="Output lineage edge TSV")
    parser.add_argument("--tree-out", required=True, help="Output lineage tree TSV")
    parser.add_argument("--json-out", required=True, help="Output lineage JSON file")
    parser.add_argument(
        "--buildtrees-dir",
        required=False,
        help="Output directory for BuildTrees artifacts (optional)",
    )
    parser.add_argument("--min-clone-size", type=int, default=2)
    parser.add_argument("--nproc", type=int, default=1)
    parser.add_argument(
        "--igphyml-path",
        default="",
        help="Optional absolute path to igphyml executable when not on PATH.",
    )
    args = parser.parse_args()

    clone_path = Path(args.input)
    edge_out = Path(args.edge_out)
    tree_out = Path(args.tree_out)
    json_out = Path(args.json_out)
    edge_out.parent.mkdir(parents=True, exist_ok=True)
    tree_out.parent.mkdir(parents=True, exist_ok=True)
    json_out.parent.mkdir(parents=True, exist_ok=True)

    buildtrees_dir = (
        Path(args.buildtrees_dir)
        if args.buildtrees_dir
        else json_out.parent / f"{args.study}_{args.sample}_buildtrees"
    )
    buildtrees_method, buildtrees_message = maybe_run_buildtrees(
        clone_tsv=clone_path,
        outdir=buildtrees_dir,
        outname=f"{args.study}_{args.sample}",
        min_seq=args.min_clone_size,
        nproc=args.nproc,
        igphyml_path=args.igphyml_path or None,
    )

    df = pd.read_csv(clone_path, sep="\t", dtype={"clone_id": "object"})
    if "clone_id" not in df.columns or "sequence_id" not in df.columns:
        raise ValueError("Input clone table must contain clone_id and sequence_id columns.")
    if "sequence_alignment" not in df.columns:
        raise ValueError("Input clone table must contain sequence_alignment column.")

    clone_sizes = df["clone_id"].astype(str).value_counts()
    valid_clones = clone_sizes[clone_sizes >= max(1, args.min_clone_size)].index.tolist()
    sub = df[df["clone_id"].astype(str).isin(valid_clones)].copy()

    all_edges: list[dict] = []
    all_clone_rows: list[dict] = []
    clone_payloads: list[dict] = []

    for clone_id, clone_df in sub.groupby("clone_id", sort=False):
        clone_id = str(clone_id)
        nodes_json, edge_json, edge_rows, _ = build_clone_tree(clone_df)

        for row in edge_rows:
            all_edges.append(
                {
                    "study": args.study,
                    "sample_id": args.sample,
                    "clone_id": clone_id,
                    "parent_node": row["parent_node"],
                    "child_node": row["child_node"],
                    "parent_type": row["parent_type"],
                    "child_type": row["child_type"],
                    "distance": row["distance"],
                    "edge_weight": row["edge_weight"],
                    "is_germline_edge": row["is_germline_edge"],
                    "build_method": buildtrees_method,
                }
            )

        clone_json_obj = {
            "clone_id": clone_id,
            "n_nodes": len(nodes_json),
            "n_edges": len(edge_json),
            "layout": {"name": "breadthfirst", "roots": [f"{clone_id}__germline"]},
            "elements": {"nodes": nodes_json, "edges": edge_json},
        }
        all_clone_rows.append(
            {
                "study": args.study,
                "sample_id": args.sample,
                "clone_id": clone_id,
                "n_nodes": len(nodes_json),
                "n_edges": len(edge_json),
                "build_method": buildtrees_method,
                "tree_json": json.dumps(clone_json_obj, ensure_ascii=False, separators=(",", ":")),
            }
        )
        clone_payloads.append(clone_json_obj)

    edge_df = pd.DataFrame(
        all_edges,
        columns=[
            "study",
            "sample_id",
            "clone_id",
            "parent_node",
            "child_node",
            "parent_type",
            "child_type",
            "distance",
            "edge_weight",
            "is_germline_edge",
            "build_method",
        ],
    )
    edge_df.to_csv(edge_out, sep="\t", index=False, quoting=csv.QUOTE_NONE, escapechar="\\")

    tree_df = pd.DataFrame(
        all_clone_rows,
        columns=[
            "study",
            "sample_id",
            "clone_id",
            "n_nodes",
            "n_edges",
            "build_method",
            "tree_json",
        ],
    )
    tree_df.to_csv(tree_out, sep="\t", index=False, quoting=csv.QUOTE_NONE, escapechar="\\")

    payload = {
        "study": args.study,
        "sample_id": args.sample,
        "build_method": buildtrees_method,
        "buildtrees_message": buildtrees_message,
        "n_clones": len(clone_payloads),
        "clones": clone_payloads,
    }
    json_out.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

import re
from pathlib import Path

import numpy as np
import pandas as pd

input_files = list(snakemake.input)
output_file = snakemake.output[0]


def load_sample_study_map():
    p = Path("metadata/sample.tsv")
    if not p.exists():
        return {}
    meta = pd.read_csv(p, sep="\t", dtype=str)
    if "sample_id" not in meta.columns or "study" not in meta.columns:
        return {}
    return dict(zip(meta["sample_id"], meta["study"]))


def sample_from_clone_path(path, sample_to_study):
    name = Path(path).stem
    if name.endswith("_clone"):
        name = name[:-len("_clone")]
    matches = [sid for sid in sample_to_study if name == sid or name.endswith(f"_{sid}")]
    if matches:
        sid = max(matches, key=len)
        return sample_to_study.get(sid, ""), sid
    m = re.match(r"^(.*?)_(.*?)$", name)
    return (m.group(1), m.group(2)) if m else ("", name)


def shannon(p):
    p = p[p > 0]
    return float(-np.sum(p * np.log(p))) if len(p) else 0.0


def simpson(p):
    return float(np.sum(p**2)) if len(p) else 0.0


def gini(x):
    x = np.asarray(x, dtype=float)
    if len(x) == 0 or np.sum(x) == 0:
        return 0.0
    x = np.sort(x)
    n = len(x)
    cumulative = np.cumsum(x)
    return float((n + 1 - 2 * np.sum(cumulative) / cumulative[-1]) / n)


def d50(clone_counts):
    if len(clone_counts) == 0 or clone_counts.sum() == 0:
        return 0.0
    total = clone_counts.sum()
    sorted_counts = clone_counts.sort_values(ascending=False)
    cumulative = sorted_counts.cumsum()
    n = (cumulative <= total * 0.5).sum()
    return float(n / len(clone_counts)) if len(clone_counts) else 0.0


sample_to_study = load_sample_study_map()
results = []

for file in input_files:
    study, sample = sample_from_clone_path(file, sample_to_study)
    df = pd.read_csv(file, sep="\t")

    if "clone_id" not in df.columns or df.empty:
        results.append({
            "study": study,
            "sample_id": sample,
            "n_cells": 0,
            "n_clones": 0,
            "shannon": 0.0,
            "simpson": 0.0,
            "inv_simpson": 0.0,
            "evenness": 0.0,
            "gini": 0.0,
            "d50": 0.0,
            "top10_fraction": 0.0,
        })
        continue

    clone_counts = df["clone_id"].dropna().astype(str).value_counts()
    total_cells = int(clone_counts.sum())
    clone_freq = clone_counts / total_cells if total_cells else clone_counts
    richness = int(len(clone_counts))
    shannon_index = shannon(clone_freq.values)
    simpson_index = simpson(clone_freq.values)
    inv_simpson = float(1 / simpson_index) if simpson_index else 0.0
    evenness = float(shannon_index / np.log(richness)) if richness > 1 else 0.0
    sorted_counts = clone_counts.sort_values(ascending=False)
    top10_fraction = float(sorted_counts.head(10).sum() / total_cells) if total_cells else 0.0

    results.append({
        "study": study,
        "sample_id": sample,
        "n_cells": total_cells,
        "n_clones": richness,
        "shannon": shannon_index,
        "simpson": simpson_index,
        "inv_simpson": inv_simpson,
        "evenness": evenness,
        "gini": gini(clone_counts.values),
        "d50": d50(clone_counts),
        "top10_fraction": top10_fraction,
    })

pd.DataFrame(results).to_csv(output_file, index=False, sep="\t")

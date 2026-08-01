from collections import defaultdict
from itertools import combinations
from pathlib import Path
import math
import re

import pandas as pd


PUBLIC_MIN_SUBJECTS = 2
PUBLIC_IDENTITY = 0.85


def hamming_distance(a, b):
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def v_family(v):
    if pd.isna(v):
        return "NA"
    return str(v).split("-")[0]


def disease_slug(value):
    slug = re.sub(r"[^A-Za-z0-9]+", "_", str(value).strip())
    slug = re.sub(r"_+", "_", slug).strip("_")
    return slug or "Disease"


def is_healthy_disease(value):
    text = str(value or "").strip().lower()
    return text in {"healthy", "health", "normal", "control", "healthy control"}


def max_mismatches_for_identity(length, identity):
    return max(0, math.floor((1 - identity) * length + 1e-9))


def load_sample_metadata():
    p = Path("metadata/sample.tsv")
    if not p.exists():
        return pd.DataFrame(columns=["sample_id", "study", "subject_id", "disease"])
    meta = pd.read_csv(p, sep="\t", dtype=str).fillna("")
    required = {"sample_id", "study", "subject_id", "disease"}
    missing = required - set(meta.columns)
    if missing:
        raise ValueError(f"metadata/sample.tsv missing columns: {','.join(sorted(missing))}")
    meta["subject_id"] = meta["subject_id"].where(meta["subject_id"].str.strip() != "", meta["sample_id"])
    return meta[["sample_id", "study", "subject_id", "disease"]].drop_duplicates("sample_id")


def sample_from_clone_path(path, sample_ids):
    name = Path(path).stem
    if name.endswith("_clone"):
        name = name[: -len("_clone")]
    matches = [sid for sid in sample_ids if name == sid or name.endswith(f"_{sid}")]
    if matches:
        return max(matches, key=len)
    return name.split("_", 1)[1] if "_" in name else name


def write_empty(out_public, out_members, out_network):
    pd.DataFrame(columns=PUBLIC_COLS).to_csv(out_public, sep="\t", index=False)
    pd.DataFrame(columns=MEMBER_COLS).to_csv(out_members, sep="\t", index=False)
    pd.DataFrame(columns=NETWORK_COLS).to_csv(out_network, sep="\t", index=False)


def segment_ranges(length, parts):
    parts = max(1, min(parts, length))
    base = length // parts
    remainder = length % parts
    ranges = []
    start = 0
    for i in range(parts):
        end = start + base + (1 if i < remainder else 0)
        ranges.append((start, end))
        start = end
    return ranges


def representative_clusters(unique_seqs, max_dist, seq_to_subjects):
    if max_dist <= 0:
        for seq in unique_seqs:
            if len(seq_to_subjects[seq]) >= PUBLIC_MIN_SUBJECTS:
                yield [seq]
        return

    ranges = segment_ranges(len(unique_seqs[0]), max_dist + 1)
    buckets = defaultdict(list)
    for seq in unique_seqs:
        for i, (start, end) in enumerate(ranges):
            buckets[(i, seq[start:end])].append(seq)

    assigned = set()
    for seed in unique_seqs:
        if seed in assigned:
            continue

        candidates = set()
        for i, (start, end) in enumerate(ranges):
            candidates.update(buckets.get((i, seed[start:end]), []))

        neighbors = [
            seq for seq in candidates
            if seq not in assigned and hamming_distance(seed, seq) <= max_dist
        ]
        if not neighbors:
            continue

        subjects = set().union(*(seq_to_subjects[seq] for seq in neighbors))
        if len(subjects) < PUBLIC_MIN_SUBJECTS:
            continue

        assigned.update(neighbors)
        yield sorted(neighbors, key=lambda s: (-len(seq_to_subjects[s]), s))


PUBLIC_COLS = ["clonotype_id", "n_cells", "n_samples", "samples"]
MEMBER_COLS = ["clonotype_id", "sequence_id", "sample_id", "cdr3", "v_gene", "j_gene"]
NETWORK_COLS = ["sample1", "sample2", "clonotype_id"]

files = list(snakemake.input)
out_public = snakemake.output.public
out_network = snakemake.output.network
out_members = snakemake.output.members

sample_meta = load_sample_metadata()
sample_ids = sample_meta["sample_id"].tolist()
sample_to_disease = dict(zip(sample_meta["sample_id"], sample_meta["disease"]))
sample_to_subject = dict(zip(sample_meta["sample_id"], sample_meta["subject_id"]))

records = []
member_rows = []
clonotype_counter = 0

dfs = []
for f in files:
    df = pd.read_csv(f, sep="\t")
    if df.empty:
        continue
    sample_id = sample_from_clone_path(f, sample_ids)
    disease = sample_to_disease.get(sample_id, "")
    if is_healthy_disease(disease):
        continue
    df["sample"] = sample_id
    df["subject"] = sample_to_subject.get(sample_id, sample_id) or sample_id
    df["disease"] = disease
    dfs.append(df)

if not dfs:
    write_empty(out_public, out_members, out_network)
    raise SystemExit(0)

df = pd.concat(dfs, ignore_index=True)
required = ["junction_aa", "v_call", "j_call", "sequence_id", "sample", "subject", "disease"]
if any(col not in df.columns for col in required):
    write_empty(out_public, out_members, out_network)
    raise SystemExit(0)

df = df.dropna(subset=["junction_aa", "v_call", "j_call", "sample", "subject", "disease"])
df = df[df["junction_aa"].astype(str).str.strip() != ""]
df = df[df["subject"].astype(str).str.strip() != ""]
df = df[df["disease"].astype(str).str.strip() != ""]
df = df[~df["disease"].map(is_healthy_disease)]

if df.empty:
    write_empty(out_public, out_members, out_network)
    raise SystemExit(0)

df["junction_aa"] = df["junction_aa"].astype(str)
df["v_family"] = df["v_call"].apply(v_family)
df["cdr3_len"] = df["junction_aa"].str.len()
print(
    f"Loaded {len(df):,} disease-state BCR records from "
    f"{df['sample'].nunique():,} samples, {df['subject'].nunique():,} subjects, "
    f"across {df['disease'].nunique():,} diseases.",
    flush=True,
)

for disease, disease_df in df.groupby("disease", sort=True):
    disease_df = disease_df.copy()
    if disease_df["subject"].nunique() < PUBLIC_MIN_SUBJECTS:
        continue

    print(
        f"Scanning disease={disease!r}: {len(disease_df):,} records, "
        f"{disease_df['sample'].nunique():,} samples, {disease_df['subject'].nunique():,} subjects.",
        flush=True,
    )

    groups = disease_df.groupby(["v_family", "j_call", "cdr3_len"], dropna=False, sort=True)
    for (_v_family, _j_call, cdr3_len), sub in groups:
        if sub["subject"].nunique() < PUBLIC_MIN_SUBJECTS:
            continue

        max_dist = max_mismatches_for_identity(int(cdr3_len), PUBLIC_IDENTITY)
        seq_to_idx = defaultdict(list)
        seq_to_subjects = defaultdict(set)
        seq_to_samples = defaultdict(set)
        for idx, row in sub.iterrows():
            seq = row["junction_aa"]
            seq_to_idx[seq].append(idx)
            seq_to_subjects[seq].add(row["subject"])
            seq_to_samples[seq].add(row["sample"])

        unique_seqs = sorted(seq_to_idx, key=lambda s: (-len(seq_to_idx[s]), s))
        if len(unique_seqs) < 2:
            continue

        for component in representative_clusters(unique_seqs, max_dist, seq_to_subjects):
            samples = set().union(*(seq_to_samples[seq] for seq in component))
            clonotype_counter += 1
            cid = f"{disease_slug(disease)}__clono_{clonotype_counter}"

            n_cells = 0
            for seq in component:
                for idx in seq_to_idx[seq]:
                    row = sub.loc[idx]
                    n_cells += 1
                    member_rows.append(
                        {
                            "clonotype_id": cid,
                            "sequence_id": row["sequence_id"],
                            "sample_id": row["sample"],
                            "cdr3": row["junction_aa"],
                            "v_gene": row["v_call"],
                            "j_gene": row["j_call"],
                        }
                    )

            records.append(
                {
                    "clonotype_id": cid,
                    "n_cells": n_cells,
                    "n_samples": len(samples),
                    "samples": ",".join(sorted(samples)),
                }
            )

public = pd.DataFrame(records, columns=PUBLIC_COLS)
public.to_csv(out_public, sep="\t", index=False)

members_df = pd.DataFrame(member_rows, columns=MEMBER_COLS)
members_df.to_csv(out_members, sep="\t", index=False)

edges = []
if not public.empty and not members_df.empty:
    for cid, sub in members_df.groupby("clonotype_id", sort=False):
        samples = sorted(set(sub["sample_id"]))
        for s1, s2 in combinations(samples, 2):
            edges.append((s1, s2, cid))

pd.DataFrame(edges, columns=NETWORK_COLS).to_csv(out_network, sep="\t", index=False)

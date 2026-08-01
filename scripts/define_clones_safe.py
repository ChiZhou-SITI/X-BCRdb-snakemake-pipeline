#!/usr/bin/env python3
"""Run DefineClones, but allow samples with no heavy-chain rows to continue."""

import subprocess
import sys
from pathlib import Path

import pandas as pd

input_file = Path(snakemake.input[0])
output_file = Path(snakemake.output.clone_tsv if hasattr(snakemake.output, "clone_tsv") else snakemake.output[0])

df = pd.read_csv(input_file, sep="\t")
if df.empty:
    if "clone_id" not in df.columns:
        df["clone_id"] = pd.Series(dtype="object")
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, sep="\t", index=False)
    sys.exit(0)

subprocess.run(
    [
        "DefineClones.py",
        "-d",
        str(input_file),
        "--act",
        "set",
        "--model",
        "ham",
        "--norm",
        "len",
        "-o",
        str(output_file),
    ],
    check=True,
)
subprocess.run(
    [sys.executable, "scripts/relabel_clone_ids.py", "--input", str(output_file)],
    check=True,
)

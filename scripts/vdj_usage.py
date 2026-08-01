import pandas as pd

study = snakemake.wildcards.study
sample = snakemake.wildcards.sample
input_file = snakemake.input[0]
output_file = snakemake.output[0]

out_cols = ["study", "sample_id", "vdj_type", "vdj_gene", "count", "frequency"]
df = pd.read_csv(input_file, sep="\t")

if "locus" in df.columns:
    df = df[df["locus"] == "IGH"].copy()

required = ["v_call", "d_call", "j_call"]
if df.empty or any(col not in df.columns for col in required):
    pd.DataFrame(columns=out_cols).to_csv(output_file, index=False, sep="\t")
    raise SystemExit(0)

df = df.dropna(subset=required)
df = df[(df["v_call"].astype(str).str.strip() != "") & (df["d_call"].astype(str).str.strip() != "") & (df["j_call"].astype(str).str.strip() != "")]
if df.empty:
    pd.DataFrame(columns=out_cols).to_csv(output_file, index=False, sep="\t")
    raise SystemExit(0)

v_usage = df["v_call"].value_counts().reset_index()
v_usage.columns = ["vdj_gene", "count"]
v_usage["vdj_type"] = "V"

j_usage = df["j_call"].value_counts().reset_index()
j_usage.columns = ["vdj_gene", "count"]
j_usage["vdj_type"] = "J"

vdj = df.groupby(["v_call", "d_call", "j_call"]).size().reset_index(name="count")
vdj["vdj_gene"] = vdj["v_call"].astype(str) + "_" + vdj["d_call"].astype(str) + "_" + vdj["j_call"].astype(str)
vdj["vdj_type"] = "VDJ"
vdj = vdj[["vdj_gene", "vdj_type", "count"]]

usage = pd.concat([
    v_usage[["vdj_gene", "vdj_type", "count"]],
    j_usage[["vdj_gene", "vdj_type", "count"]],
    vdj,
], ignore_index=True)
usage["study"] = study
usage["sample_id"] = sample
total = usage["count"].sum()
usage["frequency"] = usage["count"] / total if total else 0
usage = usage[out_cols]
usage.to_csv(output_file, index=False, sep="\t")

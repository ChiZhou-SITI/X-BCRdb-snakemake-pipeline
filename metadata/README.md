# Metadata files

The workflow expects four tab-separated metadata files in this directory:

| File | Required by | Purpose |
| --- | --- | --- |
| `pipeline_samples.tsv` | Snakemake | Maps each sample to Cell Ranger GEX and V(D)J input paths. |
| `study.tsv` | Database loader | Study-level metadata loaded into PostgreSQL. |
| `subject.tsv` | Database loader | Subject-level metadata loaded into PostgreSQL. |
| `sample.tsv` | Database loader | Sample-level metadata loaded into PostgreSQL. |

Example files are provided as:

- `pipeline_samples.example.tsv`
- `study.example.tsv`
- `subject.example.tsv`
- `sample.example.tsv`

Before running the workflow, copy the examples to the expected filenames and edit
them for your own data:

```bash
cp metadata/pipeline_samples.example.tsv metadata/pipeline_samples.tsv
cp metadata/study.example.tsv metadata/study.tsv
cp metadata/subject.example.tsv metadata/subject.tsv
cp metadata/sample.example.tsv metadata/sample.tsv
```

## Input path requirements

`pipeline_samples.tsv` must contain absolute or project-relative paths to:

- `count_dir`: a Cell Ranger filtered gene-expression matrix directory, for example
  `raw_data/STUDY/SAMPLE/count/sample_feature_bc_matrix`.
- `vdj_fasta`: a Cell Ranger BCR FASTA file, typically
  `raw_data/STUDY/SAMPLE/vdj_b/filtered_contig.fasta`.

The workflow automatically expects the matching Cell Ranger annotation file in
the same directory as `vdj_fasta`:

```text
filtered_contig_annotations.csv
```

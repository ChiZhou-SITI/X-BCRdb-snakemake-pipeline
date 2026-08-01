# X-BCRdb Snakemake Pipeline

Reproducible Snakemake workflow for building the processed transcriptome, BCR repertoire, lineage, public-clonotype and PostgreSQL-ready tables used by **X-BCRdb**, a cross-disease BCR-paired single-cell database.

The pipeline starts from Cell Ranger-style single-cell gene-expression matrices and BCR V(D)J contig files, identifies and annotates B cells, defines BCR clonotypes, calculates SHM and isotype features, builds lineage assets, precomputes expression caches for web visualization, and optionally loads all processed outputs into PostgreSQL.

## What this repository contains

```text
X-BCRdb-snakemake-pipeline/
├── Snakefile                         # Full workflow, including PostgreSQL loading
├── Snakefile.no_db                   # Processing workflow without PostgreSQL loading
├── config.yaml                       # Editable configuration template
├── config.example.yaml               # Example configuration
├── envs/
│   └── x-bcrdb-pipeline.yaml         # Conda environment template
├── metadata/
│   ├── README.md
│   ├── pipeline_samples.example.tsv  # Sample-to-raw-data mapping example
│   ├── sample.example.tsv            # Sample metadata example
│   ├── subject.example.tsv           # Subject metadata example
│   └── study.example.tsv             # Study metadata example
├── scripts/
│   ├── gex_pipeline.py               # Scanpy/CellTypist B-cell processing
│   ├── bcr_pairing.py                # Heavy/light chain pairing
│   ├── define_clones_safe.py         # Change-O DefineClones wrapper
│   ├── shazam.r                      # SHM calculation with ShazaM/Alakazam
│   ├── build_lineage_assets.py       # Dowser/IgPhyML lineage assets
│   ├── public_clonotype.py           # Disease/study public clonotype discovery
│   ├── load_database.py              # PostgreSQL loader
│   └── ...
├── sql/
│   ├── schema.sql                    # PostgreSQL schema
│   └── indexes.sql                   # Additional indexes
└── docs/
    └── output_files.md               # Output file descriptions
```

Large local data directories such as `raw_data/`, `results/`, `cache/` and `.snakemake/` are intentionally excluded from this repository.

## Pipeline overview

The workflow contains four major layers.

1. **Transcriptome processing**
   - Read Cell Ranger filtered gene-expression matrix.
   - Perform quality control, normalization, log transformation, highly variable gene selection, PCA, neighbor graph construction, UMAP and Leiden clustering.
   - Annotate broad immune cell types with CellTypist.
   - Extract B-cell-lineage cells and rerun B-cell-focused preprocessing.
   - Annotate B-cell subtypes with the high-resolution CellTypist immune model.
   - Export B-cell h5ad files, UMAP coordinates, Cell Subtype labels, marker-expression cache and Zarr stores.

2. **BCR repertoire processing**
   - Run Change-O `AssignGenes.py igblast` and `MakeDb.py`.
   - Reconstruct germlines with `CreateGermlines.py`.
   - Calculate SHM frequency with ShazaM/Alakazam.
   - Supplement missing constant-region calls from Cell Ranger `filtered_contig_annotations.csv`.
   - Pair heavy and light chains at the cell-barcode level.
   - Define sample-specific clonotypes from productive heavy-chain records.
   - Relabel clone IDs by clone size (`clone_1`, `clone_2`, ...).
   - Export per-sample BCR JSON files.

3. **BCR-specific downstream analysis**
   - V/J/VDJ usage summaries.
   - Clone diversity metrics.
   - CDR3 amino acid similarity networks.
   - Lineage tree assets using Dowser/IgPhyML for clones with at least two cells.
   - Public Clonotype discovery within study or disease context.
   - Clone-level expression summaries for BCR-transcriptome analysis.

4. **Database loading**
   - Load study, subject and sample metadata.
   - Load BCR records, clone membership, VDJ usage, diversity, lineage edges/nodes, B-cell UMAP embeddings, marker-expression caches and Zarr manifests.
   - Create PostgreSQL-ready tables used by the X-BCRdb web application.

## Input requirements

The workflow assumes each sample has Cell Ranger-style GEX and BCR outputs.

### Gene-expression input

`pipeline_samples.tsv` column `count_dir` should point to a filtered feature-barcode matrix directory, for example:

```text
raw_data/STUDY_ID/SAMPLE_ID/count/sample_feature_bc_matrix/
├── matrix.mtx.gz
├── features.tsv.gz
└── barcodes.tsv.gz
```

The pipeline also supports several non-standard matrix naming variants through fallback parsing in `scripts/gex_pipeline.py`.

### BCR V(D)J input

`pipeline_samples.tsv` column `vdj_fasta` should point to:

```text
raw_data/STUDY_ID/SAMPLE_ID/vdj_b/filtered_contig.fasta
```

The matching Cell Ranger annotation file must be in the same directory:

```text
raw_data/STUDY_ID/SAMPLE_ID/vdj_b/filtered_contig_annotations.csv
```

## Metadata files

Before running, create the expected metadata files:

```bash
cp metadata/pipeline_samples.example.tsv metadata/pipeline_samples.tsv
cp metadata/study.example.tsv metadata/study.tsv
cp metadata/subject.example.tsv metadata/subject.tsv
cp metadata/sample.example.tsv metadata/sample.tsv
```

Then edit them for your data.

### `metadata/pipeline_samples.tsv`

Required columns:

| Column | Description |
| --- | --- |
| `study` | Harmonized study ID used by Snakemake and database outputs. |
| `sample_id` | Harmonized sample ID. |
| `source_sample_id` | Original source sample ID. |
| `raw_study` | Raw-data study folder or original study label. |
| `raw_count_sample` | Raw-data sample folder for GEX matrix. |
| `raw_vdj_sample` | Raw-data sample folder for BCR V(D)J output. |
| `count_dir` | Path to Cell Ranger filtered matrix directory. |
| `vdj_fasta` | Path to Cell Ranger `filtered_contig.fasta`. |

### `metadata/study.tsv`

Required columns:

| Column | Description |
| --- | --- |
| `Study` | Study ID. |
| `Disease` | Harmonized disease label. |
| `Title` | Publication title. |
| `Publication` | Citation text. |
| `Patient count` | Number of subjects/patients reported for the study. |
| `Sample count` | Number of samples included. |
| `Accession` | GEO/SRA/GSA/Zenodo accession. |

### `metadata/subject.tsv`

Required columns:

| Column | Description |
| --- | --- |
| `subject_id` | Harmonized subject ID. |
| `species` | Species, usually `Human`. |
| `disease` | Subject-level disease/condition. |

### `metadata/sample.tsv`

Required columns:

| Column | Description |
| --- | --- |
| `sample_id` | Harmonized sample ID. |
| `study` | Study ID matching `study.tsv`. |
| `subject_id` | Subject ID matching `subject.tsv`. |
| `tissue` | Harmonized tissue label. |
| `disease` | Harmonized sample-level disease label. |
| `platform` | Sequencing platform, for example `10X 5'`. |
| `n_cells` | B-cell count after processing, if known before loading. |
| `paired_BCR` | Whether the sample has paired IGH and IGK/IGL recovery. |

## Installation

Create the conda environment:

```bash
conda env create -f envs/x-bcrdb-pipeline.yaml
conda activate xbcrdb-pipeline
```

Install or configure external immunoglobulin tools:

- IgBLAST database for human immunoglobulin sequences.
- IMGT human germline references for Change-O.
- Change-O command-line tools: `AssignGenes.py`, `MakeDb.py`, `CreateGermlines.py`, `DefineClones.py`.
- ShazaM/Alakazam R packages for SHM calculation.
- Dowser and IgPhyML for lineage tree reconstruction.
- PostgreSQL, if database loading is required.

Edit `config.yaml`:

```yaml
igblast_db: "/path/to/igblast/database"
vdj_db: "/path/to/imgt/human/vdj"
igphyml_path: "/path/to/igphyml"

pg:
  db: "bcrdb"
  user: "postgres_user"
  password: "postgres_password"
  host: "localhost"
```

## Running the workflow

### Dry run

```bash
snakemake -n --cores 1
```

### Run all processing and load PostgreSQL

```bash
snakemake --cores 16 --rerun-incomplete
```

### Run processing only, without PostgreSQL loading

```bash
snakemake -s Snakefile.no_db all_no_db --cores 16 --rerun-incomplete
```

### Run one specific sample-level target

```bash
snakemake results/bcr/STUDY_SAMPLE_clone.tsv --cores 8
```

Replace `STUDY_SAMPLE` with the output prefix generated from `{study}_{sample}`.

## PostgreSQL setup

Create a database:

```bash
createdb bcrdb
```

Load schema manually if desired:

```bash
psql -d bcrdb -f sql/schema.sql
psql -d bcrdb -f sql/indexes.sql
```

The default `rule all` runs `scripts/load_database.py` and writes:

```text
results/database_loaded.txt
```

If you only want to generate files and inspect them before database loading, use `Snakefile.no_db`.

## Public Clonotype definition

Public Clonotypes are identified within study or disease context using a biologically constrained CDR3 amino acid similarity rule. Candidate records are grouped by:

- the same V family,
- the same J gene,
- identical CDR3 amino acid length,
- CDR3 amino acid identity of at least 85%,
- recurrence in at least two subjects.

Healthy/control samples are excluded from disease-specific Public Clonotype discovery in the current implementation.

## Clonotype definition

Sample-specific clonotypes are defined from productive heavy-chain records using Change-O `DefineClones.py` with:

- `--act set`,
- a Hamming-distance model,
- length-normalized junction distance,
- distance threshold controlled by `config.yaml` where applicable.

After clone inference, clone IDs are relabeled by decreasing clone size within each sample:

```text
clone_1, clone_2, clone_3, ...
```

## SHM calculation

SHM frequency is calculated from germline-aware Change-O tables. The R script `scripts/shazam.r` calls ShazaM/Alakazam `observedMutations` using:

- observed sequence column: `sequence_alignment`,
- inferred germline column: `germline_alignment_d_mask`,
- IMGT V(D)J region definition,
- mutation frequency output enabled,
- combined mutation counts across evaluated immunoglobulin regions.

Constant-region calls are used for Isotype analysis but are not used to calculate SHM frequency.

## Output files

See [docs/output_files.md](docs/output_files.md) for a complete list of generated files.

Main output groups:

```text
results/gex/           # B-cell h5ad, UMAP metadata, marker expression cache
results/gex_de/        # clone-level expression summaries
results/gex_zarr/      # B-cell zarr expression stores
results/changeo/       # IgBLAST, MakeDb, germline and SHM intermediates
results/bcr/           # paired BCR tables, heavy-chain clone tables, BCR JSON
results/lineage/       # lineage tree nodes, edges and JSON
results/repertoire/    # diversity, VDJ usage, public clonotypes
results/cdr3_network/  # sample-level CDR3 similarity network files
```

## Notes for GitHub users

This repository is a code and workflow release. It does **not** include:

- raw FASTQ files,
- Cell Ranger raw output directories,
- processed X-BCRdb h5ad files,
- PostgreSQL dumps,
- local cache files,
- Zenodo upload packages.

For reproducible use, provide your own Cell Ranger GEX/BCR inputs and metadata files following the templates in `metadata/`.

## Citation

If you use this workflow, please cite the X-BCRdb manuscript after publication.

## License

Please add a license file before public release, according to your institutional and collaborator requirements.

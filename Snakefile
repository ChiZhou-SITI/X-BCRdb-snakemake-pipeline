configfile: "config.yaml"

import pandas as pd
import re
from pathlib import Path

########################################
# sample sheet
########################################
# metadata/sample.tsv is the database-facing metadata table.
# metadata/pipeline_samples.tsv carries the raw-data path mapping used by Snakemake.
PIPELINE_SAMPLES = pd.read_csv("metadata/pipeline_samples.tsv", sep="	")

STUDIES = PIPELINE_SAMPLES["study"].tolist()
SAMPLE_NAMES = PIPELINE_SAMPLES["sample_id"].tolist()
SAMPLE_KEYS = list(zip(STUDIES, SAMPLE_NAMES))
RAW_COUNT_DIRS = dict(zip(SAMPLE_KEYS, PIPELINE_SAMPLES["count_dir"].tolist()))
RAW_VDJ_FASTA = dict(zip(SAMPLE_KEYS, PIPELINE_SAMPLES["vdj_fasta"].tolist()))
RAW_VDJ_ANNOTATIONS = {key: str(Path(path).with_name("filtered_contig_annotations.csv")) for key, path in RAW_VDJ_FASTA.items()}
STUDY_PATTERN = "|".join(re.escape(study) for study in sorted(set(STUDIES), key=len, reverse=True))
IGPHYML_PATH = str(config.get("igphyml_path", "") or "").strip()

########################################
# wildcard constraints
########################################
wildcard_constraints:
    study=STUDY_PATTERN,
    sample=".+"

########################################
# final targets
########################################
rule all:
    input:
        expand("results/gex/{study}_{sample}_all_cells.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/gex/{study}_{sample}_B_cells.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/gex/{study}_{sample}_all.h5ad",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/gex/{study}_{sample}_b.h5ad",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/gex/{study}_{sample}_b_cell_embedding.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/gex/{study}_{sample}_b_marker_expression.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/gex_zarr/{study}_{sample}_b_zarr_manifest.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/bcr/{study}_{sample}_SHM_paired.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/bcr/{study}_{sample}_bcr.json",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/bcr/{study}_{sample}_clone.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/lineage/{study}_{sample}_lineage_edge.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/lineage/{study}_{sample}_lineage_tree.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/lineage/{study}_{sample}_lineage.json",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/repertoire/{study}_{sample}_vdj_usage.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/cdr3_network/{study}_{sample}_cdr3_edges.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/cdr3_network/{study}_{sample}_cdr3_nodes.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/cdr3_network/{study}_{sample}_cdr3_network.json",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        expand("results/cdr3_network/{study}_{sample}_cdr3_network_stats.json",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
        +
        ["results/repertoire/diversity.tsv"]
        +
        ["results/repertoire/public_clonotypes.tsv"]
        +
        ["results/repertoire/public_network.tsv"]
        +
        ["results/database_loaded.txt"]
########################################
# scRNA pipeline
########################################
rule gex:
    input:
        lambda wildcards: RAW_COUNT_DIRS[(wildcards.study, wildcards.sample)]
    output:
        all="results/gex/{study}_{sample}_all_cells.tsv",
        b="results/gex/{study}_{sample}_B_cells.tsv",
        all_h5ad="results/gex/{study}_{sample}_all.h5ad",
        b_h5ad="results/gex/{study}_{sample}_b.h5ad"
    threads: 10
    script:
        "scripts/gex_pipeline.py"


########################################
# B cell UMAP/celltype extraction
########################################
rule extract_b_cell_embedding:
    input:
        h5ad="results/gex/{study}_{sample}_b.h5ad"
    output:
        "results/gex/{study}_{sample}_b_cell_embedding.tsv"
    params:
        study=lambda wildcards: wildcards.study,
        sample=lambda wildcards: wildcards.sample,
        embedding_key="X_umap",
        cell_scope="B_cells"
    script:
        "scripts/extract_cell_embedding.py"


########################################
# B cell marker expression cache
########################################
rule extract_b_cell_marker_expression:
    input:
        h5ad="results/gex/{study}_{sample}_b.h5ad"
    output:
        "results/gex/{study}_{sample}_b_marker_expression.tsv"
    params:
        study=lambda wildcards: wildcards.study,
        sample=lambda wildcards: wildcards.sample,
        markers=None,
        cell_scope="B_cells",
        prefer_raw=True
    script:
        "scripts/extract_marker_expression.py"


########################################
# Clone-level gene expression summary for BCR-transcriptome analysis
########################################
rule precompute_clone_gene_expression:
    input:
        h5ad="results/gex/{study}_{sample}_b.h5ad",
        clone="results/bcr/{study}_{sample}_clone.tsv"
    output:
        "results/gex_de/{study}_{sample}_clone_gene_expression_summary.tsv"
    params:
        study=lambda wildcards: wildcards.study,
        sample=lambda wildcards: wildcards.sample,
        min_cells=6,
        max_clones=40,
        max_genes=1500,
        min_detected_cells=10,
        prefer_raw=True
    script:
        "scripts/precompute_clone_gene_expression.py"


########################################
# B cell h5ad to zarr
########################################
rule b_cell_h5ad_to_zarr:
    input:
        h5ad="results/gex/{study}_{sample}_b.h5ad"
    output:
        zarr=directory("results/gex_zarr/{study}_{sample}_b.zarr"),
        manifest="results/gex_zarr/{study}_{sample}_b_zarr_manifest.tsv"
    params:
        study=lambda wildcards: wildcards.study,
        sample=lambda wildcards: wildcards.sample,
        cell_scope="B_cells"
    script:
        "scripts/h5ad_to_zarr.py"


########################################
# ChangeO AssignGenes
########################################
rule assign_genes:
    input:
        fasta=lambda wildcards: RAW_VDJ_FASTA[(wildcards.study, wildcards.sample)]
    output:
        "results/changeo/{study}_{sample}_igblast.tsv"
    params:
        igblast_db=config["igblast_db"],
    threads: 20
    shell:
        """
        AssignGenes.py igblast \
            -s {input.fasta} \
            -b {params.igblast_db} \
            --organism human \
            --loci ig \
            --format blast \
            --nproc {threads} \
            -o {output}
        """

########################################
# MakeDb
########################################
rule make_db:
    input:
        igblast="results/changeo/{study}_{sample}_igblast.tsv",
        fasta=lambda wildcards: RAW_VDJ_FASTA[(wildcards.study, wildcards.sample)]
    output:
        "results/changeo/{study}_{sample}_db.tsv"
    params:
        vdj_db=config["vdj_db"]
    shell:
        """
        MakeDb.py igblast \
            -i {input.igblast} \
            -s {input.fasta} \
            -r {params.vdj_db} \
            --extended \
            -o {output}
        """
########################################
# Create Germlines
########################################
rule create_germlines:
    input:
        "results/changeo/{study}_{sample}_db.tsv"
    output:
        "results/changeo/{study}_{sample}_db_germ.tsv"
    params:
        vdj_db=config["vdj_db"]
    shell:
        """
        CreateGermlines.py \
        -d {input} \
	    -g dmask \
        -r {params.vdj_db} \
        -o {output}
        """

########################################
# SHM
########################################
rule shm:
    input:
        "results/changeo/{study}_{sample}_db_germ.tsv"
    output:
        "results/changeo/{study}_{sample}_SHM.tsv"
    shell:
        """
        Rscript scripts/shazam.r \
        {input} {output}
        """


########################################
# Fill c_call from Cell Ranger VDJ annotations
########################################
rule fill_c_call:
    input:
        shm="results/changeo/{study}_{sample}_SHM.tsv",
        annotations=lambda wildcards: RAW_VDJ_ANNOTATIONS[(wildcards.study, wildcards.sample)]
    output:
        "results/changeo/{study}_{sample}_SHM_c_call.tsv"
    shell:
        """
        python scripts/fill_c_call_from_cellranger.py \
            --input {input.shm} \
            --annotations {input.annotations} \
            --output {output}
        """

########################################
# BCR pairing
########################################
rule bcr_pair:
    input:
        "results/changeo/{study}_{sample}_SHM_c_call.tsv"
    output:
        bcr_paired="results/bcr/{study}_{sample}_SHM_paired.tsv",
        bcr_IGH="results/bcr/{study}_{sample}_SHM_IGH.tsv"
    script:
        "scripts/bcr_pairing.py"


########################################
# tsv to json
########################################
rule tsv_to_json:
    input:
        paired="results/bcr/{study}_{sample}_SHM_paired.tsv",
        heavy="results/bcr/{study}_{sample}_SHM_IGH.tsv",
        clone="results/bcr/{study}_{sample}_clone.tsv"
    output:
        "results/bcr/{study}_{sample}_bcr.json"
    shell:
        """
        python scripts/export_bcr_json.py \
            --paired {input.paired} \
            --heavy {input.heavy} \
            --clone {input.clone} \
            --sample-id {wildcards.sample} \
            --output {output}
        """

########################################
# DefineClones
########################################
rule defineclones:
    input:
        "results/bcr/{study}_{sample}_SHM_IGH.tsv"
    output:
        clone_tsv="results/bcr/{study}_{sample}_clone.tsv"
    script:
        "scripts/define_clones_safe.py"

########################################
# lineage tree assets (BuildTrees + JSON)
########################################
rule lineage_tree_assets:
    input:
        clone_tsv="results/bcr/{study}_{sample}_clone.tsv"
    output:
        edge_tsv="results/lineage/{study}_{sample}_lineage_edge.tsv",
        tree_tsv="results/lineage/{study}_{sample}_lineage_tree.tsv",
        lineage_json="results/lineage/{study}_{sample}_lineage.json"
    params:
        igphyml_arg=(lambda wildcards: f"--igphyml-path {IGPHYML_PATH}" if IGPHYML_PATH else "")
    shell:
        """
        python scripts/build_lineage_assets.py \
            --input {input.clone_tsv} \
            --study {wildcards.study} \
            --sample {wildcards.sample} \
            --edge-out {output.edge_tsv} \
            --tree-out {output.tree_tsv} \
            --json-out {output.lineage_json} \
            --buildtrees-dir results/lineage/{wildcards.study}_{wildcards.sample}_buildtrees \
            --min-clone-size 2 \
            --nproc 1 \
            {params.igphyml_arg}
        """

########################################
# VDJ usage
########################################

rule vdj_usage:
    input:
        "results/bcr/{study}_{sample}_clone.tsv"
    output:
        "results/repertoire/{study}_{sample}_vdj_usage.tsv"
    script:
        "scripts/vdj_usage.py"

########################################
# diversity
########################################
rule diversity:
    input:
        expand("results/bcr/{study}_{sample}_clone.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
    output:
        "results/repertoire/diversity.tsv"
    script:
        "scripts/diversity.py"

########################################
# cdr3_network
########################################
rule cdr3_network:
    input:
        "results/bcr/{study}_{sample}_clone.tsv"
    output:
        edges="results/cdr3_network/{study}_{sample}_cdr3_edges.tsv",
        nodes="results/cdr3_network/{study}_{sample}_cdr3_nodes.tsv",
        json="results/cdr3_network/{study}_{sample}_cdr3_network.json",
        stats="results/cdr3_network/{study}_{sample}_cdr3_network_stats.json"
    shell:
        """
        python scripts/cdr3_network.py \
            --input {input} \
            --edges {output.edges} \
            --nodes {output.nodes} \
            --json {output.json} \
            --max-dist 5 \
            --stats {output.stats}
        """

########################################
# public clonotype
########################################
rule public_clonotype:
    input:
        expand("results/bcr/{study}_{sample}_clone.tsv",
               zip, study=STUDIES, sample=SAMPLE_NAMES)
    output:
        public="results/repertoire/public_clonotypes.tsv",
        network="results/repertoire/public_network.tsv",
        members="results/repertoire/clone_members.tsv"
    script:
        "scripts/public_clonotype.py"

########################################
# load_database
########################################
rule load_database:
    input:
        study="metadata/study.tsv",
        subject="metadata/subject.tsv",
        sample="metadata/sample.tsv",
        clones=expand(
            "results/bcr/{study}_{sample}_clone.tsv",
            zip,
            study=STUDIES,
            sample=SAMPLE_NAMES
        ),
        bcr_paired=expand(
            "results/bcr/{study}_{sample}_SHM_paired.tsv",
            zip,
            study=STUDIES,
            sample=SAMPLE_NAMES
        ),
        vdj_usage=expand(
            "results/repertoire/{study}_{sample}_vdj_usage.tsv",
            zip,
            study=STUDIES,
            sample=SAMPLE_NAMES
        ),
        cell_embeddings=expand(
            "results/gex/{study}_{sample}_b_cell_embedding.tsv",
            zip,
            study=STUDIES,
            sample=SAMPLE_NAMES
        ),
        marker_expression=expand(
            "results/gex/{study}_{sample}_b_marker_expression.tsv",
            zip,
            study=STUDIES,
            sample=SAMPLE_NAMES
        ),
        clone_gene_expression=expand(
            "results/gex_de/{study}_{sample}_clone_gene_expression_summary.tsv",
            zip,
            study=STUDIES,
            sample=SAMPLE_NAMES
        ),
        expression_stores=expand(
            "results/gex_zarr/{study}_{sample}_b_zarr_manifest.tsv",
            zip,
            study=STUDIES,
            sample=SAMPLE_NAMES
        ),
        lineage_edges=expand(
            "results/lineage/{study}_{sample}_lineage_edge.tsv",
            zip,
            study=STUDIES,
            sample=SAMPLE_NAMES
        ),
        lineage_trees=expand(
            "results/lineage/{study}_{sample}_lineage_tree.tsv",
            zip,
            study=STUDIES,
            sample=SAMPLE_NAMES
        ),
        member="results/repertoire/clone_members.tsv",
        public="results/repertoire/public_clonotypes.tsv",
        network="results/repertoire/public_network.tsv",
        diversity = "results/repertoire/diversity.tsv"
    output:
        "results/database_loaded.txt"
    run:
        # 将文件列表转换为空格分隔的字符串
        clone_files = " ".join(input.clones)
        paired_files = " ".join(input.bcr_paired)
        vdj_files = " ".join(input.vdj_usage)
        cell_embedding_files = " ".join(input.cell_embeddings)
        marker_expression_files = " ".join(input.marker_expression)
        clone_gene_expression_files = " ".join(input.clone_gene_expression)
        expression_store_files = " ".join(input.expression_stores)
        lineage_edge_files = " ".join(input.lineage_edges)
        lineage_tree_files = " ".join(input.lineage_trees)
        
        shell("""
            python scripts/load_database.py \
              --study {input.study} \
              --subject {input.subject} \
                --sample {input.sample} \
                --clones {clone_files} \
                --paired-bcr {paired_files} \
                --vdj-usage {vdj_files} \
                --cell-embeddings {cell_embedding_files} \
                --marker-expression {marker_expression_files} \
                --clone-gene-expression {clone_gene_expression_files} \
                --expression-stores {expression_store_files} \
                --lineage-edges {lineage_edge_files} \
                --lineage-trees {lineage_tree_files} \
                --member {input.member} \
                --public {input.public} \
                --network {input.network} \
                --diversity {input.diversity} \
                --reload
            
            touch {output}
        """)

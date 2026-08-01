#!/usr/bin/env python3
"""
BCR数据库加载脚本
将分析结果导入PostgreSQL数据库
支持按样本导入克隆数据和V(D)J使用数据
"""

import psycopg2
from psycopg2 import sql
from pathlib import Path
import sys
import argparse
import logging
from datetime import datetime
import yaml
import os
import csv
import tempfile
import re

# 配置日志
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# 读取配置 - 使用相对于脚本位置的路径
script_dir = Path(__file__).parent
config_path = script_dir.parent / 'config.yaml'

if not config_path.exists():
    logger.error(f"配置文件不存在: {config_path}")
    sys.exit(1)

with open(config_path, 'r') as f:
    config = yaml.safe_load(f)

############################################
# 数据库连接
############################################

def get_db_connection():
    """获取数据库连接，带重试机制"""
    try:
        conn = psycopg2.connect(
            dbname=config['pg']['db'],
            user=config['pg']['user'],
            password=config['pg']['password'],
            host=config['pg']['host'],
            connect_timeout=10
        )
        return conn
    except psycopg2.Error as e:
        logger.error(f"数据库连接失败: {e}")
        sys.exit(1)


############################################
# 辅助函数
############################################


def table_exists(conn, table_name):
    """检查数据库表是否存在"""
    with conn.cursor() as cur:
        cur.execute("SELECT to_regclass(%s)", (table_name,))
        return cur.fetchone()[0] is not None


DISEASE_NORMALIZATION = {
    "follicular_lymphoma": "Lymphoma",
    "follicular lymphoma": "Lymphoma",
    "primary central nervous system lymphoma": "Lymphoma",
    "b-cell acute lymphoblastic leukaemia": "B-cell acute lymphoblastic leukemia",
    "acute lymphocytic leukemia": "B-cell acute lymphoblastic leukemia",
    "high-grade serous ovarian cancer": "Ovarian Cancer",
    "pancreatic ductal adenocarcinoma": "Pancreatic Cancer",
    "intrahepatic cholangiocarcinoma": "Liver Cancer",
    "people living with hiv": "HIV infection",
    "kawasaki disease": "Kawasaki Disease",
    "solid organ transplant recipients": "Solid Organ Transplant Recipient",
    "atrip deficient": "ATRIP deficiency",
    "parkinson disease": "Parkinson's disease",
    "parkinson's disease": "Parkinson's disease",
    "parkinson’s disease": "Parkinson's disease",
}


def normalize_disease_column(cur, temp_table):
    """Normalize disease labels in a temporary metadata import table."""
    for source, target in DISEASE_NORMALIZATION.items():
        cur.execute(
            sql.SQL("UPDATE {} SET disease = %s WHERE lower(trim(disease)) = %s").format(
                sql.Identifier(temp_table)
            ),
            (target, source),
        )


def ensure_sample_paired_bcr_column(conn):
    """确保 sample 表包含 paired_bcr 标记列。"""
    if not table_exists(conn, "sample"):
        return
    try:
        with conn.cursor() as cur:
            cur.execute("""
                ALTER TABLE sample
                ADD COLUMN IF NOT EXISTS paired_bcr BOOLEAN DEFAULT TRUE
            """)
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"确保 sample.paired_bcr 失败: {e}")
        raise

def ensure_bcr_cdr3_aa_column(conn):
    """确保 bcr_sequences 中存在 cdr3_aa 列及检索索引。"""
    if not table_exists(conn, "bcr_sequences"):
        return
    try:
        with conn.cursor() as cur:
            cur.execute("""
                ALTER TABLE bcr_sequences
                ADD COLUMN IF NOT EXISTS cdr3_aa TEXT
            """)
            cur.execute("""
                UPDATE bcr_sequences
                SET cdr3_aa = COALESCE(NULLIF(cdr3_aa, ''), NULLIF(junction_aa, ''))
                WHERE cdr3_aa IS NULL OR cdr3_aa = ''
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_bcr_sequences_cdr3_aa
                ON bcr_sequences(cdr3_aa)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_bcr_sequences_cdr3_aa_upper
                ON bcr_sequences(UPPER(COALESCE(NULLIF(cdr3_aa, ''), NULLIF(junction_aa, ''))))
                WHERE COALESCE(NULLIF(cdr3_aa, ''), NULLIF(junction_aa, '')) IS NOT NULL
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_bcr_sequences_cdr3_aa_length
                ON bcr_sequences(length(UPPER(COALESCE(NULLIF(cdr3_aa, ''), NULLIF(junction_aa, '')))))
                WHERE COALESCE(NULLIF(cdr3_aa, ''), NULLIF(junction_aa, '')) IS NOT NULL
                  AND COALESCE(cell_barcode, '') <> ''
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_bcr_sequences_cdr3_search_cell
                ON bcr_sequences(sample_id, cell_barcode, locus)
                WHERE COALESCE(cell_barcode, '') <> ''
            """)
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"确保 bcr_sequences.cdr3_aa 失败: {e}")
        raise


def ensure_analysis_indexes(conn):
    """创建Analysis页面高频聚合查询所需的组合索引。"""
    index_statements = [
        """
        CREATE INDEX IF NOT EXISTS idx_sample_study_sample
        ON sample(study, sample_id)
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_sample_disease_sample
        ON sample((COALESCE(NULLIF(disease, ''), 'Unknown')), sample_id)
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_sample_subject_sample
        ON sample(subject_id, sample_id)
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_clone_barcode
        ON bcr_sequences(sample_id, locus, clone_id, cell_barcode)
        WHERE COALESCE(clone_id, '') <> '' AND COALESCE(cell_barcode, '') <> ''
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_cdr3
        ON bcr_sequences(sample_id, locus, cdr3_aa)
        WHERE COALESCE(cdr3_aa, '') <> ''
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_junction_aa
        ON bcr_sequences(sample_id, locus, junction_aa)
        WHERE COALESCE(junction_aa, '') <> ''
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_mu
        ON bcr_sequences(sample_id, locus, mu_freq)
        WHERE mu_freq IS NOT NULL
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_c_call
        ON bcr_sequences(sample_id, locus, c_call)
        """,
        """
        CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_cell
        ON bcr_sequences(sample_id, locus, cell_barcode)
        WHERE COALESCE(cell_barcode, '') <> ''
        """,
    ]
    try:
        with conn.cursor() as cur:
            for statement in index_statements:
                cur.execute(statement)
            if table_exists(conn, "gex_cell_embedding"):
                cur.execute("""
                    CREATE INDEX IF NOT EXISTS idx_gex_embedding_adv_sample_scope_key_barcode
                    ON gex_cell_embedding(sample_id, cell_scope, embedding_key, cell_barcode, cell_type)
                """)
            cur.execute("""
                CREATE TABLE IF NOT EXISTS bcr_cdr3_vocab (
                    cdr3_key TEXT PRIMARY KEY,
                    cdr3_length INTEGER NOT NULL
                )
            """)
            cur.execute("TRUNCATE bcr_cdr3_vocab")
            cur.execute("""
                INSERT INTO bcr_cdr3_vocab (cdr3_key, cdr3_length)
                SELECT cdr3_key, length(cdr3_key) AS cdr3_length
                FROM (
                    SELECT DISTINCT UPPER(COALESCE(NULLIF(cdr3_aa, ''), NULLIF(junction_aa, ''))) AS cdr3_key
                    FROM bcr_sequences
                    WHERE COALESCE(NULLIF(cdr3_aa, ''), NULLIF(junction_aa, '')) IS NOT NULL
                ) AS vocab
                WHERE cdr3_key <> ''
                ON CONFLICT (cdr3_key) DO UPDATE SET cdr3_length = EXCLUDED.cdr3_length
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_bcr_cdr3_vocab_length
                ON bcr_cdr3_vocab(cdr3_length)
            """)
            cur.execute("CREATE EXTENSION IF NOT EXISTS pg_trgm")
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_bcr_cdr3_vocab_key_length
                ON bcr_cdr3_vocab(cdr3_key, cdr3_length)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_bcr_cdr3_vocab_trgm
                ON bcr_cdr3_vocab USING GIN (cdr3_key gin_trgm_ops)
            """)
            cur.execute("ANALYZE bcr_cdr3_vocab")
            cur.execute("""
                CREATE TABLE IF NOT EXISTS bcr_cdr3_cell_lookup (
                    cdr3_key TEXT NOT NULL,
                    sample_id TEXT NOT NULL,
                    cell_barcode TEXT NOT NULL,
                    PRIMARY KEY (cdr3_key, sample_id, cell_barcode)
                )
            """)
            cur.execute("TRUNCATE bcr_cdr3_cell_lookup")
            cur.execute("""
                INSERT INTO bcr_cdr3_cell_lookup (cdr3_key, sample_id, cell_barcode)
                SELECT DISTINCT
                    UPPER(COALESCE(NULLIF(cdr3_aa, ''), NULLIF(junction_aa, ''))) AS cdr3_key,
                    sample_id,
                    cell_barcode
                FROM bcr_sequences
                WHERE COALESCE(NULLIF(cdr3_aa, ''), NULLIF(junction_aa, '')) IS NOT NULL
                  AND COALESCE(sample_id, '') <> ''
                  AND COALESCE(cell_barcode, '') <> ''
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_bcr_cdr3_cell_lookup_sample_cell
                ON bcr_cdr3_cell_lookup(sample_id, cell_barcode)
            """)
            cur.execute("ANALYZE bcr_cdr3_cell_lookup")
            cur.execute("ANALYZE sample")
            cur.execute("ANALYZE bcr_sequences")
            if table_exists(conn, "gex_cell_embedding"):
                cur.execute("ANALYZE gex_cell_embedding")
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"创建Analysis组合索引失败: {e}")
        raise


def ensure_gex_cell_embedding_table(conn):
    """创建前端UMAP展示所需的单细胞embedding表"""
    try:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS gex_cell_embedding (
                    study TEXT REFERENCES study(study) ON DELETE CASCADE,
                    sample_id TEXT REFERENCES sample(sample_id) ON DELETE CASCADE,
                    cell_barcode TEXT NOT NULL,
                    cell_scope TEXT NOT NULL DEFAULT 'B_cells',
                    embedding_key TEXT NOT NULL DEFAULT 'X_umap',
                    umap_1 DOUBLE PRECISION NOT NULL,
                    umap_2 DOUBLE PRECISION NOT NULL,
                    cell_type TEXT,
                    cell_subtype TEXT,
                    cluster TEXT,
                    PRIMARY KEY (sample_id, cell_scope, embedding_key, cell_barcode)
                )
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_cell_embedding_study
                ON gex_cell_embedding(study)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_cell_embedding_sample
                ON gex_cell_embedding(sample_id)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_cell_embedding_cell_type
                ON gex_cell_embedding(cell_type)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_cell_embedding_cell_subtype
                ON gex_cell_embedding(cell_subtype)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_cell_embedding_cluster
                ON gex_cell_embedding(sample_id, cluster)
            """)
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"创建 gex_cell_embedding 表失败: {e}")
        raise


def ensure_gex_marker_expression_table(conn):
    """创建高频featureplot查询所需的B细胞marker表达缓存表"""
    try:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS gex_marker_expression (
                    study TEXT REFERENCES study(study) ON DELETE CASCADE,
                    sample_id TEXT REFERENCES sample(sample_id) ON DELETE CASCADE,
                    cell_barcode TEXT NOT NULL,
                    cell_scope TEXT NOT NULL DEFAULT 'B_cells',
                    gene TEXT NOT NULL,
                    expression_layer TEXT NOT NULL DEFAULT 'raw_log1p',
                    expression DOUBLE PRECISION NOT NULL,
                    detected BOOLEAN NOT NULL DEFAULT false,
                    PRIMARY KEY (sample_id, cell_scope, expression_layer, gene, cell_barcode)
                )
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_marker_expression_study
                ON gex_marker_expression(study)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_marker_expression_sample_gene
                ON gex_marker_expression(sample_id, gene)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_marker_expression_gene
                ON gex_marker_expression(gene)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_marker_expression_detected
                ON gex_marker_expression(sample_id, gene, detected)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_marker_sample_scope_barcode_gene
                ON gex_marker_expression(sample_id, cell_scope, cell_barcode, gene)
                INCLUDE (expression, detected)
            """)
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"创建 gex_marker_expression 表失败: {e}")
        raise


def ensure_clone_gene_expression_summary_table(conn):
    """创建clone-level表达摘要表，用于BCR-transcriptome快速差异分析。"""
    try:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS clone_gene_expression_summary (
                    study TEXT REFERENCES study(study) ON DELETE CASCADE,
                    sample_id TEXT REFERENCES sample(sample_id) ON DELETE CASCADE,
                    clone_id TEXT NOT NULL,
                    gene TEXT NOT NULL,
                    expression_layer TEXT NOT NULL DEFAULT 'raw_log1p',
                    n_cells INT NOT NULL,
                    mean_expression DOUBLE PRECISION NOT NULL,
                    variance_expression DOUBLE PRECISION NOT NULL DEFAULT 0,
                    detected_cells INT NOT NULL DEFAULT 0,
                    detection_fraction DOUBLE PRECISION NOT NULL DEFAULT 0,
                    PRIMARY KEY (sample_id, expression_layer, clone_id, gene)
                )
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_clone_gene_summary_study
                ON clone_gene_expression_summary(study)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_clone_gene_summary_sample_clone
                ON clone_gene_expression_summary(sample_id, clone_id)
                INCLUDE (n_cells)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_clone_gene_summary_sample_gene
                ON clone_gene_expression_summary(sample_id, gene)
            """)
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"创建 clone_gene_expression_summary 表失败: {e}")
        raise


def ensure_gex_expression_store_table(conn):
    """创建完整表达矩阵zarr存储位置表"""
    try:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS gex_expression_store (
                    study TEXT REFERENCES study(study) ON DELETE CASCADE,
                    sample_id TEXT REFERENCES sample(sample_id) ON DELETE CASCADE,
                    cell_scope TEXT NOT NULL DEFAULT 'B_cells',
                    store_format TEXT NOT NULL DEFAULT 'zarr',
                    zarr_path TEXT NOT NULL,
                    expression_layer TEXT NOT NULL DEFAULT 'raw_log1p',
                    n_cells INT,
                    n_genes INT,
                    source_h5ad_path TEXT,
                    created_at TIMESTAMP,
                    PRIMARY KEY (sample_id, cell_scope, store_format)
                )
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_expression_store_study
                ON gex_expression_store(study)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_gex_expression_store_sample
                ON gex_expression_store(sample_id)
            """)
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"创建 gex_expression_store 表失败: {e}")
        raise


def ensure_lineage_edge_table(conn):
    """创建clone lineage边表"""
    try:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS lineage_edge (
                    study TEXT REFERENCES study(study) ON DELETE CASCADE,
                    sample_id TEXT REFERENCES sample(sample_id) ON DELETE CASCADE,
                    clone_id TEXT NOT NULL,
                    parent_node TEXT NOT NULL,
                    child_node TEXT NOT NULL,
                    parent_type TEXT,
                    child_type TEXT,
                    distance INT,
                    edge_weight DOUBLE PRECISION,
                    is_germline_edge BOOLEAN DEFAULT false,
                    build_method TEXT,
                    PRIMARY KEY (sample_id, clone_id, parent_node, child_node)
                )
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_lineage_edge_sample
                ON lineage_edge(sample_id)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_lineage_edge_clone
                ON lineage_edge(sample_id, clone_id)
            """)
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"创建 lineage_edge 表失败: {e}")
        raise


def ensure_lineage_tree_table(conn):
    """创建clone lineage JSON表"""
    try:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS lineage_tree_json (
                    study TEXT REFERENCES study(study) ON DELETE CASCADE,
                    sample_id TEXT REFERENCES sample(sample_id) ON DELETE CASCADE,
                    clone_id TEXT NOT NULL,
                    n_nodes INT,
                    n_edges INT,
                    build_method TEXT,
                    tree_json JSONB NOT NULL,
                    PRIMARY KEY (sample_id, clone_id)
                )
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_lineage_tree_sample
                ON lineage_tree_json(sample_id)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_lineage_tree_clone
                ON lineage_tree_json(sample_id, clone_id)
            """)
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"创建 lineage_tree_json 表失败: {e}")
        raise


def ensure_download_file_table(conn):
    """创建前端Download页面使用的Zenodo/本地文件索引表。"""
    try:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS download_file (
                    study TEXT REFERENCES study(study) ON DELETE CASCADE,
                    sample_id TEXT REFERENCES sample(sample_id) ON DELETE CASCADE,
                    file_kind TEXT NOT NULL,
                    file_name TEXT NOT NULL,
                    upload_folder TEXT,
                    zenodo_part TEXT,
                    part_index INT,
                    part_total INT,
                    relative_upload_path TEXT,
                    source_path TEXT,
                    file_size_bytes BIGINT,
                    md5 TEXT,
                    zenodo_record_id TEXT,
                    zenodo_doi TEXT,
                    file_url TEXT,
                    created_at TIMESTAMPTZ NOT NULL DEFAULT now(),
                    updated_at TIMESTAMPTZ NOT NULL DEFAULT now(),
                    PRIMARY KEY (sample_id, file_kind)
                )
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_download_file_study
                ON download_file(study)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_download_file_record
                ON download_file(zenodo_record_id)
            """)
            cur.execute("""
                CREATE INDEX IF NOT EXISTS idx_download_file_kind
                ON download_file(file_kind)
            """)
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"创建 download_file 表失败: {e}")
        raise


def ensure_study_accession_column(conn):
    """Add the study accession field for existing databases."""
    try:
        with conn.cursor() as cur:
            cur.execute("ALTER TABLE study ADD COLUMN IF NOT EXISTS accession TEXT")
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"添加 study.accession 字段失败: {e}")
        raise


def ensure_frontend_required_tables(conn):
    """确保前端所有可选模块依赖表存在。

    这些表可能没有在某次入库命令中提供对应结果文件，但前端API仍会检查或查询它们。
    统一创建空表可以避免数据库更新后出现缺表错误。
    """
    ensure_gex_cell_embedding_table(conn)
    ensure_gex_marker_expression_table(conn)
    ensure_clone_gene_expression_summary_table(conn)
    ensure_gex_expression_store_table(conn)
    ensure_lineage_edge_table(conn)
    ensure_lineage_tree_table(conn)
    ensure_download_file_table(conn)


def copy_table(conn, file_path, table_name, sample_id=None):
    """使用COPY命令导入数据，带错误处理"""
    if not Path(file_path).exists():
        logger.warning(f"文件不存在，跳过: {file_path}")
        return False

    logger.info(f"正在导入 {file_path} -> {table_name}" +
                (f" (样本: {sample_id})" if sample_id else ""))

    try:
        with conn.cursor() as cur:
            # 如果指定了sample_id，先删除该样本的旧数据（可选）
            if sample_id and table_name in [
                'bcr_sequences',
                'vdj_usage',
                'gex_cell_embedding',
                'gex_marker_expression',
                'clone_gene_expression_summary',
                'gex_expression_store',
                'lineage_edge',
                'lineage_tree_json'
            ]:
                cur.execute(
                    sql.SQL("DELETE FROM {} WHERE sample_id = %s").format(
                        sql.Identifier(table_name)
                    ),
                    (sample_id,)
                )
                logger.info(f"已清除样本 {sample_id} 在 {table_name} 中的旧数据")

            # 创建临时表
            temp_table = f"temp_{table_name}_{datetime.now().strftime('%Y%m%d%H%M%S')}"
            cur.execute(sql.SQL("CREATE TEMPORARY TABLE {} (LIKE {} INCLUDING DEFAULTS) ON COMMIT DROP").format(
                sql.Identifier(temp_table), sql.Identifier(table_name)
            ))

            # 使用COPY命令
            with open(file_path, 'r', encoding='utf-8') as f:
                if table_name == 'bcr_sequences':
                    header_line = f.readline().strip()
                    if not header_line:
                        raise ValueError(f"{file_path} 文件表头为空")
                    input_columns = [col.strip() for col in header_line.split('\t') if col.strip()]
                    if "sequence_id" not in input_columns or "sample_id" not in input_columns:
                        raise ValueError(f"{file_path} 缺少 sequence_id/sample_id 列，无法导入 bcr_sequences")
                    f.seek(0)
                    copy_sql = sql.SQL("""
                        COPY {} ({}) FROM STDIN
                        WITH CSV HEADER DELIMITER E'\\t'
                        NULL ''
                        QUOTE E'\\b'
                    """).format(
                        sql.Identifier(temp_table),
                        sql.SQL(", ").join(sql.Identifier(col) for col in input_columns)
                    )
                elif table_name == 'clonotype_member':
                    # clone_members.tsv 的列顺序是 clonotype_id, sequence_id, sample_id, cdr3, v_gene, j_gene
                    # 显式指定目标列名，避免 COPY 按临时表定义顺序错误映射。
                    copy_sql = sql.SQL("""
                        COPY {} (clonotype_key, sequence_id, sample_id, cdr3_aa, v_call, j_call) FROM STDIN
                        WITH CSV HEADER DELIMITER E'\\t'
                        NULL ''
                        QUOTE E'\\b'
                    """).format(sql.Identifier(temp_table))
                elif table_name == 'sample':
                    header_line = f.readline().strip()
                    if not header_line:
                        raise ValueError(f"{file_path} 文件表头为空")
                    input_columns = [col.strip() for col in header_line.split('	') if col.strip()]
                    normalized_columns = [
                        'paired_bcr' if col.lower() == 'paired_bcr' else col
                        for col in input_columns
                    ]
                    f.seek(0)
                    copy_sql = sql.SQL("""
                        COPY {} ({}) FROM STDIN
                        WITH CSV HEADER DELIMITER E'\t'
                        NULL ''
                        QUOTE E'\b'
                    """).format(
                        sql.Identifier(temp_table),
                        sql.SQL(", ").join(sql.Identifier(col) for col in normalized_columns)
                    )
                elif table_name == 'study':
                    copy_sql = sql.SQL("""
                        COPY {} (study, title, year, journal, accession) FROM STDIN
                        WITH CSV HEADER DELIMITER E'\\t'
                        NULL ''
                        QUOTE E'\\b'
                    """).format(sql.Identifier(temp_table))
                else:
                    copy_sql = sql.SQL("""
                        COPY {} FROM STDIN
                        WITH CSV HEADER DELIMITER E'\\t'
                        NULL ''
                        QUOTE E'\\b'
                    """).format(sql.Identifier(temp_table))

                cur.copy_expert(copy_sql, f)

            # 获取处理的行数
            cur.execute(sql.SQL("SELECT COUNT(*) FROM {}").format(sql.Identifier(temp_table)))
            row_count = cur.fetchone()[0]

            # 根据表类型进行UPSERT
            if table_name == 'study':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (study, title, year, journal, accession)
                    SELECT study, title, year, journal, accession FROM {}
                    ON CONFLICT (study) DO UPDATE SET
                        title = EXCLUDED.title,
                        year = EXCLUDED.year,
                        journal = EXCLUDED.journal,
                        accession = EXCLUDED.accession
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'subject':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (subject_id, species, disease)
                    SELECT subject_id, species, disease FROM {}
                    ON CONFLICT (subject_id) DO UPDATE SET
                        species = EXCLUDED.species,
                        disease = EXCLUDED.disease
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'sample':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (sample_id, study, subject_id, tissue, disease, platform, n_cells, paired_bcr)
                    SELECT sample_id, study, subject_id, tissue, disease, platform, n_cells, COALESCE(paired_bcr, TRUE) FROM {}
                    ON CONFLICT (sample_id) DO UPDATE SET
                        study = EXCLUDED.study,
                        subject_id = EXCLUDED.subject_id,
                        tissue = EXCLUDED.tissue,
                        disease = EXCLUDED.disease,
                        platform = EXCLUDED.platform,
                        n_cells = EXCLUDED.n_cells,
                        paired_bcr = EXCLUDED.paired_bcr
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'bcr_sequences':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (
                        sequence_id, sequence, rev_comp, productive, v_call, d_call, j_call,
                        sequence_alignment, germline_alignment, junction, junction_aa,
                        v_cigar, d_cigar, j_cigar, stop_codon, vj_in_frame, locus, c_call,
                        junction_length, np1_length, np2_length,
                        v_sequence_start, v_sequence_end, v_germline_start, v_germline_end,
                        d_sequence_start, d_sequence_end, d_germline_start, d_germline_end,
                        j_sequence_start, j_sequence_end, j_germline_start, j_germline_end,
                        v_score, v_identity, v_support, d_score, d_identity, d_support,
                        j_score, j_identity, j_support, fwr1, fwr2, fwr3, fwr4,
                        cdr1, cdr2, cdr3, cdr3_aa, clone_id, germline_alignment_d_mask, mu_freq,
                        cell_barcode, sample_id
                    )
                    SELECT
                        sequence_id, sequence, rev_comp, productive, v_call, d_call, j_call,
                        sequence_alignment, germline_alignment, junction, junction_aa,
                        v_cigar, d_cigar, j_cigar, stop_codon, vj_in_frame, locus, c_call,
                        junction_length, np1_length, np2_length,
                        v_sequence_start, v_sequence_end, v_germline_start, v_germline_end,
                        d_sequence_start, d_sequence_end, d_germline_start, d_germline_end,
                        j_sequence_start, j_sequence_end, j_germline_start, j_germline_end,
                        v_score, v_identity, v_support, d_score, d_identity, d_support,
                        j_score, j_identity, j_support, fwr1, fwr2, fwr3, fwr4,
                        cdr1, cdr2, cdr3,
                        COALESCE(NULLIF(cdr3_aa, ''), NULLIF(junction_aa, '')) AS cdr3_aa,
                        clone_id, germline_alignment_d_mask, mu_freq,
                        cell_barcode, sample_id
                    FROM {}
                    ON CONFLICT (sequence_id, sample_id) DO UPDATE SET
                        sequence = EXCLUDED.sequence,
                        rev_comp = EXCLUDED.rev_comp,
                        productive = EXCLUDED.productive,
                        v_call = EXCLUDED.v_call,
                        d_call = EXCLUDED.d_call,
                        j_call = EXCLUDED.j_call,
                        sequence_alignment = EXCLUDED.sequence_alignment,
                        germline_alignment = EXCLUDED.germline_alignment,
                        junction = EXCLUDED.junction,
                        junction_aa = EXCLUDED.junction_aa,
                        v_cigar = EXCLUDED.v_cigar,
                        d_cigar = EXCLUDED.d_cigar,
                        j_cigar = EXCLUDED.j_cigar,
                        stop_codon = EXCLUDED.stop_codon,
                        vj_in_frame = EXCLUDED.vj_in_frame,
                        locus = EXCLUDED.locus,
                        c_call = EXCLUDED.c_call,
                        junction_length = EXCLUDED.junction_length,
                        np1_length = EXCLUDED.np1_length,
                        np2_length = EXCLUDED.np2_length,
                        v_sequence_start = EXCLUDED.v_sequence_start,
                        v_sequence_end = EXCLUDED.v_sequence_end,
                        v_germline_start = EXCLUDED.v_germline_start,
                        v_germline_end = EXCLUDED.v_germline_end,
                        d_sequence_start = EXCLUDED.d_sequence_start,
                        d_sequence_end = EXCLUDED.d_sequence_end,
                        d_germline_start = EXCLUDED.d_germline_start,
                        d_germline_end = EXCLUDED.d_germline_end,
                        j_sequence_start = EXCLUDED.j_sequence_start,
                        j_sequence_end = EXCLUDED.j_sequence_end,
                        j_germline_start = EXCLUDED.j_germline_start,
                        j_germline_end = EXCLUDED.j_germline_end,
                        v_score = EXCLUDED.v_score,
                        v_identity = EXCLUDED.v_identity,
                        v_support = EXCLUDED.v_support,
                        d_score = EXCLUDED.d_score,
                        d_identity = EXCLUDED.d_identity,
                        d_support = EXCLUDED.d_support,
                        j_score = EXCLUDED.j_score,
                        j_identity = EXCLUDED.j_identity,
                        j_support = EXCLUDED.j_support,
                        fwr1 = EXCLUDED.fwr1,
                        fwr2 = EXCLUDED.fwr2,
                        fwr3 = EXCLUDED.fwr3,
                        fwr4 = EXCLUDED.fwr4,
                        cdr1 = EXCLUDED.cdr1,
                        cdr2 = EXCLUDED.cdr2,
                        cdr3 = EXCLUDED.cdr3,
                        cdr3_aa = EXCLUDED.cdr3_aa,
                        clone_id = EXCLUDED.clone_id,
                        germline_alignment_d_mask = EXCLUDED.germline_alignment_d_mask,
                        mu_freq = EXCLUDED.mu_freq,
                        cell_barcode = EXCLUDED.cell_barcode
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'clonotype_member':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (clonotype_key, sample_id, sequence_id, cdr3_aa, v_call, j_call)
                    SELECT clonotype_key, sample_id, sequence_id, cdr3_aa, v_call, j_call FROM {}
                    ON CONFLICT (clonotype_key, sample_id, sequence_id) DO UPDATE SET
                        cdr3_aa = EXCLUDED.cdr3_aa,
                        v_call = EXCLUDED.v_call,
                        j_call = EXCLUDED.j_call
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'public_clonotype':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (clonotype_key, n_cells, n_samples, samples)
                    SELECT clonotype_key, n_cells, n_samples, samples FROM {}
                    ON CONFLICT (clonotype_key) DO UPDATE SET
                        n_cells = EXCLUDED.n_cells,
                        n_samples = EXCLUDED.n_samples,
                        samples = EXCLUDED.samples
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'public_network':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (sample1, sample2, clonotype_key)
                    SELECT sample1, sample2, clonotype_key FROM {}
                    ON CONFLICT DO NOTHING
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'vdj_usage':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (study, sample_id, vdj_type, vdj_gene, count, frequency)
                    SELECT
                        COALESCE((SELECT s.study FROM "sample" s WHERE s.sample_id = t.sample_id), t.study) AS study,
                        t.sample_id, t.vdj_type, t.vdj_gene, t.count, t.frequency
                    FROM {} t
                    ON CONFLICT (sample_id, vdj_type, vdj_gene) DO UPDATE SET
                        study = EXCLUDED.study,
                        count = EXCLUDED.count,
                        frequency = EXCLUDED.frequency
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'diversity':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (study, sample_id, n_cells, n_clones, shannon, simpson, inv_simpson, evenness, gini, d50, top10_fraction)
                    SELECT
                        COALESCE((SELECT s.study FROM "sample" s WHERE s.sample_id = t.sample_id), t.study) AS study,
                        t.sample_id, t.n_cells, t.n_clones, t.shannon, t.simpson, t.inv_simpson, t.evenness, t.gini, t.d50, t.top10_fraction
                    FROM {} t
                    ON CONFLICT (sample_id) DO UPDATE SET
                        study = EXCLUDED.study,
                        n_cells = EXCLUDED.n_cells,
                        n_clones = EXCLUDED.n_clones,
                        shannon = EXCLUDED.shannon,
                        simpson = EXCLUDED.simpson,
                        inv_simpson = EXCLUDED.inv_simpson,
                        evenness = EXCLUDED.evenness,
                        gini = EXCLUDED.gini,
                        d50 = EXCLUDED.d50,
                        top10_fraction = EXCLUDED.top10_fraction
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'gex_cell_embedding':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (
                        study, sample_id, cell_barcode, cell_scope, embedding_key,
                        umap_1, umap_2, cell_type, cell_subtype, cluster
                    )
                    SELECT
                        COALESCE((SELECT s.study FROM "sample" s WHERE s.sample_id = t.sample_id), t.study) AS study,
                        t.sample_id, t.cell_barcode, t.cell_scope, t.embedding_key,
                        t.umap_1, t.umap_2, t.cell_type, t.cell_subtype, t.cluster
                    FROM {} t
                    ON CONFLICT (sample_id, cell_scope, embedding_key, cell_barcode) DO UPDATE SET
                        study = EXCLUDED.study,
                        umap_1 = EXCLUDED.umap_1,
                        umap_2 = EXCLUDED.umap_2,
                        cell_type = EXCLUDED.cell_type,
                        cell_subtype = EXCLUDED.cell_subtype,
                        cluster = EXCLUDED.cluster
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'gex_marker_expression':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (
                        study, sample_id, cell_barcode, cell_scope, gene,
                        expression_layer, expression, detected
                    )
                    SELECT
                        COALESCE((SELECT s.study FROM "sample" s WHERE s.sample_id = t.sample_id), t.study) AS study,
                        t.sample_id, t.cell_barcode, t.cell_scope, t.gene,
                        t.expression_layer, t.expression, t.detected
                    FROM {} t
                    ON CONFLICT (sample_id, cell_scope, expression_layer, gene, cell_barcode) DO UPDATE SET
                        study = EXCLUDED.study,
                        expression = EXCLUDED.expression,
                        detected = EXCLUDED.detected
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'clone_gene_expression_summary':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (
                        study, sample_id, clone_id, gene, expression_layer,
                        n_cells, mean_expression, variance_expression,
                        detected_cells, detection_fraction
                    )
                    SELECT
                        COALESCE((SELECT s.study FROM "sample" s WHERE s.sample_id = t.sample_id), t.study) AS study,
                        t.sample_id, t.clone_id, t.gene, t.expression_layer,
                        t.n_cells, t.mean_expression, t.variance_expression,
                        t.detected_cells, t.detection_fraction
                    FROM {} t
                    ON CONFLICT (sample_id, expression_layer, clone_id, gene) DO UPDATE SET
                        study = EXCLUDED.study,
                        n_cells = EXCLUDED.n_cells,
                        mean_expression = EXCLUDED.mean_expression,
                        variance_expression = EXCLUDED.variance_expression,
                        detected_cells = EXCLUDED.detected_cells,
                        detection_fraction = EXCLUDED.detection_fraction
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'gex_expression_store':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (
                        study, sample_id, cell_scope, store_format, zarr_path,
                        expression_layer, n_cells, n_genes, source_h5ad_path, created_at
                    )
                    SELECT
                        COALESCE((SELECT s.study FROM "sample" s WHERE s.sample_id = t.sample_id), t.study) AS study,
                        t.sample_id, t.cell_scope, t.store_format, t.zarr_path,
                        t.expression_layer, t.n_cells, t.n_genes, t.source_h5ad_path, t.created_at
                    FROM {} t
                    ON CONFLICT (sample_id, cell_scope, store_format) DO UPDATE SET
                        study = EXCLUDED.study,
                        zarr_path = EXCLUDED.zarr_path,
                        expression_layer = EXCLUDED.expression_layer,
                        n_cells = EXCLUDED.n_cells,
                        n_genes = EXCLUDED.n_genes,
                        source_h5ad_path = EXCLUDED.source_h5ad_path,
                        created_at = EXCLUDED.created_at
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'lineage_edge':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (
                        study, sample_id, clone_id, parent_node, child_node,
                        parent_type, child_type, distance, edge_weight,
                        is_germline_edge, build_method
                    )
                    SELECT
                        COALESCE((SELECT s.study FROM "sample" s WHERE s.sample_id = t.sample_id), t.study) AS study,
                        t.sample_id, t.clone_id, t.parent_node, t.child_node,
                        t.parent_type, t.child_type, t.distance, t.edge_weight,
                        t.is_germline_edge, t.build_method
                    FROM {} t
                    ON CONFLICT (sample_id, clone_id, parent_node, child_node) DO UPDATE SET
                        study = EXCLUDED.study,
                        parent_type = EXCLUDED.parent_type,
                        child_type = EXCLUDED.child_type,
                        distance = EXCLUDED.distance,
                        edge_weight = EXCLUDED.edge_weight,
                        is_germline_edge = EXCLUDED.is_germline_edge,
                        build_method = EXCLUDED.build_method
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))
            elif table_name == 'lineage_tree_json':
                cur.execute(sql.SQL("""
                    INSERT INTO {} (
                        study, sample_id, clone_id, n_nodes, n_edges,
                        build_method, tree_json
                    )
                    SELECT
                        COALESCE((SELECT s.study FROM "sample" s WHERE s.sample_id = t.sample_id), t.study) AS study,
                        t.sample_id, t.clone_id, t.n_nodes, t.n_edges,
                        t.build_method, t.tree_json
                    FROM {} t
                    ON CONFLICT (sample_id, clone_id) DO UPDATE SET
                        study = EXCLUDED.study,
                        n_nodes = EXCLUDED.n_nodes,
                        n_edges = EXCLUDED.n_edges,
                        build_method = EXCLUDED.build_method,
                        tree_json = EXCLUDED.tree_json
                """).format(sql.Identifier(table_name), sql.Identifier(temp_table)))

        conn.commit()
        logger.info(f"成功导入 {file_path} -> {table_name}，处理 {row_count} 行")
        return True
    except Exception as e:
        conn.rollback()
        logger.error(f"导入 {file_path} 到 {table_name} 失败: {e}")
        return False

_SAMPLE_IDS_CACHE = None


def load_sample_ids_from_metadata():
    """读取 metadata/sample.tsv 中的真实 sample_id，并做进程内缓存。"""
    global _SAMPLE_IDS_CACHE
    if _SAMPLE_IDS_CACHE is not None:
        return _SAMPLE_IDS_CACHE

    sample_metadata = script_dir.parent / 'metadata' / 'sample.tsv'
    if not sample_metadata.exists():
        _SAMPLE_IDS_CACHE = []
        return _SAMPLE_IDS_CACHE

    try:
        with open(sample_metadata, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='	')
            _SAMPLE_IDS_CACHE = [row['sample_id'] for row in reader if row.get('sample_id')]
    except Exception as e:
        logger.warning(f"无法读取metadata/sample.tsv中的sample_id: {e}")
        _SAMPLE_IDS_CACHE = []
    return _SAMPLE_IDS_CACHE


def match_sample_id_from_metadata(name):
    """用metadata/sample.tsv从 {study}_{sample} 风格文件名中匹配真实sample_id。"""
    sample_ids = load_sample_ids_from_metadata()
    try:
        matches = [sample_id for sample_id in sample_ids if name == sample_id or name.endswith(f"_{sample_id}")]
        if matches:
            return max(matches, key=len)

        suffix_matches = [sample_id for sample_id in sample_ids if sample_id.endswith(f"_{name}")]
        if len(suffix_matches) == 1:
            return suffix_matches[0]
        if len(suffix_matches) > 1:
            logger.warning(f"{name} 可匹配多个metadata sample_id，跳过短名匹配: {suffix_matches[:5]}")
    except Exception as e:
        logger.warning(f"无法从metadata匹配 {name} 的sample_id: {e}")

    return None

def get_sample_id_from_filename(file_path):
    """从文件名提取样本ID，支持 {study}_{sample}_*.tsv。"""
    name = Path(file_path).stem

    for suffix in (
        "_clone",
        "_lineage_edge",
        "_lineage_tree",
        "_vdj_usage",
        "_vdj",
        "_SHM_paired",
        "_SHM_IGH",
    ):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break

    matched = match_sample_id_from_metadata(name)
    if matched:
        return matched

    return name


def get_sample_id_from_embedding_filename(file_path):
    """从文件名提取样本ID，支持 {study}_{sample}_*.tsv。"""
    name = Path(file_path).stem
    for suffix in (
        "_b_cell_embedding",
        "_cell_embedding",
        "_embedding",
        "_lineage_edge",
        "_lineage_tree",
        "_b_marker_expression",
        "_b_zarr_manifest",
    ):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
            break

    matched = match_sample_id_from_metadata(name)
    if matched:
        return matched

    return get_sample_id_from_filename(file_path)

def get_sample_id_from_tsv(file_path):
    """从TSV第一行读取sample_id，并校正为metadata里的真实sample_id。"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            reader = csv.DictReader(f, delimiter='	')
            row = next(reader, None)
            if row and row.get('sample_id'):
                raw_sample_id = row['sample_id'].strip()
                sample_ids = set(load_sample_ids_from_metadata())
                if raw_sample_id in sample_ids:
                    return raw_sample_id

                filename_match = get_sample_id_from_embedding_filename(file_path)
                if filename_match in sample_ids:
                    return filename_match

                raw_match = match_sample_id_from_metadata(raw_sample_id)
                if raw_match:
                    return raw_match

                return raw_sample_id
    except Exception as e:
        logger.warning(f"无法从 {file_path} 读取 sample_id，回退到文件名解析: {e}")

    return get_sample_id_from_embedding_filename(file_path)

def is_sample_data_imported(conn, sample_id, data_type):
    """检查样本的特定数据类型是否已导入"""
    try:
        with conn.cursor() as cur:
            cur.execute(
                "SELECT 1 FROM import_log WHERE sample_id = %s AND data_type = %s AND status = 'success'",
                (sample_id, data_type)
            )
            return cur.fetchone() is not None
    except psycopg2.Error as e:
        logger.error(f"检查导入状态失败: {e}")
        return False

def mark_sample_data_imported(conn, sample_id, data_type, file_path=None, status='success'):
    """标记样本的特定数据类型已导入"""
    try:
        with conn.cursor() as cur:
            cur.execute("""
                INSERT INTO import_log (sample_id, data_type, file_path, status, import_time)
                VALUES (%s, %s, %s, %s, %s)
                ON CONFLICT (sample_id, data_type) 
                DO UPDATE SET 
                    import_time = EXCLUDED.import_time,
                    status = EXCLUDED.status,
                    file_path = EXCLUDED.file_path
            """, (sample_id, data_type, str(file_path) if file_path else None, status, datetime.now()))
        conn.commit()
    except psycopg2.Error as e:
        conn.rollback()
        logger.error(f"标记导入状态失败: {e}")

############################################
# BCR轻链拼接辅助
############################################

def infer_paired_file_from_clone(clone_file):
    """根据 clone 文件路径推断对应 SHM_paired 文件。"""
    clone_path = Path(clone_file)
    name = clone_path.name
    if name.endswith("_clone.tsv"):
        candidate = clone_path.with_name(name.replace("_clone.tsv", "_SHM_paired.tsv"))
        if candidate.exists():
            return str(candidate)
    return None


def build_bcr_import_with_light(clone_file, paired_file):
    """
    生成用于 bcr_sequences 入库的临时TSV：
    - 保留 clone_file 中原有行（通常为IGH）
    - 从 paired_file 中补充轻链行（IGK/IGL），并尽量映射 clone_id
    """
    clone_path = Path(clone_file)
    paired_path = Path(paired_file)
    if not paired_path.exists():
        raise FileNotFoundError(f"paired文件不存在: {paired_file}")

    with open(clone_path, 'r', encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f, delimiter='	')
        headers = reader.fieldnames or []

    if not headers:
        raise ValueError(f"{clone_file} 表头为空")

    required_cols = {"sequence_id", "sample_id", "clone_id", "cell_barcode"}
    missing = sorted(required_cols - set(headers))
    if missing:
        raise ValueError(f"{clone_file} 缺少关键列: {','.join(missing)}")

    temp = tempfile.NamedTemporaryFile(
        mode='w',
        encoding='utf-8',
        newline='',
        suffix='_bcr_with_light.tsv',
        prefix='load_db_',
        delete=False,
    )

    heavy_count = 0
    light_count = 0
    skipped_light = 0
    sample_id_fallback = ""
    existing_sequence_ids = set()
    seq_to_clone = {}
    cell_to_clone = {}

    try:
        writer = csv.DictWriter(temp, fieldnames=headers, delimiter='	', extrasaction='ignore')
        writer.writeheader()

        # 1) 先写入原始clone（IGH）行
        with open(clone_path, 'r', encoding='utf-8', newline='') as f_clone:
            clone_reader = csv.DictReader(f_clone, delimiter='	')
            for row in clone_reader:
                seq_id = (row.get('sequence_id') or '').strip()
                barcode = (row.get('cell_barcode') or '').strip()
                clone_id = (row.get('clone_id') or '').strip()
                sample_id = (row.get('sample_id') or '').strip()

                if sample_id and not sample_id_fallback:
                    sample_id_fallback = sample_id
                if seq_id:
                    existing_sequence_ids.add(seq_id)
                    if clone_id:
                        seq_to_clone[seq_id] = clone_id
                if barcode and clone_id:
                    cell_to_clone[barcode] = clone_id

                writer.writerow({h: row.get(h, '') for h in headers})
                heavy_count += 1

        # 2) 从paired补充轻链行（_L列）
        with open(paired_path, 'r', encoding='utf-8', newline='') as f_paired:
            paired_reader = csv.DictReader(f_paired, delimiter='	')
            for prow in paired_reader:
                light_seq_id = (prow.get('sequence_id_L') or '').strip()
                light_seq = prow.get('sequence_L')
                if not light_seq_id or light_seq_id in existing_sequence_ids:
                    skipped_light += 1
                    continue
                if light_seq is None or str(light_seq).strip() == '':
                    skipped_light += 1
                    continue

                row = {h: '' for h in headers}
                for h in headers:
                    key_l = f"{h}_L"
                    if key_l in prow and prow.get(key_l) is not None:
                        row[h] = prow.get(key_l, '')

                # 覆盖关键字段
                row['sequence_id'] = light_seq_id
                row['sequence'] = prow.get('sequence_L', row.get('sequence', ''))
                row['sample_id'] = (prow.get('sample_id') or sample_id_fallback or row.get('sample_id') or '').strip()
                row['cell_barcode'] = (prow.get('cell_barcode') or row.get('cell_barcode') or '').strip()

                if not str(row.get('cdr3_aa', '')).strip():
                    row['cdr3_aa'] = (prow.get('junction_aa_L') or '').strip()

                heavy_seq_id = (prow.get('sequence_id_H') or '').strip()
                mapped_clone = seq_to_clone.get(heavy_seq_id) or cell_to_clone.get(row['cell_barcode']) or ''
                row['clone_id'] = mapped_clone

                writer.writerow(row)
                existing_sequence_ids.add(light_seq_id)
                light_count += 1

        temp.flush()
        temp.close()

        logger.info(
            f"构建BCR入库临时文件: {temp.name} (IGH={heavy_count}, light_added={light_count}, light_skipped={skipped_light})"
        )
        return temp.name
    except Exception:
        temp.close()
        try:
            Path(temp.name).unlink(missing_ok=True)
        except Exception:
            pass
        raise


def extract_year_from_text(*values):
    """Extract a four-digit year from free text values."""
    for value in values:
        text = str(value or "")
        match = re.search(r"(19|20)\d{2}", text)
        if match:
            return int(match.group(0))
    return None


def normalize_study_metadata_file(file_path):
    """Accept either database schema or rich study metadata schema.

    Supported database schema: study, title, year, journal, accession.
    Supported rich schema: Study, Disease, Title, Publication, Patient count,
    Sample count, Accession. Rich Study IDs are mapped to database-facing study
    IDs through metadata/pipeline_samples.tsv when available.
    """
    path = Path(file_path)
    with open(path, 'r', encoding='utf-8', newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames or []
        rows = list(reader)

    if {'study', 'title', 'year', 'journal', 'accession'}.issubset(set(fieldnames)):
        return str(path)

    if not {'Study', 'Title', 'Publication'}.issubset(set(fieldnames)):
        raise ValueError(
            f"{file_path} study metadata columns are not supported: {','.join(fieldnames)}"
        )

    raw_to_study = {}
    pipeline_file = script_dir.parent / 'metadata' / 'pipeline_samples.tsv'
    if pipeline_file.exists():
        with open(pipeline_file, 'r', encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                raw_study = (row.get('raw_study') or '').strip()
                study = (row.get('study') or '').strip()
                if raw_study and study:
                    raw_to_study[raw_study] = study

    temp = tempfile.NamedTemporaryFile(
        mode='w',
        encoding='utf-8',
        newline='',
        suffix='_study_normalized.tsv',
        prefix='load_db_',
        delete=False,
    )
    writer = csv.DictWriter(
        temp,
        fieldnames=['study', 'title', 'year', 'journal', 'accession'],
        delimiter='\t',
    )
    writer.writeheader()
    for row in rows:
        raw_study = (row.get('Study') or '').strip()
        if not raw_study:
            continue
        writer.writerow({
            'study': raw_study,
            'title': (row.get('Title') or '').strip(),
            'year': extract_year_from_text(row.get('Publication'), raw_study) or '',
            'journal': (row.get('Publication') or '').strip(),
            'accession': (row.get('Accession') or '').strip(),
        })
    temp.flush()
    temp.close()
    logger.info(f"已将富信息study元数据转换为数据库入库格式: {temp.name}")
    return temp.name


############################################
# 主函数
############################################

def main():
    parser = argparse.ArgumentParser(description='加载BCR分析结果到数据库')
    parser.add_argument('--study', required=True, help='study元数据文件路径')
    parser.add_argument('--subject', required=True, help='subject元数据文件路径')
    parser.add_argument('--sample', required=True, help='sample元数据文件路径')
    parser.add_argument('--clones', nargs='+', required=True, help='BCR克隆文件列表（每个样本一个）')
    parser.add_argument('--paired-bcr', nargs='+', help='BCR重轻链配对文件列表（*_SHM_paired.tsv）')
    parser.add_argument('--vdj-usage', nargs='+', help='V(D)J使用文件列表（每个样本一个）')
    parser.add_argument('--cell-embeddings', nargs='+', help='单细胞UMAP/celltype文件列表（每个样本一个）')
    parser.add_argument('--marker-expression', nargs='+', help='B细胞marker gene表达缓存文件列表（每个样本一个）')
    parser.add_argument('--clone-gene-expression', nargs='+', help='clone-level表达摘要文件列表（每个样本一个）')
    parser.add_argument('--expression-stores', nargs='+', help='B细胞完整表达矩阵zarr manifest文件列表（每个样本一个）')
    parser.add_argument('--lineage-edges', nargs='+', help='clone lineage边文件列表（每个样本一个）')
    parser.add_argument('--lineage-trees', nargs='+', help='clone lineage JSON文件列表（每个样本一个）')
    parser.add_argument('--member', help='克隆型成员文件（全局）')
    parser.add_argument('--public', help='公共克隆型文件（全局）')
    parser.add_argument('--network', help='公共网络文件（全局）')
    parser.add_argument('--diversity', help='多样性文件（全局）')
    parser.add_argument('--reload', action='store_true', help='重新导入已存在的样本')
    parser.add_argument('--clear-all', action='store_true', help='清空所有表重新导入（危险操作）')

    args = parser.parse_args()

    # 验证输入文件存在
    required_files = [args.study, args.subject, args.sample] + args.clones
    for file_path in required_files:
        if not Path(file_path).exists():
            logger.error(f"必需文件不存在: {file_path}")
            sys.exit(1)

    if args.paired_bcr:
        for file_path in args.paired_bcr:
            if not Path(file_path).exists():
                logger.warning(f"paired文件不存在，将在导入时跳过轻链补充: {file_path}")
    
    # 连接数据库
    conn = get_db_connection()
    
    try:
        ensure_sample_paired_bcr_column(conn)
        ensure_bcr_cdr3_aa_column(conn)
        ensure_study_accession_column(conn)

        # 清空旧数据：--clear-all 需要人工确认，--reload 自动进行全量重载。
        if args.clear_all or args.reload:
            if args.clear_all:
                logger.warning("⚠️ 正在清空所有表！")
                confirm = input("确认清空所有数据？(输入 'yes' 确认): ")
                if confirm.lower() != 'yes':
                    logger.info("操作已取消")
                    return
            else:
                logger.info("--reload 已启用：清空旧 metadata 和结果表后重新导入")

            tables = [
                "import_log",
                "public_network",
                "clonotype_member",
                "public_clonotype",
                "diversity",
                "bcr_sequences",
                "vdj_usage",
                "sample",
                "subject",
                "study",
            ]
            optional_tables = [
                "gex_cell_embedding",
                "gex_marker_expression",
                "clone_gene_expression_summary",
                "gex_expression_store",
                "lineage_edge",
                "lineage_tree_json",
                "download_file",
            ]
            tables = [table for table in optional_tables + tables if table_exists(conn, table)]
            with conn.cursor() as cur:
                cur.execute(
                    sql.SQL("TRUNCATE TABLE {} RESTART IDENTITY CASCADE").format(
                        sql.SQL(", ").join(sql.Identifier(table) for table in tables)
                    )
                )
            conn.commit()
            logger.info("旧数据库数据已清空")
        
        # 1. 导入元数据
        logger.info("开始导入元数据...")
        logger.info("开始导入全局分析结果...")
        study_import_file = normalize_study_metadata_file(args.study)
        metadata_files = [
            (study_import_file, "study", "study"),
            (args.subject, "subject", "subject"),
            (args.sample, "sample", "sample")
        ]
        
        for file_path, table_name, data_type in metadata_files:
            if file_path and Path(file_path).exists():

                success = copy_table(conn, file_path, table_name)
                if success:
                    mark_sample_data_imported(conn, 'global', data_type, file_path, 'success')

        logger.info("确保前端可选模块依赖表存在...")
        ensure_frontend_required_tables(conn)

        # 2. 导入单细胞UMAP/celltype数据（按样本）
        if args.cell_embeddings:
            ensure_gex_cell_embedding_table(conn)
            logger.info("开始导入单细胞UMAP/celltype数据...")
            total_embeddings = len(args.cell_embeddings)
            for i, embedding_file in enumerate(args.cell_embeddings, 1):
                sample_id = get_sample_id_from_tsv(embedding_file)
                logger.info(f"[{i}/{total_embeddings}] 导入样本 {sample_id} 的UMAP/celltype数据")

                if not args.reload and is_sample_data_imported(conn, sample_id, 'gex_cell_embedding'):
                    logger.info(f"样本 {sample_id} 的UMAP/celltype数据已存在，跳过")
                    continue

                success = copy_table(conn, embedding_file, "gex_cell_embedding", sample_id)

                if success:
                    mark_sample_data_imported(conn, sample_id, 'gex_cell_embedding', embedding_file, 'success')
                else:
                    mark_sample_data_imported(conn, sample_id, 'gex_cell_embedding', embedding_file, 'failed')
                    logger.error(f"样本 {sample_id} UMAP/celltype数据导入失败")

        # 3. 导入B细胞marker gene表达缓存（按样本）
        if args.marker_expression:
            ensure_gex_marker_expression_table(conn)
            logger.info("开始导入B细胞marker gene表达缓存...")
            total_marker_files = len(args.marker_expression)
            for i, marker_file in enumerate(args.marker_expression, 1):
                sample_id = get_sample_id_from_tsv(marker_file)
                logger.info(f"[{i}/{total_marker_files}] 导入样本 {sample_id} 的marker表达缓存")

                if not args.reload and is_sample_data_imported(conn, sample_id, 'gex_marker_expression'):
                    logger.info(f"样本 {sample_id} 的marker表达缓存已存在，跳过")
                    continue

                success = copy_table(conn, marker_file, "gex_marker_expression", sample_id)

                if success:
                    mark_sample_data_imported(conn, sample_id, 'gex_marker_expression', marker_file, 'success')
                else:
                    mark_sample_data_imported(conn, sample_id, 'gex_marker_expression', marker_file, 'failed')
                    logger.error(f"样本 {sample_id} marker表达缓存导入失败")

        # 4. 导入clone-level表达摘要（按样本）
        if args.clone_gene_expression:
            ensure_clone_gene_expression_summary_table(conn)
            logger.info("开始导入clone-level表达摘要...")
            total_clone_gene_files = len(args.clone_gene_expression)
            for i, clone_gene_file in enumerate(args.clone_gene_expression, 1):
                sample_id = get_sample_id_from_tsv(clone_gene_file)
                logger.info(f"[{i}/{total_clone_gene_files}] 导入样本 {sample_id} 的clone-level表达摘要")

                if not args.reload and is_sample_data_imported(conn, sample_id, 'clone_gene_expression_summary'):
                    logger.info(f"样本 {sample_id} 的clone-level表达摘要已存在，跳过")
                    continue

                success = copy_table(conn, clone_gene_file, "clone_gene_expression_summary", sample_id)

                if success:
                    mark_sample_data_imported(conn, sample_id, 'clone_gene_expression_summary', clone_gene_file, 'success')
                else:
                    mark_sample_data_imported(conn, sample_id, 'clone_gene_expression_summary', clone_gene_file, 'failed')
                    logger.error(f"样本 {sample_id} clone-level表达摘要导入失败")

        # 5. 导入完整表达矩阵zarr存储manifest（按样本）
        if args.expression_stores:
            ensure_gex_expression_store_table(conn)
            logger.info("开始导入B细胞完整表达矩阵zarr manifest...")
            total_store_files = len(args.expression_stores)
            for i, store_file in enumerate(args.expression_stores, 1):
                sample_id = get_sample_id_from_tsv(store_file)
                logger.info(f"[{i}/{total_store_files}] 导入样本 {sample_id} 的zarr manifest")

                if not args.reload and is_sample_data_imported(conn, sample_id, 'gex_expression_store'):
                    logger.info(f"样本 {sample_id} 的zarr manifest已存在，跳过")
                    continue

                success = copy_table(conn, store_file, "gex_expression_store", sample_id)

                if success:
                    mark_sample_data_imported(conn, sample_id, 'gex_expression_store', store_file, 'success')
                else:
                    mark_sample_data_imported(conn, sample_id, 'gex_expression_store', store_file, 'failed')
                    logger.error(f"样本 {sample_id} zarr manifest导入失败")

        # 6. 导入lineage edge（按样本）
        if args.lineage_edges:
            ensure_lineage_edge_table(conn)
            logger.info("开始导入clone lineage边数据...")
            total_lineage_edges = len(args.lineage_edges)
            for i, edge_file in enumerate(args.lineage_edges, 1):
                sample_id = get_sample_id_from_tsv(edge_file)
                logger.info(f"[{i}/{total_lineage_edges}] 导入样本 {sample_id} 的lineage边")

                if not args.reload and is_sample_data_imported(conn, sample_id, 'lineage_edge'):
                    logger.info(f"样本 {sample_id} 的lineage边已存在，跳过")
                    continue

                success = copy_table(conn, edge_file, "lineage_edge", sample_id)
                if success:
                    mark_sample_data_imported(conn, sample_id, 'lineage_edge', edge_file, 'success')
                else:
                    mark_sample_data_imported(conn, sample_id, 'lineage_edge', edge_file, 'failed')
                    logger.error(f"样本 {sample_id} lineage边导入失败")

        # 7. 导入lineage tree JSON（按样本）
        if args.lineage_trees:
            ensure_lineage_tree_table(conn)
            logger.info("开始导入clone lineage JSON数据...")
            total_lineage_trees = len(args.lineage_trees)
            for i, tree_file in enumerate(args.lineage_trees, 1):
                sample_id = get_sample_id_from_tsv(tree_file)
                logger.info(f"[{i}/{total_lineage_trees}] 导入样本 {sample_id} 的lineage JSON")

                if not args.reload and is_sample_data_imported(conn, sample_id, 'lineage_tree_json'):
                    logger.info(f"样本 {sample_id} 的lineage JSON已存在，跳过")
                    continue

                success = copy_table(conn, tree_file, "lineage_tree_json", sample_id)
                if success:
                    mark_sample_data_imported(conn, sample_id, 'lineage_tree_json', tree_file, 'success')
                else:
                    mark_sample_data_imported(conn, sample_id, 'lineage_tree_json', tree_file, 'failed')
                    logger.error(f"样本 {sample_id} lineage JSON导入失败")
        
        # 7. 导入BCR序列（按样本，IGH + 轻链）
        logger.info("开始导入BCR序列（含轻链补充）...")
        paired_map = {}
        if args.paired_bcr:
            for paired_file in args.paired_bcr:
                sid = get_sample_id_from_tsv(paired_file)
                paired_map[sid] = paired_file

        total_samples = len(args.clones)
        for i, clone_file in enumerate(args.clones, 1):
            sample_id = get_sample_id_from_tsv(clone_file)
            logger.info(f"[{i}/{total_samples}] 导入样本 {sample_id} 的BCR数据")

            if not args.reload and is_sample_data_imported(conn, sample_id, 'bcr'):
                logger.info(f"样本 {sample_id} 的BCR数据已存在，跳过")
                continue

            paired_file = paired_map.get(sample_id) or infer_paired_file_from_clone(clone_file)
            import_file = clone_file
            import_log_path = clone_file
            temp_file = None

            try:
                if paired_file and Path(paired_file).exists():
                    logger.info(f"样本 {sample_id} 检测到paired文件，补充轻链: {paired_file}")
                    temp_file = build_bcr_import_with_light(clone_file, paired_file)
                    import_file = temp_file
                    import_log_path = f"{clone_file} | {paired_file}"
                else:
                    logger.warning(f"样本 {sample_id} 未找到paired文件，仅导入clone文件（通常仅IGH）")

                success = copy_table(conn, import_file, "bcr_sequences", sample_id)

                if success:
                    mark_sample_data_imported(conn, sample_id, 'bcr', import_log_path, 'success')
                else:
                    mark_sample_data_imported(conn, sample_id, 'bcr', import_log_path, 'failed')
                    logger.error(f"样本 {sample_id} BCR数据导入失败")
            finally:
                if temp_file:
                    try:
                        Path(temp_file).unlink(missing_ok=True)
                    except Exception:
                        logger.warning(f"临时文件删除失败: {temp_file}")
        
        # 8. 导入V(D)J使用数据（按样本）
        if args.vdj_usage:
            logger.info("开始导入V(D)J使用数据...")
            for vdj_file in args.vdj_usage:
                sample_id = get_sample_id_from_embedding_filename(vdj_file)

                
                logger.info(f"导入样本 {sample_id} 的V(D)J数据")
                success = copy_table(conn, vdj_file, "vdj_usage", sample_id)
                
                if success:
                    mark_sample_data_imported(conn, sample_id, 'vdj_usage', vdj_file, 'success')
                else:
                    mark_sample_data_imported(conn, sample_id, 'vdj_usage', vdj_file, 'failed')
        
        # 9. 导入全局分析结果
        logger.info("开始导入全局分析结果...")
        global_files = [
            (args.member, "clonotype_member", "clone_members"),
            (args.public, "public_clonotype", "public_clonotypes"),
            (args.network, "public_network", "public_network"),
            (args.diversity, "diversity", "diversity")
        ]
        
        public_global_tables = {"clonotype_member", "public_clonotype", "public_network"}
        if any(file_path and Path(file_path).exists() and table_name in public_global_tables
               for file_path, table_name, _data_type in global_files):
            with conn.cursor() as cur:
                # Public clonotype results are regenerated as one cohort-level snapshot.
                # Clear the previous snapshot first so obsolete global clonotypes do not remain.
                for table_name in ("public_network", "clonotype_member", "public_clonotype"):
                    cur.execute(sql.SQL("DELETE FROM {}").format(sql.Identifier(table_name)))
            conn.commit()
            logger.info("已清空旧的公共克隆全局结果")

        for file_path, table_name, data_type in global_files:
            if file_path and Path(file_path).exists():
                
                success = copy_table(conn, file_path, table_name)
                if success:
                    mark_sample_data_imported(conn, 'global', data_type, file_path, 'success')
        
        # 10. 创建Analysis页面高频查询组合索引并更新统计信息
        logger.info("创建Analysis页面高频查询组合索引...")
        ensure_analysis_indexes(conn)

        # 创建完成标记文件
        with open("results/database_loaded.txt", "w") as f:
            f.write(f"Database loaded successfully at {datetime.now()}\n")
            f.write(f"BCR samples: {len(args.clones)}\n")
            if args.vdj_usage:
                f.write(f"VDJ usage samples: {len(args.vdj_usage)}\n")
            if args.cell_embeddings:
                f.write(f"Cell embedding samples: {len(args.cell_embeddings)}\n")
            if args.marker_expression:
                f.write(f"Marker expression samples: {len(args.marker_expression)}\n")
            if args.expression_stores:
                f.write(f"Expression zarr stores: {len(args.expression_stores)}\n")
            if args.lineage_edges:
                f.write(f"Lineage edge samples: {len(args.lineage_edges)}\n")
            if args.lineage_trees:
                f.write(f"Lineage tree samples: {len(args.lineage_trees)}\n")
        
        # 打印统计信息
        with conn.cursor() as cur:
            cur.execute("SELECT COUNT(*) FROM bcr_sequences")
            bcr_count = cur.fetchone()[0]
            cur.execute("SELECT COUNT(*) FROM vdj_usage")
            vdj_count = cur.fetchone()[0]
            cur.execute("SELECT COUNT(DISTINCT sample_id) FROM bcr_sequences")
            sample_count = cur.fetchone()[0]
            embedding_count = 0
            embedding_sample_count = 0
            marker_expression_count = 0
            marker_gene_count = 0
            expression_store_count = 0
            lineage_edge_count = 0
            lineage_clone_count = 0
            if table_exists(conn, "gex_cell_embedding"):
                cur.execute("SELECT COUNT(*) FROM gex_cell_embedding")
                embedding_count = cur.fetchone()[0]
                cur.execute("SELECT COUNT(DISTINCT sample_id) FROM gex_cell_embedding")
                embedding_sample_count = cur.fetchone()[0]
            if table_exists(conn, "gex_marker_expression"):
                cur.execute("SELECT COUNT(*) FROM gex_marker_expression")
                marker_expression_count = cur.fetchone()[0]
                cur.execute("SELECT COUNT(DISTINCT gene) FROM gex_marker_expression")
                marker_gene_count = cur.fetchone()[0]
            if table_exists(conn, "gex_expression_store"):
                cur.execute("SELECT COUNT(*) FROM gex_expression_store")
                expression_store_count = cur.fetchone()[0]
            if table_exists(conn, "lineage_edge"):
                cur.execute("SELECT COUNT(*) FROM lineage_edge")
                lineage_edge_count = cur.fetchone()[0]
            if table_exists(conn, "lineage_tree_json"):
                cur.execute("SELECT COUNT(*) FROM lineage_tree_json")
                lineage_clone_count = cur.fetchone()[0]
        
        logger.info(f"✅ 所有数据导入完成！")
        logger.info(f"   - 样本数: {sample_count}")
        logger.info(f"   - BCR序列数: {bcr_count}")
        logger.info(f"   - V(D)J记录数: {vdj_count}")
        logger.info(f"   - UMAP/celltype样本数: {embedding_sample_count}")
        logger.info(f"   - UMAP/celltype细胞数: {embedding_count}")
        logger.info(f"   - marker表达缓存记录数: {marker_expression_count}")
        logger.info(f"   - marker gene数: {marker_gene_count}")
        logger.info(f"   - zarr表达矩阵store数: {expression_store_count}")
        logger.info(f"   - lineage边数: {lineage_edge_count}")
        logger.info(f"   - lineage克隆树数: {lineage_clone_count}")
        
    except Exception as e:
        logger.error(f"导入过程中出现错误: {e}")
        conn.rollback()
        sys.exit(1)
    finally:
        conn.close()

if __name__ == "__main__":
    main()

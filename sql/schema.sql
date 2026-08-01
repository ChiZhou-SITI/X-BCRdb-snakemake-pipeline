
CREATE EXTENSION IF NOT EXISTS pg_trgm;

------------------------------------------
-- study
------------------------------------------
DROP TABLE IF EXISTS study CASCADE;
CREATE TABLE IF NOT EXISTS study (
    study TEXT PRIMARY KEY,
    title TEXT,
    year INT,
    journal TEXT,
    accession TEXT
);

------------------------------------------
-- subject
------------------------------------------
DROP TABLE IF EXISTS subject CASCADE;
CREATE TABLE IF NOT EXISTS subject (
    subject_id TEXT PRIMARY KEY,
    species TEXT,
    disease TEXT
);

------------------------------------------
-- sample
------------------------------------------
DROP TABLE IF EXISTS sample CASCADE;
CREATE TABLE IF NOT EXISTS sample (
    sample_id TEXT PRIMARY KEY,
    study TEXT REFERENCES study(study) ON DELETE CASCADE,
    subject_id TEXT REFERENCES subject(subject_id) ON DELETE CASCADE,
    tissue TEXT,
    disease TEXT,
    platform TEXT,
    n_cells INT,
    paired_bcr BOOLEAN DEFAULT TRUE
);

------------------------------------------
-- downloadable h5ad / BCR JSON files
------------------------------------------
DROP TABLE IF EXISTS download_file CASCADE;
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
);
CREATE INDEX IF NOT EXISTS idx_download_file_study ON download_file(study);
CREATE INDEX IF NOT EXISTS idx_download_file_record ON download_file(zenodo_record_id);
CREATE INDEX IF NOT EXISTS idx_download_file_kind ON download_file(file_kind);

------------------------------------------
-- bcr sequence
------------------------------------------
-- 如果表已存在，先删除（谨慎操作！）

DROP TABLE IF EXISTS bcr_sequences CASCADE;
CREATE TABLE bcr_sequences (
    -- 基本标识
    sequence_id VARCHAR(500) NOT NULL,
    sequence TEXT NOT NULL,
    rev_comp TEXT,
    
    -- 功能注释
    productive BOOLEAN DEFAULT false,
    
    -- V(D)J基因调用
    v_call VARCHAR(50),
    d_call VARCHAR(50),
    j_call VARCHAR(50),
    
    -- 序列比对
    sequence_alignment TEXT,
    germline_alignment TEXT,
    junction TEXT,
    junction_aa VARCHAR(100),
    
    -- 比对格式(CIGAR)
    v_cigar VARCHAR(500),
    d_cigar VARCHAR(500),
    j_cigar VARCHAR(500),
    
    -- 功能标志
    stop_codon BOOLEAN DEFAULT false,
    vj_in_frame BOOLEAN DEFAULT false,
    
    -- 基因座
    locus VARCHAR(20),
    c_call VARCHAR(50),
    
    -- 长度信息
    junction_length INTEGER,
    np1_length NUMERIC DEFAULT 0,
    np2_length NUMERIC DEFAULT 0,
    
    -- V基因位置
    v_sequence_start INTEGER,
    v_sequence_end INTEGER,
    v_germline_start INTEGER,
    v_germline_end INTEGER,
    
    -- D基因位置
    d_sequence_start NUMERIC,
    d_sequence_end NUMERIC,
    d_germline_start NUMERIC,
    d_germline_end NUMERIC,
    
    -- J基因位置
    j_sequence_start NUMERIC,
    j_sequence_end NUMERIC,
    j_germline_start NUMERIC,
    j_germline_end NUMERIC,
    
    -- 比对分数
    v_score NUMERIC(10,3),
    v_identity NUMERIC(5,4),
    v_support NUMERIC,
    
    d_score NUMERIC(10,3),
    d_identity NUMERIC(5,4),
    d_support NUMERIC,
    
    j_score NUMERIC(10,3),
    j_identity NUMERIC(5,4),
    j_support NUMERIC,
    
    -- 框架区和CDR区
    fwr1 TEXT,
    fwr2 TEXT,
    fwr3 TEXT,
    fwr4 TEXT,
    cdr1 TEXT,
    cdr2 TEXT,
    cdr3 TEXT,
    cdr3_aa TEXT,
    
    -- 高级特征
    clone_id TEXT,
    germline_alignment_d_mask TEXT,
    mu_freq NUMERIC(10,6),
    
    -- 样本信息
    cell_barcode VARCHAR(100),
    sample_id VARCHAR(100) REFERENCES sample(sample_id) ON DELETE CASCADE,

    
    -- 复合主键（如果sequence_id可能重复，可以加上sample_id）
    PRIMARY KEY (sequence_id, sample_id)
);

-- 创建更新时间触发器
CREATE OR REPLACE FUNCTION update_updated_at_column()
RETURNS TRIGGER AS $$
BEGIN
    NEW.updated_at = CURRENT_TIMESTAMP;
    RETURN NEW;
END;
$$ language 'plpgsql';

CREATE TRIGGER update_bcr_sequences_updated_at
    BEFORE UPDATE ON bcr_sequences
    FOR EACH ROW
    EXECUTE FUNCTION update_updated_at_column();

-- 创建索引（优化版）
-- 1. 样本相关索引
CREATE INDEX idx_bcr_sequences_sample ON bcr_sequences(sample_id);
CREATE INDEX idx_bcr_sequences_barcode ON bcr_sequences(cell_barcode);

-- 2. 基因调用索引
CREATE INDEX idx_bcr_sequences_v_call ON bcr_sequences(v_call);
CREATE INDEX idx_bcr_sequences_j_call ON bcr_sequences(j_call);
CREATE INDEX idx_bcr_sequences_locus ON bcr_sequences(locus);

-- 3. 功能特征索引
CREATE INDEX idx_bcr_sequences_productive ON bcr_sequences(productive);
CREATE INDEX idx_bcr_sequences_vj_in_frame ON bcr_sequences(vj_in_frame);
CREATE INDEX idx_bcr_sequences_junction_aa ON bcr_sequences(junction_aa);
CREATE INDEX idx_bcr_sequences_cdr3_aa ON bcr_sequences(cdr3_aa);

-- 4. 复合索引（用于常见查询模式）
CREATE INDEX idx_bcr_sequences_sample_productive ON bcr_sequences(sample_id, productive);
CREATE INDEX idx_bcr_sequences_locus_v_call ON bcr_sequences(locus, v_call);
CREATE INDEX idx_bcr_sequences_v_call_productive ON bcr_sequences(v_call, productive);

-- 4b. Analysis页面高频聚合查询组合索引
CREATE INDEX IF NOT EXISTS idx_sample_study_sample ON sample(study, sample_id);
CREATE INDEX IF NOT EXISTS idx_sample_disease_sample ON sample((COALESCE(NULLIF(disease, ''), 'Unknown')), sample_id);
CREATE INDEX IF NOT EXISTS idx_sample_subject_sample ON sample(subject_id, sample_id);
CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_clone_barcode
  ON bcr_sequences(sample_id, locus, clone_id, cell_barcode)
  WHERE COALESCE(clone_id, '') <> '' AND COALESCE(cell_barcode, '') <> '';
CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_cdr3
  ON bcr_sequences(sample_id, locus, cdr3_aa)
  WHERE COALESCE(cdr3_aa, '') <> '';
CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_junction_aa
  ON bcr_sequences(sample_id, locus, junction_aa)
  WHERE COALESCE(junction_aa, '') <> '';
CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_mu
  ON bcr_sequences(sample_id, locus, mu_freq)
  WHERE mu_freq IS NOT NULL;
CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_c_call ON bcr_sequences(sample_id, locus, c_call);
CREATE INDEX IF NOT EXISTS idx_bcr_adv_sample_locus_cell
  ON bcr_sequences(sample_id, locus, cell_barcode)
  WHERE COALESCE(cell_barcode, '') <> '';

-- 5. 全文搜索索引（如果需要搜索序列）
CREATE INDEX idx_bcr_sequences_cdr3_gin ON bcr_sequences USING GIN (to_tsvector('simple', COALESCE(cdr3, '')));
CREATE INDEX idx_bcr_sequences_cdr3_aa_trgm ON bcr_sequences USING GIN (cdr3_aa gin_trgm_ops);

-- 添加注释
COMMENT ON TABLE bcr_sequences IS 'BCR/TCR免疫组库测序数据表';
COMMENT ON COLUMN bcr_sequences.sequence_id IS '序列唯一标识';
COMMENT ON COLUMN bcr_sequences.productive IS '是否为功能性序列';
COMMENT ON COLUMN bcr_sequences.v_call IS 'V基因调用';
COMMENT ON COLUMN bcr_sequences.d_call IS 'D基因调用';
COMMENT ON COLUMN bcr_sequences.j_call IS 'J基因调用';
COMMENT ON COLUMN bcr_sequences.junction_aa IS '连接区氨基酸序列';
COMMENT ON COLUMN bcr_sequences.locus IS '基因座(IGH/IGK/IGL等)';
COMMENT ON COLUMN bcr_sequences.cdr3 IS 'CDR3核酸序列';
COMMENT ON COLUMN bcr_sequences.cdr3_aa IS 'CDR3氨基酸序列';
COMMENT ON COLUMN bcr_sequences.mu_freq IS '突变频率';
COMMENT ON COLUMN bcr_sequences.sample_id IS '样本ID';
COMMENT ON COLUMN bcr_sequences.cell_barcode IS '细胞条形码';


------------------------------------------
-- clonotype member
------------------------------------------
DROP TABLE IF EXISTS clonotype_member CASCADE;

CREATE TABLE IF NOT EXISTS clonotype_member (
    clonotype_key TEXT,
    sample_id TEXT REFERENCES sample(sample_id) ON DELETE CASCADE,
    sequence_id TEXT,
    cdr3_aa TEXT,
    v_call TEXT,
    j_call TEXT,
    PRIMARY KEY (clonotype_key, sample_id, sequence_id)
);

------------------------------------------
-- public clonotype
------------------------------------------
DROP TABLE IF EXISTS public_clonotype CASCADE;
CREATE TABLE IF NOT EXISTS public_clonotype (

    clonotype_key TEXT PRIMARY KEY,
    n_cells INT,
    n_samples INT,
    samples TEXT
);

------------------------------------------
-- public network
------------------------------------------
DROP TABLE IF EXISTS public_network CASCADE;
CREATE TABLE IF NOT EXISTS public_network (

    sample1 TEXT,
    sample2 TEXT,
    clonotype_key TEXT
);

------------------------------------------
-- lineage edge
------------------------------------------
DROP TABLE IF EXISTS lineage_edge CASCADE;
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
);

CREATE INDEX idx_lineage_edge_sample ON lineage_edge(sample_id);
CREATE INDEX idx_lineage_edge_clone ON lineage_edge(sample_id, clone_id);

------------------------------------------
-- lineage tree json (per clone)
------------------------------------------
DROP TABLE IF EXISTS lineage_tree_json CASCADE;
CREATE TABLE IF NOT EXISTS lineage_tree_json (
    study TEXT REFERENCES study(study) ON DELETE CASCADE,
    sample_id TEXT REFERENCES sample(sample_id) ON DELETE CASCADE,
    clone_id TEXT NOT NULL,
    n_nodes INT,
    n_edges INT,
    build_method TEXT,
    tree_json JSONB NOT NULL,
    PRIMARY KEY (sample_id, clone_id)
);

CREATE INDEX idx_lineage_tree_sample ON lineage_tree_json(sample_id);
CREATE INDEX idx_lineage_tree_clone ON lineage_tree_json(sample_id, clone_id);

------------------------------------------
-- repertoire
------------------------------------------
DROP TABLE IF EXISTS vdj_usage CASCADE;
CREATE TABLE IF NOT EXISTS vdj_usage (
    study TEXT REFERENCES study(study) ON DELETE CASCADE,
    sample_id TEXT REFERENCES sample(sample_id) ON DELETE CASCADE,
    vdj_type TEXT,
    vdj_gene TEXT,
    count INT,
    frequency FLOAT,
    PRIMARY KEY (sample_id, vdj_type, vdj_gene)
);


------------------------------------------
-- diversity
------------------------------------------
DROP TABLE IF EXISTS diversity CASCADE;
CREATE TABLE IF NOT EXISTS diversity (
    study TEXT REFERENCES study(study) ON DELETE CASCADE,
    sample_id TEXT REFERENCES sample(sample_id) ON DELETE CASCADE,
    n_cells INT,
    n_clones INT,
    shannon FLOAT,
    simpson FLOAT,
    inv_simpson FLOAT,
    evenness FLOAT,
    gini FLOAT,
    d50 FLOAT,
    top10_fraction FLOAT,
    PRIMARY KEY (sample_id)
);

------------------------------------------
-- GEX cell embedding
------------------------------------------
DROP TABLE IF EXISTS gex_cell_embedding CASCADE;
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
);

CREATE INDEX idx_gex_cell_embedding_study ON gex_cell_embedding(study);
CREATE INDEX idx_gex_cell_embedding_sample ON gex_cell_embedding(sample_id);
CREATE INDEX idx_gex_cell_embedding_cell_type ON gex_cell_embedding(cell_type);
CREATE INDEX idx_gex_cell_embedding_cell_subtype ON gex_cell_embedding(cell_subtype);
CREATE INDEX idx_gex_cell_embedding_cluster ON gex_cell_embedding(sample_id, cluster);

------------------------------------------
-- GEX marker expression cache
------------------------------------------
DROP TABLE IF EXISTS gex_marker_expression CASCADE;
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
);

CREATE INDEX idx_gex_marker_expression_study ON gex_marker_expression(study);
CREATE INDEX idx_gex_marker_expression_sample_gene ON gex_marker_expression(sample_id, gene);
CREATE INDEX idx_gex_marker_expression_gene ON gex_marker_expression(gene);
CREATE INDEX idx_gex_marker_expression_detected ON gex_marker_expression(sample_id, gene, detected);

------------------------------------------
-- GEX full expression zarr store
------------------------------------------
DROP TABLE IF EXISTS gex_expression_store CASCADE;
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
);

CREATE INDEX idx_gex_expression_store_study ON gex_expression_store(study);
CREATE INDEX idx_gex_expression_store_sample ON gex_expression_store(sample_id);






DROP TABLE IF EXISTS import_log CASCADE;
-- 导入日志表
CREATE TABLE IF NOT EXISTS import_log (
    -- 主键
    id SERIAL PRIMARY KEY,
    
    -- 样本标识（必须）
    sample_id VARCHAR(100) NOT NULL,
    
    -- 数据类型：'bcr', 'vdj_usage', 'gex_cell_embedding', 'gex_marker_expression', 'gex_expression_store', 'clone_members', 'public_clonotypes', 'network', 'diversity'
    data_type VARCHAR(50) NOT NULL,
    
    -- 文件路径（可选）
    file_path TEXT,
    
    -- 导入状态：'success', 'failed', 'processing'
    status VARCHAR(20) NOT NULL DEFAULT 'success',
    
    -- 导入时间
    import_time TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
    
    -- 确保每个样本的每种数据类型只有一条记录
    CONSTRAINT unique_sample_data_type UNIQUE (sample_id, data_type)
);

-- 创建索引加速查询
CREATE INDEX idx_import_log_sample ON import_log(sample_id);
CREATE INDEX idx_import_log_data_type ON import_log(data_type);
CREATE INDEX idx_import_log_status ON import_log(status);
CREATE INDEX idx_import_log_time ON import_log(import_time);

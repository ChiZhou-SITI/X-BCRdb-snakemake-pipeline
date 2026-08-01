CREATE INDEX idx_cdr3
ON clonotype_member(cdr3_aa);

CREATE INDEX idx_vgene
ON clonotype_member(v_call);

CREATE INDEX idx_jgene
ON clonotype_member(j_call);

CREATE INDEX idx_sample
ON clonotype_member(sample_id);

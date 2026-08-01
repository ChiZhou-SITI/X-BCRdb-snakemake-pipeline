#!/usr/bin/env python3
"""
Load downloadable h5ad/BCR JSON file metadata into PostgreSQL.

This script is intentionally separate from the full Snakemake database loader so
Zenodo records can be added study-by-study after file publication.
"""

import argparse
import csv
import logging
import sys
from pathlib import Path

from load_database import get_db_connection


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def ensure_download_file_table(conn):
    with conn.cursor() as cur:
        cur.execute(
            """
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
            """
        )
        cur.execute(
            """
            CREATE INDEX IF NOT EXISTS idx_download_file_study
            ON download_file(study)
            """
        )
        cur.execute(
            """
            CREATE INDEX IF NOT EXISTS idx_download_file_record
            ON download_file(zenodo_record_id)
            """
        )
        cur.execute(
            """
            CREATE INDEX IF NOT EXISTS idx_download_file_kind
            ON download_file(file_kind)
            """
        )
    conn.commit()


def parse_int(value):
    value = str(value or "").strip()
    return int(value) if value else None


def load_manifest(conn, manifest_path, study=None, require_url=True):
    manifest_path = Path(manifest_path)
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")

    loaded = 0
    skipped = 0
    missing_samples = []

    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "study",
            "sample_id",
            "file_kind",
            "file_name",
            "zenodo_record_id",
            "zenodo_doi",
            "zenodo_file_url",
        }
        missing_cols = sorted(required - set(reader.fieldnames or []))
        if missing_cols:
            raise ValueError(f"Manifest missing columns: {', '.join(missing_cols)}")

        with conn.cursor() as cur:
            for row in reader:
                row_study = (row.get("study") or "").strip()
                if study and row_study != study:
                    continue

                file_url = (row.get("zenodo_file_url") or "").strip()
                if require_url and not file_url:
                    skipped += 1
                    continue

                sample_id = (row.get("sample_id") or "").strip()
                cur.execute(
                    "SELECT 1 FROM sample WHERE sample_id = %s AND study = %s",
                    (sample_id, row_study),
                )
                if cur.fetchone() is None:
                    missing_samples.append((row_study, sample_id))
                    skipped += 1
                    continue

                cur.execute(
                    """
                    INSERT INTO download_file (
                        study, sample_id, file_kind, file_name, upload_folder,
                        zenodo_part, part_index, part_total, relative_upload_path,
                        source_path, file_size_bytes, md5, zenodo_record_id,
                        zenodo_doi, file_url, updated_at
                    )
                    VALUES (
                        %(study)s, %(sample_id)s, %(file_kind)s, %(file_name)s,
                        %(upload_folder)s, %(zenodo_part)s, %(part_index)s,
                        %(part_total)s, %(relative_upload_path)s, %(source_path)s,
                        %(file_size_bytes)s, %(md5)s, %(zenodo_record_id)s,
                        %(zenodo_doi)s, %(file_url)s, now()
                    )
                    ON CONFLICT (sample_id, file_kind) DO UPDATE SET
                        study = EXCLUDED.study,
                        file_name = EXCLUDED.file_name,
                        upload_folder = EXCLUDED.upload_folder,
                        zenodo_part = EXCLUDED.zenodo_part,
                        part_index = EXCLUDED.part_index,
                        part_total = EXCLUDED.part_total,
                        relative_upload_path = EXCLUDED.relative_upload_path,
                        source_path = EXCLUDED.source_path,
                        file_size_bytes = EXCLUDED.file_size_bytes,
                        md5 = EXCLUDED.md5,
                        zenodo_record_id = EXCLUDED.zenodo_record_id,
                        zenodo_doi = EXCLUDED.zenodo_doi,
                        file_url = EXCLUDED.file_url,
                        updated_at = now()
                    """,
                    {
                        "study": row_study,
                        "sample_id": sample_id,
                        "file_kind": (row.get("file_kind") or "").strip(),
                        "file_name": (row.get("file_name") or "").strip(),
                        "upload_folder": (row.get("upload_folder") or "").strip() or None,
                        "zenodo_part": (row.get("zenodo_part") or "").strip() or None,
                        "part_index": parse_int(row.get("part_index")),
                        "part_total": parse_int(row.get("part_total")),
                        "relative_upload_path": (row.get("relative_upload_path") or "").strip() or None,
                        "source_path": (row.get("source_path") or "").strip() or None,
                        "file_size_bytes": parse_int(row.get("file_size_bytes")),
                        "md5": (row.get("md5") or "").strip() or None,
                        "zenodo_record_id": (row.get("zenodo_record_id") or "").strip() or None,
                        "zenodo_doi": (row.get("zenodo_doi") or "").strip() or None,
                        "file_url": file_url or None,
                    },
                )
                loaded += 1

    conn.commit()
    return loaded, skipped, missing_samples


def main():
    parser = argparse.ArgumentParser(description="Load Zenodo download manifest into PostgreSQL")
    parser.add_argument("--manifest", required=True, help="Manifest TSV generated for Zenodo uploads")
    parser.add_argument("--study", help="Only import rows for this study")
    parser.add_argument(
        "--allow-empty-url",
        action="store_true",
        help="Import rows even if zenodo_file_url is empty",
    )
    args = parser.parse_args()

    conn = get_db_connection()
    try:
        ensure_download_file_table(conn)
        loaded, skipped, missing_samples = load_manifest(
            conn,
            manifest_path=args.manifest,
            study=args.study,
            require_url=not args.allow_empty_url,
        )
        logger.info("Loaded %s download rows; skipped %s rows", loaded, skipped)
        if missing_samples:
            logger.warning("Skipped %s rows because sample/study was not found", len(missing_samples))
            for row_study, sample_id in missing_samples[:10]:
                logger.warning("Missing sample in DB: study=%s sample_id=%s", row_study, sample_id)
        if loaded == 0:
            sys.exit(2)
    finally:
        conn.close()


if __name__ == "__main__":
    main()

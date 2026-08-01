#!/usr/bin/env python3
"""Populate download_file checksums from published Zenodo record metadata."""

import argparse
import json
import logging
import time
from urllib.request import Request, urlopen

from load_database import get_db_connection

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def fetch_record(record_id):
    request = Request(
        f"https://zenodo.org/api/records/{record_id}",
        headers={"Accept": "application/json", "User-Agent": "X-BCRdb-checksum-sync/1.0"},
    )
    with urlopen(request, timeout=60) as response:
        return json.load(response)


def sync_checksums(conn, force=False, pause=0.1):
    condition = "" if force else "AND NULLIF(TRIM(md5), '') IS NULL"
    with conn.cursor() as cur:
        cur.execute(
            f"""
            SELECT DISTINCT zenodo_record_id
            FROM download_file
            WHERE NULLIF(TRIM(zenodo_record_id), '') IS NOT NULL
              {condition}
            ORDER BY zenodo_record_id
            """
        )
        record_ids = [row[0] for row in cur.fetchall()]

    updated = 0
    unmatched = []
    for index, record_id in enumerate(record_ids, 1):
        payload = fetch_record(record_id)
        files = {
            item.get("key"): item
            for item in payload.get("files", [])
            if item.get("key")
        }
        with conn.cursor() as cur:
            cur.execute(
                "SELECT sample_id, file_kind, file_name FROM download_file WHERE zenodo_record_id = %s",
                (record_id,),
            )
            rows = cur.fetchall()
            for sample_id, file_kind, file_name in rows:
                item = files.get(file_name)
                if not item:
                    unmatched.append((record_id, file_name))
                    continue
                checksum = str(item.get("checksum") or "").strip()
                if checksum.lower().startswith("md5:"):
                    checksum = checksum[4:]
                cur.execute(
                    """
                    UPDATE download_file
                    SET md5 = %s,
                        file_size_bytes = COALESCE(%s, file_size_bytes),
                        updated_at = now()
                    WHERE sample_id = %s AND file_kind = %s
                    """,
                    (checksum or None, item.get("size"), sample_id, file_kind),
                )
                updated += cur.rowcount
        conn.commit()
        logger.info("[%s/%s] Synced Zenodo record %s", index, len(record_ids), record_id)
        if pause:
            time.sleep(pause)

    return len(record_ids), updated, unmatched


def main():
    parser = argparse.ArgumentParser(description="Sync Zenodo MD5 checksums into download_file")
    parser.add_argument("--force", action="store_true", help="Refresh checksums that are already populated")
    parser.add_argument("--pause", type=float, default=0.1, help="Delay between Zenodo API requests")
    args = parser.parse_args()

    conn = get_db_connection()
    try:
        records, updated, unmatched = sync_checksums(conn, force=args.force, pause=max(0, args.pause))
        logger.info("Synced %s files from %s Zenodo records", updated, records)
        if unmatched:
            logger.warning("Could not match %s database files to Zenodo metadata", len(unmatched))
            for record_id, file_name in unmatched[:20]:
                logger.warning("Unmatched: record=%s file=%s", record_id, file_name)
    finally:
        conn.close()


if __name__ == "__main__":
    main()

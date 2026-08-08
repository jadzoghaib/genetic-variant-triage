"""Database access, schema initialisation, and provenance helpers."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import duckdb

ROOT = Path(__file__).parent
DB_PATH = ROOT / "locus.duckdb"
SCHEMA = ROOT / "schema.sql"

ACTOR_INGEST = "ingest"


def connect(path: Path | str = DB_PATH, read_only: bool = False) -> duckdb.DuckDBPyConnection:
    """Open the store.

    `read_only=True` lets the UI attach while an ingest holds the file, and
    guarantees a rendering surface can never mutate the ontology.
    """
    return duckdb.connect(str(path), read_only=read_only)


def init_schema(con: duckdb.DuckDBPyConnection) -> None:
    """Apply schema.sql.

    CREATE TYPE has no IF NOT EXISTS in DuckDB, so statements are applied one at
    a time and 'already exists' is tolerated. Any other error propagates.
    """
    raw = SCHEMA.read_text(encoding="utf-8")
    # Strip line comments before splitting on ';' — prose in the comments
    # contains semicolons and would otherwise cut statements in half.
    sql = "\n".join(line.split("--", 1)[0] for line in raw.splitlines())
    for stmt in (s.strip() for s in sql.split(";")):
        if not stmt:
            continue
        try:
            con.execute(stmt)
        except duckdb.CatalogException as exc:
            if "already exists" not in str(exc):
                raise


# ────────────────────────────────────────────────────────────── provenance ──

def record_evidence(
    con: duckdb.DuckDBPyConnection,
    source_system: str,
    resource_url: str,
    payload: bytes,
    source_version: str | None = None,
    record_count: int | None = None,
) -> str:
    """Register one retrieval and return its evidence_id.

    The id is derived from the source, the URL and the retrieval instant, so a
    re-fetch produces a new row rather than overwriting history — evidence is a
    log, not current state.
    """
    retrieved_at = datetime.now(timezone.utc).replace(tzinfo=None)
    evidence_id = hashlib.sha256(
        f"{source_system}|{resource_url}|{retrieved_at.isoformat()}".encode()
    ).hexdigest()[:16]
    con.execute(
        "INSERT INTO evidence VALUES (?, ?, ?, ?, ?, ?, ?)",
        [
            evidence_id,
            source_system,
            resource_url,
            source_version,
            retrieved_at,
            hashlib.sha256(payload).hexdigest(),
            record_count,
        ],
    )
    return evidence_id


def emit_event(
    con: duckdb.DuckDBPyConnection,
    event_type: str,
    object_type: str,
    object_key: str,
    actor: str = ACTOR_INGEST,
    payload: dict | None = None,
) -> None:
    con.execute(
        "INSERT INTO event VALUES (nextval('event_seq'), ?, ?, ?, ?, ?, ?)",
        [
            event_type,
            object_type,
            object_key,
            datetime.now(timezone.utc).replace(tzinfo=None),
            actor,
            json.dumps(payload) if payload else None,
        ],
    )


# ───────────────────────────────────────────────────────────────── loading ──

def clear_target(con: duckdb.DuckDBPyConnection, uniprot_acc: str) -> None:
    """Remove a target's current state so it can be reloaded.

    Deletion runs child-first to satisfy foreign keys. Evidence and events are
    deliberately untouched: they are append-only history of what was retrieved
    and what happened, not current state.

    Variants are shared objects — a variant is only removed once no target still
    claims it, so reloading one gene cannot silently delete another's data.
    """
    con.execute("""
        DELETE FROM triage_decision WHERE variant_id IN (
            SELECT variant_id FROM target_variant WHERE uniprot_acc = ?
        )
    """, [uniprot_acc])

    # dossier layer, child-first
    for sql in (
        "DELETE FROM target_drug_report        WHERE uniprot_acc = ?",
        "DELETE FROM target_drug               WHERE uniprot_acc = ?",
        "DELETE FROM association_datatype_score WHERE uniprot_acc = ?",
        "DELETE FROM association               WHERE uniprot_acc = ?",
        "DELETE FROM tractability              WHERE uniprot_acc = ?",
        "DELETE FROM target_prioritisation     WHERE uniprot_acc = ?",
        "DELETE FROM structure_coverage        WHERE uniprot_acc = ?",
    ):
        con.execute(sql, [uniprot_acc])

    con.execute("DELETE FROM variant_effect WHERE uniprot_acc = ?", [uniprot_acc])
    con.execute("DELETE FROM target_variant WHERE uniprot_acc = ?", [uniprot_acc])
    con.execute("""
        DELETE FROM clinical_assertion
        WHERE variant_id NOT IN (SELECT variant_id FROM target_variant)
    """)
    con.execute("""
        DELETE FROM variant
        WHERE variant_id NOT IN (SELECT variant_id FROM target_variant)
    """)
    con.execute("DELETE FROM residue   WHERE uniprot_acc = ?", [uniprot_acc])
    con.execute("DELETE FROM structure WHERE uniprot_acc = ?", [uniprot_acc])
    con.execute("DELETE FROM target    WHERE uniprot_acc = ?", [uniprot_acc])

    # Shared objects: drop only once no target still references them, so
    # reloading one gene cannot delete another's diseases or drugs.
    for sql in (
        "DELETE FROM drug_mechanism WHERE chembl_id NOT IN (SELECT chembl_id FROM target_drug)",
        "DELETE FROM drug           WHERE chembl_id NOT IN (SELECT chembl_id FROM target_drug)",
        "DELETE FROM clinical_report WHERE report_id NOT IN (SELECT report_id FROM target_drug_report)",
        "DELETE FROM disease_therapeutic_area WHERE efo_id NOT IN (SELECT efo_id FROM association)",
        "DELETE FROM disease        WHERE efo_id NOT IN (SELECT efo_id FROM association)",
        "DELETE FROM therapeutic_area WHERE ta_id NOT IN (SELECT ta_id FROM disease_therapeutic_area)",
    ):
        con.execute(sql)


def insert_df(
    con: duckdb.DuckDBPyConnection,
    table: str,
    df,
    columns: list[str],
    or_ignore: bool = False,
) -> int:
    """Insert a dataframe into `table`, selecting `columns` in schema order.

    `or_ignore` is for objects shared between targets — diseases, drugs, trial
    reports — where a second target legitimately re-supplies a row that already
    exists.
    """
    if df is None or len(df) == 0:
        return 0
    con.register("_staging", df[columns])
    verb = "INSERT OR IGNORE INTO" if or_ignore else "INSERT INTO"
    con.execute(f"{verb} {table} SELECT {', '.join(columns)} FROM _staging")
    con.unregister("_staging")
    return len(df)


def insert_records(
    con: duckdb.DuckDBPyConnection,
    table: str,
    records: list[dict],
    columns: list[str],
    key: list[str],
    evidence_id: str | None = None,
    or_ignore: bool = False,
) -> int:
    """Deduplicate flat records on `key`, stamp evidence, and insert.

    Connectors emit records positionally as they walk a nested payload, so the
    same disease or drug appears many times. Deduplication belongs here rather
    than in every connector.
    """
    import pandas as pd

    if not records:
        return 0
    df = pd.DataFrame(records).drop_duplicates(subset=key)
    if evidence_id is not None:
        df = df.assign(evidence_id=evidence_id)
    return insert_df(con, table, df, columns, or_ignore=or_ignore)

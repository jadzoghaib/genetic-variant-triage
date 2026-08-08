"""Load an exported decision log back into the ontology.

The console holds decisions in browser storage, which is fine for working but
is not a record. This closes the loop: an exported log becomes rows in
`triage_decision` (empty since Phase 1) plus target status changes and events,
so analyst judgements live in the same store, with the same provenance
discipline, as the facts they were made about.

    uv run python import_decisions.py locus-decisions-2026-08-08.json
    uv run python import_decisions.py --dry-run <file>

Import is idempotent: decisions are keyed by the id the browser assigned, so
re-importing the same file changes nothing.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import locus_db as db

TARGET_STATUS = {"shortlisted": "shortlisted", "deprioritized": "deprioritized"}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("path", type=Path)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    payload = json.loads(args.path.read_text(encoding="utf-8"))
    if payload.get("kind") != "locus.decisions":
        sys.exit(f"{args.path} is not a Locus decision log")

    decisions = payload.get("decisions", [])
    print(f"{len(decisions)} decision(s), exported {payload.get('exported_at')}, "
          f"against data build {payload.get('data_build')}")

    con = db.connect(read_only=args.dry_run)
    db.init_schema(con) if not args.dry_run else None

    known = {r[0] for r in con.execute(
        "SELECT decision_id FROM triage_decision").fetchall()}
    build = con.execute(
        "SELECT MAX(retrieved_at) FROM evidence").fetchone()[0]

    stats = {"variant_new": 0, "variant_skipped": 0, "target": 0,
             "target_skipped": 0, "orphan": 0, "stale": 0}

    for d in sorted(decisions, key=lambda x: x["decided_at"]):
        if payload.get("data_build") and d.get("data_build") \
                and d["data_build"] != payload["data_build"]:
            stats["stale"] += 1

        if d["object_type"] == "variant":
            if d["id"] in known:
                stats["variant_skipped"] += 1
                continue
            # A decision about a variant the store has never seen is a real
            # inconsistency, not something to insert quietly.
            exists = con.execute("SELECT 1 FROM variant WHERE variant_id = ?",
                                 [d["object_key"]]).fetchone()
            if not exists:
                print(f"  ORPHAN: {d['object_key']} is not in this store")
                stats["orphan"] += 1
                continue
            if args.dry_run:
                stats["variant_new"] += 1
                continue
            con.execute(
                """INSERT INTO triage_decision
                   (decision_id, variant_id, analyst, prior_class, new_class,
                    rationale, decided_at, superseded_by)
                   VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
                [d["id"], d["object_key"], d["analyst"], d.get("prior_class"),
                 d["outcome"], d["rationale"], d["decided_at"],
                 d.get("superseded_by")])
            # The event's own timestamp is the load time, so the moment the
            # analyst actually decided travels in the payload — otherwise a
            # re-import silently restates an old judgement as a fresh one.
            db.emit_event(con, "variant.reviewed", "variant", d["object_key"],
                          actor=d["analyst"],
                          payload={"outcome": d["outcome"],
                                   "decision_id": d["id"],
                                   "decided_at": d["decided_at"],
                                   "data_build": d.get("data_build")})
            stats["variant_new"] += 1

        elif d["object_type"] == "target":
            status = TARGET_STATUS.get(d["outcome"])
            if not status or d.get("superseded_by"):
                continue          # only the live decision sets current status
            # Target decisions have no table of their own, so idempotency has to
            # be judged from the event log AND the effect. Checking the event
            # alone is not enough: `ingest.py` rebuilds target rows and resets
            # status, while events are append-only and survive. A re-import
            # would then see its own old event, conclude "already applied", and
            # silently leave the target un-shortlisted — a governed decision
            # dropped with no error. Skip only when the intended state is
            # actually in place.
            logged = con.execute(
                "SELECT COUNT(*) FROM event WHERE object_type = 'target' "
                "AND payload LIKE ?", [f'%"decision_id": "{d["id"]}"%']
            ).fetchone()[0]
            row = con.execute("SELECT status FROM target WHERE uniprot_acc = ?",
                              [d["object_key"]]).fetchone()
            if not row:
                print(f"  ORPHAN: target {d['object_key']} is not in this store")
                stats["orphan"] += 1
                continue
            if logged and row[0] == status:
                stats["target_skipped"] += 1
                continue
            if args.dry_run:
                stats["target"] += 1
                continue
            con.execute("UPDATE target SET status = ? WHERE uniprot_acc = ?",
                        [status, d["object_key"]])
            db.emit_event(con, f"target.{status}", "target", d["object_key"],
                          actor=d["analyst"],
                          payload={"rationale": d["rationale"],
                                   "decision_id": d["id"]})
            stats["target"] += 1

    if not args.dry_run:
        con.execute("CHECKPOINT")

    print(f"\n  variant decisions inserted : {stats['variant_new']}")
    print(f"  already present            : {stats['variant_skipped']}")
    print(f"  target status changes      : {stats['target']} "
          f"({stats['target_skipped']} already applied)")
    print(f"  orphaned (unknown variant) : {stats['orphan']}")
    if stats["stale"]:
        print(f"  taken against another build: {stats['stale']}")
    print(f"  store last refreshed       : {build}")
    if args.dry_run:
        print("\n  dry run — nothing written")

    con.close()


if __name__ == "__main__":
    main()

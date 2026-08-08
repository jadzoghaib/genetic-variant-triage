"""RCSB PDB connector — experimental structure coverage for a UniProt accession.

Two calls: the Search API resolves an accession to polymer entity ids, then the
Data API's GraphQL endpoint returns method, resolution and — the part that
matters — the UniProt-coordinate spans each entity actually covers.

Those spans are what let the confidence model distinguish "experimentally
solved" from "confidently predicted". EGFR shows why it matters: 1M14 covers
residues 695-1022 (the kinase domain) while 1IVO covers 25-646 (the
ectodomain). Coverage is domain-specific, never whole-protein.
"""

from __future__ import annotations

import httpx

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL = "https://data.rcsb.org/graphql"
SOURCE = "rcsb_pdb"

# Bounds the work for heavily-studied targets (EGFR has 392 entities, TP53 321).
# The true total is stored on the target so the dossier can report it.
MAX_ENTITIES = 250
BATCH = 50

ACCESSION_ATTR = ("rcsb_polymer_entity_container_identifiers"
                  ".reference_sequence_identifiers.database_accession")

DETAIL_QUERY = """
query E($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    rcsb_polymer_entity_align {
      reference_database_accession
      aligned_regions { ref_beg_seq_id length }
    }
    entry {
      rcsb_id
      struct { title }
      exptl { method }
      rcsb_entry_info { resolution_combined }
    }
  }
}
"""


def search_entities(acc: str) -> tuple[list[str], int]:
    """Polymer entity ids referencing this accession, plus the true total."""
    query = {
        "query": {"type": "terminal", "service": "text", "parameters": {
            "attribute": ACCESSION_ATTR, "operator": "exact_match", "value": acc}},
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": 0, "rows": MAX_ENTITIES},
            "results_content_type": ["experimental"],
        },
    }
    r = httpx.post(SEARCH, json=query, timeout=180)
    if r.status_code == 204:          # no hits at all
        return [], 0
    r.raise_for_status()
    payload = r.json()
    ids = [h["identifier"] for h in payload.get("result_set", [])]
    return ids, int(payload.get("total_count", len(ids)))


def fetch_details(entity_ids: list[str], acc: str) -> tuple[list[dict], list[dict], bytes]:
    """Return (structure rows, coverage rows, raw payload) for the given entities."""
    structures: list[dict] = []
    coverage: list[dict] = []
    raw = bytearray()

    for i in range(0, len(entity_ids), BATCH):
        chunk = entity_ids[i:i + BATCH]
        r = httpx.post(GRAPHQL, json={"query": DETAIL_QUERY, "variables": {"ids": chunk}},
                       timeout=300)
        r.raise_for_status()
        raw += r.content
        for e in (r.json().get("data") or {}).get("polymer_entities") or []:
            entry = e.get("entry") or {}
            info = entry.get("rcsb_entry_info") or {}
            methods = [m["method"] for m in (entry.get("exptl") or []) if m.get("method")]
            res = info.get("resolution_combined") or []
            structures.append({
                "structure_id": e["rcsb_id"],
                "uniprot_acc": acc,
                "source": "pdb",
                "model_version": None,
                "method": methods[0] if methods else None,
                "resolution": float(res[0]) if res else None,
                "title": ((entry.get("struct") or {}).get("title") or "")[:200],
            })
            # Keep only alignments against the accession we asked for: an entity
            # can align to several reference sequences.
            for al in e.get("rcsb_polymer_entity_align") or []:
                if al.get("reference_database_accession") != acc:
                    continue
                for reg in al.get("aligned_regions") or []:
                    beg, length = reg.get("ref_beg_seq_id"), reg.get("length")
                    if beg is None or not length:
                        continue
                    coverage.append({
                        "structure_id": e["rcsb_id"], "uniprot_acc": acc,
                        "ref_beg": int(beg), "ref_end": int(beg) + int(length) - 1,
                    })

    return structures, coverage, bytes(raw)

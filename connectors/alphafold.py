"""AlphaFold DB connector — structure metadata, per-residue pLDDT, and the
AlphaMissense predictions that AlphaFold DB serves alongside them.

Every fetch returns (parsed, raw_bytes, url) so the caller can stamp evidence.
Responses are cached on disk so a demo runs without network.
"""

from __future__ import annotations

import json
from pathlib import Path

import httpx
import pandas as pd

API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"
CACHE = Path(__file__).parent.parent / "data"

SOURCE = "alphafold_db"

# AlphaMissense ships two files per protein with different class vocabularies:
# aa-substitutions.csv uses LPath/LBen/Amb, hg38.csv spells them out. A filter
# written for one silently matches nothing in the other.
AM_CLASS_CANON = {
    "likely_pathogenic": "LPath",
    "likely_benign": "LBen",
    "ambiguous": "Amb",
    "LPath": "LPath",
    "LBen": "LBen",
    "Amb": "Amb",
}


def _cached(url: str, name: str) -> bytes:
    CACHE.mkdir(exist_ok=True)
    path = CACHE / name
    if not path.exists():
        r = httpx.get(url, timeout=180)
        r.raise_for_status()
        path.write_bytes(r.content)
    return path.read_bytes()


def entry(acc: str) -> tuple[dict, bytes, str]:
    """Canonical AlphaFold entry for an accession.

    The API returns one entry per isoform — BRCA1 returns eight. The canonical
    model is the one whose modelEntityId is exactly AF-{acc}-F1. File URLs must
    be read from the response, never constructed: the DB moved v4 -> v6 and
    hand-built v4 URLs now 404.
    """
    url = API.format(acc=acc)
    raw = _cached(url, f"af_entry_{acc}.json")
    entries = json.loads(raw)
    want = f"AF-{acc}-F1"
    for e in entries:
        if e.get("modelEntityId") == want:
            return e, raw, url
    raise RuntimeError(
        f"no canonical entry {want}; got {[e.get('modelEntityId') for e in entries]}"
    )


def sequence(e: dict) -> str:
    seq = e.get("uniprotSequence") or e.get("sequence")
    if not seq:
        raise RuntimeError(f"no sequence in entry {e.get('modelEntityId')}")
    return seq


def plddt(e: dict, symbol: str) -> tuple[pd.DataFrame, bytes, str]:
    url = e["plddtDocUrl"]
    raw = _cached(url, f"plddt_{symbol}.json")
    d = json.loads(raw)
    df = pd.DataFrame({"position": d["residueNumber"], "plddt": d["confidenceScore"]})
    return df, raw, url


def alphamissense(e: dict, symbol: str) -> tuple[pd.DataFrame, bytes, str]:
    """Per-protein AlphaMissense predictions in hg38 genomic coordinates.

    The hg38 file carries CHROM/POS/REF/ALT, which is the join key to ClinVar,
    and transcript_id, which is the grain of the prediction.
    """
    url = e["amAnnotationsHg38Url"]
    raw = _cached(url, f"am_{symbol}_hg38.csv")
    df = pd.read_csv(CACHE / f"am_{symbol}_hg38.csv")

    unknown = set(df["am_class"].unique()) - set(AM_CLASS_CANON)
    if unknown:
        raise RuntimeError(f"unrecognised am_class values for {symbol}: {unknown}")
    df["am_class"] = df["am_class"].map(AM_CLASS_CANON)

    df["chrom"] = df["CHROM"].str.removeprefix("chr")
    df = df.rename(columns={"POS": "pos", "REF": "ref", "ALT": "alt"})

    parts = df["protein_variant"].str.extract(r"^(?P<ref_aa>[A-Z])(?P<aa_pos>\d+)(?P<alt_aa>[A-Z])$")
    if parts["aa_pos"].isna().any():
        bad = df.loc[parts["aa_pos"].isna(), "protein_variant"].unique()[:5]
        raise RuntimeError(f"unparseable protein_variant values for {symbol}: {bad}")
    df["ref_aa"] = parts["ref_aa"]
    df["alt_aa"] = parts["alt_aa"]
    df["aa_pos"] = parts["aa_pos"].astype(int)

    return df, raw, url

"""Build the static site.

The console and dossier are static: nothing in this product changes at runtime,
so a server earns nothing. This script runs the same `core/` rules the tests
cover, serialises the result, and the browser only renders — which is exactly
what the pure-core architecture was built to allow.

Payloads are columnar with dictionary-encoded categoricals. An array of 12,463
row objects is mostly repeated key names; columns plus small integer codes cut
it by roughly 5x before the server ever gzips it.

    uv run python export_site.py
    uv run python -m http.server 8080 -d site
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import httpx
import pandas as pd

import queries
from core import confidence as C
from core import dossier as D
from core import triage as T
from ui import structure as S

ROOT = Path(__file__).parent
SITE = ROOT / "site"
DATA = SITE / "data"
VENDOR = SITE / "vendor"

# Vendored rather than loaded from a CDN so the site is self-contained: it runs
# from a local directory, from GitHub Pages, or offline, with no third-party
# request at view time.
MOL3D_URL = "https://unpkg.com/3dmol@2.4.2/build/3Dmol-min.js"


def _round(series, places: int):
    return [None if pd.isna(v) else round(float(v), places) for v in series]


def records(df: pd.DataFrame) -> list[dict]:
    """DataFrame -> JSON-safe records.

    `df.where(df.notna(), None)` looks like it nulls missing values but does not:
    assigning None into a float column stores NaN, and `NaN` is not valid JSON —
    the browser rejects the whole file. Casting to object first makes the None
    survive.
    """
    return df.astype(object).where(df.notna(), None).to_dict("records")


def encode(values, none_code: int = -1) -> tuple[list, list[int]]:
    """Dictionary-encode a categorical column -> (levels, codes)."""
    levels: list = []
    index: dict = {}
    codes: list[int] = []
    for v in values:
        if v is None or (isinstance(v, float) and v != v):
            codes.append(none_code)
            continue
        if v not in index:
            index[v] = len(levels)
            levels.append(v)
        codes.append(index[v])
    return levels, codes


def build_variants(df: pd.DataFrame) -> dict:
    sig_levels, sig_codes = encode(df["significance"])
    raw_levels, raw_codes = encode(df["significance_raw"])
    tri_levels, tri_codes = encode(df["triage_class"])
    pri_levels, pri_codes = encode(df["review_priority"])
    tier_levels, tier_codes = encode(df["evidence_tier"])
    amc_levels, amc_codes = encode(df["am_class"])
    # Measured: the rationale strings were 34% of the whole payload (374 KB for
    # BRCA1) across just 45 distinct values, because the same four factors
    # recombine. Encoding them is the single largest saving available.
    rea_levels, rea_codes = encode(
        [r if r else None for r in df["priority_reasons"]])

    return {
        "n": int(len(df)),
        "levels": {"sig": sig_levels, "raw": raw_levels, "triage": tri_levels,
                   "priority": pri_levels, "tier": tier_levels, "amc": amc_levels,
                   "reasons": rea_levels},
        "cols": {
            "pv": df["protein_variant"].tolist(),
            "vid": df["variant_id"].tolist(),
            "pos": [int(p) for p in df["aa_pos"]],
            "am": _round(df["am_pathogenicity"], 4),
            "plddt": _round(df["plddt"], 1),
            "stars": [-1 if pd.isna(s) else int(s) for s in df["stars"]],
            "solved": [1 if bool(s) else 0 for s in df["is_solved"]],
            "sig": sig_codes, "raw": raw_codes, "triage": tri_codes,
            "priority": pri_codes, "tier": tier_codes, "amc": amc_codes,
            "reasons": rea_codes,
        },
        "transcript": df["transcript_id"].iloc[0] if len(df) else None,
    }


def build_profile(prof: pd.DataFrame) -> dict:
    """Per-residue profile, with the tier and class already decided.

    The browser receives codes, never thresholds: `evidence_tier` comes from
    core.confidence and the class from AlphaMissense itself, so the 70 pLDDT
    cut and the 0.34/0.564 cuts exist in exactly one place each. Logic
    duplicated across a language boundary is logic that will eventually
    disagree with itself.
    """
    tiers = [C.evidence_tier(bool(s), float(p))
             for s, p in zip(prof["is_solved"], prof["plddt"])]
    tier_levels, tier_codes = encode(tiers)
    amc_levels, amc_codes = encode(prof["max_am_class"])

    return {
        "n": int(len(prof)),
        "wt": "".join(prof["wt_aa"].tolist()),
        "plddt": _round(prof["plddt"], 1),
        "maxam": _round(prof["max_am"], 3),
        "levels": {"tier": tier_levels, "amc": amc_levels},
        "tier": tier_codes,
        "amc": amc_codes,
        "runs": S.tier_runs(prof),
    }


def build_dossier(con, symbol: str, acc: str, facts: dict, triaged: pd.DataFrame) -> dict:
    facts = dict(facts)
    # `symbol` is the index of target_facts, not a column, so it is absent from
    # the row dict and the card would render "None: ...".
    facts["symbol"] = symbol
    for k in ("max_genetic_score", "pct_solved", "global_plddt"):
        if facts.get(k) is not None and not pd.isna(facts[k]):
            facts[k] = round(float(facts[k]), 3)
    facts["drug_stages"] = queries.drug_stages(con, acc)
    facts["n_upgrade_candidates"] = int((triaged.triage_class == T.UPGRADE).sum())
    facts["n_discordant"] = int((triaged.triage_class == T.DISCORDANT).sum())
    card = D.build(facts)

    prio = queries.prioritisation(con, acc)
    return {
        "card": card,
        "summary": D.summarise(card),
        "associations": records(queries.associations(con, acc)),
        "drugs": records(queries.drugs(con, acc)),
        "structures": records(queries.experimental_structures(con, acc)),
        "tractability": records(queries.tractability(con, acc)),
        "prioritisation": {r.metric_key: (None if pd.isna(r.metric_value)
                                          else round(float(r.metric_value), 3))
                           for r in prio.itertuples()},
    }


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    VENDOR.mkdir(parents=True, exist_ok=True)

    con = queries.connect()
    facts_df = queries.target_facts(con).set_index("symbol")

    manifest = {
        "generated_at": datetime.now(timezone.utc).replace(microsecond=0).isoformat(),
        "targets": [],
    }

    for symbol, facts in facts_df.iterrows():
        acc = facts["uniprot_acc"]
        print(f"  {symbol} ({acc})")

        triaged = T.assign(C.assign(queries.worklist(con, symbol)))
        prof = queries.residue_profile(con, acc)
        struct = queries.predicted_structure(con, acc)

        # allow_nan=False: Python happily emits bare `NaN`, which is not valid
        # JSON and which the browser rejects for the whole file. Fail here
        # instead of shipping a payload that cannot be parsed.
        dump = dict(separators=(",", ":"), allow_nan=False)
        (DATA / f"variants_{symbol}.json").write_text(
            json.dumps(build_variants(triaged), **dump))
        (DATA / f"profile_{symbol}.json").write_text(
            json.dumps(build_profile(prof), **dump))
        (DATA / f"dossier_{symbol}.json").write_text(
            json.dumps(build_dossier(con, symbol, acc, facts.to_dict(), triaged),
                       default=str, **dump))

        # Backbone only: a cartoon needs nothing else, and it is ~35% smaller.
        pdb = S.backbone_only(S.fetch_structure(struct["file_url"],
                                                struct["structure_id"]))
        (DATA / f"{struct['structure_id']}.pdb").write_text(pdb)

        actionable = triaged[triaged.triage_class.isin(T.ACTIONABLE)]
        manifest["targets"].append({
            "symbol": symbol,
            "acc": acc,
            "name": facts["approved_name"],
            "structure_id": struct["structure_id"],
            "kpi": {
                "predictions": int(len(triaged)),
                "actionable": int(len(actionable)),
                "high": int((actionable.review_priority == "high").sum()),
                "upgrades": int((triaged.triage_class == T.UPGRADE).sum()),
                "pct_solved": round(100.0 * prof["is_solved"].mean(), 1),
                "global_plddt": round(float(facts["global_plddt"]), 1),
            },
        })

    (DATA / "manifest.json").write_text(json.dumps(manifest, indent=1))

    vendor_js = VENDOR / "3Dmol-min.js"
    if not vendor_js.exists():
        print("  fetching 3Dmol...")
        vendor_js.write_bytes(httpx.get(MOL3D_URL, timeout=300).content)

    con.close()

    total = sum(f.stat().st_size for f in SITE.rglob("*") if f.is_file())
    print(f"\n  site: {total / 1048576:.1f} MB")
    for f in sorted(SITE.rglob("*")):
        if f.is_file() and f.stat().st_size > 120_000:
            print(f"    {f.relative_to(SITE)}  {f.stat().st_size / 1024:.0f} KB")


if __name__ == "__main__":
    main()

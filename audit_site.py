"""Audit the exported site against the store it was built from.

An export is a second copy of the truth, and a second copy can drift. This
re-derives every headline number straight from DuckDB and compares it to what
the site will actually display, so a stale or half-written export fails loudly
instead of showing plausible wrong figures.

    uv run python audit_site.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import queries
from core import confidence as C
from core import triage as T

SITE = Path(__file__).parent / "site"
DATA = SITE / "data"

problems: list[str] = []
checks = 0


def check(label: str, ok: bool, detail: str = "") -> None:
    global checks
    checks += 1
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{(' — ' + detail) if detail else ''}")
    if not ok:
        problems.append(label)


def main() -> None:
    con = queries.connect()
    manifest = json.loads((DATA / "manifest.json").read_text())
    symbols = [t["symbol"] for t in manifest["targets"]]

    print("=" * 78)
    print("1. STRUCTURAL COMPLETENESS")
    print("=" * 78)
    db_symbols = [r[0] for r in con.execute(
        "SELECT symbol FROM target ORDER BY symbol").fetchall()]
    check("manifest lists every target in the store", sorted(symbols) == db_symbols,
          f"{sorted(symbols)} vs {db_symbols}")

    for name in ["index.html", "app.js", "view-model.js", "decisions.js",
                 "styles.css", "vendor/3Dmol-min.js"]:
        check(f"{name} present", (SITE / name).exists())

    for t in manifest["targets"]:
        s = t["symbol"]
        for f in [f"variants_{s}.json", f"profile_{s}.json", f"dossier_{s}.json",
                  f"{t['structure_id']}.pdb"]:
            check(f"{f} present", (DATA / f).exists())

    print()
    print("=" * 78)
    print("2. EVERY DISPLAYED NUMBER RE-DERIVED FROM THE STORE")
    print("=" * 78)
    for t in manifest["targets"]:
        s, acc = t["symbol"], t["acc"]
        df = T.assign(C.assign(queries.worklist(con, s)))
        prof = queries.residue_profile(con, acc)
        act = df[df.triage_class.isin(T.ACTIONABLE)]
        k = t["kpi"]
        expect = {
            "predictions": len(df),
            "actionable": len(act),
            "high": int((act.review_priority == "high").sum()),
            "upgrades": int((df.triage_class == T.UPGRADE).sum()),
            "pct_solved": round(100.0 * prof["is_solved"].mean(), 1),
        }
        for key, want in expect.items():
            check(f"{s} KPI {key}", k[key] == want, f"site={k[key]} store={want}")

    print()
    print("=" * 78)
    print("3. PAYLOAD INTEGRITY")
    print("=" * 78)
    for t in manifest["targets"]:
        s, acc = t["symbol"], t["acc"]
        v = json.loads((DATA / f"variants_{s}.json").read_text())
        p = json.loads((DATA / f"profile_{s}.json").read_text())
        d = json.loads((DATA / f"dossier_{s}.json").read_text())

        n_db = con.execute(
            "SELECT COUNT(*) FROM variant_effect ve JOIN target t "
            "ON t.uniprot_acc = ve.uniprot_acc WHERE t.symbol = ?", [s]).fetchone()[0]
        check(f"{s} variant count", v["n"] == n_db, f"site={v['n']} store={n_db}")
        check(f"{s} every column is full length",
              all(len(c) == v["n"] for c in v["cols"].values()))
        check(f"{s} no categorical code escapes its dictionary",
              all(max(v["cols"][k], default=-1) < len(v["levels"][k])
                  for k in v["levels"]))

        n_res = con.execute("SELECT COUNT(*) FROM residue WHERE uniprot_acc = ?",
                            [acc]).fetchone()[0]
        check(f"{s} residue count", p["n"] == n_res, f"site={p['n']} store={n_res}")
        check(f"{s} sequence length matches residue count", len(p["wt"]) == p["n"])
        check(f"{s} profile arrays aligned",
              all(len(p[k]) == p["n"] for k in ("plddt", "maxam", "tier", "amc")))
        check(f"{s} tier runs span the whole protein",
              p["runs"][0]["start"] == 1 and p["runs"][-1]["end"] >= p["n"])

        card = d["card"]
        check(f"{s} dossier names its target", card["symbol"] == s)
        check(f"{s} dossier has all five dimensions",
              len(card["dimensions"]) == 5)

    print()
    print("=" * 78)
    print("4. CLAIMS THE DEMO MAKES OUT LOUD")
    print("=" * 78)
    arche = {t["symbol"]: json.loads(
        (DATA / f"dossier_{t['symbol']}.json").read_text())["card"]["archetype"]
        for t in manifest["targets"]}
    check("EGFR reads as a validated druggable target",
          arche.get("EGFR") == "validated druggable target", arche.get("EGFR", "?"))
    check("BRCA1 reads as genetically validated but undrugged",
          arche.get("BRCA1") == "genetically validated, chemically unexplored",
          arche.get("BRCA1", "?"))

    egfr = json.loads((DATA / "dossier_EGFR.json").read_text())
    brca = json.loads((DATA / "dossier_BRCA1.json").read_text())
    check("EGFR has chemical matter", len(egfr["drugs"]) > 0, f"{len(egfr['drugs'])} drugs")
    check("BRCA1 has none — the empty state is real, not a rendering bug",
          len(brca["drugs"]) == 0)
    check("BRCA1's experimental structures all cover the BRCT domain only",
          all(x["coverage_start"] and x["coverage_start"] > 1500
              for x in brca["structures"]),
          f"{len(brca['structures'])} structures")

    # The L858R story the spec leads with.
    row = con.execute("""
        SELECT ve.am_pathogenicity, r.plddt,
               EXISTS (SELECT 1 FROM structure_coverage sc
                       WHERE sc.uniprot_acc = ve.uniprot_acc
                         AND ve.aa_pos BETWEEN sc.ref_beg AND sc.ref_end) AS solved
        FROM variant_effect ve
        JOIN residue r ON r.uniprot_acc = ve.uniprot_acc AND r.position = ve.aa_pos
        JOIN target t ON t.uniprot_acc = ve.uniprot_acc AND t.symbol = 'EGFR'
        WHERE ve.protein_variant = 'L858R'""").fetchone()
    check("L858R is still below the pLDDT confidence cut", row[1] < 70, f"pLDDT {row[1]}")
    check("L858R is still rescued by experimental coverage", bool(row[2]))
    check("L858R is still called pathogenic", row[0] > 0.564, f"{row[0]:.4f}")

    print()
    print("=" * 78)
    print(f"{checks} checks · " + ("ALL PASS" if not problems
                                   else f"{len(problems)} FAILED: {problems}"))
    print("=" * 78)
    con.close()
    sys.exit(1 if problems else 0)


if __name__ == "__main__":
    main()

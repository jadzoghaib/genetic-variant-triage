# /// script
# requires-python = ">=3.11"
# dependencies = ["httpx", "pandas", "duckdb"]
# ///
"""
Locus Phase 0 — the join spike.

Question this answers, before any build investment:
  1. Does AlphaMissense (via AlphaFold DB, hg38 coords) actually join to ClinVar
     on (CHROM, POS, REF, ALT)? At what coverage rate?
  2. Does AlphaMissense AGREE with ClinVar where ClinVar is confident? If not,
     the reclassification premise is dead and the product changes shape.
  3. How many VUS does AlphaMissense confidently call pathogenic? (the headline)
  4. Do those reclassification candidates sit in regions where AlphaFold is
     confident enough for the 3D layer to mean anything?

If (2) fails, stop. Everything downstream assumes AlphaMissense is calibrated.
"""

import json
from pathlib import Path

import duckdb
import httpx
import pandas as pd

import clinvar_regions as cr

DATA = Path(__file__).parent / "data"
DB = Path(__file__).parent / "phase0_scratch.duckdb"

# symbol -> canonical UniProt accession
GENES = {
    "BRCA1": "P38398",   # global pLDDT 41.6 — largely disordered
    "TP53": "P04637",    # global pLDDT 75.1 — mixed
    "PTEN": "P60484",    # global pLDDT 83.0 — well folded
}

AF_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"


# ---------------------------------------------------------------- AlphaFold

def alphafold_entry(acc: str) -> dict:
    """Canonical AlphaFold entry for an accession.

    The API returns one entry per isoform (BRCA1 returns 8+). The canonical
    model is the one whose modelEntityId is exactly AF-{acc}-F1.
    """
    r = httpx.get(AF_API.format(acc=acc), timeout=60)
    r.raise_for_status()
    entries = r.json()
    want = f"AF-{acc}-F1"
    for e in entries:
        if e.get("modelEntityId") == want:
            return e
    raise RuntimeError(f"no canonical entry {want}; got {[e.get('modelEntityId') for e in entries]}")


def fetch_alphamissense(entry: dict, symbol: str) -> pd.DataFrame:
    """Per-protein AlphaMissense predictions in hg38 genomic coordinates."""
    url = entry["amAnnotationsHg38Url"]
    path = DATA / f"am_{symbol}_hg38.csv"
    if not path.exists():
        path.write_bytes(httpx.get(url, timeout=120).content)
    df = pd.read_csv(path)
    # AlphaMissense uses 'chr17'; ClinVar VCF uses '17'
    df["chrom"] = df["CHROM"].str.removeprefix("chr")
    df = df.rename(columns={"POS": "pos", "REF": "ref", "ALT": "alt"})
    # The hg38 file spells classes out ('likely_pathogenic'); the sibling
    # aa-substitutions file uses short codes ('LPath'). Normalise to the short
    # form so downstream filters work against either source.
    canon = {
        "likely_pathogenic": "LPath", "likely_benign": "LBen", "ambiguous": "Amb",
        "LPath": "LPath", "LBen": "LBen", "Amb": "Amb",
    }
    unknown = set(df["am_class"].unique()) - set(canon)
    if unknown:
        raise RuntimeError(f"unrecognised am_class values for {symbol}: {unknown}")
    df["am_class"] = df["am_class"].map(canon)
    df["symbol"] = symbol
    df["aa_pos"] = df["protein_variant"].str.extract(r"^[A-Z](\d+)[A-Z]$").astype("Int64")
    return df[["symbol", "chrom", "pos", "ref", "alt", "transcript_id",
               "protein_variant", "aa_pos", "am_pathogenicity", "am_class"]]


def fetch_plddt(entry: dict, symbol: str) -> pd.DataFrame:
    """Per-residue pLDDT."""
    path = DATA / f"plddt_{symbol}.json"
    if not path.exists():
        path.write_bytes(httpx.get(entry["plddtDocUrl"], timeout=120).content)
    d = json.loads(path.read_text())
    return pd.DataFrame({
        "symbol": symbol,
        "aa_pos": d["residueNumber"],
        "plddt": d["confidenceScore"],
    })


# ------------------------------------------------------------------ ClinVar

def parse_info(info: str) -> dict:
    out = {}
    for field in info.split(";"):
        if "=" in field:
            k, _, v = field.partition("=")
            out[k] = v
    return out


STARS = {
    "practice_guideline": 4,
    "reviewed_by_expert_panel": 3,
    "criteria_provided,_multiple_submitters,_no_conflicts": 2,
    "criteria_provided,_single_submitter": 1,
    "criteria_provided,_conflicting_classifications": 1,
    "criteria_provided,_conflicting_interpretations": 1,
}


def collapse_clnsig(s: str) -> str:
    """Collapse ClinVar significance into decision-relevant buckets.

    Order matters: 'Conflicting' and 'Uncertain' are checked before the
    Pathogenic/Benign membership tests so compound labels land correctly.
    Anything unrecognised falls to OTHER rather than being silently dropped.
    """
    if not s:
        return "OTHER"
    if "Conflicting" in s:
        return "CONFLICTING"
    if "Uncertain" in s or "no_classification" in s:
        return "VUS"
    base = s.split("|")[0]
    if base in {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}:
        return "PATH"
    if base in {"Benign", "Likely_benign", "Benign/Likely_benign"}:
        return "BENIGN"
    return "OTHER"


def fetch_clinvar(index: cr.TabixIndex, symbol: str, am: pd.DataFrame) -> pd.DataFrame:
    """ClinVar records over the gene's coding span, via tabix byte ranges.

    The region is taken from the AlphaMissense file itself rather than
    hardcoded gene coordinates, so the two sides of the join are guaranteed to
    describe the same span. The region can overlap neighbouring genes (NBR2
    overlaps BRCA1), so GENEINFO is still checked per record.
    """
    chrom = am["chrom"].iloc[0]
    beg, end = int(am["pos"].min()), int(am["pos"].max())
    lines = cr.fetch_region(index, chrom, beg, end)

    rows = []
    for line in lines:
        f = line.split("\t")
        if len(f) < 8:
            continue
        c, pos, vid, ref, alt, _q, _fl, info = f[:8]
        d = parse_info(info)
        genes = {g.split(":")[0] for g in d.get("GENEINFO", "").split("|") if g}
        if symbol not in genes:
            continue
        for a in alt.split(","):
            rows.append({
                "symbol": symbol,
                "chrom": c,
                "pos": int(pos),
                "ref": ref,
                "alt": a,
                "clnsig_raw": d.get("CLNSIG", ""),
                "clnsig": collapse_clnsig(d.get("CLNSIG", "")),
                "revstat": d.get("CLNREVSTAT", ""),
                "stars": STARS.get(d.get("CLNREVSTAT", ""), 0),
                "mc": d.get("MC", ""),
                "is_missense": "missense_variant" in d.get("MC", ""),
                "is_snv": len(ref) == 1 and len(a) == 1,
                "clinvar_id": vid,
            })
    print(f"  {symbol:6s} {chrom}:{beg:,}-{end:,}  {len(lines):>6,} region records "
          f"-> {len(rows):>6,} for this gene")
    return pd.DataFrame(rows)


# --------------------------------------------------------------------- main

def main() -> None:
    DATA.mkdir(exist_ok=True)

    print("[1/4] AlphaFold + AlphaMissense")
    am_parts, plddt_parts, meta = [], [], []
    for symbol, acc in GENES.items():
        entry = alphafold_entry(acc)
        am = fetch_alphamissense(entry, symbol)
        pl = fetch_plddt(entry, symbol)
        am_parts.append(am)
        plddt_parts.append(pl)
        meta.append({
            "symbol": symbol, "acc": acc,
            "global_plddt": entry["globalMetricValue"],
            "frac_very_low": entry["fractionPlddtVeryLow"],
            "n_am": len(am), "n_residues": len(pl),
            "transcripts": am["transcript_id"].nunique(),
        })
        print(f"  {symbol:6s} {acc}  AM={len(am):>6,}  residues={len(pl):>5,}  "
              f"pLDDT={entry['globalMetricValue']:.1f}")
    am_df = pd.concat(am_parts, ignore_index=True)
    plddt_df = pd.concat(plddt_parts, ignore_index=True)
    meta_df = pd.DataFrame(meta)

    print("\n[2/4] ClinVar via tabix byte ranges")
    index = cr.load_index()
    cv_df = pd.concat(
        [fetch_clinvar(index, sym, am_df[am_df.symbol == sym]) for sym in GENES],
        ignore_index=True,
    )

    print("\n[3/4] Join in DuckDB")
    con = duckdb.connect(str(DB))
    con.register("am", am_df)
    con.register("cv", cv_df)
    con.register("plddt", plddt_df)
    con.register("meta", meta_df)

    con.execute("""
        CREATE OR REPLACE TABLE joined AS
        SELECT
            cv.symbol, cv.chrom, cv.pos, cv.ref, cv.alt,
            cv.clinvar_id, cv.clnsig, cv.clnsig_raw, cv.stars, cv.is_missense, cv.is_snv,
            am.protein_variant, am.aa_pos, am.am_pathogenicity, am.am_class,
            p.plddt,
            (am.protein_variant IS NOT NULL) AS matched
        FROM cv
        LEFT JOIN am
          ON  cv.chrom = am.chrom AND cv.pos = am.pos
          AND cv.ref = am.ref     AND cv.alt = am.alt
        LEFT JOIN plddt p
          ON  p.symbol = am.symbol AND p.aa_pos = am.aa_pos
    """)

    def q(sql):
        return con.execute(sql).df()

    print("\n" + "=" * 78)
    print("RESULT 1 — JOIN COVERAGE (ClinVar missense SNVs -> AlphaMissense)")
    print("=" * 78)
    print(q("""
        SELECT symbol,
               COUNT(*) AS clinvar_missense_snv,
               SUM(matched::INT) AS matched,
               ROUND(100.0 * SUM(matched::INT) / COUNT(*), 1) AS pct
        FROM joined WHERE is_missense AND is_snv
        GROUP BY symbol ORDER BY symbol
    """).to_string(index=False))

    print("\n" + "=" * 78)
    print("RESULT 2 — CALIBRATION: does AlphaMissense agree where ClinVar is confident?")
    print("           (>=1 star, missense SNVs, matched only)  *** GO / NO-GO ***")
    print("=" * 78)
    print(q("""
        SELECT clnsig, am_class, COUNT(*) AS n
        FROM joined
        WHERE matched AND is_missense AND is_snv AND stars >= 1
          AND clnsig IN ('PATH','BENIGN')
        GROUP BY clnsig, am_class ORDER BY clnsig, am_class
    """).to_string(index=False))
    print("\n  agreement rates:")
    print(q("""
        SELECT clnsig,
               COUNT(*) AS n,
               ROUND(100.0 * SUM(CASE WHEN clnsig='PATH'   AND am_class='LPath' THEN 1
                                      WHEN clnsig='BENIGN' AND am_class='LBen'  THEN 1
                                      ELSE 0 END) / COUNT(*), 1) AS pct_agree,
               ROUND(100.0 * SUM(CASE WHEN clnsig='PATH'   AND am_class='LBen'  THEN 1
                                      WHEN clnsig='BENIGN' AND am_class='LPath' THEN 1
                                      ELSE 0 END) / COUNT(*), 1) AS pct_contradict,
               ROUND(100.0 * SUM((am_class='Amb')::INT) / COUNT(*), 1) AS pct_ambiguous
        FROM joined
        WHERE matched AND is_missense AND is_snv AND stars >= 1
          AND clnsig IN ('PATH','BENIGN')
        GROUP BY clnsig ORDER BY clnsig
    """).to_string(index=False))

    print("\n" + "=" * 78)
    print("RESULT 3 — THE HEADLINE: VUS / conflicting that AlphaMissense calls pathogenic")
    print("=" * 78)
    print(q("""
        SELECT symbol, clnsig,
               COUNT(*) AS n_uncertain,
               SUM((am_class='LPath')::INT) AS reclass_upgrade,
               SUM((am_class='LBen')::INT)  AS reclass_downgrade,
               SUM((am_class='Amb')::INT)   AS remains_uncertain
        FROM joined
        WHERE matched AND is_missense AND is_snv AND clnsig IN ('VUS','CONFLICTING')
        GROUP BY symbol, clnsig ORDER BY symbol, clnsig
    """).to_string(index=False))

    print("\n" + "=" * 78)
    print("RESULT 4 — CONFIDENCE GATE: can the 3D layer support these calls?")
    print("           pLDDT band of reclassification-upgrade candidates")
    print("=" * 78)
    print(q("""
        SELECT symbol,
               CASE WHEN plddt < 50 THEN 'a_very_low (<50)'
                    WHEN plddt < 70 THEN 'b_low (50-70)'
                    WHEN plddt < 90 THEN 'c_confident (70-90)'
                    ELSE 'd_very_high (>90)' END AS band,
               COUNT(*) AS n
        FROM joined
        WHERE matched AND is_missense AND is_snv
          AND clnsig IN ('VUS','CONFLICTING') AND am_class='LPath'
        GROUP BY symbol, band ORDER BY symbol, band
    """).to_string(index=False))

    print("\n" + "=" * 78)
    print("RESULT 5 — UNMATCHED DIAGNOSIS (why did ClinVar missense SNVs miss?)")
    print("=" * 78)
    print(q("""
        SELECT symbol, COUNT(*) AS unmatched
        FROM joined WHERE is_missense AND is_snv AND NOT matched
        GROUP BY symbol ORDER BY symbol
    """).to_string(index=False))
    unmatched = q("""
        SELECT symbol, chrom, pos, ref, alt, clnsig, clinvar_id
        FROM joined WHERE is_missense AND is_snv AND NOT matched
        LIMIT 8
    """)
    if len(unmatched):
        print("\n  sample:")
        print(unmatched.to_string(index=False))

    print("\n" + "=" * 78)
    print("RESULT 6 — DISCORDANT (model contradicts a confident clinical call)")
    print("=" * 78)
    print(q("""
        SELECT symbol, clnsig, am_class, stars, COUNT(*) AS n
        FROM joined
        WHERE matched AND is_missense AND is_snv AND stars >= 2
          AND ((clnsig='PATH' AND am_class='LBen') OR (clnsig='BENIGN' AND am_class='LPath'))
        GROUP BY symbol, clnsig, am_class, stars ORDER BY n DESC
    """).to_string(index=False))

    con.execute("CREATE OR REPLACE TABLE gene_meta AS SELECT * FROM meta")
    con.close()
    print(f"\nPersisted to {DB}")


if __name__ == "__main__":
    main()

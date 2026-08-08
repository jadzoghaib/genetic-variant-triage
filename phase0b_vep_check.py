# /// script
# requires-python = ">=3.11"
# dependencies = ["duckdb", "pandas", "httpx"]
# ///
"""
Phase 0b, part 2 — settle the mechanism behind the 74 unmatched records.

Section C of the first diagnostic suggested these are not coordinate defects:
at BRCA1 chr17:43047679 ClinVar asks for G>A while AlphaMissense carries G>T
(Q1811K) and G>C (Q1811E) at the same position with the same REF. On the minus
strand that makes the missing G>A a coding C>T, turning CAA (Gln) into TAA
(stop) — a nonsense variant, which AlphaMissense correctly excludes.

Ensembl VEP gives per-transcript consequences, so it can confirm directly
whether each unmatched position is (a) non-missense on MANE Select, or
(b) missense only on a non-MANE transcript.
"""

from pathlib import Path

import duckdb
import httpx
import pandas as pd

DB = Path(__file__).parent / "phase0_scratch.duckdb"
VEP = "https://rest.ensembl.org/vep/human/region"

con = duckdb.connect(str(DB))
unmatched = con.execute("""
    SELECT symbol, chrom, pos, ref, alt, clnsig, clinvar_id
    FROM joined WHERE is_missense AND is_snv AND NOT matched
    ORDER BY symbol, pos
""").df()
con.close()

# Sample across both buckets and both genes rather than the first N.
sample = pd.concat([
    unmatched[unmatched.symbol == "BRCA1"].head(6),
    unmatched[unmatched.symbol == "TP53"].head(4),
])

print(f"querying Ensembl VEP for {len(sample)} of {len(unmatched)} unmatched records\n")

rows = []
with httpx.Client(timeout=60) as client:
    for _, v in sample.iterrows():
        region = f"{v.chrom}:{v.pos}-{v.pos}:1/{v.alt}"
        r = client.get(f"{VEP}/{region}",
                       params={"content-type": "application/json", "mane": 1, "canonical": 1})
        if r.status_code != 200:
            rows.append({"symbol": v.symbol, "pos": v.pos, "sub": f"{v.ref}>{v.alt}",
                         "mane_consequence": f"HTTP {r.status_code}", "n_missense_tx": None})
            continue
        data = r.json()
        tcs = data[0].get("transcript_consequences", []) if data else []
        mane = [t for t in tcs if t.get("mane_select")]
        n_missense = sum(1 for t in tcs if "missense_variant" in t.get("consequence_terms", []))
        rows.append({
            "symbol": v.symbol,
            "pos": v.pos,
            "sub": f"{v.ref}>{v.alt}",
            "mane_tx": mane[0].get("transcript_id") if mane else "(none)",
            "mane_consequence": ",".join(mane[0].get("consequence_terms", [])) if mane else "(no MANE tx)",
            "n_tx_total": len(tcs),
            "n_tx_missense": n_missense,
        })

df = pd.DataFrame(rows)
print(df.to_string(index=False))

print("\n" + "=" * 78)
print("VERDICT")
print("=" * 78)
non_missense_on_mane = df[~df["mane_consequence"].str.contains("missense", na=False)]
print(f"  records whose MANE Select consequence is NOT missense: "
      f"{len(non_missense_on_mane)} / {len(df)}")
print(f"  of those, how many are missense on some OTHER transcript: "
      f"{(non_missense_on_mane['n_tx_missense'] > 0).sum()}")
print("\n  consequence breakdown on MANE Select:")
print(df["mane_consequence"].value_counts().to_string())
